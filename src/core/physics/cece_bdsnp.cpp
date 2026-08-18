/**
 * @file cece_bdsnp.cpp
 * @brief Standalone BDSNP soil NO physics module implementation.
 *
 * Implements the Berkeley-Dalhousie Soil NO Parameterization (BDSNP) with
 * Yienger & Levy (1995) fallback. Replaces the existing SoilNoxScheme
 * ("soil_nox") registration with a more comprehensive soil NO model.
 *
 * Two algorithms are supported, selectable via `soil_no_method` YAML key:
 *   - "yl95": Yienger & Levy (1995) — soil temperature response, soil moisture
 *     pulse, canopy reduction factor (identical to existing SoilNoxScheme)
 *   - "bdsnp" (default): biome-specific base emission factors, soil moisture
 *     dependence, nitrogen deposition fertilization, canopy reduction
 *
 * Both modes set soil NO to zero when soil temperature < 0°C (273.15 K).
 * Writes to export field "soil_nox_emissions" for consumption by MEGAN3.
 *
 * References:
 * - Yienger, J.J. and H. Levy II (1995), JGR, 100(D6), 11447-11464.
 * - Hudman et al. (2012), BDSNP parameterization.
 *
 * @author CECE Team
 * @date 2024
 */

#include "cece/physics/cece_bdsnp.hpp"

#include <Kokkos_Array.hpp>
#include <Kokkos_Core.hpp>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "cece/cece_physics_factory.hpp"

namespace cece {

/// @brief Self-registration for the BDSNP soil NO emission scheme.
static PhysicsRegistration<BdsnpScheme> reg("bdsnp");

// ============================================================================
// YL95 inline helper functions (matching existing SoilNoxScheme exactly)
// ============================================================================

/**
 * @brief Calculate temperature-dependent soil NOx emission factor (YL95).
 *
 * Computes the exponential temperature response of soil microbial activity.
 * Returns 0 when soil temperature is at or below freezing.
 *
 * @param tc Soil temperature [°C]
 * @param tc_max Maximum temperature for emission calculation [°C]
 * @param exp_coeff Exponential temperature coefficient [1/°C]
 * @return Temperature-dependent emission factor (dimensionless)
 */
KOKKOS_INLINE_FUNCTION
double bdsnp_soil_temp_term(double tc, double tc_max, double exp_coeff) {
    if (tc <= 0.0) {
        return 0.0;  // No emission below freezing
    }
    return std::exp(exp_coeff * std::min(tc_max, tc));
}

/**
 * @brief Calculate soil moisture-dependent NOx emission factor (YL95).
 *
 * Computes the water-filled pore space (WFPS) effect on soil NOx emissions
 * using a Poisson-like response function.
 *
 * @param gw Water-filled pore space fraction [0-1]
 * @param wet_c1 Moisture response coefficient 1
 * @param wet_c2 Moisture response coefficient 2 [1/WFPS²]
 * @return Moisture-dependent emission factor (dimensionless)
 */
KOKKOS_INLINE_FUNCTION
double bdsnp_soil_wet_term(double gw, double wet_c1, double wet_c2) {
    return wet_c1 * gw * std::exp(wet_c2 * gw * gw);
}

// ============================================================================
// BDSNP-specific inline helper functions
// ============================================================================

/**
 * @brief Compute BDSNP soil moisture dependence factor.
 *
 * Piecewise linear response to soil moisture for the BDSNP algorithm.
 *
 * @param soil_moisture Soil moisture fraction [0-1]
 * @return Soil moisture dependence factor (dimensionless, non-negative)
 */
KOKKOS_INLINE_FUNCTION
double bdsnp_moisture_factor(double soil_moisture) {
    // Piecewise linear: ramp up from 0 at SM=0 to 1 at SM=0.3,
    // then decrease linearly to 0.5 at SM=1.0
    if (soil_moisture <= 0.0) {
        return 0.0;
    }
    if (soil_moisture <= 0.3) {
        return soil_moisture / 0.3;
    }
    // Linear decrease from 1.0 at SM=0.3 to 0.5 at SM=1.0
    return 1.0 - 0.5 * (soil_moisture - 0.3) / 0.7;
}

/**
 * @brief Compute nitrogen deposition fertilization factor.
 *
 * Enhances soil NO emissions based on nitrogen deposition rates.
 *
 * @param ndep Nitrogen deposition rate [kg N/m²/s]
 * @param fert_emission_factor Fertilizer emission factor scaling
 * @param wet_dep_scaling Wet deposition scaling factor
 * @param dry_dep_scaling Dry deposition scaling factor
 * @return Fertilization enhancement factor (dimensionless, >= 1.0)
 */
KOKKOS_INLINE_FUNCTION
double bdsnp_ndep_factor(double ndep, double fert_emission_factor, double wet_dep_scaling, double dry_dep_scaling) {
    // Combined wet + dry deposition contribution
    double dep_contribution = ndep * (wet_dep_scaling + dry_dep_scaling);
    return 1.0 + fert_emission_factor * dep_contribution;
}

/**
 * @brief Compute canopy reduction factor for soil NO.
 *
 * Reduces soil NO emissions based on canopy uptake, parameterized by LAI.
 *
 * @param lai Leaf area index [m²/m²]
 * @return Canopy reduction factor [0-1]
 */
KOKKOS_INLINE_FUNCTION
double bdsnp_canopy_reduction(double lai) {
    if (lai <= 0.0) {
        return 1.0;
    }
    // Exponential reduction with LAI (Beer-Lambert-like)
    return std::exp(-0.24 * lai);
}

// ============================================================================
// Initialize
// ============================================================================

void BdsnpScheme::Initialize(const YAML::Node& config, CeceDiagnosticManager* diag_manager) {
    // Call base class to parse input_mapping, output_mapping, diagnostics
    BasePhysicsScheme::Initialize(config, diag_manager);

    // ---- Read soil_no_method (default "bdsnp", fallback "yl95") ----
    soil_no_method_ = "bdsnp";
    if (config["soil_no_method"]) {
        soil_no_method_ = config["soil_no_method"].as<std::string>();
    }
    if (soil_no_method_ != "bdsnp" && soil_no_method_ != "yl95" && soil_no_method_ != "hemco_3_12_1") {
        std::cout << "BdsnpScheme: WARNING - Unknown soil_no_method '" << soil_no_method_ << "', falling back to 'bdsnp'\n";
        soil_no_method_ = "bdsnp";
    }

    use_soil_temperature_ = false;
    if (config["use_soil_temperature"]) {
        use_soil_temperature_ = config["use_soil_temperature"].as<bool>();
    }

    // ---- Read YL95 parameters ----
    a_biome_wet_ = 0.5;
    tc_max_ = 30.0;
    exp_coeff_ = 0.103;
    wet_c1_ = 5.5;
    wet_c2_ = -5.55;

    if (config["biome_coefficient_wet"]) {
        a_biome_wet_ = config["biome_coefficient_wet"].as<double>();
    } else if (config["a_biome_wet"]) {
        a_biome_wet_ = config["a_biome_wet"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default a_biome_wet = " << a_biome_wet_ << "\n";
    }

    if (config["temp_limit"]) {
        tc_max_ = config["temp_limit"].as<double>();
    } else if (config["tc_max"]) {
        tc_max_ = config["tc_max"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default tc_max = " << tc_max_ << "\n";
    }

    if (config["temp_exp_coeff"]) {
        exp_coeff_ = config["temp_exp_coeff"].as<double>();
    } else if (config["exp_coeff"]) {
        exp_coeff_ = config["exp_coeff"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default exp_coeff = " << exp_coeff_ << "\n";
    }

    if (config["wet_coeff_1"]) {
        wet_c1_ = config["wet_coeff_1"].as<double>();
    } else if (config["wet_c1"]) {
        wet_c1_ = config["wet_c1"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default wet_c1 = " << wet_c1_ << "\n";
    }

    if (config["wet_coeff_2"]) {
        wet_c2_ = config["wet_coeff_2"].as<double>();
    } else if (config["wet_c2"]) {
        wet_c2_ = config["wet_c2"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default wet_c2 = " << wet_c2_ << "\n";
    }

    // ---- Read BDSNP parameters ----
    fert_emission_factor_ = 1.0;
    wet_dep_scaling_ = 1.0;
    dry_dep_scaling_ = 1.0;
    pulse_decay_constant_ = 0.5;

    if (config["fert_emission_factor"]) {
        fert_emission_factor_ = config["fert_emission_factor"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default fert_emission_factor = " << fert_emission_factor_ << "\n";
    }

    if (config["wet_dep_scaling"]) {
        wet_dep_scaling_ = config["wet_dep_scaling"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default wet_dep_scaling = " << wet_dep_scaling_ << "\n";
    }

    if (config["dry_dep_scaling"]) {
        dry_dep_scaling_ = config["dry_dep_scaling"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default dry_dep_scaling = " << dry_dep_scaling_ << "\n";
    }

    if (config["pulse_decay_constant"]) {
        pulse_decay_constant_ = config["pulse_decay_constant"].as<double>();
    } else {
        std::cout << "BdsnpScheme: Using default pulse_decay_constant = " << pulse_decay_constant_ << "\n";
    }

    std::cout << "BdsnpScheme: Initialized with soil_no_method='" << soil_no_method_ << "'\n";
}

// ============================================================================
// Run
// ============================================================================

void BdsnpScheme::Run(CeceImportState& import_state, CeceExportState& export_state) {
    if (soil_no_method_ == "hemco_3_12_1") {
        // Exact HEMCO 3.12.1 stateless cell arithmetic. Stateful quantities
        // and upstream canopy/deposited-N calculations are explicit inputs so
        // they can first be checked against a frozen HEMCO global reference.
        auto temperature = ResolveImport(use_soil_temperature_ ? "soil_temperature" : "surface_temperature", import_state);
        auto soil_moisture = ResolveImport("soil_moisture", import_state);
        auto land_fractions = ResolveImport("soilnox_land_fractions", import_state);
        auto arid_fraction = ResolveImport("soilnox_arid_fraction", import_state);
        auto nonarid_fraction = ResolveImport("soilnox_nonarid_fraction", import_state);
        auto lai = ResolveImport("leaf_area_index", import_state);
        auto canopy_nox = ResolveImport("soilnox_canopy_nox", import_state);
        auto wind_speed_squared = ResolveImport("wind_speed_squared", import_state);
        auto solar_zenith_cosine = ResolveImport("solar_zenith_cosine", import_state);
        auto soil_fertilizer = ResolveImport("soil_fertilizer", import_state);
        auto deposited_nitrogen = ResolveImport("deposited_nitrogen", import_state);
        auto pulse_factor = ResolveImport("soilnox_pulse_factor", import_state);
        auto soil_nox = ResolveExport("soil_nox_emissions", export_state);
        auto fertilizer_nox = ResolveExport("soil_nox_fertilizer_emissions", export_state);

        if (temperature.data() == nullptr || soil_moisture.data() == nullptr || land_fractions.data() == nullptr || arid_fraction.data() == nullptr ||
            nonarid_fraction.data() == nullptr || lai.data() == nullptr || canopy_nox.data() == nullptr || wind_speed_squared.data() == nullptr ||
            solar_zenith_cosine.data() == nullptr || soil_fertilizer.data() == nullptr || deposited_nitrogen.data() == nullptr ||
            pulse_factor.data() == nullptr || soil_nox.data() == nullptr || fertilizer_nox.data() == nullptr) {
            throw std::runtime_error(
                "BdsnpScheme hemco_3_12_1 mode requires surface/soil temperature, soil_moisture, "
                "soilnox_land_fractions[24], soilnox_arid_fraction, soilnox_nonarid_fraction, "
                "leaf_area_index, soilnox_canopy_nox[24], wind_speed_squared, solar_zenith_cosine, "
                "soil_fertilizer, deposited_nitrogen, soilnox_pulse_factor, and both SoilNOx exports");
        }

        const int nx = static_cast<int>(soil_nox.extent(0));
        const int ny = static_cast<int>(soil_nox.extent(1));
        auto require_scalar_shape = [nx, ny](const auto& view, const char* name) {
            if (view.extent(0) != nx || view.extent(1) != ny || view.extent(2) < 1) {
                throw std::runtime_error(std::string("BdsnpScheme hemco_3_12_1 field has incompatible shape: ") + name);
            }
        };
        require_scalar_shape(temperature, "temperature");
        require_scalar_shape(soil_moisture, "soil_moisture");
        require_scalar_shape(arid_fraction, "soilnox_arid_fraction");
        require_scalar_shape(nonarid_fraction, "soilnox_nonarid_fraction");
        require_scalar_shape(lai, "leaf_area_index");
        require_scalar_shape(wind_speed_squared, "wind_speed_squared");
        require_scalar_shape(solar_zenith_cosine, "solar_zenith_cosine");
        require_scalar_shape(soil_fertilizer, "soil_fertilizer");
        require_scalar_shape(deposited_nitrogen, "deposited_nitrogen");
        require_scalar_shape(pulse_factor, "soilnox_pulse_factor");
        require_scalar_shape(fertilizer_nox, "soil_nox_fertilizer_emissions");
        if (land_fractions.extent(0) != nx || land_fractions.extent(1) != ny || land_fractions.extent(2) != 24 || canopy_nox.extent(0) != nx ||
            canopy_nox.extent(1) != ny || canopy_nox.extent(2) != 24) {
            throw std::runtime_error("BdsnpScheme hemco_3_12_1 requires exactly 24 land-fraction and canopy-NOx layers");
        }

        const Kokkos::Array<double, 24> a_biome = {
            0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 0.84, 0.84, 0.24,
            0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 1.66, 0.08, 0.44, 0.57, 0.57, 0.57};
        const Kokkos::Array<double, 24> soil_ta = {
            0.00, 0.92, 0.00, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66,
            0.66, 0.66, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 1.03, 1.03, 1.03};
        const Kokkos::Array<double, 24> soil_tb = {
            0.00, 4.40, 0.00, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80,
            8.80, 8.80, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 2.90, 2.90, 2.90};
        const Kokkos::Array<double, 24> soil_exc = {
            0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00, 1.00, 1.00,
            1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 2.00, 0.10, 2.00};
        const bool use_soil_temperature = use_soil_temperature_;
        constexpr double unit_conversion = 1.0e-12 / 14.0 * 30.0;
        constexpr double fertilizer_scale = 0.0068;
        constexpr double seconds_per_year = 3.1536e7;
        constexpr double soil_exp_coefficient = static_cast<double>(0.103F);
        constexpr double cubic_3 = static_cast<double>(-0.009F);
        constexpr double cubic_2 = static_cast<double>(0.837F);
        constexpr double cubic_1 = static_cast<double>(-22.527F);
        constexpr double cubic_0 = static_cast<double>(196.149F);

        Kokkos::parallel_for(
            "BdsnpKernel_HEMCO_3_12_1", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
            KOKKOS_LAMBDA(int i, int j) {
                soil_nox(i, j, 0) = 0.0;
                fertilizer_nox(i, j, 0) = 0.0;

                // Match HEMCO's exact no-soil gate before updating any state.
                if (land_fractions(i, j, 0) == 1.0) {
                    return;
                }

                const double temperature_c = temperature(i, j, 0) - 273.15;
                const double gwet = soil_moisture(i, j, 0);
                const double arid = arid_fraction(i, j, 0);
                const double nonarid = nonarid_fraction(i, j, 0);
                const double wetness = (arid >= nonarid && arid > 0.0) ? 8.24 * gwet * std::exp(-12.5 * gwet * gwet)
                                                                          : 5.5 * gwet * std::exp(-5.55 * gwet * gwet);
                const double fertilizer = (soil_fertilizer(i, j, 0) + deposited_nitrogen(i, j, 0)) / seconds_per_year * fertilizer_scale;
                const double leaf_area_index = lai(i, j, 0);
                const double wind_squared = wind_speed_squared(i, j, 0);
                const double sun_cosine = solar_zenith_cosine(i, j, 0);
                const double pulse = pulse_factor(i, j, 0);

                double total = 0.0;
                double added_n_total = 0.0;
                for (int k = 0; k < 24; ++k) {
                    double adjusted_temperature = temperature_c;
                    if (!use_soil_temperature) {
                        adjusted_temperature = gwet < 0.3 ? adjusted_temperature + 5.0 : soil_ta[k] * adjusted_temperature + soil_tb[k];
                    }

                    double temperature_term = 0.0;
                    if (adjusted_temperature > 0.0) {
                        if (!use_soil_temperature) {
                            adjusted_temperature = adjusted_temperature >= 30.0 ? 30.0 : adjusted_temperature;
                            temperature_term = std::exp(0.103 * adjusted_temperature);
                        } else {
                            adjusted_temperature = adjusted_temperature >= 40.0 ? 40.0 : adjusted_temperature;
                            temperature_term = adjusted_temperature <= 20.0
                                                   ? std::exp(soil_exp_coefficient * adjusted_temperature)
                                                   : cubic_3 * adjusted_temperature * adjusted_temperature * adjusted_temperature +
                                                         cubic_2 * adjusted_temperature * adjusted_temperature + cubic_1 * adjusted_temperature +
                                                         cubic_0;
                        }
                    }

                    double canopy_reduction = 0.0;
                    const double canopy = canopy_nox(i, j, k);
                    if (leaf_area_index > 0.0 && canopy > 0.0) {
                        double ventilation = sun_cosine > 0.0 ? 1.0e-2 : 0.2e-2;
                        ventilation *= std::sqrt(wind_squared / 9.0 * 7.0 / leaf_area_index) * (soil_exc[20] / soil_exc[k]);
                        canopy_reduction = canopy / (canopy + ventilation);
                    }

                    const double common = temperature_term * wetness * pulse * land_fractions(i, j, k) * (1.0 - canopy_reduction);
                    total += (a_biome[k] * unit_conversion + fertilizer) * common;
                    added_n_total += fertilizer * common;
                }

                soil_nox(i, j, 0) = total > 0.0 ? total : 0.0;
                fertilizer_nox(i, j, 0) = added_n_total;
            });

        Kokkos::fence();
        MarkModified("soil_nox_emissions", export_state);
        MarkModified("soil_nox_fertilizer_emissions", export_state);
        return;
    }

    // ---- Resolve import fields ----
    auto soil_temp = ResolveImport("soil_temperature", import_state);
    auto soil_moisture = ResolveImport("soil_moisture", import_state);
    auto soil_nox = ResolveExport("soil_nox_emissions", export_state);

    // Fail loudly: a silent return creates a valid-looking all-zero output.
    if (soil_temp.data() == nullptr || soil_moisture.data() == nullptr || soil_nox.data() == nullptr) {
        throw std::runtime_error(
            "BdsnpScheme requires imports 'soil_temperature' and "
            "'soil_moisture' plus export 'soil_nox_emissions'");
    }

    int nx = static_cast<int>(soil_nox.extent(0));
    int ny = static_cast<int>(soil_nox.extent(1));

    if (soil_no_method_ == "yl95") {
        // ================================================================
        // YL95 mode: identical algorithm to existing SoilNoxScheme
        // ================================================================
        const double MW_NO = 30.0;
        const double UNITCONV = 1.0e-12 / 14.0 * MW_NO;  // ng N -> kg NO
        double a_biome_wet = a_biome_wet_;
        double tc_max = tc_max_;
        double exp_coeff = exp_coeff_;
        double wet_c1 = wet_c1_;
        double wet_c2 = wet_c2_;

        Kokkos::parallel_for(
            "BdsnpKernel_YL95", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
                double tc = soil_temp(i, j, 0) - 273.15;
                double gw = soil_moisture(i, j, 0);

                // Set to zero if soil temp < 0°C
                if (tc <= 0.0) {
                    soil_nox(i, j, 0) = 0.0;
                    return;
                }

                double t_term = bdsnp_soil_temp_term(tc, tc_max, exp_coeff);
                double w_term = bdsnp_soil_wet_term(gw, wet_c1, wet_c2);

                // Pulse factor placeholder
                double pulse = 1.0;

                // Total emission [kg NO/m2/s]
                double emiss = a_biome_wet * UNITCONV * t_term * w_term * pulse;
                soil_nox(i, j, 0) = emiss;
            });
    } else {
        // ================================================================
        // BDSNP mode: biome-specific emission with N-dep fertilization
        // ================================================================

        // Resolve additional BDSNP-specific import fields
        auto ndep = ResolveImport("nitrogen_deposition", import_state);
        auto land_use = ResolveImport("land_use_type", import_state);
        auto lai = ResolveImport("leaf_area_index", import_state);
        auto biome_ef = ResolveImport("biome_emission_factors", import_state);

        // BDSNP can proceed with partial fields — use defaults where missing
        bool has_ndep = (ndep.data() != nullptr);
        bool has_land_use = (land_use.data() != nullptr);
        bool has_lai = (lai.data() != nullptr);
        bool has_biome_ef = (biome_ef.data() != nullptr);

        double fert_ef = fert_emission_factor_;
        double wet_dep_s = wet_dep_scaling_;
        double dry_dep_s = dry_dep_scaling_;
        double pulse_decay = pulse_decay_constant_;

        const double MW_NO = 30.0;
        const double UNITCONV = 1.0e-12 / 14.0 * MW_NO;  // ng N -> kg NO

        Kokkos::parallel_for(
            "BdsnpKernel_BDSNP", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
            KOKKOS_LAMBDA(int i, int j) {
                double tc = soil_temp(i, j, 0) - 273.15;

                // Set to zero if soil temp < 0°C
                if (tc <= 0.0) {
                    soil_nox(i, j, 0) = 0.0;
                    return;
                }

                double sm = soil_moisture(i, j, 0);

                // Biome-specific base emission factor
                double base_ef = has_biome_ef ? biome_ef(i, j, 0) : 1.0;

                // Temperature response (exponential, same form as YL95)
                double t_response = std::exp(0.103 * std::min(30.0, tc));

                // Soil moisture dependence (BDSNP piecewise)
                double sm_factor = bdsnp_moisture_factor(sm);

                // Nitrogen deposition fertilization
                double ndep_val = has_ndep ? ndep(i, j, 0) : 0.0;
                double fert_factor = bdsnp_ndep_factor(ndep_val, fert_ef, wet_dep_s, dry_dep_s);

                // Canopy reduction
                double lai_val = has_lai ? lai(i, j, 0) : 0.0;
                double canopy_red = bdsnp_canopy_reduction(lai_val);

                // Pulse decay (simplified — stateless placeholder)
                double pulse = std::exp(-pulse_decay * 0.0);  // No prior rain info

                // Total BDSNP emission [kg NO/m2/s]
                double emiss = base_ef * UNITCONV * t_response * sm_factor * fert_factor * canopy_red * pulse;
                soil_nox(i, j, 0) = emiss;
            });
    }

    Kokkos::fence();
    MarkModified("soil_nox_emissions", export_state);
}

}  // namespace cece
