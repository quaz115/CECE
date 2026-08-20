#ifndef CECE_HEMCO_SOILNOX_STATELESS_HPP
#define CECE_HEMCO_SOILNOX_STATELESS_HPP

/**
 * @file hemco_soilnox_stateless.hpp
 * @brief Host-scalar HEMCO 3.12.1 SoilNOx equations for reference parity.
 *
 * This layer contains the host-scalar oracle used by reference tests. The
 * production Kokkos integration is implemented by BdsnpScheme. Persistent
 * pulse and deposited-N reservoir evolution remain outside this stateless layer.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cece::hemco_soilnox::v3_12_1 {

/** Number of MODIS/Koppen biome classes in HEMCO SoilNOx. */
inline constexpr std::size_t kBiomeCount = 24;

/** HEMCO 3.12.1 wet-biome background coefficients [ng N m-2 s-1]. */
inline constexpr std::array<double, kBiomeCount> kABiome = {
    0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 0.84, 0.84, 0.24, 0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 1.66, 0.08, 0.44, 0.57, 0.57, 0.57,
};

/** HEMCO 3.12.1 air-to-soil temperature slope by biome. */
inline constexpr std::array<double, kBiomeCount> kSoilTa = {
    0.00, 0.92, 0.00, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84, 1.03, 1.03, 1.03,
};

/** HEMCO 3.12.1 air-to-soil temperature intercept by biome [deg C]. */
inline constexpr std::array<double, kBiomeCount> kSoilTb = {
    0.00, 4.40, 0.00, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60, 2.90, 2.90, 2.90,
};

/** HEMCO 3.12.1 canopy wind-extinction coefficients by biome. */
inline constexpr std::array<double, kBiomeCount> kSoilExc = {
    0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00, 1.00, 1.00, 1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 2.00, 0.10, 2.00,
};

inline constexpr double kReferenceFertilizerScale = 0.0068;
inline constexpr double kSecondsPerYear = 3.1536e7;

/** Inputs for one stateless, one-biome HEMCO SoilNOx contribution. */
struct StatelessContributionInput {
    bool use_soil_temperature = false;
    int biome = 1;  // HEMCO uses one-based biome indices.
    double temperature_c = 0.0;
    double gwet = 0.0;
    double arid_fraction = 0.0;
    double nonarid_fraction = 0.0;
    double leaf_area_index = 0.0;
    double canopy_nox = 0.0;
    double wind_speed_squared = 0.0;
    double solar_zenith_cosine = 0.0;
    double soil_fertilizer = 0.0;
    double deposited_nitrogen = 0.0;
    double molecular_weight_no_g = 30.0;
};

/** Scalar terms emitted by the pinned HEMCO reference-vector generator. */
struct StatelessContributionResult {
    double temperature_term = 0.0;
    double wetness_term = 0.0;
    double canopy_reduction_fraction = 0.0;
    double fertilizer_unscaled = 0.0;
    double fertilizer_scaled = 0.0;
    double unit_conversion = 0.0;
    double cell_flux = 0.0;
};

inline std::size_t BiomeIndex(int biome) {
    if (biome < 1 || biome > static_cast<int>(kBiomeCount)) {
        throw std::out_of_range("HEMCO SoilNOx biome index must be in [1, 24]");
    }
    return static_cast<std::size_t>(biome - 1);
}

/** Convert ng N to kg NO using the configured NO molecular weight. */
inline double UnitConversion(double molecular_weight_no_g) {
    return 1.0e-12 / 14.0 * molecular_weight_no_g;
}

/** Exact HEMCO 3.12.1 SoilTemp scalar response. */
inline double SoilTemperatureTerm(bool use_soil_temperature, int biome, double temperature_c, double gwet) {
    const std::size_t index = BiomeIndex(biome);
    double temperature = temperature_c;

    if (!use_soil_temperature) {
        if (gwet < 0.3) {
            temperature += 5.0;
        } else {
            temperature = kSoilTa[index] * temperature + kSoilTb[index];
        }
    }

    if (temperature <= 0.0) {
        return 0.0;
    }

    if (!use_soil_temperature) {
        if (temperature >= 30.0) {
            temperature = 30.0;
        }
        return std::exp(0.103 * temperature);
    }

    if (temperature >= 40.0) {
        temperature = 40.0;
    }

    // These HEMCO source literals have default-real precision. Preserve the
    // float32 values promoted to double for numerical parity with USE_REAL8.
    constexpr double kExpCoefficient = static_cast<double>(0.103F);
    constexpr double kCubic3 = static_cast<double>(-0.009F);
    constexpr double kCubic2 = static_cast<double>(0.837F);
    constexpr double kCubic1 = static_cast<double>(-22.527F);
    constexpr double kCubic0 = static_cast<double>(196.149F);

    if (temperature <= 20.0) {
        return std::exp(kExpCoefficient * temperature);
    }

    return kCubic3 * temperature * temperature * temperature + kCubic2 * temperature * temperature + kCubic1 * temperature + kCubic0;
}

/** Exact HEMCO 3.12.1 SoilWet scalar response. */
inline double SoilWetnessTerm(double gwet, double arid_fraction, double nonarid_fraction) {
    if (arid_fraction >= nonarid_fraction && arid_fraction > 0.0) {
        return 8.24 * gwet * std::exp(-12.5 * gwet * gwet);
    }
    return 5.5 * gwet * std::exp(-5.55 * gwet * gwet);
}

/** Exact HEMCO 3.12.1 SoilCrf scalar canopy-reduction fraction. */
inline double CanopyReductionFraction(int biome, double leaf_area_index, double canopy_nox, double wind_speed_squared, double solar_zenith_cosine) {
    const std::size_t index = BiomeIndex(biome);
    double ventilation_velocity = solar_zenith_cosine > 0.0 ? 1.0e-2 : 0.2e-2;

    if (leaf_area_index <= 0.0 || canopy_nox <= 0.0) {
        return 0.0;
    }

    constexpr std::size_t kTropicalRainforestIndex = 20;  // HEMCO biome 21.
    ventilation_velocity *= std::sqrt(wind_speed_squared / 9.0 * 7.0 / leaf_area_index) * (kSoilExc[kTropicalRainforestIndex] / kSoilExc[index]);
    return canopy_nox / (canopy_nox + ventilation_velocity);
}

/** HEMCO FertAdd before applying the pinned reference-build scale. */
inline double FertilizerUnscaled(double soil_fertilizer, double deposited_nitrogen) {
    return (soil_fertilizer + deposited_nitrogen) / kSecondsPerYear;
}

/** Evaluate one stateless, one-biome contribution with neutral pulse/land. */
inline StatelessContributionResult EvaluateStatelessContribution(const StatelessContributionInput& input) {
    const std::size_t index = BiomeIndex(input.biome);
    StatelessContributionResult result;

    result.temperature_term = SoilTemperatureTerm(input.use_soil_temperature, input.biome, input.temperature_c, input.gwet);
    result.wetness_term = SoilWetnessTerm(input.gwet, input.arid_fraction, input.nonarid_fraction);
    result.canopy_reduction_fraction =
        CanopyReductionFraction(input.biome, input.leaf_area_index, input.canopy_nox, input.wind_speed_squared, input.solar_zenith_cosine);
    result.fertilizer_unscaled = FertilizerUnscaled(input.soil_fertilizer, input.deposited_nitrogen);
    result.fertilizer_scaled = result.fertilizer_unscaled * kReferenceFertilizerScale;
    result.unit_conversion = UnitConversion(input.molecular_weight_no_g);
    result.cell_flux = (kABiome[index] * result.unit_conversion + result.fertilizer_scaled) * result.temperature_term * result.wetness_term *
                       (1.0 - result.canopy_reduction_fraction);
    return result;
}

}  // namespace cece::hemco_soilnox::v3_12_1

#endif  // CECE_HEMCO_SOILNOX_STATELESS_HPP
