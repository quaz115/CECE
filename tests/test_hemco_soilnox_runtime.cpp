/**
 * @file test_hemco_soilnox_runtime.cpp
 * @brief Production-kernel checks against the pinned HEMCO 3.12.1 oracle.
 */

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <string>

#include "cece/cece_state.hpp"
#include "cece/physics/cece_bdsnp.hpp"
#include "cece/physics/hemco_soilnox_stateless.hpp"

namespace cece {
namespace {

using hemco_soilnox::v3_12_1::EvaluateStatelessContribution;
using hemco_soilnox::v3_12_1::StatelessContributionInput;

void ExpectHemcoNear(double actual, double expected) {
    ASSERT_TRUE(std::isfinite(actual));
    ASSERT_TRUE(std::isfinite(expected));
    if (expected == 0.0) {
        EXPECT_DOUBLE_EQ(actual, 0.0);
        return;
    }
    EXPECT_LE(std::abs(actual - expected), 2.0e-12 * std::max(std::abs(expected), 1.0e-30));
}

class HemcoSoilNoxRuntimeTest : public ::testing::Test {
   protected:
    static constexpr int kNx = 2;
    static constexpr int kNy = 2;

    CeceImportState import_state;
    CeceExportState export_state;

    DualView3D MakeField(const std::string& name, int levels, double value) {
        DualView3D field(name, kNx, kNy, levels);
        Kokkos::deep_copy(field.view_host(), value);
        field.modify<Kokkos::HostSpace>();
        field.sync<Kokkos::DefaultExecutionSpace>();
        return field;
    }

    void AddImport(const std::string& name, int levels, double value) {
        import_state.fields[name] = MakeField(name, levels, value);
    }

    void SetImportLayer(const std::string& name, int layer, double value) {
        auto& field = import_state.fields.at(name);
        field.sync<Kokkos::HostSpace>();
        auto host = field.view_host();
        for (int i = 0; i < kNx; ++i) {
            for (int j = 0; j < kNy; ++j) {
                host(i, j, layer) = value;
            }
        }
        field.modify<Kokkos::HostSpace>();
        field.sync<Kokkos::DefaultExecutionSpace>();
    }

    void SetScalarImport(const std::string& name, double value) {
        SetImportLayer(name, 0, value);
    }

    void AddOutputs() {
        export_state.fields["soil_nox_emissions"] = MakeField("soil_nox_emissions", 1, -1.0);
        export_state.fields["soil_nox_fertilizer_emissions"] = MakeField("soil_nox_fertilizer_emissions", 1, -1.0);
    }

    void AddNeutralInputs() {
        AddImport("surface_temperature", 1, 293.15);
        AddImport("soil_moisture", 1, 0.2);
        AddImport("soilnox_land_fractions", 24, 0.0);
        AddImport("soilnox_arid_fraction", 1, 0.0);
        AddImport("soilnox_nonarid_fraction", 1, 1.0);
        AddImport("leaf_area_index", 1, 0.0);
        AddImport("soilnox_canopy_nox", 24, 0.0);
        AddImport("wind_speed_squared", 1, 9.0);
        AddImport("solar_zenith_cosine", 1, 1.0);
        AddImport("soil_fertilizer", 1, 0.0);
        AddImport("deposited_nitrogen", 1, 0.0);
        AddImport("soilnox_pulse_factor", 1, 1.0);
        AddOutputs();
    }

    static YAML::Node HemcoConfig() {
        YAML::Node config;
        config["soil_no_method"] = "hemco_3_12_1";
        config["use_soil_temperature"] = false;
        return config;
    }

    double OutputAt(const std::string& name) {
        auto& field = export_state.fields.at(name);
        field.sync<Kokkos::HostSpace>();
        return field.view_host()(0, 0, 0);
    }
};

TEST_F(HemcoSoilNoxRuntimeTest, OneHotBiomeMatchesPinnedScalarOracle) {
    AddNeutralInputs();
    SetImportLayer("soilnox_land_fractions", 5, 1.0);  // HEMCO biome 6.

    BdsnpScheme scheme;
    scheme.Initialize(HemcoConfig(), nullptr);
    scheme.Run(import_state, export_state);

    StatelessContributionInput input;
    input.biome = 6;
    input.temperature_c = 20.0;
    input.gwet = 0.2;
    input.nonarid_fraction = 1.0;
    const auto expected = EvaluateStatelessContribution(input);

    ExpectHemcoNear(OutputAt("soil_nox_emissions"), expected.cell_flux);
    EXPECT_DOUBLE_EQ(OutputAt("soil_nox_fertilizer_emissions"), 0.0);
}

TEST_F(HemcoSoilNoxRuntimeTest, MixedBiomesAndAddedNitrogenAreWeightedExactly) {
    AddNeutralInputs();
    SetImportLayer("soilnox_land_fractions", 5, 0.25);  // HEMCO biome 6.
    SetImportLayer("soilnox_land_fractions", 9, 0.75);  // HEMCO biome 10.
    SetScalarImport("soil_fertilizer", 1.0e-4);
    SetScalarImport("deposited_nitrogen", 2.0e-4);

    BdsnpScheme scheme;
    scheme.Initialize(HemcoConfig(), nullptr);
    scheme.Run(import_state, export_state);

    StatelessContributionInput biome6;
    biome6.biome = 6;
    biome6.temperature_c = 20.0;
    biome6.gwet = 0.2;
    biome6.nonarid_fraction = 1.0;
    biome6.soil_fertilizer = 1.0e-4;
    biome6.deposited_nitrogen = 2.0e-4;
    auto biome10 = biome6;
    biome10.biome = 10;

    const auto result6 = EvaluateStatelessContribution(biome6);
    const auto result10 = EvaluateStatelessContribution(biome10);
    const double expected_total = 0.25 * result6.cell_flux + 0.75 * result10.cell_flux;
    const double expected_added_n =
        0.25 * result6.fertilizer_scaled * result6.temperature_term * result6.wetness_term +
        0.75 * result10.fertilizer_scaled * result10.temperature_term * result10.wetness_term;

    ExpectHemcoNear(OutputAt("soil_nox_emissions"), expected_total);
    ExpectHemcoNear(OutputAt("soil_nox_fertilizer_emissions"), expected_added_n);
}

TEST_F(HemcoSoilNoxRuntimeTest, ExactNoSoilClassClearsBothOutputs) {
    AddNeutralInputs();
    SetImportLayer("soilnox_land_fractions", 0, 1.0);
    SetScalarImport("soil_fertilizer", 1.0);
    SetScalarImport("deposited_nitrogen", 1.0);

    BdsnpScheme scheme;
    scheme.Initialize(HemcoConfig(), nullptr);
    scheme.Run(import_state, export_state);

    EXPECT_DOUBLE_EQ(OutputAt("soil_nox_emissions"), 0.0);
    EXPECT_DOUBLE_EQ(OutputAt("soil_nox_fertilizer_emissions"), 0.0);
}

TEST_F(HemcoSoilNoxRuntimeTest, MissingOrMalformedHemcoContractFailsLoudly) {
    AddOutputs();
    BdsnpScheme missing_scheme;
    missing_scheme.Initialize(HemcoConfig(), nullptr);
    EXPECT_THROW(missing_scheme.Run(import_state, export_state), std::runtime_error);

    missing_scheme.ClearPhysicsCache();
    AddNeutralInputs();
    import_state.fields["soilnox_land_fractions"] = MakeField("bad_land_fractions", 1, 1.0);
    EXPECT_THROW(missing_scheme.Run(import_state, export_state), std::runtime_error);
}

}  // namespace
}  // namespace cece

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (!Kokkos::is_initialized()) {
        Kokkos::initialize(argc, argv);
    }
    const int result = RUN_ALL_TESTS();
    if (Kokkos::is_initialized()) {
        Kokkos::finalize();
    }
    return result;
}
