/**
 * @file test_bdsnp.cpp
 * @brief Generic registration and YL95 compatibility checks for BdsnpScheme.
 */

#include <gtest/gtest.h>
#include <rapidcheck.h>
#include <rapidcheck/gtest.h>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "cece/cece_physics_factory.hpp"
#include "cece/cece_state.hpp"
#include "cece/physics/cece_bdsnp.hpp"
#include "cece/physics/cece_soil_nox.hpp"

#ifdef CECE_HAS_FORTRAN
#include "cece/physics/cece_bdsnp_fortran.hpp"
#endif

namespace cece {
namespace {

DualView3D MakeField(const std::string& name, double value) {
    DualView3D field(name, 1, 1, 1);
    Kokkos::deep_copy(field.view_host(), value);
    field.modify_host();
    field.sync_device();
    return field;
}

double RunBdsnpYl95(double temperature_k, double soil_moisture) {
    CeceImportState import_state;
    CeceExportState export_state;
    import_state.fields["soil_temperature"] = MakeField("soil_temperature", temperature_k);
    import_state.fields["soil_moisture"] = MakeField("soil_moisture", soil_moisture);
    export_state.fields["soil_nox_emissions"] = MakeField("soil_nox_emissions", -1.0);

    YAML::Node config;
    config["soil_no_method"] = "yl95";
    BdsnpScheme scheme;
    scheme.Initialize(config, nullptr);
    scheme.Run(import_state, export_state);

    auto& output = export_state.fields.at("soil_nox_emissions");
    output.sync_host();
    return output.view_host()(0, 0, 0);
}

double RunSoilNoxYl95(double temperature_k, double soil_moisture) {
    CeceImportState import_state;
    CeceExportState export_state;
    import_state.fields["temperature"] = MakeField("temperature", temperature_k);
    import_state.fields["soil_moisture"] = MakeField("soil_moisture", soil_moisture);
    // SoilNoxScheme accumulates into its export; zero isolates its one-run YL95 flux.
    export_state.fields["soil_nox_emissions"] = MakeField("soil_nox_emissions", 0.0);

    SoilNoxScheme scheme;
    scheme.Initialize(YAML::Node{}, nullptr);
    scheme.Run(import_state, export_state);

    auto& output = export_state.fields.at("soil_nox_emissions");
    output.sync_host();
    return output.view_host()(0, 0, 0);
}

TEST(BdsnpSchemeTest, FactoryCreatesBdsnpScheme) {
    PhysicsSchemeConfig config;
    config.name = "bdsnp";
    auto scheme = PhysicsFactory::CreateScheme(config);
    EXPECT_NE(scheme, nullptr);
}

TEST(BdsnpSchemeTest, SupportedMethodsInitializeAndUnknownMethodsFail) {
    BdsnpScheme canonical;
    YAML::Node bdsnp;
    bdsnp["soil_no_method"] = "bdsnp";
    EXPECT_NO_THROW(canonical.Initialize(bdsnp, nullptr));

    BdsnpScheme fallback;
    YAML::Node yl95;
    yl95["soil_no_method"] = "yl95";
    EXPECT_NO_THROW(fallback.Initialize(yl95, nullptr));

    BdsnpScheme removed_selector;
    YAML::Node old_name;
    old_name["soil_no_method"] = "hemco_3_12_1";
    EXPECT_THROW(removed_selector.Initialize(old_name, nullptr), std::invalid_argument);

    BdsnpScheme unknown;
    YAML::Node typo;
    typo["soil_no_method"] = "not-a-method";
    EXPECT_THROW(unknown.Initialize(typo, nullptr), std::invalid_argument);

    BdsnpScheme obsolete;
    YAML::Node removed_option;
    removed_option["fert_emission_factor"] = 1.0;
    EXPECT_THROW(obsolete.Initialize(removed_option, nullptr), std::invalid_argument);
}

TEST(BdsnpSchemeTest, DiagnosticFieldsRegisterWhenEnabled) {
    YAML::Node config;
    config["soil_no_method"] = "bdsnp";
    config["diagnostics"].push_back("soil_no_emission_rate");
    config["nx"] = 1;
    config["ny"] = 1;
    config["nz"] = 1;

    CeceDiagnosticManager diagnostics;
    BdsnpScheme scheme;
    EXPECT_NO_THROW(scheme.Initialize(config, &diagnostics));
}

RC_GTEST_PROP(BdsnpProperty, Yl95FreezingProducesExactlyZero, ()) {
    const double temperature_k = 200.0 + (*rc::gen::inRange(0, 7315)) / 100.0;
    RC_PRE(temperature_k < 273.15);
    const double soil_moisture = (*rc::gen::inRange(0, 10001)) / 10000.0;
    RC_ASSERT(RunBdsnpYl95(temperature_k, soil_moisture) == 0.0);
}

RC_GTEST_PROP(BdsnpProperty, Yl95MatchesExistingSoilNoxScheme, ()) {
    const double temperature_k = 274.0 + (*rc::gen::inRange(0, 5601)) / 100.0;
    const double soil_moisture = 0.01 + (*rc::gen::inRange(0, 9900)) / 10000.0;
    const double candidate = RunBdsnpYl95(temperature_k, soil_moisture);
    const double reference = RunSoilNoxYl95(temperature_k, soil_moisture);
    const double tolerance = std::max(std::abs(reference) * 1.0e-6, 1.0e-15);
    RC_ASSERT(std::abs(candidate - reference) <= tolerance);
}

#ifdef CECE_HAS_FORTRAN
RC_GTEST_PROP(BdsnpProperty, Yl95MatchesLegacyFortranInterface, ()) {
    const double temperature_k = 274.0 + (*rc::gen::inRange(0, 5601)) / 100.0;
    const double soil_moisture = 0.01 + (*rc::gen::inRange(0, 9900)) / 10000.0;

    CeceImportState import_state;
    CeceExportState export_state;
    import_state.fields["soil_temperature"] = MakeField("soil_temperature", temperature_k);
    import_state.fields["soil_moisture"] = MakeField("soil_moisture", soil_moisture);
    export_state.fields["soil_nox_emissions"] = MakeField("soil_nox_emissions", -1.0);

    YAML::Node config;
    config["soil_no_method"] = "yl95";
    BdsnpFortranScheme scheme;
    scheme.Initialize(config, nullptr);
    scheme.Run(import_state, export_state);

    auto& output = export_state.fields.at("soil_nox_emissions");
    output.sync_host();
    const double reference = output.view_host()(0, 0, 0);
    const double candidate = RunBdsnpYl95(temperature_k, soil_moisture);
    const double tolerance = std::max(std::abs(reference) * 1.0e-6, 1.0e-15);
    RC_ASSERT(std::abs(candidate - reference) <= tolerance);
}
#endif

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
