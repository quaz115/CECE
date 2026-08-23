#ifndef CECE_BDSNP_HPP
#define CECE_BDSNP_HPP

#include "cece/physics_scheme.hpp"

namespace cece {

/**
 * @class BdsnpScheme
 * @brief Standalone soil NO physics module implementing the Berkeley-Dalhousie
 * Soil NO Parameterization (BDSNP) with YL95 fallback.
 *
 * Provides a comprehensive soil-NO model alongside the legacy SoilNoxScheme
 * ("soil_nox"). Supports two algorithms selectable via the `soil_no_method`
 * YAML configuration key:
 *   - "bdsnp" (default): validated effective-input BDSNP calculation with
 *     24-biome weighting, temperature and moisture responses, canopy reduction,
 *     pulse scaling, fertilizer, and deposited nitrogen
 *   - "yl95": Yienger & Levy (1995) soil temperature response, soil moisture
 *     pulse, canopy reduction factor
 *
 * BDSNP pulse, canopy, fertilizer, and deposited-nitrogen terms are supplied as
 * effective input fields. Persistent pulse and deposited-nitrogen reservoir
 * evolution are outside this stateless scheme.
 *
 * Writes computed soil NO emissions to the export state field
 * "soil_nox_emissions" for consumption by MEGAN3 or other schemes.
 */
class BdsnpScheme : public BasePhysicsScheme {
   public:
    BdsnpScheme() = default;
    ~BdsnpScheme() override = default;

    void Initialize(const YAML::Node& config, CeceDiagnosticManager* diag_manager) override;
    void Run(CeceImportState& import_state, CeceExportState& export_state) override;

   private:
    std::string soil_no_method_ = "bdsnp";  // "bdsnp" or "yl95"
    bool use_soil_temperature_ = false;

    // YL95 parameters (reused from existing SoilNoxScheme)
    double a_biome_wet_ = 0.5;
    double tc_max_ = 30.0;
    double exp_coeff_ = 0.103;
    double wet_c1_ = 5.5;
    double wet_c2_ = -5.55;
};

}  // namespace cece

#endif  // CECE_BDSNP_HPP
