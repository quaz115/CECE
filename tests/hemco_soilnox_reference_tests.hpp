#ifndef CECE_TESTS_HEMCO_SOILNOX_REFERENCE_TESTS_HPP
#define CECE_TESTS_HEMCO_SOILNOX_REFERENCE_TESTS_HPP

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cece/physics/hemco_soilnox_stateless.hpp"

namespace cece::hemco_soilnox::v3_12_1 {
namespace {

struct ReferenceCase {
    std::string name;
    StatelessContributionInput input;
    StatelessContributionResult expected;
};

struct ReferenceData {
    std::unordered_map<std::string, std::string> metadata;
    std::array<double, kBiomeCount> a_biome{};
    std::array<double, kBiomeCount> soil_ta{};
    std::array<double, kBiomeCount> soil_tb{};
    std::array<double, kBiomeCount> soil_exc{};
    std::array<bool, kBiomeCount> saw_biome{};
    std::vector<ReferenceCase> cases;
};

std::unordered_set<std::string> ExpectedCaseNames() {
    std::unordered_set<std::string> names = {
        "air_dry_freeze",      "air_dry_at_zero",  "air_dry_cap_25",
        "air_dry_cap_40",      "air_wet_biome_2",  "air_wet_biome_6",
        "air_wet_biome_15",    "air_wet_biome_22", "soil_negative",
        "soil_zero",           "soil_10",          "soil_20",
        "soil_20_plus",        "soil_30",          "soil_40",
        "soil_45_cap",         "wet_arid_0",       "wet_arid_02",
        "wet_arid_03",         "wet_arid_1",       "wet_nonarid_02",
        "wet_nonarid_03",      "wet_nonarid_1",    "wet_tie_is_arid",
        "wet_zero_is_nonarid", "crf_day_biome_21", "crf_night_biome_21",
        "crf_day_biome_1",     "crf_lai_zero",     "crf_cpynox_zero",
        "fert_zero",           "fert_soil_only",   "fert_dep_only",
        "fert_combined",
    };
    for (std::size_t biome = 1; biome <= kBiomeCount; ++biome) {
        std::ostringstream name;
        name << "onehot_" << std::setfill('0') << std::setw(2) << biome;
        names.emplace(name.str());
    }
    return names;
}

std::string Trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> SplitCsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        fields.push_back(Trim(std::move(field)));
    }
    return fields;
}

double ParseDouble(const std::string& value) {
    std::size_t parsed = 0;
    const double result = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(result)) {
        throw std::runtime_error("Invalid finite floating-point value in HEMCO fixture: " + value);
    }
    return result;
}

int ParseInt(const std::string& value) {
    std::size_t parsed = 0;
    const int result = std::stoi(value, &parsed);
    if (parsed != value.size()) {
        throw std::runtime_error("Invalid integer value in HEMCO fixture: " + value);
    }
    return result;
}

std::filesystem::path ReferenceCsvPath() {
    const auto by_source_file = std::filesystem::path(__FILE__).parent_path() / "data/hemco_bdsnp/hemco_3_12_1_soilnox_reference.csv";
    if (std::filesystem::is_regular_file(by_source_file)) {
        return by_source_file;
    }

    // Defensive fallback for toolchains that emit a relative __FILE__ while
    // CTest launches the binary from inside the build tree.
    for (auto root = std::filesystem::current_path();;) {
        const auto candidate = root / "tests/data/hemco_bdsnp/hemco_3_12_1_soilnox_reference.csv";
        if (std::filesystem::is_regular_file(candidate)) {
            return candidate;
        }
        if (root == root.root_path()) {
            break;
        }
        root = root.parent_path();
    }

    throw std::runtime_error("Unable to locate the pinned HEMCO SoilNOx reference CSV");
}

ReferenceCase ParseCase(const std::vector<std::string>& fields) {
    if (fields.size() != 21) {
        throw std::runtime_error("HEMCO reference case must contain 21 columns");
    }

    ReferenceCase result;
    result.name = fields[0];
    const int use_soil_temperature = ParseInt(fields[1]);
    if (use_soil_temperature != 0 && use_soil_temperature != 1) {
        throw std::runtime_error("HEMCO use_soil_temperature must be zero or one");
    }
    result.input.use_soil_temperature = use_soil_temperature == 1;
    result.input.biome = ParseInt(fields[2]);
    BiomeIndex(result.input.biome);
    result.input.temperature_c = ParseDouble(fields[3]);
    result.input.gwet = ParseDouble(fields[4]);
    result.input.arid_fraction = ParseDouble(fields[5]);
    result.input.nonarid_fraction = ParseDouble(fields[6]);
    result.input.leaf_area_index = ParseDouble(fields[7]);
    result.input.canopy_nox = ParseDouble(fields[8]);
    result.input.wind_speed_squared = ParseDouble(fields[9]);
    result.input.solar_zenith_cosine = ParseDouble(fields[10]);
    result.input.soil_fertilizer = ParseDouble(fields[11]);
    result.input.deposited_nitrogen = ParseDouble(fields[12]);
    result.input.molecular_weight_no_g = ParseDouble(fields[13]);
    result.expected.temperature_term = ParseDouble(fields[14]);
    result.expected.wetness_term = ParseDouble(fields[15]);
    result.expected.canopy_reduction_fraction = ParseDouble(fields[16]);
    result.expected.fertilizer_unscaled = ParseDouble(fields[17]);
    result.expected.fertilizer_scaled = ParseDouble(fields[18]);
    result.expected.unit_conversion = ParseDouble(fields[19]);
    result.expected.cell_flux = ParseDouble(fields[20]);
    return result;
}

ReferenceData LoadReferenceData() {
    std::ifstream input(ReferenceCsvPath());
    if (!input) {
        throw std::runtime_error("Unable to open the pinned HEMCO SoilNOx reference CSV");
    }

    const std::vector<std::string> expected_header = {
        "case",          "use_soil_temperature",
        "biome",         "temp_c",
        "gwet",          "arid",
        "nonarid",       "lai",
        "cpynox",        "windsqr",
        "suncos",        "soilfrt_internal",
        "depn_internal", "mw_no_g",
        "temp_term",     "wet_term",
        "crf",           "fert_unscaled",
        "fert_scaled",   "unitconv",
        "cell_flux",
    };

    ReferenceData result;
    std::unordered_set<std::string> case_names;
    std::string line;
    bool saw_header = false;

    while (std::getline(input, line)) {
        if (line.rfind("# biome_constant,", 0) == 0) {
            const auto fields = SplitCsv(line);
            if (fields.size() != 6) {
                throw std::runtime_error("HEMCO biome constant row must contain six columns");
            }
            const int biome = ParseInt(fields[1]);
            const std::size_t index = BiomeIndex(biome);
            if (result.saw_biome[index]) {
                throw std::runtime_error("Duplicate HEMCO biome constant row");
            }
            result.saw_biome[index] = true;
            result.a_biome[index] = ParseDouble(fields[2]);
            result.soil_ta[index] = ParseDouble(fields[3]);
            result.soil_tb[index] = ParseDouble(fields[4]);
            result.soil_exc[index] = ParseDouble(fields[5]);
            continue;
        }

        if (line.rfind("# ", 0) == 0) {
            const auto separator = line.find('=');
            if (separator != std::string::npos) {
                const std::string key = Trim(line.substr(2, separator - 2));
                const std::string value = Trim(line.substr(separator + 1));
                if (key.empty() || value.empty() || !result.metadata.emplace(key, value).second) {
                    throw std::runtime_error("Invalid or duplicate HEMCO reference metadata: " + key);
                }
            }
            continue;
        }

        if (Trim(line).empty()) {
            continue;
        }

        const auto fields = SplitCsv(line);
        if (!saw_header) {
            if (fields != expected_header) {
                throw std::runtime_error("Unexpected HEMCO reference CSV header");
            }
            saw_header = true;
            continue;
        }

        ReferenceCase test_case = ParseCase(fields);
        if (!case_names.emplace(test_case.name).second) {
            throw std::runtime_error("Duplicate HEMCO reference case: " + test_case.name);
        }
        result.cases.push_back(std::move(test_case));
    }

    if (!saw_header) {
        throw std::runtime_error("HEMCO reference CSV header is missing");
    }
    return result;
}

const ReferenceCase& FindCase(const ReferenceData& reference, const std::string& name) {
    const auto found =
        std::find_if(reference.cases.begin(), reference.cases.end(), [&name](const ReferenceCase& test_case) { return test_case.name == name; });
    if (found == reference.cases.end()) {
        throw std::runtime_error("Missing HEMCO reference case: " + name);
    }
    return *found;
}

void ExpectHemcoNear(double actual, double expected) {
    ASSERT_TRUE(std::isfinite(actual));
    ASSERT_TRUE(std::isfinite(expected));
    if (expected == 0.0) {
        EXPECT_DOUBLE_EQ(actual, 0.0);
        return;
    }
    constexpr double kRelativeTolerance = 2.0e-12;
    constexpr double kScaleFloor = 1.0e-30;
    EXPECT_LE(std::abs(actual - expected), kRelativeTolerance * std::max(std::abs(expected), kScaleFloor));
}

TEST(HemcoSoilNox3121ReferenceTest, FixtureContractAndBiomeConstantsArePinned) {
    const ReferenceData reference = LoadReferenceData();

    ASSERT_EQ(reference.metadata.at("reference_repo"), "https://github.com/geoschem/HEMCO");
    ASSERT_EQ(reference.metadata.at("reference_tag"), "3.12.1");
    ASSERT_EQ(reference.metadata.at("reference_commit"), "07da3c29fd85abc3824cb6288578b0b68c2395a3");
    ASSERT_EQ(reference.metadata.at("reference_file"), "src/Extensions/hcox_soilnox_mod.F90");
    ASSERT_EQ(reference.metadata.at("reference_blob"), "a8712e3e89fc0034ccd9761b476ce07608ac4506");
    ASSERT_EQ(reference.metadata.at("no_spec_blob"), "afbaf799e66773ef95f682ad984cde86954bc922");
    ASSERT_EQ(reference.metadata.at("hemco_hp_bits"), "64");
    ASSERT_EQ(reference.metadata.at("neutral_state"), "pulse:1,land_fraction:1");
    ASSERT_TRUE(reference.metadata.contains("compiler_version"));
    ASSERT_TRUE(reference.metadata.contains("compiler_options"));
    ASSERT_FALSE(reference.metadata.at("compiler_version").empty());
    ASSERT_FALSE(reference.metadata.at("compiler_options").empty());
    ASSERT_EQ(reference.cases.size(), 58U);

    std::unordered_set<std::string> actual_case_names;
    for (const auto& test_case : reference.cases) {
        actual_case_names.emplace(test_case.name);
    }
    EXPECT_EQ(actual_case_names, ExpectedCaseNames());

    for (std::size_t biome = 1; biome <= kBiomeCount; ++biome) {
        std::ostringstream name;
        name << "onehot_" << std::setfill('0') << std::setw(2) << biome;
        EXPECT_EQ(FindCase(reference, name.str()).input.biome, static_cast<int>(biome));
    }

    for (std::size_t index = 0; index < kBiomeCount; ++index) {
        SCOPED_TRACE("biome=" + std::to_string(index + 1));
        ASSERT_TRUE(reference.saw_biome[index]);
        EXPECT_DOUBLE_EQ(reference.a_biome[index], kABiome[index]);
        EXPECT_DOUBLE_EQ(reference.soil_ta[index], kSoilTa[index]);
        EXPECT_DOUBLE_EQ(reference.soil_tb[index], kSoilTb[index]);
        EXPECT_DOUBLE_EQ(reference.soil_exc[index], kSoilExc[index]);
    }
}

TEST(HemcoSoilNox3121ReferenceTest, StatelessFunctionsMatchAllGeneratedVectors) {
    const ReferenceData reference = LoadReferenceData();

    for (const auto& test_case : reference.cases) {
        SCOPED_TRACE(test_case.name);
        const auto actual = EvaluateStatelessContribution(test_case.input);

        ExpectHemcoNear(actual.temperature_term, test_case.expected.temperature_term);
        ExpectHemcoNear(actual.wetness_term, test_case.expected.wetness_term);
        ExpectHemcoNear(actual.canopy_reduction_fraction, test_case.expected.canopy_reduction_fraction);
        ExpectHemcoNear(actual.fertilizer_unscaled, test_case.expected.fertilizer_unscaled);
        ExpectHemcoNear(actual.fertilizer_scaled, test_case.expected.fertilizer_scaled);
        ExpectHemcoNear(actual.unit_conversion, test_case.expected.unit_conversion);
        ExpectHemcoNear(actual.cell_flux, test_case.expected.cell_flux);
    }
}

TEST(HemcoSoilNox3121ReferenceTest, BoundaryAndBranchSemanticsAreExplicit) {
    const ReferenceData reference = LoadReferenceData();

    EXPECT_DOUBLE_EQ(FindCase(reference, "air_dry_freeze").expected.temperature_term, 0.0);
    EXPECT_DOUBLE_EQ(FindCase(reference, "soil_negative").expected.temperature_term, 0.0);
    EXPECT_DOUBLE_EQ(FindCase(reference, "soil_zero").expected.temperature_term, 0.0);
    EXPECT_DOUBLE_EQ(FindCase(reference, "crf_lai_zero").expected.canopy_reduction_fraction, 0.0);
    EXPECT_DOUBLE_EQ(FindCase(reference, "crf_cpynox_zero").expected.canopy_reduction_fraction, 0.0);

    ExpectHemcoNear(FindCase(reference, "air_dry_cap_25").expected.temperature_term, FindCase(reference, "air_dry_cap_40").expected.temperature_term);
    ExpectHemcoNear(FindCase(reference, "soil_40").expected.temperature_term, FindCase(reference, "soil_45_cap").expected.temperature_term);
    EXPECT_NE(FindCase(reference, "soil_20").expected.temperature_term, FindCase(reference, "soil_20_plus").expected.temperature_term);
    ExpectHemcoNear(FindCase(reference, "wet_tie_is_arid").expected.wetness_term, FindCase(reference, "wet_arid_02").expected.wetness_term);
    ExpectHemcoNear(FindCase(reference, "wet_zero_is_nonarid").expected.wetness_term, FindCase(reference, "wet_nonarid_03").expected.wetness_term);
}

TEST(HemcoSoilNox3121ReferenceTest, RejectsOutOfRangeBiomeIndices) {
    EXPECT_THROW(BiomeIndex(0), std::out_of_range);
    EXPECT_THROW(BiomeIndex(25), std::out_of_range);
}

}  // namespace
}  // namespace cece::hemco_soilnox::v3_12_1

#endif  // CECE_TESTS_HEMCO_SOILNOX_REFERENCE_TESTS_HPP
