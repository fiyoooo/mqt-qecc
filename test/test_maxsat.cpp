//
// Created by Fiona Froehler on 25.07.24.
//
// to keep 0/1 in boolean areas without clang-tidy warnings:
// NOLINTBEGIN(readability-implicit-bool-conversion,modernize-use-bool-literals)

#include "LightsOut.hpp"
#include "MaxSATDecoder.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>

class MaxSATTest : public testing::TestWithParam<gf2Vec> {};
class MaxSATLightsOutComparison : public testing::TestWithParam<std::tuple<gf2Vec, std::unordered_map<int, std::vector<int>>>> {};
class UniquelyCorrectableErrTestMaxSAT : public MaxSATTest {};
class InCorrectableErrTestMaxSAT : public MaxSATTest {};
class UpToStabCorrectableErrTestMaxSAT : public MaxSATTest {};

// not all have to be the same, as different corrections can match the same syndrome
INSTANTIATE_TEST_SUITE_P(TriTri4Bit, MaxSATLightsOutComparison,
                         ::testing::Combine(
                                 ::testing::Values(
                                         gf2Vec{false, false, false, false},
                                         gf2Vec{true, false, false, false},
                                         gf2Vec{false, true, false, false},
                                         gf2Vec{false, false, true, false},
                                         gf2Vec{false, false, false, true}),
                                 ::testing::Values(
                                         std::unordered_map<int, std::vector<int>>{
                                                 {0, {0, 1, 2}},
                                                 {1, {1, 2, 4}},
                                                 {2, {1, 3, 4}},
                                                 {3, {2, 4, 5}}})));

// not all have to be the same, as different corrections can match the same syndrome
INSTANTIATE_TEST_SUITE_P(TriHex7Bit, MaxSATLightsOutComparison,
                         ::testing::Combine(
                                 ::testing::Values(
                                         gf2Vec{0, 0, 0, 0, 0, 0, 0},
                                         gf2Vec{1, 0, 0, 0, 0, 0, 0},
                                         gf2Vec{0, 1, 0, 0, 0, 0, 0},
                                         gf2Vec{0, 0, 1, 0, 0, 0, 0},
                                         gf2Vec{0, 0, 0, 1, 0, 0, 0},
                                         gf2Vec{0, 0, 0, 0, 1, 0, 0}),
                                 ::testing::Values(
                                         std::unordered_map<int, std::vector<int>>{
                                                 {0, {0, 1, 2, 3}},
                                                 {1, {1, 3, 5, 6}},
                                                 {2, {2, 3, 4, 5, 7, 8}},
                                                 {3, {4, 7, 10, 11}},
                                                 {4, {5, 6, 8, 9, 12, 13}},
                                                 {5, {7, 8, 11, 12}},
                                                 {6, {9, 13}}})));

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrsSteane, UniquelyCorrectableErrTestMaxSAT,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, InCorrectableErrTestMaxSAT,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 0, 1, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, UpToStabCorrectableErrTestMaxSAT,
                         testing::Values(
                                 std::vector<bool>{0, 0, 0, 0, 0, 0, 1},
                                 std::vector<bool>{1, 1, 0, 0, 0, 0, 0},
                                 std::vector<bool>{0, 0, 0, 0, 1, 1, 0},
                                 std::vector<bool>{1, 0, 0, 0, 0, 0, 1}));

/**
 * Tests for comparing MaxSAT without exact translation LightsOut. Not all have to pass.
 */
TEST_P(MaxSATLightsOutComparison, ComparingWithLightsOut) {
    auto                                            params           = GetParam();
    gf2Vec const                                    syndrome         = std::get<0>(params);
    std::unordered_map<int, std::vector<int>> const lightsToSwitches = std::get<1>(params);

    // transform lightsToSwitches to switchesToLights
    std::unordered_map<int, std::vector<int>> switchesToLights;
    for (const auto& lightSwitchPair : lightsToSwitches) {
        int const light = lightSwitchPair.first;
        for (int const switchIndex : lightSwitchPair.second) {
            switchesToLights[switchIndex].push_back(light);
        }
    }

    // initialize LightsOutDecoder solver
    LightsOut lightsOutSolver(lightsToSwitches, switchesToLights);
    lightsOutSolver.preconstructZ3Instance();
    auto [switchesLO, constrTimeLO, solveTimeLO] = lightsOutSolver.solve(syndrome);

    // initialize hZ
    size_t const k = lightsToSwitches.size(); // number of lights
    size_t const n = switchesToLights.size(); // number of switches
    gf2Mat       hz(k, gf2Vec(n, false));
    for (const auto& entry : lightsToSwitches) {
        auto const light = static_cast<const size_t>(entry.first);
        for (int const switchIndex : entry.second) {
            hz[light][static_cast<const size_t>(switchIndex)] = true;
        }
    }

    // initialize MaxSATDecoder solver
    Code          code(hz);
    MaxSATDecoder maxSATSolver(code);
    maxSATSolver.preconstructZ3Instance();
    maxSATSolver.decode(syndrome);
    std::vector<bool>      boolVector = maxSATSolver.result.estimBoolVector;
    std::vector<int> const switchesMs(boolVector.begin(), boolVector.end());

    // compare results
    std::cout << "Syndrome: ";
    for (auto sw : syndrome) {
        std::cout << sw << " ";
    }
    std::cout << "\nLightsOut switches: ";
    for (auto& sw : switchesLO) {
        std::cout << sw << " ";
    }
    std::cout << "\nMaxSATDecoder switches: ";
    for (auto sw : switchesMs) {
        std::cout << sw << " ";
    }
    std::cout << "\n\n";

    EXPECT_EQ(switchesLO, switchesMs);
}

/**
 * Tests for unambiguous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTestMaxSAT, SteaneCodeDecodingTestEstim) {
    auto          code = SteaneXCode();
    MaxSATDecoder decoder(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    const std::vector<bool> err = GetParam();

    auto syndr = code.getXSyndrome(err);
    std::cout << "syndrome: ";
    Utils::printGF2vector(syndr);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << std::endl;
    const gf2Vec sol = GetParam();

    std::cout << "Estim: " << std::endl;
    Utils::printGF2vector(estim);
    std::cout << "EstimIdx: " << std::endl;
    Utils::printGF2vector(estim2);
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_TRUE(sol == estim);
    EXPECT_TRUE(sol == estim2);
}

/**
 * Tests for ambiguous errors that cannot be corrected
 */
TEST_P(InCorrectableErrTestMaxSAT, SteaneCodeDecodingTestEstim) {
    auto          code = SteaneXCode();
    MaxSATDecoder decoder(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();

    auto syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }

    const gf2Vec sol = GetParam();

    std::cout << "Estim: " << std::endl;
    Utils::printGF2vector(estim);
    std::cout << "EstimIdx: " << std::endl;
    Utils::printGF2vector(estim2);
    std::cout << "Sol: " << std::endl;
    Utils::printGF2vector(sol);
    EXPECT_FALSE(sol == estim);
    EXPECT_FALSE(sol == estim2);
}

/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_P(UpToStabCorrectableErrTestMaxSAT, SteaneCodeDecodingTest) {
    auto          code = SteaneCode();
    MaxSATDecoder decoder(code);
    std::cout << "code: " << std::endl
              << code << std::endl;
    std::vector<bool> err = GetParam();
    std::cout << "err :" << std::endl;
    Utils::printGF2vector(err);
    auto syndr = code.getXSyndrome(err);
    decoder.decode(syndr);
    const auto& decodingResult = decoder.result;
    const auto& estim          = decodingResult.estimBoolVector;
    const auto& estimIdx       = decodingResult.estimNodeIdxVector;
    gf2Vec      estim2(err.size());
    std::cout << "estiIdxs: ";
    for (auto idx : estimIdx) {
        estim2.at(idx) = true;
        std::cout << idx;
    }
    std::cout << std::endl;
    std::vector<bool> residualErr(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
    }
    std::vector<bool> residualErr2(err.size());
    for (size_t i = 0; i < err.size(); i++) {
        residualErr2.at(i) = (err[i] != estim2[i]);
    }

    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethX()->pcm, residualErr));
    EXPECT_TRUE(Utils::isVectorInRowspace(*code.gethX()->pcm, residualErr2));
}

// NOLINTEND(readability-implicit-bool-conversion,modernize-use-bool-literals)
