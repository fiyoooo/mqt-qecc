//
// Created by Fiona Froehler on 29.08.24.
//

#include "Codes.hpp"
#include "MaxSATDecoder.hpp"
#include "UFDecoder.hpp"
#include "Utils.hpp"
#include "union_find.hpp"

#include <chrono>
#include <gtest/gtest.h>
#include <tuple>

class DecoderComparison : public testing::TestWithParam<gf2Vec> {};

INSTANTIATE_TEST_SUITE_P(SingleBitErrs, DecoderComparison,
                         ::testing::Values(
                                 gf2Vec{0, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{1, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 1, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 0, 1, 0, 0, 0, 0},
                                 gf2Vec{0, 0, 0, 1, 0, 0, 0},
                                 gf2Vec{0, 0, 0, 0, 1, 0, 0},
                                 gf2Vec{0, 0, 0, 0, 0, 1, 0},
                                 gf2Vec{0, 0, 0, 0, 0, 0, 1}));

gf2Vec runUFDecoder(Code& code, const gf2Vec& syndrome) {
    UFDecoder ufDecoder(code);
    ufDecoder.decode(syndrome);
    return ufDecoder.result.estimBoolVector;
}

gf2Vec runMaxSatDecoder(Code& code, const gf2Vec& syndrome) {
    MaxSATDecoder maxSATSolver(code);
    maxSATSolver.preconstructZ3Instance();
    maxSATSolver.decode(syndrome);
    return maxSATSolver.result.estimBoolVector;
}

std::vector<uint8_t> runUfdPeelDecoder(ldpc::bp::BpSparse& pcm, const std::vector<uint8_t>& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.peel_decode(syndrome);
}

std::vector<uint8_t> runUfdMatrixDecoder(ldpc::bp::BpSparse& pcm, const std::vector<uint8_t>& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.matrix_decode(syndrome);
}

std::vector<uint8_t> runUfdMaxSatDecoder(ldpc::bp::BpSparse& pcm, const std::vector<uint8_t>& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.maxsat_decode(syndrome);
}

template <typename Func, typename ReturnType = std::invoke_result_t<Func>>
std::tuple<double, ReturnType> measureRuntime(Func func) {
    auto                                start    = std::chrono::high_resolution_clock::now();
    ReturnType                          estim    = func();
    auto                                end      = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const duration = end - start;
    return std::make_tuple(duration.count(), estim);
}

void testAllDecoders(Code& code, const gf2Vec& err, const gf2Vec& syndr, bool peel) {
    std::cout << "Sol: ";
    Utils::printGF2vector(err);

    auto [runtimeUFD, estimUFD] = measureRuntime([&]() { return runUFDecoder(code, syndr); });
    std::cout << "\nEstimUFD: ";
    Utils::printGF2vector(estimUFD);

    auto [runtimeMS, estimMS] = measureRuntime([&]() { return runMaxSatDecoder(code, syndr); });
    std::cout << "\nEstimMS: ";
    Utils::printGF2vector(estimMS);

    // transform pcm and syndrome to right format
    ldpc::bp::BpSparse         pcm      = Utils::toBpSparse(code);
    std::vector<uint8_t> const syndrInt = Utils::toUint8(syndr);

    gf2Vec estimUfdPeel;
    double runtimeUfdPeel = 0;
    if (peel) {
        auto [runtimeUfdPeel, estimUfdPeelInt] = measureRuntime([&]() { return runUfdPeelDecoder(pcm, syndrInt); });
        estimUfdPeel                           = Utils::toGf2Vec(estimUfdPeelInt);
        std::cout << "\nEstimUfdPeel: ";
        Utils::printGF2vector(estimUfdPeel);
    }

    auto [runtimeUfdMatrix, estimUfdMatrixInt] = measureRuntime([&]() { return runUfdMatrixDecoder(pcm, syndrInt); });
    gf2Vec estimUfdMatrix                      = Utils::toGf2Vec(estimUfdMatrixInt);
    std::cout << "\nEstimUfdMatrix: ";
    Utils::printGF2vector(estimUfdMatrix);

    auto [runtimeUfdMaxSat, estimUfdMaxSatInt] = measureRuntime([&]() { return runUfdMaxSatDecoder(pcm, syndrInt); });
    gf2Vec estimUfdMaxSat                      = Utils::toGf2Vec(estimUfdMaxSatInt);
    std::cout << "\nEstimUfdMaxSat: ";
    Utils::printGF2vector(estimUfdMaxSat);
    std::cout << "\n\n";

    // compare results with actual error
    EXPECT_EQ(estimUFD, err);
    EXPECT_EQ(estimMS, err);
    if (peel) {
        EXPECT_EQ(estimUfdPeel, err);
    }
    EXPECT_EQ(estimUfdMatrix, err);
    EXPECT_EQ(estimUfdMaxSat, err);

    std::cout << "Runtime (UFDecoder): " << runtimeUFD << " seconds\n";
    std::cout << "Runtime (MaxSatDecoder): " << runtimeMS << " seconds\n";
    if (peel) {
        std::cout << "Runtime (ldpc::uf::UfDecoder peel_decode): " << runtimeUfdPeel << " seconds\n";
    }
    std::cout << "Runtime (ldpc::uf::UfDecoder matrix_decode): " << runtimeUfdMatrix << " seconds\n";
    std::cout << "Runtime (ldpc::uf::UfDecoder maxsat_decode): " << runtimeUfdMaxSat << " seconds\n\n";
}

TEST_P(DecoderComparison, SteaneXCodeTest) {
    auto code = SteaneXCode();
    std::cout << "Code: " << '\n'
              << code << '\n';

    gf2Vec const err   = GetParam();
    gf2Vec const syndr = code.getXSyndrome(err);
    std::cout << "Syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << '\n';

    testAllDecoders(code, err, syndr, false);
}

TEST_P(DecoderComparison, ToricCode8Test) {
    auto code = ToricCode8();
    std::cout << "Code: " << '\n'
              << code << '\n';

    gf2Vec err = GetParam();
    err.push_back('0');
    gf2Vec const syndr = code.getXSyndrome(err);
    std::cout << "Syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << '\n';

    testAllDecoders(code, err, syndr, true);
}