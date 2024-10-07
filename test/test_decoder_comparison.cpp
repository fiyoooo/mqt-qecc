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
#include <random>
#include <tuple>

class DecoderComparison : public testing::TestWithParam<gf2Vec> {};
class RandomDecoderComparison : public testing::TestWithParam<double> {};

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

INSTANTIATE_TEST_SUITE_P(firstOrderErrorRate, RandomDecoderComparison,
                         ::testing::Values(
                                 0,
                                 0.01,
                                 0.02,
                                 0.03,
                                 0.04,
                                 0.05,
                                 0.06,
                                 0.07,
                                 0.08,
                                 0.09,
                                 0.1));

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
    ldpc::bp::BpSparse         pcm      = code.toBpSparse();
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

gf2Vec generateRandomError(size_t size, double prob) {
    std::random_device          rd;
    std::mt19937                gen(rd());
    std::bernoulli_distribution dist(prob); // prob is probability for a bit to be 1 (to have flipped)

    gf2Vec err(size);
    for (size_t i = 0; i < size; ++i) {
        err[i] = dist(gen);
    }
    return err;
}

std::tuple<std::vector<double>, std::vector<double>> compareDecodersWithRandomErrors(Code& code, double prob, int numRounds, bool peel) {
    int failuresUFD       = 0;
    int failuresMS        = 0;
    int failuresUfdPeel   = 0;
    int failuresUfdMatrix = 0;
    int failuresUfdMaxSat = 0;

    double totalRuntimeUFD       = 0;
    double totalRuntimeMS        = 0;
    double totalRuntimeUfdPeel   = 0;
    double totalRuntimeUfdMatrix = 0;
    double totalRuntimeUfdMaxSat = 0;

    int zeroBitErrors = 0;
    int oneBitErrors  = 0;
    int twoBitErrors  = 0;
    int otherErrors   = 0;

    size_t const size = code.getN();

    for (int i = 0; i < numRounds; ++i) {
        gf2Vec const err   = generateRandomError(size, prob);
        gf2Vec       syndr = code.getXSyndrome(err);

        int const bitCount = std::count(err.begin(), err.end(), 1);

        if (bitCount == 0) {
            ++zeroBitErrors;
        } else if (bitCount == 1) {
            ++oneBitErrors;
        } else if (bitCount == 2) {
            ++twoBitErrors;
        } else {
            ++otherErrors;
        }

        auto [runtimeUFD, estimUFD] = measureRuntime([&]() { return runUFDecoder(code, syndr); });
        totalRuntimeUFD += runtimeUFD;
        if (estimUFD != err) {
            ++failuresUFD;
        }

        auto [runtimeMaxSat, estimMaxSat] = measureRuntime([&]() { return runMaxSatDecoder(code, syndr); });
        totalRuntimeMS += runtimeMaxSat;
        if (estimMaxSat != err) {
            ++failuresMS;
        }

        ldpc::bp::BpSparse   pcm      = code.toBpSparse();
        std::vector<uint8_t> syndrInt = Utils::toUint8(syndr);

        if (peel) {
            auto [runtimeUfdPeel, estimUfdPeelInt] = measureRuntime([&]() { return runUfdPeelDecoder(pcm, syndrInt); });
            gf2Vec const estimUfdPeel              = Utils::toGf2Vec(estimUfdPeelInt);
            totalRuntimeUfdPeel += runtimeUfdPeel;
            if (estimUfdPeel != err) {
                ++failuresUfdPeel;
            }
        }

        auto [runtimeUfdMatrix, estimUfdMatrixInt] = measureRuntime([&]() { return runUfdMatrixDecoder(pcm, syndrInt); });
        gf2Vec const estimUfdMatrix                = Utils::toGf2Vec(estimUfdMatrixInt);
        totalRuntimeUfdMatrix += runtimeUfdMatrix;
        if (estimUfdMatrix != err) {
            ++failuresUfdMatrix;
        }

        auto [runtimeUfdMaxSat, estimUfdMaxSatInt] = measureRuntime([&]() { return runUfdMaxSatDecoder(pcm, syndrInt); });
        gf2Vec const estimUfdMaxSat                = Utils::toGf2Vec(estimUfdMaxSatInt);
        totalRuntimeUfdMaxSat += runtimeUfdMaxSat;
        if (estimUfdMaxSat != err) {
            ++failuresUfdMaxSat;
        }
    }

    // Output results
    std::cout << "Total random errors tested: " << numRounds << "\n";
    std::cout << "Zero-bit, one-bit, two-bit or other errors: "
              << zeroBitErrors << ", "
              << oneBitErrors << ", "
              << twoBitErrors << ", "
              << otherErrors << "\n\n";

    double const logErrRateUFD       = static_cast<double>(failuresUFD) / numRounds;
    double const logErrRateMS        = static_cast<double>(failuresMS) / numRounds;
    double const logErrRateUfdPeel   = peel ? static_cast<double>(failuresUfdPeel) / numRounds : 0;
    double const logErrRateUfdMatrix = static_cast<double>(failuresUfdMatrix) / numRounds;
    double const logErrRateUfdMaxSat = static_cast<double>(failuresUfdMaxSat) / numRounds;

    std::cout << "UFD, MaxSatD, ldpc::uf::UfD peel, matrix, maxsat:"
              << "\n";
    std::cout << "Logical-error rate: "
              << logErrRateUFD << ", "
              << logErrRateMS << ", "
              << logErrRateUfdPeel << ", "
              << logErrRateUfdMatrix << ", "
              << logErrRateUfdMaxSat << "\n";

    double const avgRuntimeUFD       = totalRuntimeUFD / numRounds;
    double const avgRuntimeMS        = totalRuntimeMS / numRounds;
    double const avgRuntimeUfdPeel   = peel ? (totalRuntimeUfdPeel / numRounds) : 0;
    double const avgRuntimeUfdMatrix = totalRuntimeUfdMatrix / numRounds;
    double const avgRuntimeUfdMaxSat = totalRuntimeUfdMaxSat / numRounds;

    std::cout << "Average runtime: "
              << avgRuntimeUFD << ", "
              << avgRuntimeMS << ", "
              << avgRuntimeUfdPeel << ", "
              << avgRuntimeUfdMatrix << ", "
              << avgRuntimeUfdMaxSat << "\n\n";

    std::vector<double> const logErrRate  = {logErrRateUFD, logErrRateMS, logErrRateUfdPeel, logErrRateUfdMatrix, logErrRateUfdMaxSat};
    std::vector<double> const avgRuntimes = {avgRuntimeUFD, avgRuntimeMS, avgRuntimeUfdPeel, avgRuntimeUfdMatrix, avgRuntimeUfdMaxSat};

    return std::make_tuple(logErrRate, avgRuntimes);
}

void printPythonVector(std::vector<double>& vec) {
    std::cout << '[';
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ']' << '\n';
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
    err.push_back(0);
    gf2Vec const syndr = code.getXSyndrome(err);
    std::cout << "Syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << '\n';

    testAllDecoders(code, err, syndr, true);
}

TEST_P(RandomDecoderComparison, SteaneXCodeTest) {
    auto code = SteaneXCode();
    std::cout << "Code: " << '\n'
              << code << '\n';

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    compareDecodersWithRandomErrors(code, prob, 100, false);
}

TEST_P(RandomDecoderComparison, ToricCode8Test) {
    auto code = ToricCode8();
    std::cout << "Code: " << '\n'
              << code << '\n';

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    compareDecodersWithRandomErrors(code, prob, 100, true);
}