#include "DecoderComparison.hpp"

#include "MaxSATDecoder.hpp"
#include "UFDecoder.hpp"
#include "Utils.hpp"
#include "ldpc/osd.hpp"
#include "union_find.hpp"

#include <chrono>
#include <random>

std::tuple<std::vector<gf2Vec>, std::vector<double>>
DecoderComparisonHelper::testWithError(const gf2Vec& err, bool verbose) {
    gf2Vec const syndr = code.getXSyndrome(err);

    if (verbose) {
        std::cout << "Syndrome: ";
        Utils::printGF2vector(syndr);
        std::cout << '\n';

        std::cout << "Sol: ";
        Utils::printGF2vector(err);
    }

    auto [estims, runtimes] = runDecoders(syndr);

    if (verbose) {
        std::cout << "\nEstimUFD: ";
        Utils::printGF2vector(estims[0]);
        if (maxsat) {
            std::cout << "\nEstimMS: ";
            Utils::printGF2vector(estims[1]);
        }
        if (peel) {
            std::cout << "\nEstimUfdPeel: ";
            Utils::printGF2vector(estims[2]);
        }
        std::cout << "\nEstimUfdMatrix: ";
        Utils::printGF2vector(estims[3]);
        std::cout << "\nEstimUfdMaxSat: ";
        Utils::printGF2vector(estims[4]);
        std::cout << "\n\n";

        std::cout << "Runtime (UFDecoder): " << runtimes[0] << " seconds\n";
        if (maxsat) {
            std::cout << "Runtime (MaxSatDecoder): " << runtimes[1] << " seconds\n";
        }
        if (peel) {
            std::cout << "Runtime (ldpc::uf::UfDecoder peel_decode): " << runtimes[2] << " seconds\n";
        }
        std::cout << "Runtime (ldpc::uf::UfDecoder matrix_decode): " << runtimes[3] << " seconds\n";
        std::cout << "Runtime (ldpc::uf::UfDecoder maxsat_decode): " << runtimes[4] << " seconds\n\n";
    }

    return std::make_tuple(estims, runtimes);
}

std::tuple<std::vector<double>, std::vector<double>>
DecoderComparisonHelper::testWithErrorRate(double prob, int numRounds, bool verbose) {
    size_t const        numDecoders = 6; // bc we test 6 decoders
    std::vector<int>    failures(numDecoders);
    std::vector<double> totalRuntimes(numDecoders);

    std::vector<int> bitErrorCounts(4); // count 0-bit, 1-bit, 2-bit, and 3+ bit errors

    size_t const size = code.getN();

    // test & track the failures
    int const minErrors = 500; // minimal number of non-zero errors
    int       i;
    for (i = 0; (i < numRounds) || (i - bitErrorCounts[0] < minErrors); ++i) {
        gf2Vec const err   = generateRandomBitFlipError(size, prob);
        gf2Vec const syndr = code.getXSyndrome(err);

        // count 0-bit, 1-bit, 2-bit, and 3+ bit errors
        size_t const bitCount = static_cast<const size_t>(std::count(err.begin(), err.end(), 1));
        if (bitCount <= 2) {
            ++bitErrorCounts[bitCount]; // increment count for 0, 1, or 2 bit flips
        } else {
            ++bitErrorCounts[3]; // increment count for 3 or more bit flips
        }

        // run all decoders
        auto [estims, runtimes] = runDecoders(syndr, prob);
        for (size_t j = 0; j < numDecoders; ++j) {
            totalRuntimes[j] += runtimes[j];
            if ((j != 1 || maxsat) && (j != 3 || peel) && !isCorrectable(err, estims[j])) {
                ++failures[j];
            }
        }
    }

    // get averages by dividing through number of rounds
    std::vector<double> errorRates;
    std::vector<double> avgRuntimes;
    errorRates.reserve(failures.size());
    for (int const& failure : failures) {
        errorRates.push_back(static_cast<double>(failure) / i);
    }
    avgRuntimes.reserve(totalRuntimes.size());
    for (double const& runtime : totalRuntimes) {
        avgRuntimes.push_back(runtime / i);
    }

    if (verbose) {
        std::cout << "Total random errors tested: " << i << " with physical error rate " << prob << "\n";
        std::cout << "Zero-bit, one-bit, two-bit or other errors: ";
        for (const auto& bitCount : bitErrorCounts) {
            std::cout << bitCount << " ";
        }
        std::cout << "\n\n";

        std::cout << "UFD, MaxSatD, (ldpc::osd) OSDD, (ldpc::uf::UfD) peel, matrix, maxsat:\n";
        std::cout << "Failures: ";
        for (const auto& failure : failures) {
            std::cout << failure << "  ";
        }
        std::cout << "\nError rates: ";
        for (const auto& errorRate : errorRates) {
            std::cout << errorRate << "  ";
        }
        std::cout << "\nAverage runtime: ";
        for (const auto& runtime : avgRuntimes) {
            std::cout << runtime << "  ";
        }
        std::cout << "\n\n";
    }

    return std::make_tuple(errorRates, avgRuntimes);
}

gf2Vec DecoderComparisonHelper::generateRandomBitFlipError(size_t size, double prob) {
    std::random_device          rd;
    std::mt19937                gen(rd());
    std::bernoulli_distribution dist(prob); // prob is probability for a bit to be 1 (to have flipped)

    gf2Vec err(size);
    for (size_t i = 0; i < size; ++i) {
        err[i] = dist(gen);
    }
    return err;
}

std::tuple<std::vector<gf2Vec>, std::vector<double>> DecoderComparisonHelper::runDecoders(gf2Vec const& syndr, double prob) {
    auto [runtimeUFD, estimUFD] = measureRuntime([&]() { return runUFDecoder(syndr); });

    gf2Vec estimMS;
    double runtimeMS = 0;
    if (maxsat) {
        std::tie(runtimeMS, estimMS) = measureRuntime([&]() { return runMaxSatDecoder(syndr); });
    }

    // transform syndrome to right format
    std::vector<uint8_t> syndrInt = Utils::toUint8(syndr);

    auto [runtimeOsd, estimOsdInt] = measureRuntime([&]() { return runOsdDecoder(syndrInt, prob); });
    gf2Vec const estimOsd          = Utils::toGf2Vec(estimOsdInt);

    gf2Vec estimUfdPeel;
    double runtimeUfdPeel = 0;
    if (peel) {
        auto [runtimeUfdPeel_, estimUfdPeelInt] = measureRuntime([&]() { return runUfdPeelDecoder(syndrInt); });
        runtimeUfdPeel                          = runtimeUfdPeel_;
        estimUfdPeel                            = Utils::toGf2Vec(estimUfdPeelInt);
    }

    auto [runtimeUfdMatrix, estimUfdMatrixInt] = measureRuntime([&]() { return runUfdMatrixDecoder(syndrInt); });
    gf2Vec const estimUfdMatrix                = Utils::toGf2Vec(estimUfdMatrixInt);

    auto [runtimeUfdMaxSat, estimUfdMaxSatInt] = measureRuntime([&]() { return runUfdMaxSatDecoder(syndrInt); });
    gf2Vec const estimUfdMaxSat                = Utils::toGf2Vec(estimUfdMaxSatInt);

    std::vector<gf2Vec> const estims   = {estimUFD, estimMS, estimOsd, estimUfdPeel, estimUfdMatrix, estimUfdMaxSat};
    std::vector<double> const runtimes = {runtimeUFD, runtimeMS, runtimeOsd, runtimeUfdPeel, runtimeUfdMatrix, runtimeUfdMaxSat};

    return std::make_tuple(estims, runtimes);
}

gf2Vec DecoderComparisonHelper::runUFDecoder(gf2Vec const& syndrome) {
    UFDecoder ufDecoder(code);
    ufDecoder.decode(syndrome);
    return ufDecoder.result.estimBoolVector;
}

gf2Vec DecoderComparisonHelper::runMaxSatDecoder(gf2Vec const& syndrome) {
    MaxSATDecoder maxSATSolver(code);
    maxSATSolver.preconstructZ3Instance();
    maxSATSolver.decode(syndrome);
    return maxSATSolver.result.estimBoolVector;
}

std::vector<uint8_t> DecoderComparisonHelper::runOsdDecoder(std::vector<uint8_t>& syndrome, double prob) {
    auto                  errorChannel = std::vector<double>(pcm.n, prob);
    ldpc::osd::OsdDecoder osdDecoderLdpc(pcm, ldpc::osd::OSD_0, 0, errorChannel);
    std::vector<double>   lbr; // log-likelihood ratios
    lbr.reserve(pcm.n);
    for (auto err : errorChannel) {
        lbr.push_back(log((1 - err) / err));
    }
    return osdDecoderLdpc.decode(syndrome, lbr);
}

std::vector<uint8_t> DecoderComparisonHelper::runUfdPeelDecoder(std::vector<uint8_t> const& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.peel_decode(syndrome);
}

std::vector<uint8_t> DecoderComparisonHelper::runUfdMatrixDecoder(std::vector<uint8_t> const& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.matrix_decode(syndrome);
}

std::vector<uint8_t> DecoderComparisonHelper::runUfdMaxSatDecoder(std::vector<uint8_t> const& syndrome) {
    ldpc::uf::UfDecoder ufDecoderLdpc(pcm);
    return ufDecoderLdpc.maxsat_decode(syndrome);
}

template <typename Func, typename ReturnType>
std::tuple<double, ReturnType> DecoderComparisonHelper::measureRuntime(Func func) {
    auto                                start    = std::chrono::high_resolution_clock::now();
    ReturnType                          estim    = func();
    auto                                end      = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const duration = end - start;
    return std::make_tuple(duration.count(), estim);
}

bool DecoderComparisonHelper::isCorrectable(gf2Vec const& error, gf2Vec const& estimation) {
    gf2Vec residualErr(error.size());
    for (size_t i = 0; i < error.size(); i++) {
        residualErr.at(i) = (error[i] != estimation[i]);
    }
    return Utils::isVectorInRowspace(*code.gethZ()->pcm, residualErr);
}
