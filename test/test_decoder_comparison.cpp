#include "Codes.hpp"
#include "DecoderComparison.hpp"

#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::basic_json<>;

class DecoderComparison : public ::testing::TestWithParam<gf2Vec> {};
class RandomErrorDecoderComparison : public ::testing::TestWithParam<double> {};

INSTANTIATE_TEST_SUITE_P(SingleBitErrs, DecoderComparison,
                         ::testing::Values(
                                 gf2Vec{0, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{1, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 1, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(firstOrderErrorRates, RandomErrorDecoderComparison,
                         ::testing::Values(
                                 0,
                                 0.01, 0.02,
                                 0.03, 0.04,
                                 0.05, 0.06,
                                 0.07, 0.08,
                                 0.09, 0.1));

TEST_P(DecoderComparison, SteaneXCodeTest) {
    auto       code = SteaneXCode();
    bool const peel = false; // peel decoder doesn't work
    std::cout << "Code: " << '\n'
              << code << '\n';
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, peel);

    gf2Vec const err        = GetParam();
    auto [estims, runtimes] = dch.testWithError(err);

    // compare results with actual error
    EXPECT_EQ(estims[0], err);
    EXPECT_EQ(estims[1], err);
    // skipping 2 bc peel decoder doesn't work here
    EXPECT_EQ(estims[3], err);
    EXPECT_EQ(estims[4], err);
}

TEST_P(RandomErrorDecoderComparison, SteaneXCodeTest) {
    auto       code = SteaneXCode();
    bool const peel = false;
    std::cout << "Code: " << '\n'
              << code << '\n';
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, peel);

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    dch.testWithErrorRate(prob, 100);
}

TEST_P(RandomErrorDecoderComparison, ToricCode8Test) {
    auto       code = ToricCode8();
    bool const peel = true;
    std::cout << "Code: " << '\n'
              << code << '\n';
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, peel);

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    dch.testWithErrorRate(prob, 100);
}

// TODO integrate in DecoderComparison?
class DecoderErrorRateTest : public ::testing::TestWithParam<std::reference_wrapper<Code>> {
protected:
    std::vector<std::string> decoders           = {"ufDecoder", "maxSatDecoder", "ufdPeelDecoder", "ufdMatrixDecoder", "ufdMaxSatDecoder"};
    std::vector<double>      errorProbabilities = {0.01, 0.02, 0.03, 0.04, 0.05,
                                                   0.06, 0.07, 0.08, 0.09, 0.1};
    int                      numRounds          = 1000;
    bool                     peel               = true; // whether to use the peel decoder

    std::vector<std::vector<double>> errorRates = std::vector<std::vector<double>>(5);

    void runTestForDecoders(Code& code) {
        DecoderComparisonHelper dch = DecoderComparisonHelper(code, peel);
        for (auto& er : errorRates) {
            er.resize(errorProbabilities.size());
        }

        for (size_t i = 0; i < errorProbabilities.size(); ++i) {
            double const prob = errorProbabilities[i];

            auto [testErrorRates, runtimes] = dch.testWithErrorRate(prob, numRounds, false);
            for (size_t j = 0; j < errorRates.size(); ++j) {
                errorRates[j][i] = testErrorRates[j];
            }
        }
    }

    // WIP
    void storeResultsAsJson(int distance) {
        for (size_t i = 0; i < decoders.size(); ++i) {
            json results = json::array();

            for (size_t j = 0; j < errorProbabilities.size(); ++j) {
                json entry;
                entry["distance"]       = distance;
                entry["p"]              = errorProbabilities[j];
                entry["avg_total_time"] = errorRates[i][j]; // Or any other data you want

                results.push_back(entry);
            }

            // create directory for solver if not exists and save results
            std::string const dir = "./decoder_results/" + decoders[i];
            std::filesystem::create_directories(dir);

            // save the file
            std::ofstream file(dir + "/results_" + std::to_string(distance) + ".json");
            file << results.dump(4); // pretty print JSON with 4 spaces indentation
            file.close();
        }
    }

    static void printPythonVector(std::vector<double>& vec) {
        std::cout << '[';
        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << vec[i];
            if (i < vec.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ']' << '\n';
    }

    void printResults() {
        std::cout << "uf_err_rates = ";
        printPythonVector(errorRates[0]);

        std::cout << "ms_err_rates = ";
        printPythonVector(errorRates[1]);

        if (peel) {
            std::cout << "peel_err_rates = ";
            printPythonVector(errorRates[2]);
        }

        std::cout << "matrix_err_rates = ";
        printPythonVector(errorRates[3]);

        std::cout << "ufms_err_rates = ";
        printPythonVector(errorRates[4]);
    }
};

ToricCode8       code8;
ToricCode18      code18;
ToricCode32      code32;
BivarBikeCode72  code72; // maxSAT throws error
BivarBikeCode90  code90;
BivarBikeCode144 code144; // runs too long
BivarBikeCode288 code288; // runs too long

INSTANTIATE_TEST_SUITE_P(ToricCodes, DecoderErrorRateTest,
                         ::testing::Values(
                                 std::ref(code8), std::ref(code18),
                                 std::ref(code32) //,
                                 // std::ref(code72)//,
                                 // std::ref(code90)
                                 // std::ref(code144)
                                 // std::ref(code288)
                                 ));

TEST_P(DecoderErrorRateTest, RunErrorProbabilityTest) {
    Code& code = GetParam().get();
    peel       = false;
    runTestForDecoders(code);
    storeResultsAsJson(4);
    printResults();
}
