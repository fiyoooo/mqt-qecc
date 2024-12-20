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
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, true, peel);

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
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, true, peel);

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    dch.testWithErrorRate(prob, 100);
}

TEST_P(RandomErrorDecoderComparison, ToricCode8Test) {
    auto       code = ToricCode8();
    bool const peel = true;
    std::cout << "Code: " << '\n'
              << code << '\n';
    DecoderComparisonHelper dch = DecoderComparisonHelper(code, true, peel);

    double const prob = GetParam();
    std::cout << "Bit-error rate: " << prob << '\n';

    dch.testWithErrorRate(prob, 100);
}

// TODO integrate in DecoderComparison?
class DecoderErrorRateTest : public ::testing::TestWithParam<std::reference_wrapper<Code>> {
protected:
    std::vector<std::string> decoders           = {"UF", "MaxSAT", "BP", "OSD", "LSD", "UF-peel", "UF-matrix", "UF-maxSAT"};
    std::vector<double>      errorProbabilitiesSmall = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01};
    std::vector<double>      errorProbabilities = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1};
    int                      numRounds          = 1000;
    bool                     maxsat             = true; // whether to use the maxsat decoder
    bool                     peel               = true; // whether to use the peel decoder

    std::vector<std::vector<double>> errorRates    = std::vector<std::vector<double>>(decoders.size());
    std::vector<std::vector<double>> avgRuntimes   = std::vector<std::vector<double>>(decoders.size());
    std::vector<int>                 realNumRounds = std::vector<int>(errorProbabilities.size());

    void runTestForDecoders(Code& code) {
        DecoderComparisonHelper dch = DecoderComparisonHelper(code, maxsat, peel);
        for (size_t i = 0; i < decoders.size(); ++i) {
            errorRates[i].resize(errorProbabilities.size());
            avgRuntimes[i].resize(errorProbabilities.size());
        }

        // for each physical error rate...
        for (size_t i = 0; i < errorProbabilities.size(); ++i) {
            double const prob = errorProbabilities[i];

            // ... run all decoders...
            auto [testErrorRates, runtimes, realNumRound] = dch.testWithErrorRate(prob, numRounds, true);
            // ... and store the outcome
            for (size_t j = 0; j < errorRates.size(); ++j) {
                errorRates[j][i]  = testErrorRates[j];
                avgRuntimes[j][i] = runtimes[j];
            }
            realNumRounds[i] = realNumRound;
        }
    }

    void storeResultsAsJson(size_t distance, const std::string& codeParams) {
        for (size_t i = 0; i < decoders.size(); ++i) {
            json results = json::array();

            for (size_t j = 0; j < errorProbabilities.size(); ++j) {
                json entry;
                entry["code"]                   = codeParams;
                entry["d"]                      = distance;
                entry["p"]                      = errorProbabilities[j];
                entry["logical_error_rates"]    = errorRates[i][j];
                entry["logical_error_rate_ebs"] = sqrt((1 - errorRates[i][j]) * errorRates[i][j] / realNumRounds[j]); // error bars
                entry["avg_total_time"]         = avgRuntimes[i][j];
                entry["min_wts_logical_err"]    = 0; // TODO change, what does it mean?

                results.push_back(entry);
            }

            // create directory for solver if not existent and save results
            std::string const dir = std::filesystem::absolute("../../test/decoder_results/" + decoders[i]).string(); // TODO change, make prettier
            std::filesystem::create_directories(dir);

            // save file
            std::string const filename = dir + "/results_" + codeParams + ".json";
            std::ofstream     file(filename);
            if (!file) {
                std::cerr << "Error opening file: " << filename << '\n';
                continue; // skip to next decoder if file opening fails
            }

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
        for (size_t i = 0; i < decoders.size(); ++i) {
            if ((i != 1 || maxsat) && (i != 5 || peel)) {
                std::cout << decoders[i] << "_err_rates = ";
                printPythonVector(errorRates[i]);
            }
        }
    }
};

// toric codes X
ToricXCode8  toricX8;
ToricXCode18 toricX18;
ToricXCode32 toricX32;

// even toric codes
ToricCode8   toric8;
ToricCode32  toric32;
ToricCode72  toric72;
ToricCode128 toric128;
ToricCode200 toric200;

// odd toric codes
ToricCode18  toric18;
ToricCode50  toric50;
ToricCode98  toric98;
ToricCode162 toric162;

// ldpc codes
BivarBikeCode72  bb72; // pure maxSAT throws error, TODO investigate
BivarBikeCode90  bb90;
BivarBikeCode144 bb144; // pure maxsat runs too long
BivarBikeCode288 bb288; // runs too long in general

INSTANTIATE_TEST_SUITE_P(ToricAndBBCodes, DecoderErrorRateTest,
                         ::testing::Values(
                                 // std::ref(toricX8), std::ref(toricX18), std::ref(toricX32)
                                 // std::ref(toric8), std::ref(toric32) //, std::ref(toric72), std::ref(toric128), std::ref(toric200)
                                 std::ref(toric50), std::ref(toric98), std::ref(toric162)
                                 //std::ref(bb72), std::ref(bb90), std::ref(bb144)//, std::ref(bb288)
                                 ));

TEST_P(DecoderErrorRateTest, RunErrorProbabilityTest) {
    Code& code = GetParam().get();
    maxsat     = false; // set to false if code too big
    peel       = false; // set to false if not usable
    runTestForDecoders(code);
    std::string const codeParams = "[[" + std::to_string(code.n) + "," + std::to_string(code.k) + "," + std::to_string(code.d) + "]]";
    storeResultsAsJson(code.d, codeParams);
    printResults();
}
