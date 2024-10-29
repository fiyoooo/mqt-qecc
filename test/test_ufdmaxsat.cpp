//
// Created by Fiona Froehler on 12.08.24.
//
// to keep 0/1 in boolean areas without clang-tidy warnings:
// NOLINTBEGIN(readability-implicit-bool-conversion,modernize-use-bool-literals)

#include "Codes.hpp"
#include "UFDecoder.hpp"
#include "ldpc/gf2codes.hpp"
#include "union_find.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace ldpc::uf;

// for changing data types between decoder packages
class VectorTranslationTest : public testing::TestWithParam<gf2Vec> {};
class MatrixTranslationTest : public testing::TestWithParam<std::reference_wrapper<Code>> {};

class UFDMaxSATTest : public testing::TestWithParam<std::vector<uint8_t>> {};
class UFDComparison : public testing::TestWithParam<gf2Vec> {};
class UniquelyCorrectableErrTestUFDMaxSAT : public UFDMaxSATTest {};
class InCorrectableErrTestUFDMaxSAT : public UFDMaxSATTest {};
class UpToStabCorrectableErrTestUFDMaxSAT : public UFDMaxSATTest {};

// toric codes X
ToricXCode8  tX8;
ToricXCode18 tX18;
ToricXCode32 tX32;

// even toric codes
ToricCode8  t8;
ToricCode32 t32;
ToricCode72 t72;

// odd toric codes
ToricCode18 t18;
ToricCode50 t50;
ToricCode98 t98;

// ldpc codes
BivarBikeCode72  bivb72; // pure maxSAT throws error
BivarBikeCode90  bivb90;
BivarBikeCode144 bivb144; // pure maxsat runs too long
BivarBikeCode288 bivb288; // pure maxsat runs too long

INSTANTIATE_TEST_SUITE_P(ToricCodes, MatrixTranslationTest,
                         ::testing::Values(
                                 std::ref(tX8), std::ref(tX18), std::ref(tX32),
                                 std::ref(t8), std::ref(t32), std::ref(t72),
                                 std::ref(t18), std::ref(t50), std::ref(t98),
                                 std::ref(bivb72), std::ref(bivb90), std::ref(bivb144), std::ref(bivb288)));

INSTANTIATE_TEST_SUITE_P(UFDComparisonSteane, UFDComparison,
                         testing::Values(
                                 gf2Vec{0, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{1, 0, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 1, 0, 0, 0, 0, 0},
                                 gf2Vec{0, 0, 1, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(CorrectableSingleBitErrsSteane, UniquelyCorrectableErrTestUFDMaxSAT,
                         testing::Values(
                                 std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 0},
                                 std::vector<uint8_t>{1, 0, 0, 0, 0, 0, 0},
                                 std::vector<uint8_t>{0, 1, 0, 0, 0, 0, 0},
                                 std::vector<uint8_t>{0, 0, 1, 0, 0, 0, 0},
                                 std::vector<uint8_t>{0, 0, 0, 1, 0, 0, 0},
                                 std::vector<uint8_t>{0, 0, 0, 0, 1, 0, 0},
                                 std::vector<uint8_t>{0, 0, 0, 0, 0, 1, 0},
                                 std::vector<uint8_t>{0, 0, 0, 0, 0, 0, 1}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, InCorrectableErrTestUFDMaxSAT,
                         testing::Values(
                                 std::vector<uint8_t>{0, 1, 0, 0, 0, 0, 1},
                                 std::vector<uint8_t>{0, 0, 1, 0, 1, 0, 0},
                                 std::vector<uint8_t>{0, 0, 1, 1, 0, 0, 0},
                                 std::vector<uint8_t>{0, 0, 1, 0, 0, 1, 0},
                                 std::vector<uint8_t>{1, 1, 0, 0, 0, 0, 0}));

INSTANTIATE_TEST_SUITE_P(UptoStabCorrSteane, UpToStabCorrectableErrTestUFDMaxSAT,
                         testing::Values(
                                 std::vector<uint8_t>{1, 0, 1, 0, 1, 0, 0}));

gf2Vec generateRandomGf2Vec(size_t size) {
    std::random_device          rd;
    std::mt19937                gen(rd());
    std::bernoulli_distribution dist(0.5); // prob is probability for a bit to be 1 (to have flipped)

    gf2Vec err(size);
    for (size_t i = 0; i < size; ++i) {
        err[i] = dist(gen);
    }
    return err;
}

void printPcm(ldpc::bp::BpSparse& pcm) {
    for (int j = 0; j < pcm.m; j++) {
        for (int i = 0; i < pcm.n; i++) {
            // if no entry in that position of matrix
            if (pcm.get_entry(j, i).at_end()) {
                std::cout << "0";
            } else {
                std::cout << "1";
            }
        }
        std::cout << '\n';
    }
}

TEST(VectorTranslationTest, CompareVectors) {
    for (int i = 0; i < 10000; ++i) {
        gf2Vec const vec = generateRandomGf2Vec(100);
        Utils::printGF2vector(vec);
        auto vecInt = Utils::toUint8(vec);
        for (const auto& j : vecInt) {
            std::cout << j;
        }
        std::cout << '\n';
        EXPECT_TRUE(vec == Utils::toGf2Vec(vecInt));
    }
}

TEST_P(MatrixTranslationTest, CompareMatrices) {
    Code& code = GetParam();
    std::cout << "Code: " << '\n'
              << code << '\n';
    auto pcm = code.hZToBpSparse();
    printPcm(pcm);

    EXPECT_TRUE(code.getHzMat() == Code::toGf2Mat(pcm));
}

/**
 * Tests for comparing UfDecoder of mqt and ldpc.
 */
TEST_P(UFDComparison, ComparingWithOriginalUfd) {
    auto code = SteaneXCode();
    std::cout << "Code: " << '\n'
              << code << '\n';

    gf2Vec const err   = GetParam();
    gf2Vec const syndr = code.getXSyndrome(err);
    std::cout << "Syndrome: ";
    Utils::printGF2vector(syndr);
    std::cout << '\n';

    // initialize original UFDecoder
    UFDecoder ufDecoder;
    ufDecoder.setCode(code);
    ufDecoder.decode(syndr);
    const auto&   decodingResult = ufDecoder.result;
    const gf2Vec& estim          = decodingResult.estimBoolVector;
    const auto&   estimIdx       = decodingResult.estimNodeIdxVector;

    // transform pcm and syndrome to right format
    ldpc::bp::BpSparse pcm(code.getK(), code.getN());
    for (int i = 0; i < pcm.m; ++i) {
        gf2Vec const row = code.gethZ()->pcm->operator[](i);
        for (int j = 0; j < row.size(); ++j) {
            if (row[j]) {
                pcm.insert_entry(i, j);
            }
        }
    }
    std::vector<uint8_t> syndrInt;
    for (bool i : syndr) {
        syndrInt.push_back(static_cast<uint8_t>(i));
    }

    // initialize UfDecoder of ldpc
    ldpc::uf::UfDecoder         ufDecoderLdpc(pcm);
    const std::vector<uint8_t>& estim2Int = ufDecoderLdpc.maxsat_decode(syndrInt);
    gf2Vec                      estim2;
    for (auto i : estim2Int) {
        estim2.push_back(static_cast<bool>(i));
    }

    // compare results
    std::cout << "\nEstim: ";
    Utils::printGF2vector(estim);
    std::cout << '\n';

    std::cout << "\nEstim of ldpc: ";
    Utils::printGF2vector(estim2);
    std::cout << "\n\n";

    EXPECT_EQ(estim, estim2);
}

/**
 * Tests for unambiguous syndromes, estimates must be computed exactly
 */
TEST_P(UniquelyCorrectableErrTestUFDMaxSAT, SteaneCodeDecodingTestEstim) {
    auto pcm = ldpc::gf2codes::hamming_code<ldpc::bp::BpEntry>(3);

    std::cout << "Pcm: " << '\n';
    printPcm(pcm);

    std::vector<uint8_t>       err   = GetParam();
    const std::vector<uint8_t> syndr = pcm.mulvec(err); // ldpc::util::decimal_to_binary(3, 3);

    std::cout << "Syndrome: ";
    for (const auto& i : syndr) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    auto ufd   = UfDecoder(pcm);
    auto estim = ufd.maxsat_decode(syndr);

    std::cout << "Estim: ";
    for (const auto& i : estim) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    std::cout << "Sol: ";
    for (const auto& i : err) {
        std::cout << static_cast<int>(i);
    }
    std::cout << "\n\n";

    EXPECT_TRUE(estim == err);
}

/**
 * Tests for ambiguous errors that cannot be corrected
 */
TEST_P(InCorrectableErrTestUFDMaxSAT, SteaneCodeDecodingTestEstim) {
    auto pcm = ldpc::gf2codes::hamming_code<ldpc::bp::BpEntry>(3);

    std::cout << "Pcm: " << '\n';
    printPcm(pcm);

    std::vector<uint8_t>       err   = GetParam();
    const std::vector<uint8_t> syndr = pcm.mulvec(err); // ldpc::util::decimal_to_binary(3, 3);

    std::cout << "Syndrome: ";
    for (const auto& i : syndr) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    auto ufd   = UfDecoder(pcm);
    auto estim = ufd.maxsat_decode(syndr);

    std::cout << "Estim: ";
    for (const auto& i : estim) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    std::cout << "Sol: ";
    for (const auto& i : err) {
        std::cout << static_cast<int>(i);
    }
    std::cout << "\n\n";

    EXPECT_FALSE(estim == err);
}

/**
 * Tests for errors that are correctable up to stabilizer
 */
TEST_P(UpToStabCorrectableErrTestUFDMaxSAT, SteaneCodeDecodingTest) {
    auto pcm = ldpc::gf2codes::hamming_code<ldpc::bp::BpEntry>(3);

    std::cout << "Pcm: " << '\n';
    printPcm(pcm);

    std::vector<uint8_t>       err   = GetParam();
    const std::vector<uint8_t> syndr = pcm.mulvec(err); // ldpc::util::decimal_to_binary(3, 3);

    std::cout << "Syndrome: ";
    for (const auto& i : syndr) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    auto ufd   = UfDecoder(pcm);
    auto estim = ufd.maxsat_decode(syndr);

    std::cout << "Estim: ";
    for (const auto& i : estim) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    std::cout << "Sol: ";
    for (const auto& i : err) {
        std::cout << static_cast<int>(i);
    }
    std::cout << '\n';

    std::vector<uint8_t> residualErr(err.size());
    std::cout << "ResidualErr: ";
    for (size_t i = 0; i < err.size(); i++) {
        residualErr.at(i) = (err[i] != estim[i]);
        std::cout << static_cast<int>(residualErr.at(i));
    }
    std::cout << "\n\n";

    EXPECT_FALSE(estim == err);
    EXPECT_TRUE(Utils::isVectorInRowspace(Code::toGf2Mat(pcm), Utils::toGf2Vec(residualErr)));
}
// NOLINTEND(readability-implicit-bool-conversion,modernize-use-bool-literals)
