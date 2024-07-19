#pragma once

#include "QeccException.hpp"
#include "TreeNode.hpp"
#include "nlohmann/json.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <vector>

using gf2Mat = std::vector<std::vector<bool>>;
using gf2Vec = std::vector<bool>;
using intMat = std::vector<std::vector<int>>;

class Utils {
public:
    /**
     * Solve the system given by Mx=b
     * Returns x if there is a solution, or an empty vector if there is no solution
     * If there are multiple valid solutions, one is returned
     * @param inmat
     * @param vec
     * @return
     */
    static gf2Vec solveSystem(const gf2Mat& inmat, const gf2Vec& vec) {
        assertMatrixPresent(inmat);
        assertVectorPresent(vec);

        if (inmat.size() != vec.size()) {
            throw QeccException("Cannot solve system, dimensions do not match");
        }

        size_t const rowCount = inmat.size();
        size_t const colCount = inmat[0].size();

        // create augmented matrix [inmat|vec]
        gf2Mat augmentedMatrix = inmat;
        for (size_t i = 0; i < augmentedMatrix.size(); ++i) {
            augmentedMatrix[i].push_back(vec[i]);
        }

        // gaussian elimination of augmented matrix
        auto reducedMatrix = gauss(augmentedMatrix);

        // check for consistency & extract solution if it exists
        gf2Vec result(colCount, false);
        for (size_t i = 0; i < rowCount; ++i) {
            bool allZeros = true;
            for (size_t j = 0; j < colCount; ++j) {
                if (reducedMatrix[i][j] == 1) {
                    allZeros = false;
                    result[j] = reducedMatrix[i][colCount] % 2;
                    break;
                }
            }
            if (allZeros && reducedMatrix[i][colCount] == 1) {
                // inconsistent system
                return {};
            }
        }

        return result;
    }

    static intMat gauss(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);
        auto res = convertToIntMatrix(matrix);
        res      = getRREF(res); // reduced row echelon form
        return res;
    }

    static intMat convertToIntMatrix(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);

        intMat result(matrix.size(), std::vector<int>(matrix.front().size(), 0));
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.front().size(); j++) {
                result[i][j] = static_cast<int>(matrix[i][j]);
            }
        }
        return result;
    }

    static intMat getRREF(intMat& matrix) {
        size_t lead     = 0;
        size_t const rowCount = matrix.size();
        size_t const colCount = matrix[0].size();

        for (size_t r = 0; r < rowCount; ++r) {
            if (lead >= colCount) {
                return matrix;
            }
            size_t i = r;
            while (matrix[i][lead] == 0) {
                ++i;
                if (i == rowCount) {
                    i = r;
                    ++lead;
                    if (lead == colCount) {
                        return matrix;
                    }
                }
            }
            std::swap(matrix[i], matrix[r]);
            int lv = matrix[r][lead];
            for (size_t j = 0; j < colCount; ++j) {
                matrix[r][j] /= lv;
            }
            for (size_t i = 0; i < rowCount; ++i) {
                if (i != r) {
                    lv = matrix[i][lead];
                    for (size_t j = 0; j < colCount; ++j) {
                        matrix[i][j] -= lv * matrix[r][j];
                    }
                }
            }
        }
        return matrix;
    }

    /**
     * Checks if the given vector is in the rowspace of matrix M
     * @param inmat
     * @param vec
     * @return
     */
    static bool isVectorInRowspace(const gf2Mat& inmat, const gf2Vec& vec) {
        assertMatrixPresent(inmat);
        assertVectorPresent(vec);

        if (std::all_of(vec.begin(), vec.end(), [](bool val) { return !val; })) { // all zeros vector trivial
            return true;
        }

        if (vec.size() != inmat.at(0).size()) {
            throw QeccException("Cannot check if in rowspace, dimensions of matrix and vector do not match");
        }

        gf2Mat matrix = getTranspose(inmat); // v is in rowspace of M <=> v is in col space of M^T
        for (std::size_t i = 0; i < matrix.size(); i++) {
            matrix.at(i).emplace_back(vec.at(i));
        }
        auto reduced = gauss(matrix);
        //  check consistency, inconsistent <=> vec not in rowspace
        for (const auto& row : reduced) {
            if (row.back() == 1) {
                if (std::all_of(row.begin(), row.end() - 1, [](int val) { return val == 0; })) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Computes the transpose of the given matrix
     * @param matrix
     * @return
     */
    static gf2Mat getTranspose(const gf2Mat& matrix) {
        assertMatrixPresent(matrix);

        gf2Mat transp(matrix.at(0).size(), gf2Vec(matrix.size()));
        for (std::size_t i = 0; i < matrix.size(); i++) {
            for (std::size_t j = 0; j < matrix.at(i).size(); j++) {
                transp.at(j).at(i) = matrix.at(i).at(j);
            }
        }
        return transp;
    }

    /**
     * Computes matrix vector product and stores it in result vector
     * @param m1
     * @param vec
     * @param result
     */
    static void rectMatrixMultiply(const gf2Mat& m1, const gf2Vec& vec, gf2Vec& result) {
        assertMatrixPresent(m1);
        assertVectorPresent(vec);

        if (m1.front().size() != vec.size() || m1.size() > result.capacity()) {
            throw QeccException("Cannot multiply, dimensions wrong");
        }
        for (std::size_t i = 0; i < m1.size(); i++) {
            const auto& row = m1.at(i);
            for (std::size_t k = 0; k < vec.size(); k++) {
                result.at(i) = result.at(i) ^ (row.at(k) && vec.at(k));
            }
        }
    }

    static void assertMatrixPresent(const gf2Mat& matrix) {
        if (matrix.empty() || matrix.at(0).empty()) {
            throw QeccException("Matrix is empty");
        }
    }

    static void assertVectorPresent(const gf2Vec& vector) {
        if (vector.empty()) {
            throw QeccException("Vector is empty");
        }
    }

    static std::string getStringFrom(const gf2Mat& matrix) {
        if (matrix.empty()) {
            return "[]";
        }

        const auto&       nrows = matrix.size();
        const auto&       ncols = matrix.at(0).size();
        std::stringstream s;
        s << nrows << "x" << ncols << "matrix [" << '\n';
        for (std::size_t i = 0; i < nrows; i++) {
            s << "[";
            for (std::size_t j = 0; j < ncols; j++) {
                s << matrix.at(i).at(j);
                if (j != ncols - 1) {
                    s << ",";
                }
            }
            s << "]";
            if (i != nrows - 1) {
                s << ",";
            }
            s << '\n';
        }
        s << "]";
        return s.str();
    }

    static std::string getStringFrom(const gf2Vec& vector) {
        if (vector.empty()) {
            return "[]";
        }

        const auto&       nelems = vector.size();
        std::stringstream s;
        s << "[";
        for (std::size_t j = 0; j < nelems; j++) {
            s << vector.at(j);
            if (j != nelems - 1) {
                s << ",";
            }
        }
        s << "]";
        return s.str();
    }
};