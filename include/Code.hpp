//
// Created by lucas on 26/04/2022.
//

#ifndef QUNIONFIND_CODE_HPP
#define QUNIONFIND_CODE_HPP
#include "QeccException.hpp"
#include "TreeNode.hpp"
#include "Utils.hpp"

#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
using json = nlohmann::json;

/**
 * Used as return object in @Code
 */
struct CodeProperties {
    std::size_t n;
    std::size_t k;
    std::size_t d;
};

struct ParityCheckMatrix {
    std::unique_ptr<gf2Mat>                                   pcm;
    std::unordered_map<std::size_t, std::vector<std::size_t>> nbrCache;

    ParityCheckMatrix(const ParityCheckMatrix &m) = delete;
    ParityCheckMatrix & operator= (const ParityCheckMatrix &) = delete;

    explicit ParityCheckMatrix(gf2Mat& pcm):
        pcm(std::make_unique<gf2Mat>(pcm)) {}

    explicit ParityCheckMatrix(const std::string& filePath) {
        if (filePath.empty()) {
            throw QeccException("[PCM::ctor] - Cannot open pcm, filepath empty");
        }
        std::string   line;
        int           word;
        std::ifstream inFile;
        gf2Mat        result;
        //std::cout << "[PCM::ctor] - reading pcm from codefile" << std::endl;
        try {
            inFile.open(filePath);
            while (getline(inFile, line, '\n')) {
                gf2Vec             tempVec;
                std::istringstream instream(line);
                while (instream >> word) {
                    tempVec.push_back(word);
                }
                result.emplace_back(tempVec);
            }
            pcm = std::make_unique<gf2Mat>(result);
        } catch (const std::exception& e) {
            std::cerr << "[PCM::ctor] - error opening file " << filePath << std::endl;
            throw QeccException(e.what());
        }
        inFile.close();
        //std::cout << "[PCM::ctor] - importing from codefile done" << std::endl;
    }

    /**
     * If H is nxm we have n checks and m bit nodes.
     * Indices of bit nodes range from 0 to m-1, and indices of check nodes from m to n-1
     * @param nodeIdx
     * @return a list of node indices of adjacent nodes
     */
    std::vector<std::size_t> getNbrs(const std::size_t& nodeIdx) {
        if (auto it = nbrCache.find(nodeIdx); it != nbrCache.end()) {
            return it->second;
        } else {
            if (pcm->empty() || pcm->front().empty()) {
                std::cerr << "error getting nbrs for node " << nodeIdx << std::endl;
                throw QeccException("Cannot return neighbours, pcm empty");
            }
            const auto               nrChecks = pcm->size();
            const auto               nrBits   = pcm->front().size();
            std::vector<std::size_t> res;
            if (nodeIdx < nrBits) {
                for (std::size_t i = 0; i < nrChecks; i++) {
                    if (pcm->at(i).at(nodeIdx)) {
                        res.emplace_back(nrBits + i);
                    }
                }
            } else {
                for (std::size_t i = 0; i < nrBits; i++) {
                    if (pcm->at(nodeIdx - nrBits).at(i)) {
                        res.emplace_back(i);
                    }
                }
            }
            const auto [nbrIt, inserted] = nbrCache.try_emplace(nodeIdx, res);
            return nbrIt->second;
        }
    }
    [[nodiscard]] json to_json() const {
        return json{
                {"pcm", *this->pcm}};
    }

    void from_json(const json& j) {

    }


    [[nodiscard]] std::string toString() const {
        return this->to_json().dump(2U);
    }
};

/**
 * Considers X errors with Z checks only.
 * Z errors are corrected completely analogously in a symmetric way for Z and X.
 */
class Code {
private:
    std::unique_ptr<ParityCheckMatrix> Hz;
    std::unique_ptr<ParityCheckMatrix> Hx;
public:
    std::size_t                        K = 0U;
    std::size_t                        N = 0U;
    std::size_t                        D = 0U;

    Code() = default;

    [[nodiscard]] const std::unique_ptr<ParityCheckMatrix>& getHz() const {
        return Hz;
    }

    gf2Mat getHzMat() {
        return *this->Hz->pcm;
    }

    gf2Mat getHxMat() {
        return *this->Hx->pcm;
    }

    void setHz(std::vector<std::vector<bool>>& hz) {
        Hz = std::make_unique<ParityCheckMatrix>(hz);
    }

    void setHx(std::vector<std::vector<bool>>& hx) {
        Hx = std::make_unique<ParityCheckMatrix>(hx);
    }
    /*
     * Takes matrix Hz over GF(2) and constructs respective code for X errors with Z checks represented by Hz
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hz):
        Hz(std::make_unique<ParityCheckMatrix>(hz)) {
        N = Hz->pcm->front().size();
    }

    /*
     * Takes two pcms over GF(2) and constructs respective code
     * Convention: Rows in first dim, columns in second
     */
    explicit Code(std::vector<std::vector<bool>>& hz, std::vector<std::vector<bool>>& hx):
        Hz(std::make_unique<ParityCheckMatrix>(hz)), Hx(std::make_unique<ParityCheckMatrix>(hx)) {
        N = Hz->pcm->front().size();
    }


    /**
     * Constructs the X check part of a code given
     * @param pathToPcm
     */
    explicit Code(const std::string& pathToPcm):
        Hz(std::make_unique<ParityCheckMatrix>(pathToPcm)) {
        if (Hz->pcm->empty() || Hz->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, Hz empty");
        }
        N = Hz->pcm->front().size();
    }

    explicit Code(const std::string& pathToHx, const std::string& pathToHz):
        Hz(std::make_unique<ParityCheckMatrix>(pathToHz)), Hx(std::make_unique<ParityCheckMatrix>(pathToHx)) {
        if (Hz->pcm->empty() || Hz->pcm->front().empty() || Hx->pcm->empty() || Hx->pcm->front().empty()) {
            throw QeccException("[Code::ctor] - Cannot construct Code, Hx or Hz empty");
        }
        N = Hz->pcm->front().size();
    }

    [[nodiscard]] std::size_t getN() const {
        return N;
    }

    [[nodiscard]] std::size_t getK() const {
        return K;
    }

    /**
     * Returns the syndrome, given a X-error represented as binary vector
     * @param err
     * @return
     */
    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& err) const {
        if (err.empty()) {
            throw QeccException("Cannot compute syndrome, err empy");
        }
        gf2Vec syndr((*Hz->pcm).size(), false);
        Utils::rectMatrixMultiply(*Hz->pcm, err, syndr);
        return syndr;
    }

    /**
     * Returns the syndrome, given an error represented as two binary vector, one holding the X and one the Z part
     * @param err
     * @return
     */
    [[nodiscard]] gf2Vec getSyndrome(const gf2Vec& Xerr, const gf2Vec& Zerr) const {
        if (Xerr.empty() || Zerr.empty()) {
            throw QeccException("Cannot compute syndrome, err empy");
        }

        gf2Vec Xsyndr((*Hz->pcm).size(), false);
        Utils::rectMatrixMultiply(*Hz->pcm, Xerr, Xsyndr);

        gf2Vec Zsyndr((*Hx->pcm).size(), false);
        Utils::rectMatrixMultiply(*Hx->pcm, Zerr, Zsyndr);
        gf2Vec res(Xsyndr.size()+Zsyndr.size());
        std::move(Xsyndr.begin(), Xsyndr.end(), std::back_inserter(res));
        std::move(Zsyndr.begin(), Zsyndr.end(), std::back_inserter(res));
        return res;
    }

    /**
     * Checks if the given vector is a X stabilizer of the code
     * @param est
     * @return
     */
    [[nodiscard]] bool isVectorStabilizer(const gf2Vec& est) const {
        return Utils::isVectorInRowspace(*Hz->pcm, est);
    }

    /**
     * Determines if the given vector represented as two components, X and Z
     * Is a stabilizer
     * @param Xest
     * @param Zest
     * @return
     */
    bool isStabilizer(const gf2Vec& Xest, const gf2Vec& Zest) const {
        return Utils::isVectorInRowspace(*Hz->pcm, Xest) && Utils::isVectorInRowspace(*Hx->pcm, Zest);
    }

    [[nodiscard]] CodeProperties getProperties() const {
        CodeProperties res{.n = N, .k = getK(), .d = D};
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, const Code& c) {
        auto   nrChecks = c.Hz->pcm->size();
        auto   nrData   = c.Hz->pcm->front().size();
        auto   dim      = nrChecks + nrData;
        gf2Mat res(dim);
        os << "Hz: " << std::endl;
        for (size_t i = 0; i < dim; i++) {
            gf2Vec row(dim);
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    row.at(nrData + j) = c.Hz->pcm->at(j).at(i);
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    row.at(j) = c.Hz->pcm->at(i - nrData).at(j);
                }
            }
            res.at(i) = row;
        }
        os << Utils::getStringFrom(res) << "Hx: " << std::endl;
        for (size_t i = 0; i < dim; i++) {
            gf2Vec row(dim);
            if (i < dim - nrChecks) {
                for (size_t j = 0; j < nrChecks; j++) {
                    row.at(nrData + j) = c.Hx->pcm->at(j).at(i);
                }
            } else {
                for (size_t j = 0; j < nrData; j++) {
                    row.at(j) = c.Hx->pcm->at(i - nrData).at(j);
                }
            }
            res.at(i) = row;
        }
        return os;
    }
    [[nodiscard]] json to_json() const {
        return json{
                {"Hz", Hz->to_json()},
                {"Hx", Hx->to_json()}
                {"N", N},
                {"K", K},
                {"D", D}};
    }

    void from_json(const json& j) {
        j.at("N").get_to(N);
        j.at("K").get_to(K);
        j.at("D").get_to(D);
    }

    [[nodiscard]] std::string toString() const {
            return this->to_json().dump(2U);
    }
};
#endif //QUNIONFIND_CODE_HPP
