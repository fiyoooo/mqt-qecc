//
// Created by lucas on 21/04/2022.
//

#include "UFDecoder.hpp"

#include "Decoder.hpp"

#include <chrono>
#include <gf2dense.hpp>
#include <queue>
#include <random>
#include <set>

/**
 * @brief Decodes the syndrome using the union-find decoder algorithm.
 *
 * Original implementation of the generalized decoder for QLDPC codes using Gaussian elimination.
 * @param syndrome
 */
void UFDecoder::decode(const gf2Vec& syndrome) {
    // if syndrome dimensions larger than parity-check-matrix
    if (syndrome.size() > this->getCode()->gethZ()->pcm->size()) {
        // split into two parts
        auto         mid = syndrome.begin() + (static_cast<std::int64_t>(syndrome.size()) / 2U);
        gf2Vec const xSyndr(syndrome.begin(), mid);
        gf2Vec const zSyndr(mid, syndrome.end());

        // decode
        doDecode(xSyndr, this->getCode()->gethZ());
        auto xres = this->result;
        this->reset();
        doDecode(zSyndr, this->getCode()->gethX());

        // combine
        this->result.decodingTime += xres.decodingTime;
        std::move(xres.estimBoolVector.begin(), xres.estimBoolVector.end(),
                  std::back_inserter(this->result.estimBoolVector));
        std::move(xres.estimNodeIdxVector.begin(), xres.estimNodeIdxVector.end(),
                  std::back_inserter(this->result.estimNodeIdxVector));
    } else {
        doDecode(syndrome, this->getCode()->gethZ());
    }
}

void UFDecoder::doDecode(const gf2Vec& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    const auto decodingTimeBegin = std::chrono::high_resolution_clock::now();
    NodeSet    components; // used to store vertex indices in E set
    NodeSet    syndr;      // vertex indices of syndrome nodes
    for (std::size_t i = 0; i < syndrome.size(); i++) {
        if (syndrome.at(i)) {
            syndr.insert(this->getCode()->getN() + i);
        }
    }

    if (!syndr.empty()) {
        // Set set of nodes equal to syndrome E = syndrome
        components.insert(syndr.begin(), syndr.end());
        std::vector<NodeSet> invalidComponents;

        while (containsInvalidComponents(components, syndr, invalidComponents, pcm) &&
               components.size() < (pcm->pcm->size() + pcm->pcm->front().size())) {
            switch (this->growth) {
            case GrowthVariant::AllComponents:
                // grow all components (including valid ones) by 1
                standardGrowth(components);
                break;
            case GrowthVariant::InvalidComponents:
                // not implemented yet
                break;
            case GrowthVariant::SingleSmallest:
                // grow only by neighbours of single smallest cluster
                singleClusterSmallestFirstGrowth(components);
                break;
            case GrowthVariant::SingleRandom:
                // grow only by neighbours of single random cluster
                singleClusterRandomFirstGrowth(components);
                break;
            case GrowthVariant::SingleQubitRandom:
                // grow only by neighbours of single qubit
                singleQubitRandomFirstGrowth(components);
                break;
            default:
                throw std::invalid_argument("Unsupported growth variant");
            }
        }
    }

    // collect all estimates in a set
    auto    connectedComponents = getConnectedComps(components);
    NodeSet estimates;
    for (const auto& comp : connectedComponents) {
        auto compEstimate = getEstimateForComponent(comp, syndr, pcm);
        estimates.insert(compEstimate.begin(), compEstimate.end());
    }

    const auto decodingTimeEnd = std::chrono::high_resolution_clock::now();
    result.decodingTime        = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                                                           decodingTimeEnd - decodingTimeBegin)
                                                           .count());

    // store result
    result.estimBoolVector = gf2Vec(getCode()->getN());
    for (auto re : estimates) {
        result.estimBoolVector.at(re) = true;
    }
    result.estimNodeIdxVector.assign(estimates.begin(), estimates.end());
}

/**
 * @brief Checks if there is a component in the list that is not valid.
 * @param nodeSet
 * @param syndrome
 * @param invalidComps
 * @param pcm
 * @return True if there are invalid components, false otherwise.
 */
bool UFDecoder::containsInvalidComponents(const NodeSet&                            nodeSet,
                                          const NodeSet&                            syndrome,
                                          std::vector<NodeSet>&                     invalidComps,
                                          const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    std::vector<NodeSet> ccomps = getConnectedComps(nodeSet);
    return std::any_of(ccomps.begin(), ccomps.end(), [&](const auto& comp) {
        bool const isValid = isValidComponent(comp, syndrome, pcm);
        if (!isValid) {
            invalidComps.emplace_back(comp.begin(), comp.end());
        }
        return !isValid;
    });
}

/**
 * @brief Computes the connected components of a given set of nodes.
 *
 * The set of all nodes considered by the algorithm in the Tanner graph, // TODO really always all nodes?
 * compute the connected components in the Tanner graph.
 * @param nodeSet
 * @return A vector of connected components.
 */
std::vector<NodeSet> UFDecoder::getConnectedComps(const NodeSet& nodeSet) const {
    NodeSet              visited;
    std::vector<NodeSet> res;

    for (auto node : nodeSet) {
        // if node not yet visited
        if (visited.find(node) == visited.end()) {
            visited.insert(node);
            NodeSet ccomp;

            // start from 'node' to compute connected component with DFS-like algorithm
            std::queue<std::size_t> stack; // can be changed: DFS uses stack, BFS uses queue
            stack.push(node);
            while (!stack.empty()) {
                // get current element
                auto curr = stack.back();
                stack.pop();

                // if not yet in component
                if (ccomp.find(curr) == ccomp.end()) {
                    visited.insert(curr);
                    ccomp.insert(curr);

                    // check all neighbours of current element
                    auto nbrs = getCode()->gethZ()->getNbrs(curr);
                    for (auto n : nbrs) {
                        // if 'n' is not in component yet but in nodeSet
                        if (ccomp.find(n) == ccomp.end() && nodeSet.find(n) != nodeSet.end()) {
                            stack.push(n);
                        }
                    }
                }
            }
            res.emplace_back(ccomp);
        }
    }
    return res;
}

/**
 * @brief Checks if a component is valid.
 *
 * A component is valid if there is a set of (bit) nodes in its interior whose syndrome is equal to the given syndrome.
 * @param nodeSet
 * @param syndrome
 * @param pcm
 * @return True if the component is valid, false otherwise.
 */
bool UFDecoder::isValidComponent(const NodeSet&                            nodeSet,
                                 const NodeSet&                            syndrome,
                                 const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    return !getEstimateForComponent(nodeSet, syndrome, pcm).empty();
}

/**
 * @brief Computes a set of nodes s.t. for each n in the list, all neighbours of n are in the component.
 * @param nodeSet of the component
 * @return A vector of interior bit node indices.
 */
NodeVector UFDecoder::computeInteriorBitNodes(const NodeSet& nodeSet) const {
    NodeVector res;
    for (auto node : nodeSet) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(node);
        // if node is a bit node and all its neighbours are in nodeSet
        if (node < getCode()->getN() && std::includes(nodeSet.begin(), nodeSet.end(), nbrs.begin(), nbrs.end())) {
            res.emplace_back(node);
        }
    }
    return res;
}

/**
 * @brief Computes estimate vector x for a component and a syndrome.
 *
 * This is done by considering all vertices in Tanner Graph that are in the Interior
 * of the given node set and additionally the neighbours of the bit vertices in the interior.
 * Then, using Gaussian elimination, it is checked whether a solution for the local cluster
 * that is consistent with the syndrome can be found. If so, this local estimate is returned.
 * @param nodeSet
 * @param syndrome
 * @param pcm
 * @return A set of estimated node indices.
 */
NodeSet UFDecoder::getEstimateForComponent(const NodeSet&                            nodeSet,
                                           const NodeSet&                            syndrome,
                                           const std::unique_ptr<ParityCheckMatrix>& pcm) const {
    NodeSet                   res{};
    std::vector<size_t> const intNodes = computeInteriorBitNodes(nodeSet);
    if (intNodes.empty()) {
        return NodeSet{};
    }

    gf2Mat            redHz;
    gf2Vec            redSyndr;
    std::vector<bool> used(pcm->pcm->size());

    // fill up reduced pcm and syndrome
    for (auto node : nodeSet) {
        if (node >= getCode()->getN()) { // is a check node
            auto idx = node - getCode()->getN();
            if (!used.at(idx)) {
                used.at(idx) = true;
                redHz.emplace_back(pcm->pcm->at(idx));
                // if the check node is in the syndrome we need to satisfy check=1
                redSyndr.emplace_back(syndrome.find(idx) != syndrome.end());
            }
        } else { // is a bit node
            const auto nbrs = this->getCode()->gethZ()->getNbrs(node);
            // add neighbouring checks (maybe not in the interior but include to stay consistent with the syndrome)
            for (auto n : nbrs) {
                auto idx = n - getCode()->getN();
                if (!used.at(idx)) {
                    used.at(idx) = true;
                    redHz.emplace_back(pcm->pcm->at(idx));
                    redSyndr.emplace_back(syndrome.find(n) != syndrome.end()); // TODO why n and not idx
                }
            }
        }
    }
    auto                 redHzCsc = Utils::toCsc(redHz);
    std::vector<uint8_t> redSyndInt(redSyndr.size());
    for (std::size_t i = 0; i < redSyndr.size(); i++) {
        redSyndInt.at(i) = redSyndr.at(i) ? 1 : 0;
    }
    auto pluDec = ldpc::gf2dense::PluDecomposition(static_cast<int>(redHz.size()), static_cast<int>(redHz.at(0).size()),
                                                   redHzCsc);
    pluDec.rref();

    auto estim = pluDec.lu_solve(
            redSyndInt); // solves the system redHz*x=redSyndr by x to see if a solution can be found
    for (std::size_t i = 0; i < estim.size(); i++) {
        if (estim.at(i) != 0U) {
            auto inst = res.insert(static_cast<size_t>(i));
            std::cout << inst.second;
        }
    }
    return res;
}

/**
 * @brief Grows the node set by the neighbours of ALL clusters.
 * @param nodeSet
 */
void UFDecoder::standardGrowth(NodeSet& nodeSet) {
    NodeSet newNodes;
    for (const auto& node : nodeSet) {
        const auto nbrs = getCode()->gethZ()->getNbrs(node);
        newNodes.insert(nbrs.begin(), nbrs.end());
    }
    nodeSet.insert(newNodes.begin(), newNodes.end());
}

/**
 * @brief Grows the node set by the neighbours of the single smallest cluster.
 * @param nodeSet
 */
void UFDecoder::singleClusterSmallestFirstGrowth(NodeSet& nodeSet) {
    std::vector<NodeSet> ccomps = getConnectedComps(nodeSet);
    if (ccomps.empty()) { // to avoid undefined behavior
        return;
    }

    NodeSet const smallestComponent = *std::min_element(ccomps.begin(), ccomps.end(),
                                                        [](const auto& lhs, const auto& rhs) {
                                                            return lhs.size() < rhs.size();
                                                        });
    for (const auto& node : smallestComponent) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * @brief Grows the node set by the neighbours of a single random cluster.
 * @param nodeSet
 */
void UFDecoder::singleClusterRandomFirstGrowth(NodeSet& nodeSet) {
    std::vector<NodeSet> ccomps = getConnectedComps(nodeSet);
    if (ccomps.empty()) { // to avoid undefined behavior
        return;
    }

    // choose index at random
    std::random_device                         rd;
    std::mt19937                               gen(rd());
    std::uniform_int_distribution<std::size_t> dis(0, ccomps.size() - 1);
    auto                                       chosenIdx = dis(gen);

    // get component at chosen index
    auto           it              = std::next(ccomps.begin(), static_cast<int>(chosenIdx));
    const NodeSet& chosenComponent = *it;

    for (const auto& node : chosenComponent) {
        const auto& nbrs = getCode()->gethZ()->getNbrs(node);
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * @brief Grows the node set by the neighbours of a single random qubit.
 * @param nodeSet
 */
void UFDecoder::singleQubitRandomFirstGrowth(NodeSet& nodeSet) {
    std::vector<NodeSet> ccomps = getConnectedComps(nodeSet);
    if (ccomps.empty()) { // to avoid undefined behavior
        return;
    }

    // choose index at random
    std::random_device                         rd;
    std::mt19937                               gen(rd());
    std::uniform_int_distribution<std::size_t> dis(0, ccomps.size() - 1);
    auto                                       chosenIdx = dis(gen);

    // get component at chosen index
    auto           it              = std::next(ccomps.begin(), static_cast<int>(chosenIdx));
    const NodeSet& chosenComponent = *it;

    if (!chosenComponent.empty()) { // to avoid undefined behavior
        const auto& nbrs = getCode()->gethZ()->getNbrs(*chosenComponent.begin());
        nodeSet.insert(nbrs.begin(), nbrs.end());
    }
}

/**
 * @brief Resets temporarily computed data.
 */
void UFDecoder::reset() {
    this->result = {};
    this->growth = GrowthVariant::AllComponents;
}
