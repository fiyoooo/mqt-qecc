//
// Created by lucas on 09/06/22.
//
#pragma once

#include "Decoder.hpp"

using NodeSet    = std::unordered_set<std::size_t>; // TODO can also change name
using NodeVector = std::vector<std::size_t>;        // TODO can also change name

class UFDecoder : public Decoder {
public:
    using Decoder::Decoder;

    void decode(const gf2Vec& syndrome) override;

    void reset() override;

private:
    void doDecode(const gf2Vec& syndrome, const std::unique_ptr<ParityCheckMatrix>& pcm);

    bool containsInvalidComponents(const NodeSet&                            nodeSet,
                                   const NodeSet&                            syndrome,
                                   std::vector<NodeSet>&                     invalidComps,
                                   const std::unique_ptr<ParityCheckMatrix>& pcm) const;

    [[nodiscard]] std::vector<NodeSet> getConnectedComps(const NodeSet& nodeSet) const;

    [[nodiscard]] bool isValidComponent(const NodeSet& nodeSet, const NodeSet& syndrome,
                                        const std::unique_ptr<ParityCheckMatrix>& pcm) const;

    [[nodiscard]] NodeVector computeInteriorBitNodes(const NodeSet& nodeSet) const;

    [[nodiscard]] NodeSet getEstimateForComponent(const NodeSet&                            nodeSet,
                                                  const NodeSet&                            syndrome,
                                                  const std::unique_ptr<ParityCheckMatrix>& pcm) const;

    void standardGrowth(NodeSet& nodeSet);

    void singleClusterSmallestFirstGrowth(NodeSet& nodeSet);

    void singleClusterRandomFirstGrowth(NodeSet& nodeSet);

    void singleQubitRandomFirstGrowth(NodeSet& nodeSet);
};
