//
// Created by Fiona Froehler on 03.06.24.
//
#pragma once

#include "Decoder.hpp"

#include <unordered_map>
#include <z3++.h>

class MaxSATDecoder : public Decoder {
public:
    // using Decoder::Decoder;

    MaxSATDecoder(Code& c);

    void preconstructZ3Instance();

    void decode(const gf2Vec& syndrome) override;

    static std::vector<unsigned> getSwitches(std::size_t light, const std::unique_ptr<ParityCheckMatrix>& pcm);

    void reset() override;

private:
    z3::context                                      ctx_;       // manages memory allocation and other resources of z3 instance
    z3::optimize                                     optimizer_; // optimizes constraints and solves problem
    z3::expr_vector                                  switch_vars_;
    std::unordered_map<std::size_t, z3::expr_vector> helper_vars_;

    void preconstructParityConstraint(std::size_t light, const std::vector<unsigned>& switches);

    void completeParityConstraint(std::size_t light, const std::vector<unsigned>& switches, bool val);

    bool validateModel(const z3::model& model, const gf2Vec& syndrome);

    [[maybe_unused]] [[maybe_unused]] std::size_t countSwitches(const z3::model& model);
};