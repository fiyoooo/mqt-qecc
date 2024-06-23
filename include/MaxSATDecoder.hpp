//
// Created by Fiona Froehler on 03.06.24.
//
#pragma once

#include <unordered_map>
#include <z3++.h>

#include "Decoder.hpp"

// using idk = std::vector<std::size_t> // TODO

class MaxSATDecoder : public Decoder {
public:
    // using Decoder::Decoder;

    MaxSATDecoder(Code& c);

    void decode(const gf2Vec &syndrome) override;

    void reset() override;

    void preconstructZ3Instance();

private:
    z3::context ctx_; // manages memory allocation and other resources of z3 instance
    z3::optimize optimizer_; // optimizes constraints and solves problem
    z3::expr_vector switch_vars_;

    // TODO vectors or sets?
    // std::unordered_map<std::size_t, std::vector<std::size_t>> lights_to_switches_;
    // std::unordered_map<std::size_t, std::vector<std::size_t>> switches_to_lights_;
    std::unordered_map<std::size_t, z3::expr_vector> helper_vars_;

    static std::vector<std::size_t> getSwitches(std::size_t light, const std::unique_ptr<ParityCheckMatrix> &pcm);

    void preconstructParityConstraint(std::size_t light, const std::vector<std::size_t> &switches);

    void completeParityConstraint(std::size_t light, const std::vector<std::size_t> &switches, bool val);

    // TODO use gf2Vec
    bool validateModel(const z3::model &model, const std::vector<bool> &lights);

    [[maybe_unused]] std::size_t countSwitches(const z3::model &model);
};