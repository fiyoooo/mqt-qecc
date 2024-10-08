//
// Created by Fiona Froehler on 03.06.24.
//

#include "MaxSATDecoder.hpp"

#include <algorithm>
#include <iostream>
#include <random>

MaxSATDecoder::MaxSATDecoder(Code& c)
    : Decoder(c), optimizer_(ctx_), switch_vars_(ctx_) {}

/**
 * @brief Preconstruct the z3 instance for the lights-out problem.
 *
 * Creates all necessary variables, adds the known parts of the parity constraints.
 * Soft constraints are added to the optimizer with default weights.
 */
void MaxSATDecoder::preconstructZ3Instance() {
    const std::unique_ptr<ParityCheckMatrix>& pcm = this->getCode()->gethZ();
    std::size_t const                         k   = pcm->pcm->size();        // number of lights
    std::size_t const                         n   = this->getCode()->getN(); // number of switches

    // declare switches as variables
    if (switch_vars_.empty()) {
        for (std::size_t i = 0; i < n; ++i) {
            switch_vars_.push_back(ctx_.bool_const(("switch_" + std::to_string(i)).c_str()));
        }
    }

    // create helper variables and add known parity constraints
    for (std::size_t light = 0; light < k; light++) {
        const std::vector<unsigned>& switches = getSwitches(light, pcm); // TODO change to gf2Vec

        // check if helper variables already exist for light
        if (helper_vars_.find(light) == helper_vars_.end()) {
            z3::expr_vector exprVec(ctx_);
            for (std::size_t i = 0; i < switches.size() - 1; ++i) {
                exprVec.push_back(
                        ctx_.bool_const(("helper_" + std::to_string(light) + "_" + std::to_string(i)).c_str()));
            }
            helper_vars_.emplace(light, exprVec);
        }

        // if helper_vars_ empty at light then no constraints needed TODO correct?
        if (!helper_vars_.at(light).empty()) {
            preconstructParityConstraint(light, switches);
        }
    }

    // add soft constraints (switches set to False)
    for (const auto& switchVar : switch_vars_) {
        optimizer_.add_soft(!switchVar, 1); // default weight 1 explicitly needed in cpp
    }
}

/**
 * @brief Solve the lights-out problem for a given pattern.
 *
 * Assumes that the Z3 instance has already been pre-constructed.
 *
 * @param syndrome to solve the problem for
 */
void MaxSATDecoder::decode(const gf2Vec& syndrome) {
    // push a new context to the optimizer
    optimizer_.push();

    const auto constructionTimeBegin = std::chrono::high_resolution_clock::now();

    const std::unique_ptr<ParityCheckMatrix>& pcm = this->getCode()->gethZ();
    std::size_t const                         k   = pcm->pcm->size(); // number of lights

    // add the problem specific constraints
    for (std::size_t light = 0; light < k; light++) {
        const std::vector<unsigned>& switches = getSwitches(light, pcm); // TODO change to gf2Vec
        completeParityConstraint(light, switches, syndrome[light]);
    }
    const auto                  constructionTimeEnd = std::chrono::high_resolution_clock::now();
    [[maybe_unused]] const auto constructionTime    = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                                                                                    constructionTimeEnd - constructionTimeBegin)
                                                                                    .count());

    // solve the problem
    const auto decodingTimeBegin = std::chrono::high_resolution_clock::now();
    auto       result            = optimizer_.check();
    const auto decodingTimeEnd   = std::chrono::high_resolution_clock::now();
    this->result.decodingTime    = static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                                                                 decodingTimeEnd - decodingTimeBegin)
                                                                 .count());

    // set the switches
    auto                model = optimizer_.get_model();
    gf2Vec              switchesBool;
    std::vector<size_t> switchesIdx;
    for (unsigned i = 0; i < switch_vars_.size(); i++) {
        z3::expr const evalResult = model.eval(switch_vars_[i], true);
        if (evalResult.is_true()) {
            switchesBool.push_back(true);
            switchesIdx.push_back(i);
        } else if (evalResult.is_false()) {
            switchesBool.push_back(false);
        } else {
            throw std::logic_error("Expression evaluation result is not boolean");
        }
    }

    // pop the context from the optimizer, restores previous state of optimizer
    optimizer_.pop();

    // store in result
    this->result.estimBoolVector    = std::move(switchesBool);
    this->result.estimNodeIdxVector = std::move(switchesIdx);

    // check found solution
    if (result != z3::sat) {
        throw std::logic_error("No solution found");
    }
    if (!validateModel(model, syndrome)) {
        throw std::logic_error("Model is invalid");
    }
}

/**
 * @brief Returns all switches toggling a given light.
 *
 * @param light in question
 * @param pcm parity check matrix of the code
 */
std::vector<unsigned> MaxSATDecoder::getSwitches(std::size_t light, const std::unique_ptr<ParityCheckMatrix>& pcm) {
    std::vector<unsigned> switches;
    const auto&           vec = pcm->pcm->at(light);

    for (unsigned i = 0; i < vec.size(); ++i) {
        if (vec[i]) {
            switches.push_back(i);
        }
    }
    return switches;
}

/**
 * @brief Preconstructs the parity constraints for a light.
 *
 * Adds all constraints to the optimizer that are independent of the value of the light.
 *
 * @param light in question
 * @param switches that toggle the light
 */
void MaxSATDecoder::preconstructParityConstraint(std::size_t light, const std::vector<unsigned>& switches) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    // get helper variables for the given light
    auto& helperVars = helper_vars_.at(light);

    // split into smaller constraints
    for (unsigned i = 1; i < switches.size() - 1; ++i) {
        // switch_i XOR h_i == h_{i-1}
        z3::expr const constraint = (switch_vars_[switches[i]] ^ helperVars[i]) == helperVars[i - 1];
        optimizer_.add(constraint.simplify());
    }

    // add last constraint
    z3::expr const constraint = switch_vars_[switches.back()] == helperVars.back();
    optimizer_.add(constraint.simplify());
}

/**
 * @brief Completes the parity constraints for a light.
 *
 * Adds the constraint that is dependent on the value of the light.
 *
 * @param light in question
 * @param switches that toggle the light
 * @param val of the light in the syndrome
 */
void MaxSATDecoder::completeParityConstraint(std::size_t light, const std::vector<unsigned>& switches, bool val) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    // get helper variables for the given light
    auto& helperVars = helper_vars_.at(light);

    // if only one switch toggles light TODO correct?
    if (helperVars.empty()) {
        return;
    }

    // switch_0 XOR h_0 == val
    z3::expr const constraint = (switch_vars_[switches[0]] ^ helperVars[0]) == ctx_.bool_val(val);
    optimizer_.add(constraint.simplify());
}

/**
 * @brief Validates the model by checking if pressing the switches turns off all lights.
 *
 * @param model to be validated
 * @param lights booleans representing current state of the lights
 * @return true if all lights are switched off after simulating switch presses, false otherwise
 */
bool MaxSATDecoder::validateModel(const z3::model& model, const gf2Vec& syndrome) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    gf2Vec copy(syndrome);
    for (unsigned i = 0; i < switch_vars_.size(); ++i) {
        // if switch set to true
        if (model.eval(switch_vars_[i], true).is_true()) {
            // flip all lights...
            for (std::size_t light = 0; light < syndrome.size(); light++) {
                // ...that are controlled by this switch
                bool const isFlipped = this->getCode()->gethZ()->pcm->at(light).at(i); // TODO check if correct column or row
                copy[light]          = copy[light] ^ isFlipped;
            }
        }
    }

    // check if all lights are switched off now
    return std::all_of(copy.begin(), copy.end(), [](bool light) { return !light; });
}

/**
 * @brief Counts the number of switches that are set to true.
 *
 * @param model to be analyzed
 * @return number of switches that are set to true
 */
[[maybe_unused]] std::size_t MaxSATDecoder::countSwitches(const z3::model& model) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    std::size_t count = 0;
    for (const auto& var : switch_vars_) {
        if (model.eval(var) == z3::sat) {
            count++;
        }
    }
    return count;

    // this doesn't work because of iterator incompatibility
    // return std::count_if(switch_vars_.begin(), switch_vars_.end(), [&model](const z3::expr &var) { return model.eval(var) == z3::sat; });
}

/**
 * @brief Resets temporarily computed data.
 */
void MaxSATDecoder::reset() {
    this->result = {};
    this->growth = GrowthVariant::AllComponents;
    // TODO reset my attributes?
}
