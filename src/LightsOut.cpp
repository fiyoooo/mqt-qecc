//
// Created by Fiona Froehler on 08.04.24.
//

#include "LightsOut.hpp"

#include <algorithm>
#include <iostream>

LightsOut::LightsOut(const std::unordered_map<int, std::vector<int>>& lightsToSwitches,
                     const std::unordered_map<int, std::vector<int>>& switchesToLights)
    : ctx_(), optimizer_(ctx_), switch_vars_(ctx_),
      lights_to_switches_(lightsToSwitches),
      switches_to_lights_(switchesToLights) {}

void LightsOut::preconstructZ3Instance() {
    // declare switches as variables
    if (switch_vars_.empty()) {
        for (size_t i = 0; i < switches_to_lights_.size(); ++i) {
            switch_vars_.push_back(ctx_.bool_const(("switch_" + std::to_string(i)).c_str()));
        }
    }

    // create helper variables and add known parity constraints
    for (const auto& entry : lights_to_switches_) {
        int const               light    = entry.first;
        const std::vector<int>& switches = entry.second;

        // check if helper variables already exist for light
        if (helper_vars_.find(light) == helper_vars_.end()) {
            z3::expr_vector exprVec(ctx_);
            for (size_t i = 0; i < switches.size() - 1; ++i) {
                exprVec.push_back(
                        ctx_.bool_const(("helper_" + std::to_string(light) + "_" + std::to_string(i)).c_str()));
            }
            helper_vars_.emplace(light, exprVec);
        }
        preconstructParityConstraint(light, switches);
    }

    // add soft constraints (switches set to False)
    for (const auto& switchVar : switch_vars_) {
        optimizer_.add_soft(!switchVar, 1); // default weight 1 explicitly needed in cpp
    }
}

std::tuple<std::vector<int>, std::chrono::duration<double>, std::chrono::duration<double>>
LightsOut::solve(const std::vector<bool>& syndrome) {
    // push a new context to the optimizer
    optimizer_.push();

    auto startTime = std::chrono::high_resolution_clock::now();

    // add the problem specific constraints
    for (const auto& entry : lights_to_switches_) {
        completeParityConstraint(entry.first, entry.second, syndrome[entry.first]);
    }
    auto constrTime = std::chrono::high_resolution_clock::now() - startTime;

    std::vector<int>              switches;
    std::chrono::duration<double> solveTime{};

    // solve the problem
    startTime   = std::chrono::high_resolution_clock::now();
    auto result = optimizer_.check();
    solveTime   = std::chrono::high_resolution_clock::now() - startTime;
    if (result != z3::sat) {
        throw std::logic_error("No solution found");
    }

    // validate the model
    auto model = optimizer_.get_model();
    if (!validateModel(model, syndrome)) {
        throw std::logic_error("Model is invalid");
    }

    // set the switches
    for (z3::expr const switchVar : switch_vars_) {
        z3::expr const evalResult = model.eval(switchVar, true);
        if (evalResult.is_true()) {
            switches.push_back(1);
        } else if (evalResult.is_false()) {
            switches.push_back(0);
        } else {
            throw std::logic_error("Expression evaluation result is not boolean");
        }
    }

    // pop the context from the optimizer, restores previous state of optimizer
    optimizer_.pop();

    return std::make_tuple(switches, constrTime, solveTime);
}

void LightsOut::preconstructParityConstraint(int light, const std::vector<int>& switches) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    // get helper variables for the given light
    auto& helperVars = helper_vars_.at(light);

    // split into smaller constraints
    for (size_t i = 1; i < switches.size() - 1; ++i) {
        // switch_i XOR h_i == h_{i-1}
        z3::expr const constraint = (switch_vars_[switches[i]] ^ helperVars[i]) == helperVars[i - 1];
        optimizer_.add(constraint.simplify());
    }

    // add last constraint
    z3::expr const constraint = switch_vars_[switches.back()] == helperVars.back();
    optimizer_.add(constraint.simplify());
}

void LightsOut::completeParityConstraint(int light, const std::vector<int>& switches, bool val) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    // get helper variables for the given light
    auto& helperVars = helper_vars_.at(light);

    // switch_0 XOR h_0 == val
    z3::expr const constraint = (switch_vars_[switches[0]] ^ helperVars[0]) == ctx_.bool_val(val);
    optimizer_.add(constraint.simplify());
}

bool LightsOut::validateModel(const z3::model& model, const std::vector<bool>& lights) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    std::vector<bool> copy(lights);
    for (size_t i = 0; i < switch_vars_.size(); ++i) {
        // if switch set to true
        if (model.eval(switch_vars_[i], true).is_true()) {
            // flip all lights that are controlled by this switch
            for (int const light : switches_to_lights_[i]) {
                copy[light] = !copy[light];
            }
        }
    }

    // check if all lights are switched off now
    return std::all_of(copy.begin(), copy.end(), [](bool light) { return !light; });
}

[[maybe_unused]] int LightsOut::countSwitches(const z3::model& model) {
    // ensure switch variables are initialized
    if (switch_vars_.empty()) {
        throw std::logic_error("Switch variables not initialized");
    }

    int count = 0;
    for (const auto& var : switch_vars_) {
        if (model.eval(var) == z3::sat) {
            count++;
        }
    }
    return count;

    // this doesn't work because of iterator incompatibility
    // return std::count_if(switch_vars_.begin(), switch_vars_.end(), [&model](const z3::expr &var) { return model.eval(var) == z3::sat; });
}
