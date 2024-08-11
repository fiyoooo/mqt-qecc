//
// Created by Fiona Froehler on 08.04.24.
//
#pragma once

#include <chrono>
#include <unordered_map>
#include <vector>
#include <z3++.h>

class LightsOut {
public:
    LightsOut(const std::unordered_map<int, std::vector<int>>& lightsToSwitches,
              const std::unordered_map<int, std::vector<int>>& switchesToLights);

    /**
     * @brief Preconstruct the z3 instance for the lights-out problem.
     *
     * Creates all necessary variables, adds the known parts of the parity constraints.
     * Soft constraints are added to the optimizer with default weights.
     */
    void preconstructZ3Instance();

    /**
     * @brief Solve the lights-out problem for a given pattern.
     *
     * Assumes that the Z3 instance has already been pre-constructed.
     *
     * @param syndrome The pattern of syndrome to solve the problem for.
     * @param solver_path The path to the solver executable (default: "z3").
     * @return A tuple containing the solution switches, time taken for constraint completion, and time taken for solving the problem.
     */
    std::tuple<std::vector<int>, std::chrono::duration<double>, std::chrono::duration<double>> solve(
            const std::vector<bool>& syndrome);

private:
    z3::context     ctx_;       // manages memory allocation and other resources of z3 instance
    z3::optimize    optimizer_; // optimizes constraints and solves problem
    z3::expr_vector switch_vars_;

    std::unordered_map<int, std::vector<int>> lights_to_switches_;
    std::unordered_map<int, std::vector<int>> switches_to_lights_;
    std::unordered_map<int, z3::expr_vector>  helper_vars_;

    /**
     * @brief Preconstructs the parity constraints for a light.
     *
     * Adds all constraints to the optimizer that are independent of the value of the light.
     *
     * @param light in question
     * @param switches that toggle the light
     */
    void preconstructParityConstraint(int light, const std::vector<int>& switches);

    /**
     * @brief Completes the parity constraints for a light.
     *
     * Adds the constraint that is dependent on the value of the light.
     *
     * @param light in question
     * @param switches that toggle the light
     * @param val of the light in the syndrome
     */
    void completeParityConstraint(int light, const std::vector<int>& switches, bool val);

    /**
     * @brief Validates the model by checking if pressing the switches turns off all lights.
     *
     * @param model to be validated
     * @param lights booleans representing current state of the lights
     * @return true if all lights are switched off after simulating switch presses, false otherwise
     */
    bool validateModel(const z3::model& model, const std::vector<bool>& lights);

    /**
     * @brief Counts the number of switches that are set to true.
     *
     * @param model to be analyzed
     * @return number of switches that are set to true
     */
    [[maybe_unused]] int countSwitches(const z3::model& model);
};