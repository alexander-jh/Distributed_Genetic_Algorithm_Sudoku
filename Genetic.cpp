// Genetic.cpp
//
// Main method for Genetic algorithm.
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "OptimizationProblem.hpp"

int main(int argv, char **argc) {
    // Create new Genetic object
    auto dna = new Genetic(*(argc+1));
    // Run algorithm
    bool status = dna->genetic_run();
    // Report end state
    dna->genetic_report(status);
    return 0;
}
