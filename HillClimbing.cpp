// HillClimbing.cpp
//
// Main calling program for HillClimbing algorithm.
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "OptimizationProblem.hpp"

int main(int argv, char **argc) {
    // Create new HillClimber object
    auto hill = new HillClimber(*(argc+1));
    // Run algorithm
    bool status = hill->hill_climb();
    // Report Results
    hill->hill_report(status);
    return 0;
}
