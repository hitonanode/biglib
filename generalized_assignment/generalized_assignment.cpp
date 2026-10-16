#include "generalized_assignment.hpp"
#include <iostream>

int main() {
    // job-major: cost[job][agent], resource[job][agent]
    cpopt::gap::Solver<> solver({{1, 7}, {4, 2}, {3, 6}}, {{2, 3}, {3, 2}, {2, 2}}, {4, 4});
    solver.set_start({0, -1, -1}); // -1: unassigned; assigned jobs are not fixed.
    cpopt::gap::Options options;
    options.time_limit = 0.05;
    const auto &first = solver.solve(options);
    if (first.feasible) std::cout << "cost: " << first.objective << '\n';

    // Continue with RNG, tabu memory, penalties, incumbent, and current assignment retained.
    const auto &continued = solver.solve(options);
    if (continued.feasible) {
        for (int agent : continued.assignment) std::cout << agent << ' ';
        std::cout << '\n';
    }
    cpopt::gap::Solver<double> decimals(
        {{0.1, 0.7}, {0.4, 0.2}}, {{0.2, 0.3}, {0.3, 0.2}}, {0.4, 0.4});
    const auto &decimal_result = decimals.solve(options);
    if (decimal_result.feasible) std::cout << "decimal cost: " << decimal_result.objective << '\n';
}
