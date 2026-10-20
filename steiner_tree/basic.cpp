// SPDX-License-Identifier: CC0-1.0
#include "steiner_tree.hpp"
#include <iostream>

int main() {
    using Solver = steiner_tree::Solver<long long>;
    std::vector<Solver::Edge> edges{{0, 3, 2}, {1, 3, 2}, {2, 3, 2}, {0, 1, 5}, {1, 2, 5}};
    Solver solver(4, edges, {0, 1, 2}, {3}); // A partial edge set is repaired automatically.
    auto result = solver.solve(0.01);
    result = solver.solve(0.05); // Continue the same search; the best cost never increases.
    std::cout << "cost=" << static_cast<long long>(result.cost) << " optimal=" << result.optimal
              << '\n';
    for (int id : result.edges) std::cout << id << ' ';
    std::cout << '\n';
    auto checkpoint = solver; // Full in-memory checkpoint, including pending work and RNG.
    (void)checkpoint;
}
