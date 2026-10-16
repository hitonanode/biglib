// SPDX-License-Identifier: CC0-1.0
#include "maximum_weight_matching.hpp"
#include <iostream>

int main() {
    using Solver = maximum_weight_matching::Solver<>;
    std::vector<Solver::Edge> edges = {{0, 1, 6}, {0, 2, 10}, {1, 3, 1}};
    Solver solver(4, edges, {0}); // 初期解: 辺0。未マッチ頂点があってよい。
    auto result = solver.solve();
    std::cout << maximum_weight_matching::to_string(result.weight) << '\n'; // 11
    for (int id : result.edges) std::cout << edges[id].u << ' ' << edges[id].v << '\n';
}
