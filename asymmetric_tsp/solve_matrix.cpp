// SPDX-License-Identifier: CC0-1.0
// stdin: n, then an n*n asymmetric cost matrix. stdout: cost, then a permutation.
#include "asymmetric_tsp.hpp"
#include <iomanip>
#include <iostream>

int main() {
    int n;
    if (!(std::cin >> n) || n < 0) return 1;
    std::vector<std::vector<long long>> cost(n, std::vector<long long>(n));
    for (auto &row : cost)
        for (auto &v : row)
            if (!(std::cin >> v)) return 1;
    asymmetric_tsp::Options options;
    options.time_limit = 1.8;
    auto result = asymmetric_tsp::solve(cost, options);
    std::cout << std::setprecision(18) << result.cost << '\n';
    for (int v : result.tour) std::cout << v << ' ';
    std::cout << '\n';
}
