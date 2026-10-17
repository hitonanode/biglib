// SPDX-License-Identifier: CC0-1.0
#include "symmetric_tsp.hpp"
#include <iomanip>
#include <iostream>

int main() {
    int n;
    double seconds;
    if (!(std::cin >> n >> seconds) || n < 0) return 1;
    std::vector<std::vector<double>> cost(n, std::vector<double>(n));
    for (auto &row : cost)
        for (auto &v : row)
            if (!(std::cin >> v)) return 1;
    symmetric_tsp::Options options;
    options.time_limit = seconds;
    auto result = symmetric_tsp::solve(cost, options);
    std::cout << std::setprecision(18) << result.cost << '\n';
    for (int v : result.tour) std::cout << v << ' ';
    std::cout << '\n';
}
