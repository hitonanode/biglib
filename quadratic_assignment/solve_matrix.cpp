// SPDX-License-Identifier: CC0-1.0
#include "quadratic_assignment.hpp"
#include <iomanip>
#include <iostream>

int main() {
    // Input: n, seconds, A (n*n values), B (n*n values).
    int n;
    double seconds;
    if (!(std::cin >> n >> seconds) || n < 0) return 1;
    std::vector<std::vector<double>> a(n, std::vector<double>(n)), b = a;
    for (auto &row : a)
        for (auto &x : row)
            if (!(std::cin >> x)) return 1;
    for (auto &row : b)
        for (auto &x : row)
            if (!(std::cin >> x)) return 1;
    quadratic_assignment::Options options;
    options.time_limit = seconds;
    auto r = quadratic_assignment::solve(a, b, options);
    std::cout << std::setprecision(17) << r.cost << '\n';
    for (int i = 0; i < n; ++i) std::cout << (i ? " " : "") << r.assignment[i];
    std::cout << '\n';
}
