// SPDX-License-Identifier: CC0-1.0
#include "set_cover_partition.hpp"
#include <iomanip>
#include <iostream>
#include <string>
int main() {
    // Input: cover|partition, element count, set count, seconds; then cost size elements...
    // Elements and set IDs are zero-based. Optional final line: seed size, selected IDs.
    using namespace set_cover_partition;
    std::string mode;
    int n, m;
    double seconds;
    if (!(std::cin >> mode >> n >> m >> seconds) || (mode != "cover" && mode != "partition"))
        return 1;
    Model problem(n, mode == "cover" ? Mode::Cover : Mode::Partition);
    for (int j = 0; j < m; ++j) {
        Cost cost;
        int k;
        if (!(std::cin >> cost >> k) || k < 0) return 1;
        std::vector<int> elements(k);
        for (int &i : elements)
            if (!(std::cin >> i)) return 1;
        problem.add_set(cost, std::move(elements));
    }
    Options options;
    options.time_limit = seconds;
    int k;
    if (std::cin >> k) {
        if (k < 0) return 1;
        options.initial_sets.resize(k);
        for (int &j : options.initial_sets)
            if (!(std::cin >> j)) return 1;
    }
    auto r = solve(problem, options);
    std::cout << (r.feasible ? "feasible" : "incomplete") << ' ' << std::setprecision(18) << r.cost
              << " violation=" << r.violation << '\n';
    for (int j : r.selected) std::cout << j << ' ';
    std::cout << '\n';
}
