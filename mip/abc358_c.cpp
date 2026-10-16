// https://atcoder.jp/contests/abc358/tasks/abc358_c
#include "mip.hpp"
#include <cassert>
#include <iostream>
#include <string>
int main() {
    int n, k;
    std::cin >> n >> k;
    std::vector<std::string> s(n);
    for (auto &row : s) std::cin >> row;
    cpopt::Model m;
    cpopt::Terms cost;
    for (int i = 0; i < n; ++i) cost.emplace_back(m.add_binary(), 1);
    for (int j = 0; j < k; ++j) {
        cpopt::Terms row;
        for (int i = 0; i < n; ++i)
            if (s[i][j] == 'o') row.emplace_back(i, 1);
        m.add_ge(row, 1);
    }
    m.minimize(cost);
    auto r = m.solve();
    assert(r.status == cpopt::Status::Optimal);
    std::cout << std::llround(r.objective) << '\n';
}
