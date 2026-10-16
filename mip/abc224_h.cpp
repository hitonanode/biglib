// https://atcoder.jp/contests/abc224/tasks/abc224_h
#include "mip.hpp"
#include <cassert>
#include <iostream>
int main() {
    int l, r;
    std::cin >> l >> r;
    cpopt::Model m;
    cpopt::Terms cost;
    for (int i = 0; i < l + r; ++i) {
        int c;
        std::cin >> c;
        cost.emplace_back(m.add_variable(), c);
    }
    for (int i = 0; i < l; ++i)
        for (int j = 0; j < r; ++j) {
            int c;
            std::cin >> c;
            m.add_ge({{i, 1}, {l + j, 1}}, c);
        }
    m.minimize(cost);
    auto result = m.solve();
    assert(result.status == cpopt::Status::Optimal);
    std::cout << std::llround(result.objective) << '\n';
}
