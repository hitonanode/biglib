// https://atcoder.jp/contests/abc364/tasks/abc364_e
#include "mip.hpp"
#include <cassert>
#include <iostream>
int main() {
    int n, x, y;
    std::cin >> n >> x >> y;
    cpopt::Model m;
    cpopt::Terms a, b, cost;
    for (int i = 0; i < n; ++i) {
        int ai, bi;
        std::cin >> ai >> bi;
        int v = m.add_binary();
        a.emplace_back(v, ai);
        b.emplace_back(v, bi);
        cost.emplace_back(v, 1);
    }
    m.add_le(a, x);
    m.add_le(b, y);
    m.maximize(cost);
    auto r = m.solve();
    assert(r.status == cpopt::Status::Optimal);
    std::cout << std::min<long long>(n, std::llround(r.objective) + 1) << '\n';
}
