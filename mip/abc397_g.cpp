// https://atcoder.jp/contests/abc397/tasks/abc397_g
// 辺の重みを0/1変数、距離を連続変数とし、h[v] <= h[u] + w を課す。
#include "mip.hpp"
#include <cassert>
#include <iostream>

int main() {
    int n, m, k;
    std::cin >> n >> m >> k;
    cpopt::Model model;
    std::vector<int> h(n);
    for (int i = 0; i < n; ++i) h[i] = model.add_variable(0, i == 0 ? 0 : n - 1);
    cpopt::Terms selected;
    for (int i = 0; i < m; ++i) {
        int u, v;
        std::cin >> u >> v;
        --u;
        --v;
        int w = model.add_binary();
        selected.emplace_back(w, 1);
        model.add_le({{h[v], 1}, {h[u], -1}, {w, -1}}, 0);
    }
    model.add_eq(selected, k);
    model.maximize({{h.back(), 1}});
    auto result = model.solve();
    assert(result.status == cpopt::Status::Optimal);
    std::cout << std::llround(result.objective) << '\n';
}
