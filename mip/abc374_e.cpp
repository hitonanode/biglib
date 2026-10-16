// https://atcoder.jp/contests/abc374/tasks/abc374_e
// 各機械の購入数と生産量を整数変数にし、予算内で全工程の最小生産量を最大化する。
#include "mip.hpp"
#include <cassert>
#include <iostream>

int main() {
    int n, budget;
    std::cin >> n >> budget;
    cpopt::Model model;
    int production = model.add_integer();
    cpopt::Terms cost;
    for (int i = 0; i < n; ++i) {
        int a, p, b, q;
        std::cin >> a >> p >> b >> q;
        int x = model.add_integer(0, budget / p);
        int y = model.add_integer(0, budget / q);
        cost.emplace_back(x, p);
        cost.emplace_back(y, q);
        model.add_le({{production, 1}, {x, -double(a)}, {y, -double(b)}}, 0);
    }
    model.add_le(cost, budget);
    model.maximize({{production, 1}});
    auto result = model.solve();
    assert(result.status == cpopt::Status::Optimal);
    std::cout << std::llround(result.objective) << '\n';
}
