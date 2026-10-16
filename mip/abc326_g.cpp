// https://atcoder.jp/contests/abc326/tasks/abc326_g
// スキルのレベルは連続変数、課題を達成するかは 0/1 変数。
#include "mip.hpp"
#include <cassert>
#include <iostream>
int main() {
    int n, k;
    std::cin >> n >> k;
    std::vector<int> c(n), a(k), x(n), y(k);
    for (int &v : c) std::cin >> v;
    for (int &v : a) std::cin >> v;
    cpopt::Model m;
    cpopt::Terms cost;
    double constant = 0;
    for (int i = 0; i < n; ++i) {
        x[i] = m.add_variable(1, 5);
        cost.emplace_back(x[i], -c[i]);
        constant += c[i];
    }
    for (int j = 0; j < k; ++j) {
        y[j] = m.add_binary();
        cost.emplace_back(y[j], a[j]);
        for (int i = 0; i < n; ++i) {
            int level;
            std::cin >> level;
            m.add_ge({{x[i], 1}, {y[j], -double(level)}}, 0);
        }
    }
    m.maximize(cost, constant);
    auto r = m.solve();
    assert(r.status == cpopt::Status::Optimal);
    std::cout << std::llround(r.objective) << '\n';
}
