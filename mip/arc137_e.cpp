// https://atcoder.jp/contests/arc137/tasks/arc137_e
// 雇用量と日ごとの販売量を変数にする。区間の制約行列は全単模なので LP で解ける。
#include "mip.hpp"
#include <cassert>
#include <iostream>

int main() {
    int n, m, d;
    std::cin >> n >> m >> d;
    cpopt::Model model;
    cpopt::Terms profit;
    std::vector<cpopt::Terms> daily(n);
    for (int i = 0; i < n; ++i) {
        int a;
        std::cin >> a;
        int sell = model.add_variable(0, a);
        daily[i].emplace_back(sell, 1);
        profit.emplace_back(sell, d);
    }
    for (int i = 0; i < m; ++i) {
        int l, r, c;
        std::cin >> l >> r >> c;
        int hire = model.add_variable(0, 1);
        profit.emplace_back(hire, -c);
        for (int day = l - 1; day < r; ++day) daily[day].emplace_back(hire, -1);
    }
    for (const auto &row : daily) model.add_le(row, 0);
    model.maximize(profit);
    auto result = model.solve();
    assert(result.status == cpopt::Status::Optimal);
    std::cout << std::llround(result.objective) << '\n';
}
