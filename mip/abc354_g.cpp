// https://atcoder.jp/contests/abc354/tasks/abc354_g
// 文字列の採用を0/1変数で表し、部分文字列の関係にある2文字列の同時採用を禁止する。
#include "mip.hpp"
#include <cassert>
#include <iostream>
#include <string>

int main() {
    int n;
    std::cin >> n;
    std::vector<std::string> s(n);
    for (auto &v : s) std::cin >> v;
    cpopt::Model model;
    std::vector<int> x(n);
    cpopt::Terms value;
    for (int i = 0; i < n; ++i) {
        int a;
        std::cin >> a;
        x[i] = model.add_binary();
        value.emplace_back(x[i], a);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (s[i].find(s[j]) != std::string::npos || s[j].find(s[i]) != std::string::npos)
                model.add_le({{x[i], 1}, {x[j], 1}}, 1);
        }
    }
    model.maximize(value);
    auto result = model.solve();
    assert(result.status == cpopt::Status::Optimal);
    std::cout << std::llround(result.objective) << '\n';
}
