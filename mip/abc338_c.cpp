// https://atcoder.jp/contests/abc338/tasks/abc338_c
#include "mip.hpp"
#include <cassert>
#include <iostream>
int main() {
    int n;
    std::cin >> n;
    std::vector<int> q(n), a(n), b(n);
    for (int &v : q) std::cin >> v;
    for (int &v : a) std::cin >> v;
    for (int &v : b) std::cin >> v;
    cpopt::Model m;
    int x = m.add_integer(), y = m.add_integer();
    for (int i = 0; i < n; ++i) m.add_le({{x, double(a[i])}, {y, double(b[i])}}, q[i]);
    m.maximize({{x, 1}, {y, 1}});
    auto r = m.solve();
    assert(r.status == cpopt::Status::Optimal);
    std::cout << std::llround(r.objective) << '\n';
}
