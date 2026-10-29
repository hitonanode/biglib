// SPDX-License-Identifier: CC0-1.0
#include "constraint_search.hpp"
#include <iostream>
int main() {
    using namespace constraint_search;
    Model m;
    int a = m.add_integer(0, 30), b = m.add_integer(0, 30), c = m.add_integer(0, 30);
    m.add_no_overlap({Interval(a, 4), Interval(b, 7), Interval(c, 3)});
    m.add_ge({{c, 1}, {a, -1}}, 4);
    int finish = m.add_maximum({Affine(a, 1, 4), Affine(b, 1, 7), Affine(c, 1, 3)});
    m.minimize({{finish, 1}});
    Options o;
    o.initial_solution = {0, unknown, 4};
    o.seed = 42;
    Solver solver(m, o);
    solver.run(.01);
    const auto &r = solver.run(.01);
    if (r.feasible) {
        std::cout << "makespan " << r.objective << '\n';
        for (int v : {a, b, c}) std::cout << r.x[v] << ' ';
        std::cout << '\n';
    }
}
