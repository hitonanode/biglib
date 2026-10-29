#ifndef CONSTRAINT_SEARCH_HPP
#define CONSTRAINT_SEARCH_HPP
#include <algorithm>
#include <chrono>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace constraint_search {
using Int = std::int64_t;
using Wide = __int128_t;
using Clock = std::chrono::steady_clock;
using Terms = std::vector<std::pair<int, Int>>;
constexpr Int bound = Int(1) << 60;
constexpr Int unknown = std::numeric_limits<Int>::min();
struct Affine {
    int var = -1;
    Int coefficient = 1, offset = 0;
    Affine() = default;
    explicit Affine(int v, Int a = 1, Int b = 0) : var(v), coefficient(a), offset(b) {}
    static Affine constant(Int b) {
        Affine a;
        a.coefficient = 0;
        a.offset = b;
        return a;
    }
};
struct Literal {
    int var = -1;
    bool positive = true;
    Literal() = default;
    explicit Literal(int v, bool p = true) : var(v), positive(p) {}
};
inline Literal pos(int v) { return Literal(v); }
inline Literal neg(int v) { return Literal(v, false); }
struct Interval {
    Affine start, duration;
    Literal present;
    Interval(int s, Int d, Literal p = {}) : start(s), duration(Affine::constant(d)), present(p) {}
    Interval(Affine s, Affine d, Literal p = {}) : start(s), duration(d), present(p) {}
};
struct Evaluation {
    bool feasible = false;
    int violated = 0;
    long double violation = 0;
    Int objective = 0;
};
namespace detail {
inline Wide absw(Wide a) { return a < 0 ? -a : a; }
inline Int narrow(Wide a) {
    if (a < -Wide(bound) || a > Wide(bound))
        throw std::invalid_argument("integer expression exceeds +/-2^60");
    return Int(a);
}
inline Wide floor_div(Wide a, Wide b) {
    if (b < 0) {
        a = -a;
        b = -b;
    }
    Wide q = a / b;
    return q - (a % b < 0);
}
enum class Kind { Linear, Different, Overlap, Circuit, Xor };
struct Constraint {
    Kind kind = Kind::Linear;
    Terms terms;
    Int lower = -bound, upper = bound;
    std::vector<Affine> expressions;
    std::vector<Interval> intervals;
    std::vector<int> successors;
    std::vector<Literal> guards, literals;
    bool self_loops = false;
    long double scale = 1;
};
struct Definition {
    enum Kind { Element, Maximum } kind;
    int output, index = -1;
    std::vector<Int> table;
    std::vector<Affine> args;
};
struct Reference {
    int constraint, slot;
};
struct NoTick {
    void operator()() const {}
};
// A small, resumable CDCL engine used only by recognized Boolean encodings.
// Literal 2*v denotes v, and 2*v+1 denotes its negation.
class BooleanSearch {
    struct Interrupted {};
    struct Clause {
        std::size_t start, length;
    };
    struct Row {
        int *data;
        std::size_t length;
        int &operator[](std::size_t i) { return data[i]; }
        std::size_t size() const { return length; }
        bool empty() const { return !length; }
        int *begin() const { return data; }
        int *end() const { return data + length; }
    };
    std::vector<Clause> clauses_;
    std::vector<int> literals_;
    std::vector<std::vector<int>> watches_;
    std::vector<std::vector<int>> groups_, membership_;
    std::vector<int> group_free_, group_true_;
    std::vector<int> assignment_, level_, reason_, phase_, trail_, limits_;
    std::vector<double> activity_;
    std::priority_queue<std::pair<double, int>> heap_;
    std::size_t head_ = 0, watch_cursor_ = 0, literal_count_ = 0;
    int conflict_ = -1;
    bool heap_ready_ = false, inconsistent_ = false, solved_ = false, saturated_ = false;
    bool group_choice_ = false;
    bool sampled_groups_ = false;
    std::uint64_t conflict_limit_ = 20000;
    std::size_t group_cursor_ = 0;
    double increment_ = 1;
    std::uint64_t conflicts_ = 0, restart_index_ = 1, next_restart_ = 64;
    Clock::time_point deadline_;
    Row clause(int id) {
        const auto c = clauses_[id];
        return {literals_.data() + c.start, c.length};
    }
    void check() const {
        if (Clock::now() >= deadline_) throw Interrupted{};
    }
    int value(int literal) const {
        int a = assignment_[literal / 2];
        return a < 0 ? -1 : a ^ (literal & 1);
    }
    void change_group(int v, int before, int after) {
        for (int group : membership_[v]) {
            group_free_[group] += int(after < 0) - int(before < 0);
            group_true_[group] += int(after == 1) - int(before == 1);
        }
    }
    bool enqueue(int literal, int reason) {
        int v = literal / 2, a = 1 ^ (literal & 1);
        if (assignment_[v] >= 0) return assignment_[v] == a;
        change_group(v, -1, a);
        assignment_[v] = phase_[v] = a;
        level_[v] = int(limits_.size());
        reason_[v] = reason;
        trail_.push_back(literal);
        return true;
    }
    void backtrack(int target) {
        if (int(limits_.size()) <= target) return;
        for (std::size_t i = trail_.size(); i > std::size_t(limits_[target]);) {
            int v = trail_[--i] / 2;
            change_group(v, assignment_[v], -1);
            assignment_[v] = -1;
            reason_[v] = -1;
            // Group branching may never consult the fallback heap. Keep it
            // dirty until first use, and bound stale entries after that use.
            if (heap_ready_) {
                if (heap_.size() < assignment_.size() * 4 + 1024)
                    heap_.emplace(activity_[v], v);
                else
                    heap_ready_ = false;
            }
        }
        trail_.resize(limits_[target]);
        limits_.resize(target);
        head_ = std::min(head_, trail_.size());
        watch_cursor_ = 0;
    }
    int propagate() {
        while (head_ < trail_.size()) {
            int false_literal = trail_[head_] ^ 1;
            auto &list = watches_[false_literal];
            while (watch_cursor_ < list.size()) {
                check();
                int id = list[watch_cursor_];
                auto row = clause(id);
                if (row[0] == false_literal) std::swap(row[0], row[1]);
                assert(row[1] == false_literal);
                if (value(row[0]) == 1) {
                    ++watch_cursor_;
                    continue;
                }
                std::size_t replacement = 2;
                while (replacement < row.size() && value(row[replacement]) == 0) {
                    if ((replacement & 127) == 0) check();
                    ++replacement;
                }
                if (replacement < row.size()) {
                    std::swap(row[1], row[replacement]);
                    list[watch_cursor_] = list.back();
                    list.pop_back();
                    watches_[row[1]].push_back(id);
                } else {
                    if (!enqueue(row[0], id)) return id;
                    ++watch_cursor_;
                }
            }
            ++head_;
            watch_cursor_ = 0;
        }
        return -1;
    }
    std::pair<std::vector<int>, int> analyze() {
        std::vector<bool> seen(assignment_.size(), false);
        std::vector<int> learned(1);
        int paths = 0, pivot = -1, row = conflict_, at = int(trail_.size()) - 1;
        std::size_t scanned = 0;
        do {
            check();
            assert(row >= 0);
            for (int literal : clause(row)) {
                int v = literal / 2;
                if ((++scanned & 127) == 0) check();
                if (v == pivot / 2 && pivot >= 0) continue;
                if (seen[v] || level_[v] == 0) continue;
                seen[v] = true;
                activity_[v] += increment_;
                if (level_[v] == int(limits_.size()))
                    ++paths;
                else
                    learned.push_back(literal);
            }
            while (at >= 0 && !seen[trail_[at] / 2]) {
                if ((at & 127) == 0) check();
                --at;
            }
            assert(at >= 0);
            pivot = trail_[at--];
            seen[pivot / 2] = false;
            --paths;
            row = reason_[pivot / 2];
        } while (paths > 0);
        learned[0] = pivot ^ 1;
        int target = 0;
        for (std::size_t i = 1; i < learned.size(); ++i) {
            if ((i & 127) == 0) check();
            if (level_[learned[i] / 2] > target) {
                target = level_[learned[i] / 2];
                std::swap(learned[1], learned[i]);
            }
        }
        return {std::move(learned), target};
    }
    int pick() {
        group_choice_ = false;
        int group = -1;
        double group_activity = -1;
        for (std::size_t offset = 0; offset < groups_.size(); ++offset) {
            if ((offset & 127) == 0) check();
            if (sampled_groups_ && offset >= 64 && group >= 0) break;
            int c = int((group_cursor_ + offset) % groups_.size());
            if (group_true_[c] || !group_free_[c]) continue;
            double score = 0;
            if (sampled_groups_)
                for (int v : groups_[c])
                    if (assignment_[v] < 0) score = std::max(score, activity_[v]);
            if (group < 0 || group_free_[c] < group_free_[group] ||
                (sampled_groups_ && group_free_[c] == group_free_[group] && score > group_activity)) {
                group = c;
                group_activity = score;
            }
            if (!sampled_groups_ && group_free_[group] <= 2) break;
        }
        if (group >= 0) {
            int best = -1;
            for (int v : groups_[group])
                if (assignment_[v] < 0) {
                    if (best < 0 || activity_[v] > activity_[best] ||
                        (activity_[v] == activity_[best] && phase_[v] > phase_[best]))
                        best = v;
                }
            assert(best >= 0);
            group_choice_ = true;
            group_cursor_ = (std::size_t(group) + 1) % groups_.size();
            return best;
        }
        if (!heap_ready_ || heap_.size() > assignment_.size() * 4 + 1024) {
            std::priority_queue<std::pair<double, int>> fresh;
            for (int v = 0; v < int(assignment_.size()); ++v) {
                if ((v & 127) == 0) check();
                if (assignment_[v] < 0) fresh.emplace(activity_[v], v);
            }
            heap_ = std::move(fresh);
            heap_ready_ = true;
        }
        while (!heap_.empty()) {
            check();
            auto a = heap_.top();
            heap_.pop();
            if (assignment_[a.second] < 0 && activity_[a.second] == a.first) return a.second;
        }
        return -1;
    }
    static std::uint64_t luby(std::uint64_t index) {
        std::uint64_t size = 1, power = 1;
        while (size < index) {
            size = size * 2 + 1;
            power *= 2;
        }
        while (size != index) {
            size = (size - 1) / 2;
            power /= 2;
            if (index > size) index -= size;
        }
        return power;
    }

public:
    void configure_groups(bool sampled, std::uint64_t conflict_limit) {
        sampled_groups_ = sampled;
        conflict_limit_ = conflict_limit;
    }
    int add_variable(bool preferred = false) {
        int v = int(assignment_.size());
        assignment_.push_back(-1);
        level_.push_back(0);
        reason_.push_back(-1);
        phase_.push_back(preferred);
        activity_.push_back(0);
        watches_.resize(2 * (v + 1));
        membership_.emplace_back();
        heap_ready_ = false;
        return v;
    }
    void add_branch_group(const std::vector<int> &positive_literals) {
        int id = int(groups_.size()), free = 0, yes = 0;
        std::vector<int> variables;
        for (int literal : positive_literals) {
            assert(!(literal & 1));
            int v = literal / 2;
            variables.push_back(v);
            membership_[v].push_back(id);
            free += assignment_[v] < 0;
            yes += assignment_[v] == 1;
        }
        groups_.push_back(std::move(variables));
        group_free_.push_back(free);
        group_true_.push_back(yes);
    }
    // Input clauses must have distinct variables and contain no tautology.
    int add_clause(std::vector<int> row, bool learned = false) {
        int id = int(clauses_.size());
        literal_count_ += row.size();
        clauses_.push_back({literals_.size(), row.size()});
        literals_.insert(literals_.end(), row.begin(), row.end());
        if (!learned)
            for (int literal : row) activity_[literal / 2] += 1;
        if (row.empty())
            inconsistent_ = true;
        else if (row.size() == 1) {
            if (!enqueue(row[0], id)) inconsistent_ = true;
        } else {
            watches_[row[0]].push_back(id);
            watches_[row[1]].push_back(id);
        }
        return id;
    }
    // 1: satisfying assignment; -1: root contradiction; 0: paused.
    int run(Clock::time_point deadline) {
        if (inconsistent_) return -1;
        if (solved_) return 1;
        if (saturated_) return 0;
        deadline_ = deadline;
        try {
            for (;;) {
                check();
                if (conflict_ < 0) conflict_ = propagate();
                if (conflict_ >= 0) {
                    if (limits_.empty()) {
                        inconsistent_ = true;
                        return -1;
                    }
                    auto analysis = analyze();
                    if (literal_count_ + analysis.first.size() > 2000000 ||
                        conflicts_ >= conflict_limit_) {
                        saturated_ = true;
                        return 0;
                    }
                    backtrack(analysis.second);
                    int literal = analysis.first[0],
                        id = add_clause(std::move(analysis.first), true);
                    if (!enqueue(literal, id)) {
                        inconsistent_ = true;
                        return -1;
                    }
                    conflict_ = -1;
                    ++conflicts_;
                    increment_ *= 1.05;
                    if (increment_ > 1e100) {
                        for (double &a : activity_) a *= 1e-100;
                        increment_ *= 1e-100;
                        heap_ready_ = false;
                    }
                    if (conflicts_ >= next_restart_) {
                        backtrack(0);
                        next_restart_ = conflicts_ + 64 * luby(++restart_index_);
                    }
                } else {
                    int v = pick();
                    if (v < 0) {
                        solved_ = true;
                        return 1;
                    }
                    limits_.push_back(int(trail_.size()));
                    enqueue(2 * v + (group_choice_ ? 0 : !phase_[v]), -1);
                }
            }
        } catch (const Interrupted &) { return 0; }
    }
    const std::vector<int> &assignment() const { return assignment_; }
    std::size_t branch_heap_size() const { return heap_.size(); }
    std::uint64_t conflicts() const { return conflicts_; }
    bool saturated() const { return saturated_; }
};
template <class Get> inline Int value(Affine a, const Get &get) {
    return Int(Wide(a.offset) + (a.var < 0 ? 0 : Wide(a.coefficient) * get(a.var)));
}
template <class Get> inline bool truth(Literal a, const Get &get) {
    return a.var < 0 ? a.positive : (get(a.var) != 0) == a.positive;
}
template <class Get> inline bool active(const Constraint &c, const Get &get) {
    for (auto a : c.guards)
        if (!truth(a, get)) return false;
    return true;
}
template <class Get> inline bool collides(const Interval &a, const Interval &b, const Get &get) {
    if (!truth(a.present, get) || !truth(b.present, get)) return false;
    const Wide sa = value(a.start, get), sb = value(b.start, get);
    return sa + value(a.duration, get) > sb && sb + value(b.duration, get) > sa;
}
template <class Get, class Tick>
Wide raw_violation(const Constraint &c, const Get &get, Tick tick) {
    if (c.kind == Kind::Linear) {
        Wide a = 0;
        for (auto t : c.terms) {
            tick();
            a += Wide(t.second) * get(t.first);
        }
        return std::max({Wide(0), Wide(c.lower) - a, a - c.upper});
    }
    if (c.kind == Kind::Different) {
        std::unordered_map<Int, int> frequency;
        Wide v = 0;
        for (auto a : c.expressions) {
            tick();
            v += frequency[value(a, get)]++;
        }
        return v;
    }
    if (c.kind == Kind::Overlap) {
        // Sort zero-duration intervals first at equal starts. They conflict with
        // an interval containing their point strictly, but not with its ends.
        std::vector<std::pair<Int, Int>> ranges;
        for (const auto &a : c.intervals) {
            tick();
            if (truth(a.present, get)) {
                Int s = value(a.start, get);
                ranges.emplace_back(s, Int(Wide(s) + value(a.duration, get)));
            }
        }
        std::sort(ranges.begin(), ranges.end(), [&](const auto &a, const auto &b) {
            tick();
            return a < b;
        });
        std::priority_queue<Int, std::vector<Int>, std::greater<Int>> ends;
        Wide v = 0;
        for (auto a : ranges) {
            tick();
            while (!ends.empty() && ends.top() <= a.first) {
                tick();
                ends.pop();
            }
            v += ends.size();
            if (a.second > a.first) ends.push(a.second);
        }
        return v;
    }
    if (c.kind == Kind::Xor) {
        bool parity = false;
        for (auto a : c.literals) {
            tick();
            parity ^= truth(a, get);
        }
        return !parity;
    }
    const int n = int(c.successors.size());
    std::vector<int> degree(n, 0), color(n, -1);
    Wide v = 0;
    for (int i = 0; i < n; ++i) {
        tick();
        int j = int(get(c.successors[i]));
        ++degree[j];
        if (!c.self_loops && n > 1 && i == j) ++v;
    }
    for (int d : degree) {
        tick();
        v += std::abs(d - 1);
    }
    int cycles = 0;
    for (int i = 0; i < n; ++i) {
        tick();
        if (color[i] >= 0) continue;
        int j = i;
        while (color[j] < 0) {
            tick();
            color[j] = i;
            j = int(get(c.successors[j]));
        }
        if (color[j] == i && !(c.self_loops && get(c.successors[j]) == j)) ++cycles;
    }
    v += std::max(0, cycles - 1);
    return v;
}
} // namespace detail

class Solver;
class Model {
    friend class Solver;
    struct Variable {
        Int lower, upper;
        int definition = -1;
    };
    std::vector<Variable> variables_;
    std::vector<detail::Constraint> constraints_;
    std::vector<detail::Definition> definitions_;
    std::vector<std::vector<detail::Reference>> incidence_;
    std::vector<std::vector<int>> dependents_;
    std::vector<Int> objective_;
    Int offset_ = 0;
    int direction_ = 1;
    std::uint64_t revision_ = 0;
    void check_var(int v) const {
        if (v < 0 || v >= int(variables_.size()))
            throw std::invalid_argument("invalid variable ID");
    }
    void check_literal(Literal a) const {
        if (a.var < 0) {
            if (a.var != -1) throw std::invalid_argument("invalid literal");
            return;
        }
        check_var(a.var);
        if (variables_[a.var].lower < 0 || variables_[a.var].upper > 1)
            throw std::invalid_argument("literal requires Boolean domain");
    }
    std::pair<Int, Int> range(Affine a) const {
        detail::narrow(a.offset);
        detail::narrow(a.coefficient);
        if (a.var < 0) {
            if (a.var != -1) throw std::invalid_argument("invalid affine expression");
            return {a.offset, a.offset};
        }
        check_var(a.var);
        Int l = detail::narrow(Wide(a.offset) + Wide(a.coefficient) * variables_[a.var].lower);
        Int u = detail::narrow(Wide(a.offset) + Wide(a.coefficient) * variables_[a.var].upper);
        return {std::min(l, u), std::max(l, u)};
    }
    int insert(detail::Constraint c) {
        const int id = int(constraints_.size());
        auto ref = [&](int v, int s) {
            if (v >= 0) incidence_[v].push_back({id, s});
        };
        for (auto a : c.guards) check_literal(a);
        for (auto a : c.guards) ref(a.var, -1);
        for (int i = 0; i < int(c.terms.size()); ++i) ref(c.terms[i].first, i);
        for (int i = 0; i < int(c.expressions.size()); ++i) ref(c.expressions[i].var, i);
        for (int i = 0; i < int(c.intervals.size()); ++i) {
            auto a = c.intervals[i];
            ref(a.start.var, i);
            ref(a.duration.var, i);
            ref(a.present.var, i);
        }
        for (int i = 0; i < int(c.successors.size()); ++i) ref(c.successors[i], i);
        for (int i = 0; i < int(c.literals.size()); ++i) ref(c.literals[i].var, i);
        constraints_.push_back(std::move(c));
        ++revision_;
        return id;
    }

public:
    int num_variables() const { return int(variables_.size()); }
    int num_constraints() const { return int(constraints_.size()); }
    int add_integer(Int lower, Int upper) {
        detail::narrow(lower);
        detail::narrow(upper);
        if (lower > upper) throw std::invalid_argument("empty variable domain");
        variables_.push_back({lower, upper, -1});
        incidence_.emplace_back();
        dependents_.emplace_back();
        objective_.push_back(0);
        ++revision_;
        return num_variables() - 1;
    }
    int add_bool() { return add_integer(0, 1); }
    int add_linear(Terms terms, Int lower, Int upper, std::vector<Literal> guards = {}) {
        detail::narrow(lower);
        detail::narrow(upper);
        if (lower > upper) throw std::invalid_argument("empty linear constraint range");
        for (auto t : terms) {
            check_var(t.first);
            detail::narrow(t.second);
        }
        std::sort(terms.begin(), terms.end());
        Terms normalized;
        for (std::size_t i = 0; i < terms.size();) {
            const int v = terms[i].first;
            Wide a = 0;
            do { a += terms[i++].second; } while (i < terms.size() && terms[i].first == v);
            if (a) normalized.emplace_back(v, detail::narrow(a));
        }
        Wide magnitude = 0;
        long double scale = 1;
        for (auto t : normalized) {
            magnitude += detail::absw(t.second) * std::max(detail::absw(variables_[t.first].lower),
                                                           detail::absw(variables_[t.first].upper));
            detail::narrow(magnitude);
            scale = std::max(scale, std::abs(static_cast<long double>(t.second)));
        }
        detail::narrow(magnitude);
        detail::Constraint c;
        c.terms = std::move(normalized);
        c.lower = lower;
        c.upper = upper;
        c.guards = std::move(guards);
        c.scale = scale;
        return insert(std::move(c));
    }
    int add_le(Terms t, Int b, std::vector<Literal> g = {}) {
        return add_linear(std::move(t), -bound, b, std::move(g));
    }
    int add_ge(Terms t, Int b, std::vector<Literal> g = {}) {
        return add_linear(std::move(t), b, bound, std::move(g));
    }
    int add_eq(Terms t, Int b, std::vector<Literal> g = {}) {
        return add_linear(std::move(t), b, b, std::move(g));
    }
    int add_bool_or(std::vector<Literal> literals, std::vector<Literal> guards = {}) {
        Terms t;
        Int b = 1;
        for (auto a : literals) {
            check_literal(a);
            if (a.var < 0) {
                b -= a.positive;
            } else {
                t.emplace_back(a.var, a.positive ? 1 : -1);
                b -= !a.positive;
            }
        }
        return add_ge(std::move(t), b, std::move(guards));
    }
    int add_implication(Literal a, Literal b) {
        a.positive = !a.positive;
        return add_bool_or({a, b});
    }
    int add_bool_xor(std::vector<Literal> literals) {
        for (auto a : literals) check_literal(a);
        detail::Constraint c;
        c.kind = detail::Kind::Xor;
        c.literals = std::move(literals);
        return insert(std::move(c));
    }
    int add_all_different(std::vector<Affine> expressions, std::vector<Literal> guards = {}) {
        for (auto a : expressions) range(a);
        detail::Constraint c;
        c.kind = detail::Kind::Different;
        c.expressions = std::move(expressions);
        c.guards = std::move(guards);
        return insert(std::move(c));
    }
    int add_all_different(const std::vector<int> &vars) {
        std::vector<Affine> a;
        for (int v : vars) a.emplace_back(v);
        return add_all_different(std::move(a));
    }
    int add_all_different(std::initializer_list<int> vars) {
        return add_all_different(std::vector<int>(vars));
    }
    int add_no_overlap(std::vector<Interval> intervals, std::vector<Literal> guards = {}) {
        for (auto a : intervals) {
            const auto s = range(a.start), d = range(a.duration);
            check_literal(a.present);
            if (d.first < 0) throw std::invalid_argument("duration domain must be nonnegative");
            detail::narrow(Wide(s.second) + d.second);
        }
        detail::Constraint c;
        c.kind = detail::Kind::Overlap;
        c.intervals = std::move(intervals);
        c.guards = std::move(guards);
        return insert(std::move(c));
    }
    // Successor IDs are local to this vector. allow_self_loops means subcircuit:
    // self-loop vertices are absent, and at most one nontrivial cycle remains.
    int add_circuit(std::vector<int> successors, bool allow_self_loops = false) {
        if (successors.empty())
            throw std::invalid_argument("circuit requires at least one vertex");
        for (int v : successors) {
            check_var(v);
            if (variables_[v].lower < 0 || variables_[v].upper >= Int(successors.size()))
                throw std::invalid_argument("successor domain out of range");
        }
        detail::Constraint c;
        c.kind = detail::Kind::Circuit;
        c.successors = std::move(successors);
        c.self_loops = allow_self_loops;
        return insert(std::move(c));
    }
    // Functional variables are exact invariants, not independent search variables.
    int add_element(int index, std::vector<Int> values) {
        check_var(index);
        if (values.empty() || variables_[index].lower < 0 ||
            variables_[index].upper >= Int(values.size()))
            throw std::invalid_argument("element index domain out of range");
        for (auto a : values) detail::narrow(a);
        const auto mm = std::minmax_element(values.begin(), values.end());
        int out = add_integer(*mm.first, *mm.second), d = int(definitions_.size());
        detail::Definition f;
        f.kind = detail::Definition::Element;
        f.output = out;
        f.index = index;
        f.table = std::move(values);
        definitions_.push_back(std::move(f));
        variables_[out].definition = d;
        dependents_[index].push_back(out);
        return out;
    }
    int add_maximum(std::vector<Affine> args) {
        if (args.empty()) throw std::invalid_argument("maximum requires arguments");
        Int lo = -bound, hi = -bound;
        for (auto a : args) {
            auto r = range(a);
            lo = std::max(lo, r.first);
            hi = std::max(hi, r.second);
        }
        int out = add_integer(lo, hi), d = int(definitions_.size());
        for (auto a : args)
            if (a.var >= 0) dependents_[a.var].push_back(out);
        detail::Definition f;
        f.kind = detail::Definition::Maximum;
        f.output = out;
        f.args = std::move(args);
        definitions_.push_back(std::move(f));
        variables_[out].definition = d;
        return out;
    }
    void minimize(Terms terms, Int offset = 0) {
        detail::narrow(offset);
        std::vector<Wide> a(num_variables(), 0);
        for (auto t : terms) {
            check_var(t.first);
            detail::narrow(t.second);
            a[t.first] += t.second;
        }
        Wide magnitude = detail::absw(offset);
        for (int i = 0; i < num_variables(); ++i) {
            detail::narrow(a[i]);
            magnitude += detail::absw(a[i]) * std::max(detail::absw(variables_[i].lower),
                                                       detail::absw(variables_[i].upper));
            detail::narrow(magnitude);
        }
        detail::narrow(magnitude);
        for (int i = 0; i < num_variables(); ++i) objective_[i] = Int(a[i]);
        offset_ = offset;
        direction_ = 1;
        ++revision_;
    }
    void maximize(Terms terms, Int offset = 0) {
        minimize(std::move(terms), offset);
        direction_ = -1;
    }
    template <class Get, class Tick> Int defined_value(int v, const Get &get, Tick tick) const {
        const auto &d = definitions_[variables_[v].definition];
        if (d.kind == detail::Definition::Element) return d.table[std::size_t(get(d.index))];
        Int a = -bound;
        for (auto e : d.args) {
            tick();
            a = std::max(a, detail::value(e, get));
        }
        return a;
    }
    Evaluation evaluate(const std::vector<Int> &x) const {
        if (x.size() != variables_.size()) throw std::invalid_argument("solution size mismatch");
        Evaluation r;
        auto get = [&](int v) { return x[v]; };
        bool domains = true;
        for (int i = 0; i < num_variables(); ++i)
            if (x[i] < variables_[i].lower || x[i] > variables_[i].upper) {
                domains = false;
                ++r.violated;
            }
        if (!domains) {
            r.violation = r.violated;
            return r;
        }
        Wide objective = offset_, violation = 0;
        for (int i = 0; i < num_variables(); ++i) {
            objective += Wide(objective_[i]) * x[i];
            if (variables_[i].definition >= 0) {
                Wide v = detail::absw(Wide(x[i]) - defined_value(i, get, detail::NoTick{}));
                violation += v;
                r.violated += v > 0;
            }
        }
        for (const auto &c : constraints_)
            if (detail::active(c, get)) {
                Wide v = detail::raw_violation(c, get, detail::NoTick{});
                violation += v;
                r.violated += v > 0;
            }
        r.objective = Int(objective);
        r.violation = static_cast<long double>(violation);
        r.feasible = r.violated == 0;
        return r;
    }
};

struct Options {
    double time_limit = 1.0;
    std::uint64_t seed = 1;
    std::uint64_t iteration_limit = std::numeric_limits<std::uint64_t>::max();
    std::vector<Int> initial_solution; // prefix or full vector; unknown marks holes
    Clock::time_point deadline = Clock::time_point::max();
    bool compound_moves = true, adaptive_weights = true, structural_initialization = true;
    bool schedule_moves = true;
    bool critical_block_moves = true;
    bool cover_moves = true;
    bool domain_reduction = true, alldifferent_reduction = true, permutation_moves = true,
         packing_moves = true, clause_moves = true;
    bool clause_scores = true, clause_pair_moves = true;
    bool latin_sat = true, permutation_restarts = true;
    bool fast_schedule = true;
    bool schedule_restarts = true;
    bool packing_dual_core = true, packing_cached_capacity = true;
    bool packing_subgradient = true;
    bool clause_multiscale = true;
    int clause_pair_samples = 8;
    int candidates = 24;
};
enum class StopReason { TimeLimit, IterationLimit, Satisfied, Fixed, ObjectiveBound };
struct Statistics {
    std::uint64_t iterations = 0, evaluations = 0, moves = 0, weight_updates = 0, restarts = 0;
    std::uint64_t latin_candidates = 0, latin_cells = 0, boolean_variables = 0,
                  boolean_clauses = 0, boolean_conflicts = 0;
    double elapsed = 0, total_elapsed = 0, first_feasible = -1;
};
struct Result : Evaluation {
    bool evaluated = false;
    std::vector<Int> x;
    Statistics stats;
    StopReason reason = StopReason::TimeLimit;
};

class Solver {
    friend Result solve(const Model &, Options);
    const Model &model_;
    Options options_;
    const std::uint64_t revision_;
    std::mt19937_64 rng_;
    Result result_;
    Statistics stats_;
    struct Cache {
        Wide raw = 0, value = 0, activity = 0;
        std::unordered_map<Int, int> frequency;
        std::unordered_map<Int, int> owner;
        std::vector<int> conflict_slots;
        std::uint64_t conflict_version = std::numeric_limits<std::uint64_t>::max();
    };
    struct Change {
        int var;
        Int value;
    };
    using Move = std::vector<Change>;
    struct Timeout {};
    Clock::time_point deadline_, started_;
    std::uint64_t ticks_ = 0, epoch_ = 0;
    std::vector<Int> x_, old_, undo_;
    std::vector<std::pair<Int, Int>> domains_;
    std::vector<std::uint64_t> domain_masks_;
    std::vector<std::uint64_t> changed_stamp_, constraint_stamp_, queued_stamp_, tabu_;
    std::vector<int> changed_, touched_, violated_, position_, decisions_, objective_vars_;
    std::vector<int> permutation_;
    std::vector<bool> permutation_constraint_;
    std::vector<std::vector<int>> slots_;
    std::vector<Cache> cache_;
    std::vector<Wide> next_raw_, next_value_, next_activity_;
    std::vector<long double> weights_;
    int phase_ = -2, cursor_ = 0;
    bool allocated_ = false, pending_ = false, has_objective_ = false;
    std::size_t commit_cursor_ = 0;
    Wide raw_total_ = 0;
    Int objective_ = 0, next_objective_ = 0;
    Int objective_bound_ = 0;
    long double objective_scale_ = 1, objective_weight_ = 1;
    std::uint64_t stall_ = 0;
    struct Schedule {
        std::vector<int> vars, machine, indegree;
        std::vector<Int> duration;
        std::vector<std::vector<std::pair<int, Int>>> edges;
        std::vector<std::vector<int>> order, best_order;
        std::vector<std::pair<int, int>> critical;
        std::unordered_map<std::uint64_t, std::uint64_t> tabu_arcs;
        std::uint64_t last_restart = 0;
    } schedule_;
    bool schedule_tried_ = false, schedule_mode_ = false;
    int route_constraint_ = -1;
    bool cover_mode_ = false;
    bool permutation_mode_ = false;
    bool recolor_mode_ = false;
    std::vector<std::uint64_t> color_tabu_;
    std::uint64_t last_recolor_restart_ = 0;
    std::vector<std::pair<int, std::unordered_map<Int, int>>> permutation_owner_reset_;
    struct LatinEncoding {
        std::vector<int> variables;
        std::vector<Int> values, fixed;
        std::vector<std::vector<int>> groups;
        int cells = 0;
    } latin_encoding_;
    bool latin_encoding_tried_ = false;
    detail::BooleanSearch latin_sat_;
    bool latin_sat_tried_ = false, latin_sat_done_ = false;
    double latin_sat_elapsed_ = 0;
    struct Packing {
        std::vector<std::vector<Int>> weight;
        std::vector<Int> capacity, profit;
        std::vector<long double> price;
        std::vector<long double> reduced;
        std::vector<long double> walk_price;
        long double dual_best = std::numeric_limits<long double>::infinity();
        int rounds = 0;
    } packing_;
    bool packing_mode_ = false;
    std::vector<int> clause_slack_;
    bool clause_mode_ = false;
    struct ClauseSearch {
        std::vector<Int> base_weight;
        std::vector<Wide> gain, penalty;
        std::vector<Wide> pair_adjustment;
        std::vector<std::uint64_t> pair_stamp;
        std::uint64_t pair_epoch = 0;
        std::uint64_t last_scale = 0;
        int scale = 1;
        std::vector<int> good, good_pos, bad, bad_pos;
        std::vector<std::uint64_t> last_flip;
        bool ready = false;
        int rebuild_phase = -1, rebuild_row = 0, rebuild_slot = 0;
    } clauses_;
    void check() const {
        if (Clock::now() >= deadline_) throw Timeout{};
    }
    static double cleanup_margin(const Model &model, const Options &options) {
        double margin =
            std::min(0.01, 0.00005 + model.num_constraints() * 4e-8 + model.num_variables() * 2e-8);
        if (options.latin_sat && model.definitions_.empty()) {
            bool objective = false;
            for (Int coefficient : model.objective_) objective |= coefficient != 0;
            if (!objective) {
                // SAT is considered only after dense fixed-value reduction,
                // for balanced groups of at most 64 values. Avoid reserving its
                // memory cleanup cost for ordinary queens/coloring instances.
                int potentially_fixed = 0;
                bool balanced = false;
                for (const auto &variable : model.variables_)
                    potentially_fixed += variable.lower == variable.upper;
                for (const auto &row : model.constraints_) {
                    if (!row.guards.empty()) continue;
                    if (row.kind == detail::Kind::Linear && row.terms.size() == 1)
                        ++potentially_fixed;
                    if (row.kind != detail::Kind::Different || row.expressions.empty() ||
                        row.expressions.size() > 64)
                        continue;
                    auto a = row.expressions[0];
                    if (a.var >= 0 && a.coefficient == 1 && !a.offset) {
                        const auto d = model.variables_[a.var];
                        balanced |= d.upper - d.lower + 1 == Int(row.expressions.size());
                    }
                }
                if (balanced && 2LL * potentially_fixed >= model.num_variables())
                    margin += std::min(0.002, 0.0002 + model.num_variables() * 3e-7);
            }
        }
        return margin;
    }
    void tick() {
        if ((++ticks_ & 127) == 0) check();
    }
    int random_index(std::size_t n) { return int(rng_() % n); }
    bool allowed(int v, Wide value) const {
        if (value < domains_[v].first || value > domains_[v].second) return false;
        return !domain_masks_[v] ||
               ((domain_masks_[v] >> int(value - model_.variables_[v].lower)) & 1);
    }
    Int random_value(int v) {
        if (domain_masks_[v]) {
            auto mask = domain_masks_[v];
            int k = random_index(__builtin_popcountll(mask));
            while (k--) mask &= mask - 1;
            return model_.variables_[v].lower + __builtin_ctzll(mask);
        }
        const auto d = domains_[v];
        return Int(Wide(d.first) + rng_() % std::uint64_t(Wide(d.second) - d.first + 1));
    }
    bool hinted(int v) const {
        return v < int(options_.initial_solution.size()) && options_.initial_solution[v] != unknown;
    }
    Int old_value(int v) const { return changed_stamp_[v] == epoch_ ? old_[v] : x_[v]; }
    void set_violation(int c, bool bad) {
        if (bad && position_[c] < 0) {
            position_[c] = int(violated_.size());
            violated_.push_back(c);
        }
        if (!bad && position_[c] >= 0) {
            int p = position_[c], last = violated_.back();
            violated_[p] = last;
            position_[last] = p;
            violated_.pop_back();
            position_[c] = -1;
        }
    }
    void remember() {
        const bool feasible = violated_.empty();
        const long double raw = static_cast<long double>(raw_total_);
        bool better = !result_.evaluated;
        if (result_.evaluated) {
            if (feasible != result_.feasible)
                better = feasible;
            else if (feasible)
                better = Wide(model_.direction_) * objective_ <
                         Wide(model_.direction_) * result_.objective;
            else
                better =
                    raw < result_.violation ||
                    (raw == result_.violation && Wide(model_.direction_) * objective_ <
                                                     Wide(model_.direction_) * result_.objective);
        }
        if (!better) return;
        result_.x = x_;
        result_.evaluated = true;
        result_.feasible = feasible;
        result_.violation = raw;
        result_.violated = int(violated_.size());
        result_.objective = objective_;
        stall_ = 0;
        if (feasible && stats_.first_feasible < 0)
            stats_.first_feasible = stats_.total_elapsed +
                                    std::chrono::duration<double>(Clock::now() - started_).count();
    }
    void structural_initialize(const detail::Constraint &c, int id) {
        if (c.kind == detail::Kind::Different && !c.expressions.empty() && c.guards.empty()) {
            const int n = int(c.expressions.size());
            const auto first = c.expressions.front();
            if (first.var < 0) return;
            const Int lower = model_.variables_[first.var].lower;
            std::vector<int> ids;
            for (auto a : c.expressions) {
                tick();
                if (a.var < 0 || a.coefficient != 1 || a.offset ||
                    model_.variables_[a.var].definition >= 0 ||
                    model_.variables_[a.var].lower != lower ||
                    Wide(model_.variables_[a.var].upper) != Wide(lower) + n - 1)
                    return;
                ids.push_back(a.var);
            }
            std::sort(ids.begin(), ids.end());
            if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) return;
            for (int v : ids)
                if (permutation_[v] >= 0 && permutation_[v] != id) return;
            if (options_.structural_initialization) {
                // A matching preserves fixed clues and valid hints. Publish it
                // only when complete, so a deadline leaves initialization resumable.
                std::vector<int> owner(n, -1), assignment(n, -1), seen(n, 0), order(n);
                std::iota(order.begin(), order.end(), 0);
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                    tick();
                    return domains_[ids[a]].second - domains_[ids[a]].first <
                           domains_[ids[b]].second - domains_[ids[b]].first;
                });
                for (int i : order)
                    if (owner[std::size_t(x_[ids[i]] - lower)] < 0) {
                        assignment[i] = int(x_[ids[i]] - lower);
                        owner[assignment[i]] = i;
                    }
                int stamp = 0;
                // Iterative augmenting paths avoid recursion on large groups.
                for (int start : order)
                    if (assignment[start] < 0) {
                        ++stamp;
                        std::vector<int> queue{start}, parent(n, -1);
                        bool found = false;
                        for (std::size_t q = 0; q < queue.size() && !found; ++q) {
                            int i = queue[q], v = ids[i];
                            auto d = domains_[v];
                            int width = int(d.second - d.first + 1), offset = random_index(width);
                            for (int k = 0; k < width; ++k) {
                                tick();
                                int value = int(d.first - lower) + (offset + k) % width;
                                if (!allowed(v, lower + value)) continue;
                                if (seen[value] == stamp) continue;
                                seen[value] = stamp;
                                parent[value] = i;
                                if (owner[value] >= 0) {
                                    queue.push_back(owner[value]);
                                    continue;
                                }
                                while (value >= 0) {
                                    int a = parent[value], previous = assignment[a];
                                    assignment[a] = value;
                                    owner[value] = a;
                                    value = previous;
                                }
                                found = true;
                                break;
                            }
                        }
                        if (!found) return;
                    }
                for (int i = 0; i < n; ++i) x_[ids[i]] = lower + assignment[i];
            }
            for (int v : ids) permutation_[v] = id;
            permutation_constraint_[id] = true;
        }
        if (c.kind == detail::Kind::Circuit && options_.structural_initialization) {
            const int n = int(c.successors.size());
            for (int v : c.successors)
                if (hinted(v) || model_.variables_[v].definition >= 0) return;
            std::vector<int> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::shuffle(order.begin(), order.end(), rng_);
            for (int i = 0; i < n; ++i) {
                tick();
                int v = c.successors[order[i]], next = order[(i + 1) % n];
                if (allowed(v, next)) x_[v] = next;
            }
        }
    }
    void cover_initialize() {
        if (!model_.definitions_.empty() || model_.direction_ != 1) return;
        const int n = model_.num_variables(), m = model_.num_constraints();
        if (!m) return;
        for (int v = 0; v < n; ++v) {
            tick();
            if (model_.variables_[v].lower < 0 || model_.variables_[v].upper > 1) return;
        }
        for (const auto &c : model_.constraints_) {
            tick();
            if (c.kind != detail::Kind::Linear || !c.guards.empty() || c.lower != 1 ||
                c.upper != bound)
                return;
            for (auto t : c.terms) {
                tick();
                if (t.second != 1) return;
            }
        }
        cover_mode_ = options_.cover_moves;
        if (!options_.structural_initialization || !options_.initial_solution.empty()) return;
        std::vector<Int> values(n, 0);
        std::vector<int> covered(m, 0), gain(n, 0);
        for (int v = 0; v < n; ++v) {
            tick();
            if (model_.variables_[v].lower == 1 ||
                (model_.variables_[v].upper == 1 && model_.objective_[v] <= 0))
                values[v] = 1;
            if (values[v])
                for (auto ref : model_.incidence_[v]) {
                    tick();
                    ++covered[ref.constraint];
                }
        }
        int missing = 0;
        for (int c = 0; c < m; ++c)
            if (!covered[c]) {
                ++missing;
                for (auto t : model_.constraints_[c].terms) {
                    tick();
                    ++gain[t.first];
                }
            }
        using Entry = std::pair<long double, int>;
        std::priority_queue<Entry> heap;
        auto price = [&](int v) {
            return static_cast<long double>(gain[v]) / std::max<Int>(1, model_.objective_[v]);
        };
        for (int v = 0; v < n; ++v)
            if (!values[v] && model_.variables_[v].upper == 1) heap.emplace(price(v), v);
        while (missing && !heap.empty()) {
            tick();
            auto entry = heap.top();
            heap.pop();
            int v = entry.second;
            if (values[v] || !gain[v]) continue;
            if (entry.first != price(v)) {
                heap.emplace(price(v), v);
                continue;
            }
            values[v] = 1;
            for (auto ref : model_.incidence_[v])
                if (covered[ref.constraint]++ == 0) {
                    tick();
                    --missing;
                    for (auto t : model_.constraints_[ref.constraint].terms) {
                        tick();
                        --gain[t.first];
                    }
                }
        }
        if (missing) return;
        std::vector<int> selected;
        for (int v = 0; v < n; ++v)
            if (values[v] && model_.variables_[v].lower == 0 && model_.objective_[v] > 0)
                selected.push_back(v);
        std::sort(selected.begin(), selected.end(), [&](int a, int b) {
            tick();
            return model_.objective_[a] > model_.objective_[b];
        });
        for (int v : selected) {
            tick();
            bool redundant = true;
            for (auto ref : model_.incidence_[v])
                if (covered[ref.constraint] <= 1) {
                    redundant = false;
                    break;
                }
            if (redundant) {
                values[v] = 0;
                for (auto ref : model_.incidence_[v]) {
                    tick();
                    --covered[ref.constraint];
                }
            }
        }
        x_ = std::move(values);
    }
    void packing_initialize() {
        if (!options_.packing_moves || !model_.definitions_.empty()) return;
        const int n = model_.num_variables(), m = model_.num_constraints();
        if (!n || !m || m > 64) return;
        for (int v = 0; v < n; ++v) {
            tick();
            auto d = domains_[v];
            if (d.first < 0 || d.second > 1 || Wide(model_.direction_) * model_.objective_[v] > 0)
                return;
        }
        for (const auto &c : model_.constraints_) {
            tick();
            if (c.kind != detail::Kind::Linear || !c.guards.empty() || c.lower > 0 ||
                c.upper < 0 || c.upper == bound)
                return;
            for (auto t : c.terms) {
                tick();
                if (t.second < 0) return;
            }
        }
        Packing p;
        p.weight.assign(n, std::vector<Int>(m, 0));
        p.profit.resize(n);
        p.capacity.resize(m);
        p.price.resize(m);
        Wide sum = 0;
        std::vector<Int> values(n), remaining(m);
        for (int v = 0; v < n; ++v) {
            tick();
            values[v] = domains_[v].first;
            p.profit[v] = -model_.direction_ * model_.objective_[v];
            sum += p.profit[v];
        }
        for (int c = 0; c < m; ++c) {
            p.capacity[c] = remaining[c] = model_.constraints_[c].upper;
            for (auto t : model_.constraints_[c].terms) {
                tick();
                p.weight[t.first][c] = t.second;
                remaining[c] -= t.second * values[t.first];
            }
            if (remaining[c] < 0) return;
            p.price[c] = static_cast<long double>(sum) / (2.0L * m * std::max<Int>(1, remaining[c]));
        }
        std::vector<std::pair<long double, int>> order;
        p.reduced.resize(n);
        for (int v = 0; v < n; ++v) {
            p.reduced[v] = p.profit[v];
            for (int c = 0; c < m; ++c) {
                tick();
                p.reduced[v] -= p.price[c] * p.weight[v][c];
            }
        }
        for (int v : decisions_) {
            long double cost = 0;
            for (int c = 0; c < m; ++c) {
                tick();
                cost += p.price[c] * p.weight[v][c];
            }
            order.emplace_back(-p.profit[v] / std::max(1e-30L, cost), v);
        }
        std::sort(order.begin(), order.end(), [&](auto a, auto b) {
            tick();
            return a < b;
        });
        for (auto entry : order) {
            int v = entry.second;
            bool fits = p.profit[v] > 0;
            for (int c = 0; c < m; ++c) {
                tick();
                if (p.weight[v][c] > remaining[c]) fits = false;
            }
            if (fits) {
                values[v] = 1;
                for (int c = 0; c < m; ++c) remaining[c] -= p.weight[v][c];
            }
        }
        // Keep a complete feasible warm start if it is better than construction.
        bool valid = true;
        Wide previous = 0, constructed = 0;
        for (int c = 0; c < m; ++c) {
            Wide load = 0;
            for (auto t : model_.constraints_[c].terms) {
                tick();
                load += Wide(t.second) * x_[t.first];
            }
            if (load > p.capacity[c]) valid = false;
        }
        for (int v = 0; v < n; ++v) {
            tick();
            previous += Wide(p.profit[v]) * x_[v];
            constructed += Wide(p.profit[v]) * values[v];
        }
        if (!valid || constructed > previous) x_ = std::move(values);
        packing_ = std::move(p);
        packing_mode_ = true;
    }
    void clause_initialize() {
        if (!options_.clause_moves || !options_.compound_moves || !model_.definitions_.empty() ||
            model_.direction_ != 1)
            return;
        const int n = model_.num_variables(), m = model_.num_constraints();
        if (!m) return;
        for (int v = 0; v < n; ++v) {
            tick();
            if (domains_[v].first < 0 || domains_[v].second > 1 || model_.objective_[v] < 0)
                return;
        }
        std::vector<int> slack(m, -1);
        std::vector<bool> auxiliary(n, false);
        auto values = x_;
        for (int c = 0; c < m; ++c) {
            const auto &row = model_.constraints_[c];
            if (row.kind != detail::Kind::Linear || !row.guards.empty() || row.upper != bound)
                return;
            int negative = 0;
            Wide activity = 0;
            for (auto t : row.terms) {
                tick();
                if (t.second != 1 && t.second != -1) return;
                if (model_.objective_[t.first] > 0) {
                    if (slack[c] >= 0 || t.second != 1 ||
                        domains_[t.first] != std::make_pair(Int(0), Int(1)) ||
                        model_.incidence_[t.first].size() != 1)
                        return;
                    slack[c] = t.first;
                } else {
                    activity += t.second * values[t.first];
                    negative += t.second < 0;
                }
            }
            if (slack[c] < 0 || row.lower != 1 - negative) return;
            auxiliary[slack[c]] = true;
            values[slack[c]] = activity < row.lower;
        }
        for (int v = 0; v < n; ++v)
            if (model_.objective_[v] && !auxiliary[v]) return;
        std::vector<Int> base_weight(m);
        Int divisor = 0;
        for (int c = 0; c < m; ++c) {
            tick();
            base_weight[c] = model_.objective_[slack[c]];
            if (divisor != 1) divisor = std::gcd(divisor, base_weight[c]);
        }
        if (divisor > 1)
            for (auto &weight : base_weight) {
                tick();
                weight /= divisor;
            }
        // Integer scaling leaves every source objective comparison unchanged
        // and makes search independent of a common multiplier on all weights.
        clauses_.base_weight = std::move(base_weight);
        decisions_.erase(std::remove_if(decisions_.begin(), decisions_.end(),
                                        [&](int v) { return auxiliary[v]; }),
                         decisions_.end());
        clause_slack_ = std::move(slack);
        x_ = std::move(values);
        clause_mode_ = true;
    }
    void reduce_different() {
        if (!options_.alldifferent_reduction) return;
        const int n = model_.num_variables(), m = model_.num_constraints();
        int fixed_variables = 0;
        for (int v = 0; v < n; ++v) {
            tick();
            fixed_variables += domains_[v].first == domains_[v].second;
        }
        // In lightly restricted permutations, unrestricted row exchanges are
        // cheaper and connect the search space without alternating-cycle repair.
        if (2 * fixed_variables < n) return;
        std::vector<std::uint64_t> masks(n, 0);
        for (int v = 0; v < n; ++v) {
            tick();
            const auto d = model_.variables_[v];
            if (d.definition >= 0 || d.upper - d.lower >= 64) continue;
            int low = int(domains_[v].first - d.lower), high = int(domains_[v].second - d.lower);
            masks[v] = (~std::uint64_t(0) << low) & (~std::uint64_t(0) >> (63 - high));
        }
        std::vector<int> queue;
        std::vector<bool> queued(m, false), supported(m, false);
        for (int c = 0; c < m; ++c) {
            tick();
            const auto &row = model_.constraints_[c];
            if (row.kind != detail::Kind::Different || !row.guards.empty()) continue;
            bool ok = true;
            for (auto a : row.expressions) {
                tick();
                if (a.var < 0 || a.coefficient != 1 || a.offset || !masks[a.var]) ok = false;
            }
            if (ok) {
                supported[c] = queued[c] = true;
                queue.push_back(c);
            }
        }
        auto restrict = [&](int v, std::uint64_t mask) {
            if (!mask || mask == masks[v]) return;
            masks[v] = mask;
            for (auto ref : model_.incidence_[v]) {
                tick();
                int c = ref.constraint;
                if (supported[c] && !queued[c]) {
                    queued[c] = true;
                    queue.push_back(c);
                }
            }
        };
        for (std::size_t head = 0; head < queue.size(); ++head) {
            tick();
            int c = queue[head];
            queued[c] = false;
            const auto &row = model_.constraints_[c];
            std::vector<Int> fixed;
            for (auto a : row.expressions) {
                tick();
                auto mask = masks[a.var];
                if (!(mask & (mask - 1)))
                    fixed.push_back(model_.variables_[a.var].lower + __builtin_ctzll(mask));
            }
            for (auto a : row.expressions) {
                int v = a.var;
                auto mask = masks[v];
                if (!(mask & (mask - 1))) continue;
                for (Int value : fixed) {
                    tick();
                    Wide bit = Wide(value) - model_.variables_[v].lower;
                    if (bit >= 0 && bit < 64) mask &= ~(std::uint64_t(1) << int(bit));
                }
                restrict(v, mask);
            }
            if (row.expressions.empty() || row.expressions.size() > 64) continue;
            const auto first = model_.variables_[row.expressions.front().var];
            bool permutation = first.upper - first.lower + 1 == Int(row.expressions.size());
            std::uint64_t once = 0, multiple = 0;
            for (auto a : row.expressions) {
                tick();
                const auto d = model_.variables_[a.var];
                if (d.lower != first.lower || d.upper != first.upper) permutation = false;
                multiple |= once & masks[a.var];
                once |= masks[a.var];
            }
            if (!permutation) continue;
            const auto unique = once & ~multiple;
            for (auto a : row.expressions) {
                tick();
                auto mask = masks[a.var] & unique;
                // Two uniquely supported values in one variable imply failure;
                // retain the original row instead of manufacturing a singleton.
                if (mask && !(mask & (mask - 1))) restrict(a.var, mask);
            }
        }
        for (int v = 0; v < n; ++v)
            if (masks[v]) {
                // No clock throw between publishing domains and their masks.
                domains_[v] = {model_.variables_[v].lower + __builtin_ctzll(masks[v]),
                               model_.variables_[v].lower + 63 - __builtin_clzll(masks[v])};
            }
        domain_masks_ = std::move(masks);
    }
    void prepare() {
        const int n = model_.num_variables(), m = model_.num_constraints();
        if (!allocated_) {
            check();
            x_.resize(n);
            old_.resize(n);
            undo_.resize(n);
            domains_.resize(n);
            domain_masks_.assign(n, 0);
            changed_stamp_.assign(n, 0);
            queued_stamp_.assign(n, 0);
            tabu_.assign(n, 0);
            permutation_.assign(n, -1);
            permutation_constraint_.assign(m, false);
            constraint_stamp_.assign(m, 0);
            position_.assign(m, -1);
            cache_.resize(m);
            slots_.resize(m);
            next_raw_.resize(m);
            next_value_.resize(m);
            next_activity_.resize(m);
            weights_.assign(m, 1);
            changed_.reserve(n);
            touched_.reserve(m);
            violated_.reserve(m);
            result_.x.reserve(n);
            allocated_ = true;
        }
        auto get = [&](int v) { return x_[v]; };
        if (phase_ == -2) {
            while (cursor_ < n) {
                tick();
                auto d = model_.variables_[cursor_];
                domains_[cursor_++] = {d.lower, d.upper};
            }
            phase_ = -1;
            cursor_ = 0;
        }
        if (phase_ == -1) {
            while (cursor_ < m) {
                tick();
                const auto &c = model_.constraints_[cursor_++];
                if (!options_.domain_reduction || c.kind != detail::Kind::Linear ||
                    !c.guards.empty() || c.terms.size() != 1)
                    continue;
                auto t = c.terms[0];
                if (model_.variables_[t.first].definition >= 0) continue;
                Wide l = c.lower, u = c.upper, a = t.second;
                if (a < 0) {
                    a = -a;
                    Wide tmp = l;
                    l = -u;
                    u = -tmp;
                }
                auto &d = domains_[t.first];
                l = std::max(Wide(d.first), -detail::floor_div(-l, a));
                u = std::min(Wide(d.second), detail::floor_div(u, a));
                // An inconsistent row remains in the model; never invent a domain.
                if (l <= u) d = {Int(l), Int(u)};
            }
            reduce_different();
            phase_ = 0;
            cursor_ = 0;
        }
        if (phase_ == 0) {
            while (cursor_ < n) {
                tick();
                const int v = cursor_++;
                const auto d = model_.variables_[v];
                if (hinted(v) && d.definition < 0) {
                    x_[v] = options_.initial_solution[v];
                    if (x_[v] < d.lower || x_[v] > d.upper)
                        throw std::invalid_argument("warm start outside variable domain");
                    x_[v] = std::max(domains_[v].first, std::min(domains_[v].second, x_[v]));
                    if (!allowed(v, x_[v])) x_[v] = random_value(v);
                } else {
                    x_[v] = domains_[v].first;
                    if (d.definition < 0 && model_.objective_[v] * model_.direction_ < 0)
                        x_[v] = domains_[v].second;
                    else if (d.definition < 0 && d.upper - d.lower > 1)
                        x_[v] = random_value(v);
                }
                if (d.definition < 0 && domains_[v].first < domains_[v].second)
                    decisions_.push_back(v);
                if (model_.objective_[v]) {
                    has_objective_ = true;
                    objective_vars_.push_back(v);
                }
            }
            phase_ = 1;
            cursor_ = 0;
        }
        if (phase_ == 1) {
            while (cursor_ < m) {
                check();
                structural_initialize(model_.constraints_[cursor_], cursor_);
                ++cursor_;
            }
            cover_initialize();
            packing_initialize();
            clause_initialize();
            phase_ = 2;
            cursor_ = 0;
        }
        if (phase_ == 2) {
            while (cursor_ < int(model_.definitions_.size())) {
                check();
                const int v = model_.definitions_[cursor_].output;
                x_[v] = model_.defined_value(v, get, [&] { tick(); });
                ++cursor_;
            }
            phase_ = 3;
            cursor_ = 0;
        }
        if (phase_ == 3) {
            while (cursor_ < m) {
                check();
                const auto &c = model_.constraints_[cursor_];
                Cache a;
                if (c.kind == detail::Kind::Different) {
                    for (auto e : c.expressions) {
                        tick();
                        Int value = detail::value(e, get);
                        a.raw += a.frequency[value]++;
                        if (permutation_constraint_[cursor_]) a.owner[value] = e.var;
                    }
                } else
                    a.raw = detail::raw_violation(c, get, [&] { tick(); });
                if (c.kind == detail::Kind::Linear)
                    for (auto t : c.terms) {
                        tick();
                        a.activity += Wide(t.second) * x_[t.first];
                    }
                a.value = detail::active(c, get) ? a.raw : 0;
                raw_total_ += a.value;
                set_violation(cursor_, a.value != 0);
                cache_[cursor_] = std::move(a);
                ++cursor_;
            }
            Wide obj = model_.offset_, optimistic = model_.offset_;
            long double sum = 0;
            for (int i = 0; i < n; ++i) {
                obj += Wide(model_.objective_[i]) * x_[i];
                optimistic +=
                    Wide(model_.objective_[i]) * (model_.objective_[i] * model_.direction_ >= 0
                                                      ? model_.variables_[i].lower
                                                      : model_.variables_[i].upper);
                if (model_.objective_[i]) {
                    long double s = std::abs(static_cast<long double>(model_.objective_[i]));
                    if (model_.variables_[i].definition >= 0)
                        s *= std::max(1.0L, static_cast<long double>(model_.variables_[i].upper -
                                                                     model_.variables_[i].lower) /
                                                10);
                    sum += s;
                }
            }
            objective_ = Int(obj);
            objective_scale_ =
                std::max(1.0L, sum / std::max<std::size_t>(1, objective_vars_.size()));
            objective_bound_ = Int(optimistic);
            if (m == 1 && model_.constraints_[0].kind == detail::Kind::Circuit) {
                std::vector<bool> successor(n, false);
                for (int v : model_.constraints_[0].successors) successor[v] = true;
                bool only_route = true;
                for (int v : decisions_)
                    if (!successor[v]) only_route = false;
                if (only_route) route_constraint_ = 0;
            }
            permutation_mode_ = options_.compound_moves && options_.permutation_moves &&
                                !has_objective_ && model_.definitions_.empty();
            for (int v : decisions_)
                if (permutation_[v] < 0) permutation_mode_ = false;
            for (const auto &c : model_.constraints_) {
                tick();
                if (!c.guards.empty() || (c.kind != detail::Kind::Different &&
                                          !(c.kind == detail::Kind::Linear && c.terms.size() == 1)))
                    permutation_mode_ = false;
                for (auto a : c.expressions) {
                    tick();
                    if (a.var < 0 || a.coefficient != 1 || a.offset) permutation_mode_ = false;
                }
            }
            std::uint64_t remaining_values = 0, original_values = 0;
            bool masked_decisions = true;
            for (int v : decisions_) {
                tick();
                remaining_values += domain_masks_[v] ? __builtin_popcountll(domain_masks_[v])
                                                     : domains_[v].second - domains_[v].first + 1;
                if (!domain_masks_[v]) masked_decisions = false;
                original_values += model_.variables_[v].upper - model_.variables_[v].lower + 1;
            }
            recolor_mode_ = permutation_mode_ && masked_decisions && original_values &&
                            remaining_values < original_values * 0.6L;
            phase_ = 4;
            remember();
        }
    }
    void rollback() {
        for (int v : changed_) x_[v] = old_[v];
        changed_.clear();
    }
    void set_temporary(int v, Int value) {
        if (x_[v] == value) return;
        if (changed_stamp_[v] != epoch_) {
            changed_stamp_[v] = epoch_;
            old_[v] = x_[v];
            changed_.push_back(v);
        }
        x_[v] = value;
    }
    void evaluate_move(const Move &move) {
        ++epoch_;
        changed_.clear();
        touched_.clear();
        std::priority_queue<int, std::vector<int>, std::greater<int>> queue;
        auto enqueue = [&](int v) {
            for (int a : model_.dependents_[v]) {
                tick();
                if (queued_stamp_[a] != epoch_) {
                    queued_stamp_[a] = epoch_;
                    queue.push(a);
                }
            }
        };
        for (auto a : move) {
            set_temporary(a.var, a.value);
            enqueue(a.var);
        }
        auto get = [&](int v) { return x_[v]; };
        auto oldget = [&](int v) { return old_value(v); };
        while (!queue.empty()) {
            tick();
            int v = queue.top();
            queue.pop();
            Int value = model_.defined_value(v, get, [&] { tick(); });
            if (value != x_[v]) {
                set_temporary(v, value);
                enqueue(v);
            }
        }
        Wide obj = objective_;
        for (int v : changed_) {
            tick();
            obj += Wide(model_.objective_[v]) * (Wide(x_[v]) - old_[v]);
            for (auto ref : model_.incidence_[v]) {
                tick();
                const int c = ref.constraint;
                if (constraint_stamp_[c] != epoch_) {
                    constraint_stamp_[c] = epoch_;
                    touched_.push_back(c);
                    slots_[c].clear();
                }
                if (ref.slot >= 0) slots_[c].push_back(ref.slot);
            }
        }
        next_objective_ = Int(obj);
        for (int id : touched_) {
            tick();
            const auto &c = model_.constraints_[id];
            auto &slots = slots_[id];
            std::sort(slots.begin(), slots.end());
            slots.erase(std::unique(slots.begin(), slots.end()), slots.end());
            Wide raw = cache_[id].raw;
            if (c.kind == detail::Kind::Linear) {
                Wide a = cache_[id].activity;
                for (int s : slots) {
                    tick();
                    auto t = c.terms[s];
                    a += Wide(t.second) * (Wide(x_[t.first]) - old_value(t.first));
                }
                next_activity_[id] = a;
                raw = std::max({Wide(0), Wide(c.lower) - a, a - c.upper});
            } else if (c.kind == detail::Kind::Different) {
                std::unordered_map<Int, int> adjustment;
                auto count = [&](Int k) {
                    auto it = cache_[id].frequency.find(k);
                    return (it == cache_[id].frequency.end() ? 0 : it->second) + adjustment[k];
                };
                for (int s : slots) {
                    tick();
                    Int a = detail::value(c.expressions[s], oldget);
                    raw -= count(a) - 1;
                    --adjustment[a];
                }
                for (int s : slots) {
                    tick();
                    Int a = detail::value(c.expressions[s], get);
                    raw += count(a);
                    ++adjustment[a];
                }
            } else if (c.kind == detail::Kind::Overlap) {
                for (int i : slots) {
                    std::size_t p = 0;
                    for (int j = 0; j < int(c.intervals.size()); ++j) {
                        tick();
                        while (p < slots.size() && slots[p] < j) ++p;
                        if (i == j || (j < i && p < slots.size() && slots[p] == j)) continue;
                        raw += int(detail::collides(c.intervals[i], c.intervals[j], get)) -
                               int(detail::collides(c.intervals[i], c.intervals[j], oldget));
                    }
                }
            } else
                raw = detail::raw_violation(c, get, [&] { tick(); });
            next_raw_[id] = raw;
            next_value_[id] = detail::active(c, get) ? raw : 0;
        }
        ++stats_.evaluations;
    }
    void finish_commit() {
        auto get = [&](int v) { return x_[v]; };
        auto oldget = [&](int v) { return old_value(v); };
        while (commit_cursor_ < touched_.size()) {
            check();
            const int id = touched_[commit_cursor_];
            const auto &c = model_.constraints_[id];
            auto &a = cache_[id];
            if (c.kind == detail::Kind::Different) {
                const bool rebuild_owner =
                    permutation_constraint_[id] && a.raw > 0 && next_raw_[id] == 0;
                std::unordered_map<Int, int> owner;
                if (rebuild_owner)
                    for (auto expression : c.expressions) {
                        // Before changing this row, rebuild privately: duplicate
                        // old values may share an owner with an unchanged cell.
                        tick();
                        owner[detail::value(expression, get)] = expression.var;
                    }
                // This short per-constraint transaction cannot be interrupted.
                for (int s : slots_[id]) {
                    Int key = detail::value(c.expressions[s], oldget);
                    auto it = a.frequency.find(key);
                    if (--it->second == 0) a.frequency.erase(it);
                }
                for (int s : slots_[id]) ++a.frequency[detail::value(c.expressions[s], get)];
                if (rebuild_owner)
                    a.owner = std::move(owner);
                else if (permutation_constraint_[id]) {
                    for (int s : slots_[id])
                        a.owner.erase(detail::value(c.expressions[s], oldget));
                    for (int s : slots_[id])
                        a.owner[detail::value(c.expressions[s], get)] = c.expressions[s].var;
                }
            }
            raw_total_ += next_value_[id] - a.value;
            a.raw = next_raw_[id];
            a.value = next_value_[id];
            if (c.kind == detail::Kind::Linear) a.activity = next_activity_[id];
            set_violation(id, a.value != 0);
            ++commit_cursor_;
        }
        while (!permutation_owner_reset_.empty()) {
            check();
            auto &entry = permutation_owner_reset_.back();
            cache_[entry.first].owner = std::move(entry.second);
            permutation_owner_reset_.pop_back();
        }
        for (int v : changed_)
            if (model_.variables_[v].definition < 0) {
                undo_[v] = old_[v];
                tabu_[v] = stats_.iterations + 3 + rng_() % 7;
            }
        objective_ = next_objective_;
        pending_ = false;
        ++stats_.moves;
        ++stall_;
        remember();
        changed_.clear();
    }
    long double objective_merit(Int obj) const {
        if (!has_objective_ || !result_.feasible) return 0;
        const Wide excess = Wide(model_.direction_) * (Wide(obj) - result_.objective) + 1;
        return static_cast<long double>(std::max(Wide(0), excess)) * objective_weight_ /
               objective_scale_;
    }
    long double move_score() const {
        long double score = objective_merit(next_objective_) - objective_merit(objective_);
        for (int id : touched_)
            score += weights_[id] * static_cast<long double>(next_value_[id] - cache_[id].value) /
                     model_.constraints_[id].scale;
        return score;
    }
    bool add_change(Move &move, int v, Wide value) {
        if (v < 0 || model_.variables_[v].definition >= 0) return false;
        if (!allowed(v, value)) return false;
        for (auto &a : move)
            if (a.var == v) {
                a.value = Int(value);
                return true;
            }
        if (value != x_[v]) move.push_back({v, Int(value)});
        return true;
    }
    bool assign_affine(Move &move, Affine a, Wide value) {
        if (a.var < 0 || !a.coefficient) return false;
        Wide v = detail::floor_div(value - a.offset, a.coefficient);
        if (rng_() % 2 && (value - a.offset) % a.coefficient) ++v;
        return add_change(move, a.var, v);
    }
    int root_variable(int v) {
        while (model_.variables_[v].definition >= 0) {
            tick();
            const auto &d = model_.definitions_[model_.variables_[v].definition];
            if (d.kind == detail::Definition::Element)
                v = d.index;
            else {
                std::vector<int> candidates;
                for (auto a : d.args) {
                    tick();
                    if (a.var >= 0 && (rng_() % 3 == 0 ||
                                       detail::value(a, [&](int i) { return x_[i]; }) == x_[v]))
                        candidates.push_back(a.var);
                }
                if (candidates.empty())
                    for (auto a : d.args)
                        if (a.var >= 0) candidates.push_back(a.var);
                if (candidates.empty()) return -1;
                v = candidates[random_index(candidates.size())];
            }
        }
        return v;
    }
    void preserve_permutations(Move &move) {
        if (!options_.compound_moves) return;
        Move out;
        auto current = [&](int v) {
            for (auto a : out)
                if (a.var == v) return a.value;
            return x_[v];
        };
        auto put = [&](int v, Int a) {
            for (auto &b : out)
                if (b.var == v) {
                    b.value = a;
                    return;
                }
            out.push_back({v, a});
        };
        for (auto a : move) {
            tick();
            const int id = permutation_[a.var];
            if (id < 0 || cache_[id].raw != 0) {
                put(a.var, a.value);
                continue;
            }
            int owner = -1;
            for (auto b : out)
                if (permutation_[b.var] == id && b.value == a.value) owner = b.var;
            if (owner < 0) {
                auto it = cache_[id].owner.find(a.value);
                if (it != cache_[id].owner.end()) owner = it->second;
            }
            if (owner < 0) continue;
            Int old = current(a.var);
            put(a.var, a.value);
            put(owner, old);
        }
        out.erase(
            std::remove_if(out.begin(), out.end(), [&](Change a) { return a.value == x_[a.var]; }),
            out.end());
        for (auto a : out)
            if (!allowed(a.var, a.value)) {
                move.clear();
                return;
            }
        move = std::move(out);
    }
    Move propose(int id) {
        Move move;
        auto get = [&](int v) { return x_[v]; };
        if (id < 0) {
            int v = -1;
            if (!objective_vars_.empty() && rng_() % 4)
                v = root_variable(objective_vars_[random_index(objective_vars_.size())]);
            if (v < 0 && !decisions_.empty()) v = decisions_[random_index(decisions_.size())];
            if (v < 0) return move;
            const auto d = model_.variables_[v];
            const int mode = int(rng_() % 6);
            Int value = mode == 0   ? d.lower
                        : mode == 1 ? d.upper
                        : mode == 2 ? std::max(d.lower, x_[v] - 1)
                        : mode == 3 ? std::min(d.upper, x_[v] + 1)
                                    : random_value(v);
            add_change(move, v, value);
            return move;
        }
        const auto &c = model_.constraints_[id];
        if (!c.guards.empty() && rng_() % 8 == 0) {
            auto a = c.guards[random_index(c.guards.size())];
            add_change(move, a.var, a.positive ? 0 : 1);
            return move;
        }
        if (c.kind == detail::Kind::Linear && !c.terms.empty()) {
            auto t = c.terms[random_index(c.terms.size())];
            int v = root_variable(t.first);
            if (v < 0) return move;
            if (v != t.first) {
                add_change(move, v, random_value(v));
                return move;
            }
            Wide target = cache_[id].activity < c.lower   ? c.lower
                          : cache_[id].activity > c.upper ? c.upper
                                                          : (rng_() % 2 ? c.lower : c.upper);
            Wide value = Wide(x_[v]) + detail::floor_div(target - cache_[id].activity, t.second);
            if (rng_() % 2 && (target - cache_[id].activity) % t.second) ++value;
            value = std::max(Wide(model_.variables_[v].lower),
                             std::min(Wide(model_.variables_[v].upper), value));
            if (value == x_[v]) value = random_value(v);
            add_change(move, v, value);
            if (options_.compound_moves && c.terms.size() > 1 && rng_() % 2 == 0) {
                auto u = c.terms[random_index(c.terms.size())];
                Wide difference = t.second * (value - x_[v]);
                if (u.first != v && difference % u.second == 0)
                    add_change(move, u.first, Wide(x_[u.first]) - difference / u.second);
            }
        } else if (c.kind == detail::Kind::Different && !c.expressions.empty()) {
            int i = random_index(c.expressions.size());
            if (cache_[id].value) {
                auto &cache = cache_[id];
                if (cache.conflict_version != stats_.moves) {
                    cache.conflict_slots.clear();
                    for (int k = 0; k < int(c.expressions.size()); ++k) {
                        tick();
                        auto it = cache.frequency.find(detail::value(c.expressions[k], get));
                        if (it != cache.frequency.end() && it->second > 1)
                            cache.conflict_slots.push_back(k);
                    }
                    cache.conflict_version = stats_.moves;
                }
                if (!cache.conflict_slots.empty())
                    i = cache.conflict_slots[random_index(cache.conflict_slots.size())];
            }
            auto a = c.expressions[i];
            if (a.var < 0) return move;
            if (options_.compound_moves && c.expressions.size() > 1 && rng_() % 4 != 0) {
                int j = random_index(c.expressions.size());
                auto b = c.expressions[j];
                if (a.var != b.var && a.var >= 0 && b.var >= 0) {
                    // For plain permutations, this preserves distinctness exactly.
                    if (permutation_[a.var] >= 0 && permutation_[a.var] == permutation_[b.var]) {
                        if (!add_change(move, a.var, x_[b.var]) ||
                            !add_change(move, b.var, x_[a.var]))
                            move.clear();
                    } else if (!assign_affine(move, a, detail::value(b, get)) ||
                               !assign_affine(move, b, detail::value(a, get)))
                        move.clear();
                }
            } else
                add_change(move, a.var, random_value(a.var));
        } else if (c.kind == detail::Kind::Overlap && !c.intervals.empty()) {
            const int n = int(c.intervals.size());
            int i = random_index(n), j = random_index(n);
            if (cache_[id].value) {
                bool found = false;
                for (int k = 0; k < 24; ++k) {
                    tick();
                    i = random_index(n);
                    j = random_index(n);
                    if (i != j && detail::collides(c.intervals[i], c.intervals[j], get)) {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    for (int k = 0; k < n; ++k) {
                        tick();
                        j = (j + 1) % n;
                        if (i != j && detail::collides(c.intervals[i], c.intervals[j], get)) break;
                    }
            }
            if (i == j) return move;
            const auto &a = c.intervals[i];
            const auto &b = c.intervals[j];
            if (a.present.var >= 0 && rng_() % 12 == 0)
                add_change(move, a.present.var, 1 - x_[a.present.var]);
            else {
                Wide target =
                    rng_() % 2 ? Wide(detail::value(b.start, get)) + detail::value(b.duration, get)
                               : Wide(detail::value(b.start, get)) - detail::value(a.duration, get);
                assign_affine(move, a.start, target);
            }
        } else if (c.kind == detail::Kind::Circuit) {
            const int n = int(c.successors.size());
            if (n < 2) return move;
            int i = random_index(n), j = random_index(n);
            if (i == j) return move;
            int a = int(x_[c.successors[i]]), b = int(x_[c.successors[j]]);
            if (cache_[id].raw == 0 && options_.compound_moves) {
                if (c.self_loops && a == i && b != j) {
                    if (!add_change(move, c.successors[i], b) ||
                        !add_change(move, c.successors[j], i))
                        move.clear();
                    return move;
                }
                if (c.self_loops && a != i && rng_() % 8 == 0) {
                    int predecessor = -1;
                    for (int k = 0; k < n; ++k) {
                        tick();
                        if (x_[c.successors[k]] == i) predecessor = k;
                    }
                    if (predecessor >= 0 && (!add_change(move, c.successors[predecessor], a) ||
                                             !add_change(move, c.successors[i], i)))
                        move.clear();
                    return move;
                }
                if (a == j || a == i || b == i || b == j) return move;
                if (rng_() % 2) {
                    // Reverse the path a,...,j. The objective handles directed
                    // costs through its element invariants, without symmetry assumptions.
                    std::vector<int> path;
                    int v = a;
                    while (v != i && path.size() < 256) {
                        tick();
                        path.push_back(v);
                        if (v == j) break;
                        v = int(x_[c.successors[v]]);
                    }
                    if (!path.empty() && path.back() == j) {
                        bool ok = add_change(move, c.successors[i], j) &&
                                  add_change(move, c.successors[a], b);
                        for (std::size_t k = 1; k < path.size(); ++k)
                            ok = add_change(move, c.successors[path[k]], path[k - 1]) && ok;
                        if (!ok) move.clear();
                    }
                    return move;
                }
                int after = int(x_[c.successors[a]]);
                // Remove a from behind i and insert it behind j (three arcs).
                if (!add_change(move, c.successors[i], after) ||
                    !add_change(move, c.successors[j], a) || !add_change(move, c.successors[a], b))
                    move.clear();
            } else if (options_.compound_moves && rng_() % 2) {
                if (!add_change(move, c.successors[i], b) || !add_change(move, c.successors[j], a))
                    move.clear();
            } else
                add_change(move, c.successors[i], j);
        } else if (c.kind == detail::Kind::Xor && !c.literals.empty()) {
            auto a = c.literals[random_index(c.literals.size())];
            if (a.var >= 0) add_change(move, a.var, 1 - x_[a.var]);
        }
        return move;
    }
    void bump_weights() {
        if (!options_.adaptive_weights) return;
        ++stats_.weight_updates;
        const std::size_t n = std::min<std::size_t>(64, violated_.size());
        for (std::size_t k = 0; k < n; ++k) {
            tick();
            int c =
                violated_[n == violated_.size() ? k : std::size_t(random_index(violated_.size()))];
            weights_[c] = std::min(100000.0L, weights_[c] + 1);
        }
        if (result_.feasible) objective_weight_ = std::min(100000.0L, objective_weight_ + 0.2L);
    }
    bool project_schedule(const Schedule &s, const std::vector<std::vector<int>> &order,
                          std::vector<Int> &values, Int &objective,
                          std::vector<std::pair<int, int>> &critical, bool collect_critical = true) {
        const int n = int(s.vars.size());
        std::vector<int> next(n, -1), indegree = s.indegree, queue, parent(n, -1), place(n);
        std::vector<Int> earliest(n);
        for (int i = 0; i < n; ++i) {
            tick();
            earliest[i] = model_.variables_[s.vars[i]].lower;
        }
        for (const auto &machine : order)
            for (std::size_t i = 0; i < machine.size(); ++i) {
                tick();
                place[machine[i]] = int(i);
                if (i) {
                    next[machine[i - 1]] = machine[i];
                    ++indegree[machine[i]];
                }
            }
        for (int i = 0; i < n; ++i)
            if (!indegree[i]) queue.push_back(i);
        auto relax = [&](int i, int j, Int lag) {
            Wide v = Wide(earliest[i]) + lag;
            if (v > bound) return false;
            if (v >= earliest[j]) {
                earliest[j] = Int(v);
                parent[j] = i;
            }
            if (--indegree[j] == 0) queue.push_back(j);
            return true;
        };
        for (std::size_t k = 0; k < queue.size(); ++k) {
            tick();
            int i = queue[k];
            if (earliest[i] > model_.variables_[s.vars[i]].upper) return false;
            for (auto e : s.edges[i]) {
                tick();
                if (!relax(i, e.first, e.second)) return false;
            }
            if (next[i] >= 0 && !relax(i, next[i], s.duration[i])) return false;
        }
        if (int(queue.size()) != n) return false;
        values = x_;
        for (int i = 0; i < n; ++i) {
            tick();
            values[s.vars[i]] = earliest[i];
        }
        for (const auto &d : model_.definitions_) {
            tick();
            values[d.output] =
                model_.defined_value(d.output, [&](int v) { return values[v]; }, [&] { tick(); });
        }
        Wide obj = model_.offset_;
        for (int i = 0; i < model_.num_variables(); ++i) {
            tick();
            obj += Wide(model_.objective_[i]) * values[i];
        }
        objective = Int(obj);
        if (!collect_critical) return true;
        critical.clear();
        int last = 0;
        for (int i = 1; i < n; ++i)
            if (Wide(earliest[i]) + s.duration[i] > Wide(earliest[last]) + s.duration[last])
                last = i;
        for (int i = last; parent[i] >= 0; i = parent[i]) {
            tick();
            int p = parent[i];
            if (next[p] == i) critical.emplace_back(s.machine[i], place[p]);
        }
        if (options_.critical_block_moves) {
            // All tight machine arcs on any longest path, not just one path.
            std::vector<Wide> tail(n);
            for (auto it = queue.rbegin(); it != queue.rend(); ++it) {
                tick();
                int i = *it;
                tail[i] = s.duration[i];
                for (auto e : s.edges[i]) {
                    tick();
                    tail[i] = std::max(tail[i], Wide(e.second) + tail[e.first]);
                }
                if (next[i] >= 0) tail[i] = std::max(tail[i], Wide(s.duration[i]) + tail[next[i]]);
            }
            const Wide makespan = Wide(earliest[last]) + s.duration[last];
            critical.clear();
            for (int i = 0; i < n; ++i) {
                tick();
                if (next[i] >= 0 && Wide(earliest[i]) + s.duration[i] + tail[next[i]] == makespan)
                    critical.emplace_back(s.machine[i], place[i]);
            }
            std::sort(critical.begin(), critical.end());
        }
        return true;
    }
    void try_schedule() {
        if (schedule_tried_) return;
        if (!options_.schedule_moves || !options_.compound_moves) {
            schedule_tried_ = true;
            return;
        }
        Schedule s;
        std::vector<int> task(model_.num_variables(), -1);
        auto unsupported = [&] { schedule_tried_ = true; };
        for (const auto &c : model_.constraints_) {
            tick();
            if (!c.guards.empty()) {
                unsupported();
                return;
            }
            if (c.kind != detail::Kind::Overlap && c.kind != detail::Kind::Linear) {
                unsupported();
                return;
            }
            if (c.kind != detail::Kind::Overlap) continue;
            s.order.emplace_back();
            for (const auto &a : c.intervals) {
                tick();
                int v = a.start.var;
                if (v < 0 || a.start.coefficient != 1 || a.start.offset || a.duration.var >= 0 ||
                    a.present.var >= 0 || !a.present.positive ||
                    model_.variables_[v].definition >= 0 || task[v] >= 0) {
                    unsupported();
                    return;
                }
                int t = int(s.vars.size());
                task[v] = t;
                s.vars.push_back(v);
                s.duration.push_back(a.duration.offset);
                s.machine.push_back(int(s.order.size()) - 1);
                s.order.back().push_back(t);
            }
        }
        if (s.vars.empty()) {
            unsupported();
            return;
        }
        for (int v = 0; v < model_.num_variables(); ++v) {
            tick();
            if ((model_.variables_[v].definition < 0 && task[v] < 0) ||
                model_.objective_[v] * model_.direction_ < 0) {
                unsupported();
                return;
            }
        }
        // Earliest scheduling is valid only for monotone maximum invariants.
        for (const auto &d : model_.definitions_) {
            if (d.kind != detail::Definition::Maximum) {
                unsupported();
                return;
            }
            for (auto a : d.args)
                if (a.coefficient < 0) {
                    unsupported();
                    return;
                }
        }
        s.edges.resize(s.vars.size());
        s.indegree.assign(s.vars.size(), 0);
        for (const auto &c : model_.constraints_)
            if (c.kind == detail::Kind::Linear) {
                tick();
                if (c.terms.size() != 2) {
                    unsupported();
                    return;
                }
                int plus = -1, minus = -1;
                for (auto t : c.terms) {
                    if (task[t.first] < 0) {
                        unsupported();
                        return;
                    }
                    if (t.second == 1)
                        plus = task[t.first];
                    else if (t.second == -1)
                        minus = task[t.first];
                }
                if (plus < 0 || minus < 0 || (c.lower != -bound && c.lower < 0) ||
                    (c.upper != bound && c.upper > 0)) {
                    unsupported();
                    return;
                }
                if (c.lower != -bound) {
                    s.edges[minus].emplace_back(plus, c.lower);
                    ++s.indegree[plus];
                }
                if (c.upper != bound) {
                    s.edges[plus].emplace_back(minus, -c.upper);
                    ++s.indegree[minus];
                }
            }
        // A serial schedule generation scheme gives acyclic machine orders even
        // when the incoming start-time hint overlaps or contains holes.
        std::vector<int> indegree = s.indegree, topo;
        for (int i = 0; i < int(indegree.size()); ++i)
            if (!indegree[i]) topo.push_back(i);
        for (std::size_t k = 0; k < topo.size(); ++k)
            for (auto e : s.edges[topo[k]]) {
                tick();
                if (--indegree[e.first] == 0) topo.push_back(e.first);
            }
        if (topo.size() != s.vars.size()) {
            unsupported();
            return;
        }
        std::vector<long double> tail(s.vars.size());
        for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
            tick();
            int i = *it;
            tail[i] += s.duration[i];
            for (auto e : s.edges[i])
                tail[i] = std::max(tail[i], static_cast<long double>(e.second) + tail[e.first]);
        }
        bool any_hint = !options_.initial_solution.empty();
        std::priority_queue<std::pair<long double, int>, std::vector<std::pair<long double, int>>,
                            std::greater<std::pair<long double, int>>>
            ready;
        auto push = [&](int i) {
            ready.emplace(any_hint ? static_cast<long double>(x_[s.vars[i]]) : -tail[i], i);
        };
        indegree = s.indegree;
        for (auto &machine : s.order) machine.clear();
        for (int i = 0; i < int(indegree.size()); ++i)
            if (!indegree[i]) push(i);
        while (!ready.empty()) {
            tick();
            int i = ready.top().second;
            ready.pop();
            s.order[s.machine[i]].push_back(i);
            for (auto e : s.edges[i]) {
                tick();
                if (--indegree[e.first] == 0) push(e.first);
            }
        }
        std::vector<Int> values;
        Int obj;
        if (!project_schedule(s, s.order, values, obj, s.critical)) {
            unsupported();
            return;
        }
        schedule_ = std::move(s);
        schedule_mode_ = true;
        schedule_tried_ = true;
        x_ = std::move(values);
        objective_ = obj;
        raw_total_ = 0;
        violated_.clear();
        remember();
        // If a feasible hint was better, its order is a valid restart state.
        schedule_.best_order = schedule_.order;
        if (result_.feasible)
            for (auto &machine : schedule_.best_order)
                std::sort(machine.begin(), machine.end(), [&](int a, int b) {
                    tick();
                    return std::make_pair(result_.x[schedule_.vars[a]], schedule_.duration[a]) <
                           std::make_pair(result_.x[schedule_.vars[b]], schedule_.duration[b]);
                });
    }
    void schedule_step() {
        if (options_.critical_block_moves) {
            block_schedule_step();
            return;
        }
        check();
        ++stats_.iterations;
        auto order = schedule_.order;
        int machine = -1, position = -1;
        if (!schedule_.critical.empty() && rng_() % 5) {
            auto p = schedule_.critical[random_index(schedule_.critical.size())];
            machine = p.first;
            position = p.second;
        } else {
            machine = random_index(order.size());
            if (order[machine].size() < 2) return;
            position = random_index(order[machine].size() - 1);
        }
        std::swap(order[machine][position], order[machine][position + 1]);
        std::vector<Int> values;
        std::vector<std::pair<int, int>> critical;
        Int obj;
        if (!project_schedule(schedule_, order, values, obj, critical)) return;
        ++stats_.evaluations;
        long double delta =
            static_cast<long double>(Wide(model_.direction_) * (Wide(obj) - objective_));
        long double temperature =
            std::max(0.2L, std::abs(static_cast<long double>(result_.objective)) * 0.1L) *
            std::pow(0.02L, static_cast<long double>(stats_.iterations % 2000) / 2000);
        if (delta <= 0 ||
            std::generate_canonical<double, 53>(rng_) < std::exp(-delta / temperature)) {
            bool improve = !result_.feasible || Wide(model_.direction_) * obj <
                                                    Wide(model_.direction_) * result_.objective;
            schedule_.order = std::move(order);
            schedule_.critical = std::move(critical);
            x_ = std::move(values);
            objective_ = obj;
            ++stats_.moves;
            remember();
            if (improve) schedule_.best_order = schedule_.order;
        }
        if (stats_.iterations % 2000 == 0) {
            std::vector<std::pair<int, int>> restart_critical;
            if (project_schedule(schedule_, schedule_.best_order, values, obj, restart_critical)) {
                schedule_.order = schedule_.best_order;
                schedule_.critical = std::move(restart_critical);
                x_ = std::move(values);
                objective_ = obj;
                ++stats_.restarts;
                remember();
            }
        }
    }
    void block_schedule_step() {
        check();
        ++stats_.iterations;
        struct Insertion {
            int machine, from, to;
        };
        std::vector<Insertion> moves;
        for (std::size_t i = 0; i < schedule_.critical.size();) {
            int m = schedule_.critical[i].first, l = schedule_.critical[i].second, r = l + 1;
            ++i;
            while (i < schedule_.critical.size() && schedule_.critical[i] == std::make_pair(m, r)) {
                ++r;
                ++i;
            }
            for (int k = l + 1; k <= r; ++k) {
                tick();
                moves.push_back({m, k, l});
                moves.push_back({m, l, k});
            }
            for (int k = l; k < r; ++k) {
                tick();
                moves.push_back({m, k, r});
                moves.push_back({m, r, k});
            }
        }
        for (int k = 0; k < 8; ++k) {
            int m = random_index(schedule_.order.size());
            if (schedule_.order[m].size() < 2) continue;
            int from = random_index(schedule_.order[m].size() - 1);
            moves.push_back({m, from, from + 1});
        }
        std::shuffle(moves.begin(), moves.end(), rng_);
        std::vector<Int> best_values;
        std::vector<std::pair<int, int>> best_critical;
        std::vector<std::vector<int>> best_order;
        Wide best_objective = Wide(bound) * 2;
        int moved = -1;
        auto arc = [](int a, int b) { return (std::uint64_t(a) << 32) | std::uint32_t(b); };
        std::vector<int> current_next(schedule_.vars.size(), -1);
        for (const auto &machine : schedule_.order)
            for (std::size_t i = 1; i < machine.size(); ++i)
                current_next[machine[i - 1]] = machine[i];
        for (std::size_t k = 0; k < std::min<std::size_t>(48, moves.size()); ++k) {
            check();
            auto move = moves[k];
            auto order = schedule_.order;
            auto &machine = order[move.machine];
            int task = machine[move.from], v = schedule_.vars[task];
            machine.erase(machine.begin() + move.from);
            machine.insert(machine.begin() + move.to, task);
            bool tabu = false;
            for (std::size_t i = 1; i < machine.size(); ++i)
                if (current_next[machine[i - 1]] != machine[i]) {
                    auto it = schedule_.tabu_arcs.find(arc(machine[i - 1], machine[i]));
                    if (it != schedule_.tabu_arcs.end() && it->second > stats_.iterations)
                        tabu = true;
                }
            std::vector<Int> values;
            std::vector<std::pair<int, int>> critical;
            Int obj;
            if (!project_schedule(schedule_, order, values, obj, critical, !options_.fast_schedule))
                continue;
            ++stats_.evaluations;
            Wide score = Wide(model_.direction_) * obj;
            bool aspiration =
                !result_.feasible || score < Wide(model_.direction_) * result_.objective;
            if (tabu && !aspiration) continue;
            if (score < best_objective || (score == best_objective && rng_() % 2)) {
                best_values = std::move(values);
                best_critical = std::move(critical);
                best_order = std::move(order);
                best_objective = score;
                moved = v;
            }
        }
        ++stall_;
        if (moved >= 0) {
            if (options_.fast_schedule) {
                Int obj;
                if (!project_schedule(schedule_, best_order, best_values, obj, best_critical))
                    return;
                best_objective = Wide(model_.direction_) * obj;
            }
            bool improve =
                !result_.feasible || best_objective < Wide(model_.direction_) * result_.objective;
            std::vector<int> following(schedule_.vars.size(), -1);
            for (const auto &machine : best_order)
                for (std::size_t i = 1; i < machine.size(); ++i)
                    following[machine[i - 1]] = machine[i];
            for (std::size_t i = 0; i < current_next.size(); ++i)
                if (current_next[i] >= 0 && current_next[i] != following[i])
                    schedule_.tabu_arcs[arc(int(i), current_next[i])] =
                        stats_.iterations + 5 + rng_() % 8;
            schedule_.order = std::move(best_order);
            schedule_.critical = std::move(best_critical);
            x_ = std::move(best_values);
            objective_ = Int(best_objective * model_.direction_);
            ++stats_.moves;
            remember();
            if (improve) schedule_.best_order = schedule_.order;
        }
        if (options_.schedule_restarts && stall_ > 384 &&
            stats_.iterations - schedule_.last_restart > 384) {
            auto order = schedule_.best_order;
            std::vector<Int> values;
            std::vector<std::pair<int, int>> critical;
            Int obj;
            if (project_schedule(schedule_, order, values, obj, critical)) {
                int accepted = 0, wanted = 2 + int(stats_.restarts % 5);
                for (int attempt = 0; attempt < 4 * wanted && accepted < wanted; ++attempt) {
                    check();
                    int m = random_index(order.size());
                    if (order[m].size() < 2) continue;
                    int p = random_index(order[m].size() - 1);
                    std::swap(order[m][p], order[m][p + 1]);
                    std::vector<Int> trial_values;
                    std::vector<std::pair<int, int>> trial_critical;
                    Int trial_obj;
                    if (project_schedule(schedule_, order, trial_values, trial_obj, trial_critical)) {
                        values = std::move(trial_values);
                        critical = std::move(trial_critical);
                        obj = trial_obj;
                        ++accepted;
                    } else
                        std::swap(order[m][p], order[m][p + 1]);
                }
                schedule_.order = std::move(order);
                schedule_.critical = std::move(critical);
                x_ = std::move(values);
                objective_ = obj;
                schedule_.tabu_arcs.clear();
                if (!result_.feasible ||
                    Wide(model_.direction_) * obj < Wide(model_.direction_) * result_.objective)
                    schedule_.best_order = schedule_.order;
                ++stats_.restarts;
                schedule_.last_restart = stats_.iterations;
                remember();
            }
        } else if (!options_.schedule_restarts && stall_ > 200 && stats_.iterations % 64 == 0) {
            std::vector<Int> values;
            std::vector<std::pair<int, int>> critical;
            Int obj;
            if (project_schedule(schedule_, schedule_.best_order, values, obj, critical)) {
                schedule_.order = schedule_.best_order;
                schedule_.critical = std::move(critical);
                x_ = std::move(values);
                objective_ = obj;
                ++stats_.restarts;
                remember();
            }
        }
        if (stats_.iterations % 512 == 0)
            for (auto it = schedule_.tabu_arcs.begin(); it != schedule_.tabu_arcs.end();) {
                tick();
                if (it->second <= stats_.iterations)
                    it = schedule_.tabu_arcs.erase(it);
                else
                    ++it;
            }
    }
    void cover_step() {
        check();
        ++stats_.iterations;
        const int n = model_.num_variables(), m = model_.num_constraints();
        auto values = result_.x;
        std::vector<int> cover(m, 0), gain(n, 0), selected;
        for (int v = 0; v < n; ++v) {
            tick();
            if (model_.objective_[v] <= 0 && model_.variables_[v].upper == 1) values[v] = 1;
            if (values[v] && model_.objective_[v] > 0 && model_.variables_[v].lower == 0)
                selected.push_back(v);
        }
        std::shuffle(selected.begin(), selected.end(), rng_);
        const std::size_t remove =
            selected.empty() ? 0 : 1 + rng_() % std::min<std::size_t>(5, selected.size());
        for (std::size_t i = 0; i < remove; ++i) values[selected[i]] = 0;
        for (int v = 0; v < n; ++v)
            if (values[v])
                for (auto ref : model_.incidence_[v]) {
                    tick();
                    ++cover[ref.constraint];
                }
        int missing = 0;
        for (int c = 0; c < m; ++c)
            if (!cover[c]) {
                ++missing;
                for (auto t : model_.constraints_[c].terms) {
                    tick();
                    ++gain[t.first];
                }
            }
        std::vector<long double> price(n);
        std::priority_queue<std::pair<long double, int>> heap;
        for (int v = 0; v < n; ++v) {
            tick();
            price[v] = std::max<Int>(1, model_.objective_[v]) *
                       (0.75L + static_cast<long double>(rng_() % 1000) / 2000);
            if (!values[v] && model_.variables_[v].upper == 1 && gain[v])
                heap.emplace(gain[v] / price[v], v);
        }
        while (missing && !heap.empty()) {
            tick();
            auto entry = heap.top();
            heap.pop();
            int v = entry.second;
            if (values[v] || !gain[v]) continue;
            if (entry.first != gain[v] / price[v]) {
                heap.emplace(gain[v] / price[v], v);
                continue;
            }
            values[v] = 1;
            for (auto ref : model_.incidence_[v]) {
                tick();
                if (cover[ref.constraint]++ == 0) {
                    --missing;
                    for (auto t : model_.constraints_[ref.constraint].terms) {
                        tick();
                        --gain[t.first];
                    }
                }
            }
        }
        if (missing) return;
        selected.clear();
        for (int v = 0; v < n; ++v)
            if (values[v] && model_.variables_[v].lower == 0 && model_.objective_[v] > 0)
                selected.push_back(v);
        std::shuffle(selected.begin(), selected.end(), rng_);
        for (int v : selected) {
            tick();
            bool redundant = true;
            for (auto ref : model_.incidence_[v]) {
                tick();
                if (cover[ref.constraint] <= 1) {
                    redundant = false;
                    break;
                }
            }
            if (redundant) {
                values[v] = 0;
                for (auto ref : model_.incidence_[v]) {
                    tick();
                    --cover[ref.constraint];
                }
            }
        }
        Wide objective = model_.offset_;
        for (int v = 0; v < n; ++v) {
            tick();
            objective += Wide(model_.objective_[v]) * values[v];
        }
        if (objective > result_.objective) return;
        Move move;
        for (int v = 0; v < n; ++v)
            if (values[v] != x_[v]) move.push_back({v, values[v]});
        if (move.empty()) return;
        evaluate_move(move);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
        // Equal-cost covers are useful distinct repair bases.
        if (result_.feasible && objective == result_.objective && rng_() % 2) result_.x = x_;
    }
    bool route_kick() {
        if (route_constraint_ < 0 || !result_.feasible || !options_.compound_moves) return false;
        const auto &c = model_.constraints_[route_constraint_];
        int n = int(c.successors.size());
        if (n < 8 || c.self_loops) return false;
        std::vector<int> tour;
        int v = 0;
        for (int k = 0; k < n; ++k) {
            tick();
            tour.push_back(v);
            v = int(result_.x[c.successors[v]]);
        }
        std::vector<int> cuts;
        while (cuts.size() < 4) {
            int a = random_index(n);
            if (std::find(cuts.begin(), cuts.end(), a) == cuts.end()) cuts.push_back(a);
        }
        std::sort(cuts.begin(), cuts.end());
        std::vector<Int> values(n);
        for (int i = 0; i < n; ++i) values[i] = result_.x[c.successors[i]];
        int a = tour[cuts[0]], b = tour[cuts[1]], d = tour[cuts[2]], e = tour[cuts[3]];
        Int an = values[a], bn = values[b], dn = values[d], en = values[e];
        values[a] = dn;
        values[b] = en;
        values[d] = an;
        values[e] = bn;
        Move move;
        for (int i = 0; i < n; ++i)
            if (!add_change(move, c.successors[i], values[i])) return false;
        if (move.empty()) return false;
        evaluate_move(move);
        ++stats_.iterations;
        ++stats_.restarts;
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
        return true;
    }
    static void clause_set(std::vector<int> &list, std::vector<int> &pos, int v, bool present) {
        if (present && pos[v] < 0) {
            pos[v] = int(list.size());
            list.push_back(v);
        } else if (!present && pos[v] >= 0) {
            int last = list.back();
            list[pos[v]] = last;
            pos[last] = pos[v];
            list.pop_back();
            pos[v] = -1;
        }
    }
    Wide clause_weight(int c) const {
        return Wide(clauses_.scale) * clauses_.base_weight[c] + clauses_.penalty[c];
    }
    Wide clause_contribution(int c, int v, Int coefficient, Wide activity) const {
        const auto &row = model_.constraints_[c];
        return (int(activity < row.lower) -
                int(activity + coefficient * (1 - 2 * x_[v]) < row.lower)) *
               clause_weight(c);
    }
    void clause_rebuild_scores() {
        const int n = model_.num_variables(), m = model_.num_constraints();
        if (clauses_.rebuild_phase < 0) {
            clauses_.ready = false;
            clauses_.good.clear();
            clauses_.bad.clear();
            clauses_.gain.assign(n, 0);
            clauses_.good_pos.assign(n, -1);
            clauses_.bad_pos.assign(m, -1);
            if (clauses_.penalty.empty()) clauses_.penalty.assign(m, 0);
            if (clauses_.last_flip.empty()) clauses_.last_flip.assign(n, 0);
            clauses_.rebuild_phase = 0;
            clauses_.rebuild_row = clauses_.rebuild_slot = 0;
        }
        if (clauses_.rebuild_phase == 0) {
            while (clauses_.rebuild_row < m) {
                int c = clauses_.rebuild_row, z = clause_slack_[c];
                const auto &row = model_.constraints_[c];
                clause_set(clauses_.bad, clauses_.bad_pos, c, x_[z] != 0);
                while (clauses_.rebuild_slot < int(row.terms.size())) {
                    tick();
                    auto t = row.terms[clauses_.rebuild_slot];
                    int v = t.first;
                    if (v != z && domains_[v].first < domains_[v].second)
                        clauses_.gain[v] +=
                            clause_contribution(c, v, t.second, cache_[c].activity - x_[z]);
                    ++clauses_.rebuild_slot;
                }
                ++clauses_.rebuild_row;
                clauses_.rebuild_slot = 0;
            }
            clauses_.rebuild_phase = 1;
            clauses_.rebuild_row = 0;
        }
        while (clauses_.rebuild_row < int(decisions_.size())) {
            tick();
            int v = decisions_[clauses_.rebuild_row++];
            clause_set(clauses_.good, clauses_.good_pos, v, clauses_.gain[v] > 0);
        }
        clauses_.ready = true;
        clauses_.rebuild_phase = -1;
    }
    void clause_update_after_flip() {
        clauses_.ready = false;
        for (int c : touched_) {
            const auto &row = model_.constraints_[c];
            const int z = clause_slack_[c];
            Wide after = cache_[c].activity - x_[z], before = after;
            for (auto t : row.terms) {
                tick();
                if (t.first != z) before -= t.second * (x_[t.first] - old_value(t.first));
            }
            const Wide weight = clause_weight(c);
            for (auto t : row.terms) {
                tick();
                int v = t.first;
                if (v == z || domains_[v].first == domains_[v].second) continue;
                int old_gain = int(before < row.lower) -
                               int(before + t.second * (1 - 2 * old_value(v)) < row.lower);
                int new_gain =
                    int(after < row.lower) - int(after + t.second * (1 - 2 * x_[v]) < row.lower);
                clauses_.gain[v] += (new_gain - old_gain) * weight;
                clause_set(clauses_.good, clauses_.good_pos, v, clauses_.gain[v] > 0);
            }
            clause_set(clauses_.bad, clauses_.bad_pos, c, x_[z] != 0);
        }
        clauses_.ready = true;
    }
    Move clause_flip_move(const std::vector<int> &variables) {
        Move move;
        std::vector<int> rows;
        for (int v : variables) {
            move.push_back({v, 1 - x_[v]});
            for (auto ref : model_.incidence_[v]) {
                tick();
                rows.push_back(ref.constraint);
            }
        }
        std::sort(rows.begin(), rows.end(), [&](int a, int b) {
            tick();
            return a < b;
        });
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        for (int c : rows) {
            const auto &row = model_.constraints_[c];
            int z = clause_slack_[c];
            Wide activity = cache_[c].activity - x_[z];
            for (auto t : row.terms) {
                tick();
                if (std::find(variables.begin(), variables.end(), t.first) != variables.end())
                    activity += t.second * (1 - 2 * x_[t.first]);
            }
            Int value = activity < row.lower;
            if (value != x_[z]) move.push_back({z, value});
        }
        return move;
    }
    std::vector<int> clause_pair() {
        std::vector<int> first;
        for (int k = 0; k < options_.clause_pair_samples && !clauses_.bad.empty(); ++k) {
            int c = clauses_.bad[random_index(clauses_.bad.size())];
            for (auto t : model_.constraints_[c].terms) {
                tick();
                if (t.first != clause_slack_[c] && domains_[t.first].first < domains_[t.first].second)
                    first.push_back(t.first);
            }
        }
        std::sort(first.begin(), first.end(), [&](int a, int b) {
            tick();
            return a < b;
        });
        first.erase(std::unique(first.begin(), first.end()), first.end());
        if (first.size() > 32) std::shuffle(first.begin(), first.end(), rng_);
        if (first.size() > 32) first.resize(32);
        Wide best = 0;
        std::vector<int> chosen;
        std::vector<int> adjacent;
        if (clauses_.pair_adjustment.empty()) {
            clauses_.pair_adjustment.resize(model_.num_variables());
            clauses_.pair_stamp.resize(model_.num_variables());
        }
        int ties = 0;
        for (int v : first) {
            adjacent.clear();
            if (++clauses_.pair_epoch == 0) {
                std::fill(clauses_.pair_stamp.begin(), clauses_.pair_stamp.end(), 0);
                ++clauses_.pair_epoch;
            }
            for (auto ref : model_.incidence_[v]) {
                int c = ref.constraint;
                const auto &row = model_.constraints_[c];
                Wide activity = cache_[c].activity - x_[clause_slack_[c]];
                Wide after = activity + row.terms[ref.slot].second * (1 - 2 * x_[v]);
                for (auto t : row.terms) {
                    tick();
                    int w = t.first;
                    if (w == v || w == clause_slack_[c] || domains_[w].first == domains_[w].second)
                        continue;
                    Wide delta = clause_contribution(c, w, t.second, after) -
                                 clause_contribution(c, w, t.second, activity);
                    if (clauses_.pair_stamp[w] != clauses_.pair_epoch) {
                        clauses_.pair_stamp[w] = clauses_.pair_epoch;
                        clauses_.pair_adjustment[w] = 0;
                        adjacent.push_back(w);
                    }
                    clauses_.pair_adjustment[w] += delta;
                }
            }
            for (int w : adjacent) {
                tick();
                ++stats_.evaluations;
                Wide gain = clauses_.gain[v] + clauses_.gain[w] + clauses_.pair_adjustment[w];
                if (gain > best) {
                    best = gain;
                    chosen = {v, w};
                    ties = 1;
                } else if (gain == best && gain > 0 && rng_() % ++ties == 0)
                    chosen = {v, w};
            }
        }
#ifdef CONSTRAINT_SEARCH_DEBUG_CHECKS
        if (!chosen.empty()) {
            Wide direct = 0;
            for (int c = 0; c < model_.num_constraints(); ++c) {
                int before = 0, after = 0;
                for (auto t : model_.constraints_[c].terms) {
                    tick();
                    if (t.first == clause_slack_[c]) continue;
                    bool value = x_[t.first] != 0;
                    before += value == (t.second > 0);
                    if (t.first == chosen[0] || t.first == chosen[1]) value = !value;
                    after += value == (t.second > 0);
                }
                direct += (int(before == 0) - int(after == 0)) * clause_weight(c);
            }
            assert(direct == best);
        }
#endif
        return chosen;
    }
    void scored_clause_step() {
        check();
        ++stats_.iterations;
        if (options_.clause_multiscale && clauses_.ready && stall_ >= 8192 &&
            stats_.iterations - clauses_.last_scale >= 8192) {
            clauses_.scale = clauses_.scale < 8 ? 2 * clauses_.scale : 1;
            clauses_.last_scale = stats_.iterations;
            clauses_.ready = false;
            clauses_.rebuild_phase = -1;
            std::fill(clauses_.penalty.begin(), clauses_.penalty.end(), 0);
        }
        if (!clauses_.ready) clause_rebuild_scores();
        if (clauses_.bad.empty()) return;
        std::vector<int> selected;
        auto choose_good = [&]() {
            int best = -1;
            const auto size = clauses_.good.size();
            for (std::size_t k = 0; k < std::min<std::size_t>(64, size); ++k) {
                tick();
                int v = clauses_.good[size <= 64 ? k : random_index(size)];
                if (best < 0 || clauses_.gain[v] > clauses_.gain[best] ||
                    (clauses_.gain[v] == clauses_.gain[best] &&
                     clauses_.last_flip[v] < clauses_.last_flip[best]))
                    best = v;
            }
            if (best >= 0) selected = {best};
        };
        choose_good();
        if (selected.empty() && options_.clause_pair_moves) selected = clause_pair();
        if (selected.empty()) {
            if (options_.adaptive_weights) {
                clauses_.ready = false;
                ++stats_.weight_updates;
                auto rows = clauses_.bad;
                auto change_weight = [&](int c, Wide delta) {
                    clauses_.penalty[c] += delta;
                    const auto &row = model_.constraints_[c];
                    int z = clause_slack_[c];
                    Wide activity = cache_[c].activity - x_[z];
                    for (auto term : row.terms) {
                        tick();
                        int v = term.first;
                        if (v == z || domains_[v].first == domains_[v].second) continue;
                        int contribution = int(activity < row.lower) -
                                           int(activity + term.second * (1 - 2 * x_[v]) < row.lower);
                        clauses_.gain[v] += contribution * delta;
                        clause_set(clauses_.good, clauses_.good_pos, v, clauses_.gain[v] > 0);
                    }
                };
                // Bounded dynamic weights remain separate from the objective.
                // Periodic smoothing removes old penalties after the search moves.
                if (stats_.weight_updates % 64 == 0) {
                    rows.resize(model_.num_constraints());
                    std::iota(rows.begin(), rows.end(), 0);
                    for (int c : rows) {
                        tick();
                        change_weight(c, clauses_.penalty[c] / 2 - clauses_.penalty[c]);
                    }
                } else
                    for (int c : rows) {
                        tick();
                        change_weight(c, std::min(Wide(100000), clauses_.penalty[c] + 1) -
                                             clauses_.penalty[c]);
                    }
                clauses_.ready = true;
                choose_good();
            }
            if (selected.empty()) {
                int c = clauses_.bad[random_index(clauses_.bad.size())], best = -1;
                for (auto t : model_.constraints_[c].terms) {
                    tick();
                    int v = t.first;
                    if (v == clause_slack_[c] || domains_[v].first == domains_[v].second) continue;
                    if (best < 0 || clauses_.gain[v] > clauses_.gain[best] ||
                        (clauses_.gain[v] == clauses_.gain[best] &&
                         clauses_.last_flip[v] < clauses_.last_flip[best]))
                        best = v;
                }
                if (best >= 0) selected = {best};
            }
        }
        if (selected.empty()) {
            ++stall_;
            return;
        }
        Move move = clause_flip_move(selected);
        clauses_.ready = false;
        evaluate_move(move);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
        for (int v : selected) clauses_.last_flip[v] = stats_.iterations;
        clause_update_after_flip();
    }
    void clause_step() {
        if (options_.clause_scores) {
            scored_clause_step();
            return;
        }
        check();
        ++stats_.iterations;
        std::vector<int> unsatisfied;
        for (int c = 0; c < model_.num_constraints(); ++c) {
            tick();
            if (x_[clause_slack_[c]]) unsatisfied.push_back(c);
        }
        if (unsatisfied.empty()) return;
        int row = unsatisfied[random_index(unsatisfied.size())];
        Move best, fallback;
        long double best_score = 1e100L, fallback_score = 1e100L;
        for (auto t : model_.constraints_[row].terms) {
            int v = t.first;
            if (v == clause_slack_[row] || domains_[v].first == domains_[v].second) continue;
            Move move{{v, 1 - x_[v]}};
            long double score = 0;
            for (auto ref : model_.incidence_[v]) {
                tick();
                int c = ref.constraint, z = clause_slack_[c];
                const auto &constraint = model_.constraints_[c];
                Int coefficient = constraint.terms[ref.slot].second;
                Wide activity = cache_[c].activity - x_[z] + coefficient * (1 - 2 * x_[v]);
                Int value = activity < constraint.lower;
                if (value != x_[z]) move.push_back({z, value});
                score += (value - x_[z]) * (model_.objective_[z] + weights_[c] - 1);
            }
            if (score < fallback_score || (score == fallback_score && rng_() % 2)) {
                fallback = move;
                fallback_score = score;
            }
            if (tabu_[v] <= stats_.iterations &&
                (score < best_score || (score == best_score && rng_() % 2))) {
                best = std::move(move);
                best_score = score;
            }
        }
        if (best.empty()) {
            best = std::move(fallback);
            best_score = fallback_score;
        }
        if (best.empty()) {
            ++stall_;
            return;
        }
        if (best_score >= 0 && options_.adaptive_weights) {
            ++stats_.weight_updates;
            for (int c : unsatisfied) {
                tick();
                weights_[c] = std::min(100000.0L, weights_[c] + 1);
            }
        }
        // A random walk in a false clause remains a valid MaxSAT assignment.
        if (stall_ > 50 && rng_() % 12 == 0) {
            int v = decisions_[random_index(decisions_.size())];
            best = {{v, 1 - x_[v]}};
            for (auto ref : model_.incidence_[v]) {
                tick();
                int c = ref.constraint, z = clause_slack_[c];
                const auto &constraint = model_.constraints_[c];
                Wide activity = cache_[c].activity - x_[z] +
                                constraint.terms[ref.slot].second * (1 - 2 * x_[v]);
                Int value = activity < constraint.lower;
                if (value != x_[z]) best.push_back({z, value});
            }
        }
        evaluate_move(best);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
    }
    void packing_price_step() {
        auto &p = packing_;
        const int n = int(p.profit.size()), m = int(p.capacity.size());
        const auto &price = p.walk_price.empty() ? p.price : p.walk_price;
        std::vector<long double> reduced(n), capacity(p.capacity.begin(), p.capacity.end()),
            gradient;
        long double upper = 0;
        for (int v = 0; v < n; ++v)
            if (domains_[v].first == 1) {
                tick();
                upper += p.profit[v];
                for (int c = 0; c < m; ++c) {
                    tick();
                    capacity[c] -= p.weight[v][c];
                }
            }
        gradient = capacity;
        for (int v = 0; v < n; ++v) {
            reduced[v] = p.profit[v];
            for (int c = 0; c < m; ++c) {
                tick();
                reduced[v] -= price[c] * p.weight[v][c];
            }
            if (domains_[v].first == 0 && domains_[v].second == 1 && reduced[v] > 0) {
                upper += reduced[v];
                for (int c = 0; c < m; ++c) {
                    tick();
                    gradient[c] -= p.weight[v][c];
                }
            }
        }
        long double norm = 0;
        for (int c = 0; c < m; ++c) {
            upper += price[c] * capacity[c];
            norm += gradient[c] * gradient[c];
        }
        long double primal = static_cast<long double>(-Wide(model_.direction_) *
                                                      (Wide(result_.objective) - model_.offset_));
        long double step = norm > 0 ? 1.5L * std::max(0.0L, upper - primal) / norm : 0;
        auto next = price;
        for (int c = 0; c < m; ++c)
            next[c] = std::max(
                0.0L, std::min(static_cast<long double>(bound), price[c] - step * gradient[c]));
        // These floating estimates guide neighborhoods only. All pruning below
        // continues to use independently valid integer/rational upper bounds.
        if (upper < p.dual_best) {
            p.price = price;
            p.reduced = std::move(reduced);
            p.dual_best = upper;
        }
        p.walk_price = std::move(next);
    }
    void packing_step() {
        check();
        ++stats_.iterations;
        auto &p = packing_;
        const int n = int(p.profit.size()), m = int(p.capacity.size());
        // Coordinate minimization of the Lagrangian bound. Floating prices are
        // used only to choose a neighborhood, never to discard integer solutions.
        if (p.rounds < 4 * m) {
            int c = p.rounds % m;
            std::vector<std::pair<long double, int>> points;
            Int capacity = p.capacity[c];
            for (int v = 0; v < n; ++v) {
                tick();
                if (domains_[v].first == 1) {
                    capacity -= p.weight[v][c];
                    continue;
                }
                if (domains_[v].second == 0 || !p.weight[v][c]) continue;
                long double profit = p.reduced[v] + p.price[c] * p.weight[v][c];
                if (profit > 0) points.emplace_back(profit / p.weight[v][c], v);
            }
            std::sort(points.begin(), points.end(), [&](auto a, auto b) {
                tick();
                return a > b;
            });
            Wide load = 0;
            long double price = 0;
            for (auto a : points) {
                tick();
                load += p.weight[a.second][c];
                if (load > capacity) {
                    price = a.first;
                    break;
                }
            }
            auto reduced = p.reduced;
            for (int v = 0; v < n; ++v) {
                tick();
                reduced[v] -= (price - p.price[c]) * p.weight[v][c];
            }
            p.reduced = std::move(reduced);
            p.price[c] = price;
            ++p.rounds;
        }
        if (options_.packing_subgradient && p.rounds >= 4 * m && stats_.iterations % 8 == 0)
            packing_price_step();
        std::vector<std::pair<long double, int>> ranked;
        const auto &base = options_.packing_cached_capacity ? x_ : result_.x;
        std::vector<Int> values = base, remaining = p.capacity;
        for (int v : decisions_) {
            tick();
            long double distance = std::abs(p.reduced[v]);
            // Occasionally sample beyond the estimated core to escape its bias.
            if (rng_() % 12 == 0) distance *= 0.02L;
            ranked.emplace_back(distance * (0.5L + (rng_() % 1000) / 1000.0L), v);
        }
        std::sort(ranked.begin(), ranked.end(), [&](auto a, auto b) {
            tick();
            return a < b;
        });
        if (ranked.empty()) return;
        const int k = std::min<int>(int(ranked.size()), 24 + rng_() % 25);
        std::vector<int> core;
        std::vector<Wide> surrogate(n, 0);
        for (int i = 0; i < k; ++i) {
            int v = ranked[i].second;
            core.push_back(v);
            values[v] = 0;
        }
        if (options_.packing_cached_capacity) {
            for (int c = 0; c < m; ++c) remaining[c] -= Int(cache_[c].activity);
            for (int v : core)
                if (base[v])
                    for (int c = 0; c < m; ++c) {
                        tick();
                        remaining[c] += p.weight[v][c];
                    }
        } else
            for (int v = 0; v < n; ++v)
                if (values[v])
                    for (int c = 0; c < m; ++c) {
                        tick();
                        remaining[c] -= p.weight[v][c];
                    }
        std::vector<Int> aggregate(m, 1);
        int units = 64 - m;
        long double total_price = 0;
        for (auto price : p.price) total_price += price;
        if (options_.packing_dual_core && total_price > 0) {
            for (int c = 0; c < m; ++c) {
                int extra = std::min(units, int((64 - m) * p.price[c] / total_price));
                aggregate[c] += extra;
                units -= extra;
            }
        }
        for (int v : core)
            for (int c = 0; c < m; ++c) {
                tick();
                surrogate[v] += Wide(aggregate[c]) * p.weight[v][c];
            }
        // Sum of at most 64 rows: cross-products fit signed __int128 because
        // every original row and the whole objective are bounded by 2^60.
        std::sort(core.begin(), core.end(), [&](int a, int b) {
            if (!surrogate[a] || !surrogate[b])
                return !surrogate[a] != !surrogate[b] ? !surrogate[a] : a < b;
            tick();
            Wide left = Wide(p.profit[a]) * surrogate[b], right = Wide(p.profit[b]) * surrogate[a];
            return left != right ? left > right : a < b;
        });
        Wide best = 0;
        std::uint64_t best_mask = 0;
        for (int i = 0; i < k; ++i)
            if (base[core[i]]) {
                best += p.profit[core[i]];
                best_mask |= std::uint64_t(1) << i;
            }
        std::uint64_t nodes = 0;
        const auto local_end = std::min(deadline_ - std::chrono::microseconds(100),
                                        Clock::now() + std::chrono::milliseconds(2));
        bool stopped = false;
        auto search = [&](auto &&self, int depth, Wide profit, std::uint64_t mask) -> void {
            if (stopped) return;
            if ((++nodes & 63) == 0 && (nodes >= 32768 || Clock::now() >= local_end)) {
                stopped = true;
                return;
            }
            if (profit > best) {
                best = profit;
                best_mask = mask;
            }
            Wide capacity = 0;
            for (int c = 0; c < m; ++c) capacity += Wide(aggregate[c]) * remaining[c];
            Wide upper = profit;
            for (int i = depth; i < k; ++i) {
                tick();
                int v = core[i];
                if (surrogate[v] <= capacity) {
                    capacity -= surrogate[v];
                    upper += p.profit[v];
                } else {
                    upper += Wide(p.profit[v]) * capacity / surrogate[v];
                    break;
                }
            }
            if (upper <= best || depth == k) return;
            int v = core[depth];
            bool fits = true;
            for (int c = 0; c < m; ++c) {
                tick();
                if (p.weight[v][c] > remaining[c]) fits = false;
            }
            if (fits) {
                for (int c = 0; c < m; ++c) remaining[c] -= p.weight[v][c];
                self(self, depth + 1, profit + p.profit[v], mask | (std::uint64_t(1) << depth));
                for (int c = 0; c < m; ++c) remaining[c] += p.weight[v][c];
            }
            self(self, depth + 1, profit, mask);
        };
        search(search, 0, 0, 0);
        stats_.evaluations += nodes;
        for (int i = 0; i < k; ++i) values[core[i]] = (best_mask >> i) & 1;
        Move move;
        for (int v = 0; v < n; ++v)
            if (x_[v] != values[v]) move.push_back({v, values[v]});
        if (move.empty()) return;
        evaluate_move(move);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
    }
    Move permutation_cycle(int v, int u, bool longer) {
        const auto &group = model_.constraints_[permutation_[v]].expressions;
        const int n = int(group.size());
        std::vector<int> queue{u}, parent(model_.num_variables(), -2);
        parent[u] = -1;
        int last = -1;
        for (std::size_t i = 0; i < queue.size() && last < 0; ++i) {
            int a = queue[i];
            if ((!longer || a != u) && allowed(a, x_[v])) {
                last = a;
                break;
            }
            int offset = random_index(n);
            for (int k = 0; k < n; ++k) {
                tick();
                int b = group[(offset + k) % n].var;
                if (b == v || parent[b] != -2 || !allowed(a, x_[b])) continue;
                parent[b] = a;
                queue.push_back(b);
            }
        }
        if (last < 0) return {};
        Move move{{v, x_[u]}};
        Int value = x_[v];
        for (int a = last; a >= 0; a = parent[a]) {
            move.push_back({a, value});
            value = x_[a];
        }
        return move;
    }
    void prepare_latin_sat() {
        if (!latin_encoding_tried_) prepare_latin_encoding();
        const auto &e = latin_encoding_;
        if (e.groups.empty() || e.variables.size() > 6000) {
            latin_sat_tried_ = latin_sat_done_ = true;
            return;
        }
        std::size_t variables = e.variables.size(), clauses = 0, literals = 0;
        for (const auto &group : e.groups) {
            tick();
            std::size_t k = group.size(), pairs = k * (k ? k - 1 : 0) / 2;
            ++clauses;
            literals += k;
            if (k <= 6) {
                clauses += pairs;
                literals += 2 * pairs;
            } else {
                variables += k - 1;
                clauses += 3 * k - 4;
                literals += 2 * (3 * k - 4);
            }
        }
        stats_.latin_candidates = e.variables.size();
        stats_.latin_cells = e.cells;
        stats_.boolean_variables = variables;
        stats_.boolean_clauses = clauses;
        if (variables > 20000 || clauses > 100000 || literals > 300000) {
            latin_sat_tried_ = latin_sat_done_ = true;
            return;
        }
        detail::BooleanSearch sat;
        sat.configure_groups(true, 60000);
        for (std::size_t i = 0; i < e.variables.size(); ++i) {
            tick();
            sat.add_variable(e.values[i] == result_.x[e.variables[i]]);
        }
        for (const auto &group : e.groups) {
            tick();
            std::vector<int> candidates;
            for (int row : group) {
                tick();
                candidates.push_back(2 * row);
            }
            sat.add_branch_group(candidates);
            sat.add_clause(candidates);
            if (candidates.size() <= 6) {
                for (std::size_t i = 0; i < candidates.size(); ++i)
                    for (std::size_t j = 0; j < i; ++j) {
                        tick();
                        sat.add_clause({candidates[i] ^ 1, candidates[j] ^ 1});
                    }
            } else {
                int previous = 2 * sat.add_variable();
                sat.add_clause({candidates.front() ^ 1, previous});
                for (std::size_t i = 1; i + 1 < candidates.size(); ++i) {
                    tick();
                    int next = 2 * sat.add_variable();
                    sat.add_clause({candidates[i] ^ 1, next});
                    sat.add_clause({previous ^ 1, next});
                    sat.add_clause({candidates[i] ^ 1, previous ^ 1});
                    previous = next;
                }
                sat.add_clause({candidates.back() ^ 1, previous ^ 1});
            }
        }
        latin_sat_ = std::move(sat);
        latin_sat_tried_ = true;
    }
    bool latin_sat_step() {
        if (!latin_sat_tried_ && deadline_ - Clock::now() < std::chrono::milliseconds(3))
            return false;
        struct Elapsed {
            double &total;
            Clock::time_point begin = Clock::now();
            ~Elapsed() { total += std::chrono::duration<double>(Clock::now() - begin).count(); }
        } elapsed{latin_sat_elapsed_};
        if (!latin_sat_tried_) prepare_latin_sat();
        if (latin_sat_done_) return false;
        auto before = latin_sat_.conflicts();
        int status =
            latin_sat_.run(std::min(deadline_, Clock::now() + std::chrono::milliseconds(2)));
        stats_.evaluations += latin_sat_.conflicts() - before;
        stats_.boolean_conflicts = latin_sat_.conflicts();
        if (status == -1 || latin_sat_.saturated()) {
            latin_sat_done_ = true;
            return false;
        }
        if (!status) return false;
        const auto &e = latin_encoding_;
        auto values = e.fixed;
        for (std::size_t i = 0; i < e.variables.size(); ++i) {
            tick();
            if (latin_sat_.assignment()[i] == 1) values[e.variables[i]] = e.values[i];
        }
        Move move;
        for (int v = 0; v < model_.num_variables(); ++v) {
            tick();
            if (values[v] != x_[v]) move.push_back({v, values[v]});
        }
        evaluate_move(move);
        Wide raw = raw_total_;
        for (int c : touched_) raw += next_value_[c] - cache_[c].value;
        latin_sat_done_ = true;
        if (raw) {
            rollback();
            return false;
        }
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
        return true;
    }
    void prepare_latin_encoding() {
        if (!model_.definitions_.empty() || !options_.domain_reduction) {
            latin_encoding_tried_ = true;
            return;
        }
        const int n = model_.num_variables(), m = model_.num_constraints();
        LatinEncoding e;
        e.fixed.assign(n, unknown);
        auto masks = domain_masks_;
        std::vector<int> variable_column(n, -1), different;
        std::vector<std::vector<int>> group_column(m), membership(n);
        std::size_t free_values = 0;
        for (int v = 0; v < n; ++v) {
            tick();
            if (!masks[v]) {
                const auto &d = model_.variables_[v];
                if (d.upper - d.lower >= 64) {
                    latin_encoding_tried_ = true;
                    return;
                }
                int low = int(domains_[v].first - d.lower), high = int(domains_[v].second - d.lower);
                masks[v] = (~std::uint64_t(0) << low) & (~std::uint64_t(0) >> (63 - high));
            }
            if (domains_[v].first == domains_[v].second)
                e.fixed[v] = domains_[v].first;
            else if (__builtin_popcountll(masks[v]) == 1)
                e.fixed[v] = model_.variables_[v].lower + __builtin_ctzll(masks[v]);
            else
                free_values += __builtin_popcountll(masks[v]);
        }
        if (free_values > 6000) {
            latin_encoding_tried_ = true;
            return;
        }
        for (int c = 0; c < m; ++c)
            if (model_.constraints_[c].kind == detail::Kind::Different) {
                const auto &row = model_.constraints_[c];
                int k = int(row.expressions.size());
                if (!k || k > 64) {
                    latin_encoding_tried_ = true;
                    return;
                }
                Int lower = model_.variables_[row.expressions[0].var].lower;
                std::vector<bool> occupied(k, false);
                std::vector<int> ids;
                for (auto a : row.expressions) {
                    tick();
                    const auto &d = model_.variables_[a.var];
                    if (d.lower != lower || Wide(d.upper) - lower + 1 != k) {
                        latin_encoding_tried_ = true;
                        return;
                    }
                    ids.push_back(a.var);
                    if (e.fixed[a.var] != unknown) {
                        int bit = int(e.fixed[a.var] - lower);
                        if (occupied[bit]) {
                            latin_encoding_tried_ = true;
                            return;
                        }
                        occupied[bit] = true;
                    } else {
                        membership[a.var].push_back(c);
                        if (membership[a.var].size() > 8) {
                            latin_encoding_tried_ = true;
                            return;
                        }
                    }
                }
                std::sort(ids.begin(), ids.end());
                if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
                    latin_encoding_tried_ = true;
                    return;
                }
                group_column[c].assign(k, -1);
                for (int a = 0; a < k; ++a)
                    if (!occupied[a]) group_column[c][a] = -2;
                different.push_back(c);
            }
        auto add_column = [&]() {
            int id = int(e.groups.size());
            e.groups.emplace_back();
            return id;
        };
        for (int v = 0; v < n; ++v)
            if (e.fixed[v] == unknown) {
                tick();
                variable_column[v] = add_column();
                ++e.cells;
            }
        for (int c : different)
            for (int &column : group_column[c])
                if (column == -2) {
                    tick();
                    column = add_column();
                }
        // Incumbent-compatible candidates come first, without changing domains.
        for (int pass = 0; pass < 2; ++pass)
            for (int v = 0; v < n; ++v)
                if (e.fixed[v] == unknown) {
                    for (auto bits = masks[v]; bits; bits &= bits - 1) {
                        tick();
                        Int value = model_.variables_[v].lower + __builtin_ctzll(bits);
                        if ((value != result_.x[v]) != (pass != 0)) continue;
                        std::vector<int> columns{variable_column[v]};
                        bool possible = true;
                        for (int c : membership[v]) {
                            tick();
                            Int lower =
                                model_.variables_[model_.constraints_[c].expressions[0].var].lower;
                            int id = group_column[c][int(value - lower)];
                            if (id < 0) {
                                possible = false;
                                break;
                            }
                            columns.push_back(id);
                        }
                        if (!possible) continue;
                        int row = int(e.variables.size());
                        e.variables.push_back(v);
                        e.values.push_back(value);
                        for (int c : columns) {
                            tick();
                            e.groups[c].push_back(row);
                        }
                    }
                }
        latin_encoding_ = std::move(e);
        latin_encoding_tried_ = true;
    }
    void permutation_step() {
        check();
        if (options_.latin_sat && recolor_mode_ && !latin_sat_done_) {
            if (latin_sat_elapsed_ < .05 || stats_.iterations % 1024 == 0) {
                bool solved = latin_sat_step();
                if (solved || ((latin_sat_elapsed_ < .05) && latin_sat_tried_ && !latin_sat_done_)) {
                    ++stats_.iterations;
                    return;
                }
            }
        }
        const std::uint64_t restart_period = recolor_mode_ ? 8192 : 1024;
        if (options_.permutation_restarts && result_.violation <= 2 && stall_ >= restart_period &&
            stats_.iterations - last_recolor_restart_ >= restart_period) {
            auto values = result_.x;
            int perturbations = 2 + int(stats_.restarts % 8) * 2;
            if (recolor_mode_)
                for (int i = 0; i < perturbations; ++i) {
                    tick();
                    int v = decisions_[random_index(decisions_.size())];
                    values[v] = random_value(v);
                }
            else
                for (int attempt = 0, changed = 0;
                     attempt < 8 * perturbations && changed < perturbations; ++attempt) {
                    tick();
                    int v = decisions_[random_index(decisions_.size())];
                    const auto &group = model_.constraints_[permutation_[v]].expressions;
                    int u = group[random_index(group.size())].var;
                    if (u != v && permutation_[u] == permutation_[v] && allowed(v, values[u]) &&
                        allowed(u, values[v])) {
                        std::swap(values[v], values[u]);
                        ++changed;
                    }
                }
            Move move;
            for (int v = 0; v < model_.num_variables(); ++v) {
                tick();
                if (values[v] != x_[v]) move.push_back({v, values[v]});
            }
            evaluate_move(move);
            pending_ = true;
            commit_cursor_ = 0;
            finish_commit();
            ++stats_.iterations;
            ++stats_.restarts;
            last_recolor_restart_ = stats_.iterations;
            return;
        }
        std::vector<int> candidates;
        for (int id : violated_) {
            const auto &c = model_.constraints_[id];
            if (c.kind != detail::Kind::Different) continue;
            for (auto a : c.expressions) {
                tick();
                if (a.var < 0 || domains_[a.var].first == domains_[a.var].second) continue;
                auto it = cache_[id].frequency.find(detail::value(a, [&](int v) { return x_[v]; }));
                if (it != cache_[id].frequency.end() && it->second > 1)
                    candidates.push_back(a.var);
            }
        }
        if (candidates.empty()) {
            step();
            return;
        }
        int v = candidates[random_index(candidates.size())], id = permutation_[v];
        if (id < 0) {
            step();
            return;
        }
        Move best, fallback;
        long double score_best = 1e100L, score_fallback = 1e100L;
        if (recolor_mode_) {
            if (color_tabu_.empty())
                color_tabu_.assign(std::size_t(model_.num_variables()) * 64, 0);
            int ties = 0;
            for (int a : candidates) {
                for (auto bits = domain_masks_[a]; bits; bits &= bits - 1) {
                    tick();
                    int bit = __builtin_ctzll(bits);
                    Int value = model_.variables_[a].lower + bit;
                    if (value == x_[a]) continue;
                    Wide delta = 0;
                    long double score = 0;
                    for (auto ref : model_.incidence_[a]) {
                        tick();
                        const int c = ref.constraint;
                        if (model_.constraints_[c].kind != detail::Kind::Different) continue;
                        const auto &frequency = cache_[c].frequency;
                        auto from = frequency.find(x_[a]), to = frequency.find(value);
                        int removed = from->second - 1,
                            added = to == frequency.end() ? 0 : to->second;
                        delta += added - removed;
                        score += added - removed;
                    }
                    ++stats_.evaluations;
                    bool tabu = color_tabu_[std::size_t(a) * 64 + bit] > stats_.iterations;
                    if (!tabu || raw_total_ + delta < result_.violation) {
                        if (score < score_best) {
                            best = {{a, value}};
                            score_best = score;
                            ties = 1;
                        } else if (score == score_best && rng_() % ++ties == 0)
                            best = {{a, value}};
                    }
                    if (score < score_fallback || (score == score_fallback && rng_() % 2)) {
                        fallback = {{a, value}};
                        score_fallback = score;
                    }
                }
            }
            if (best.empty()) {
                best = std::move(fallback);
                score_best = score_fallback;
            }
            ++stats_.iterations;
            if (best.empty()) {
                ++stall_;
                return;
            }
            int a = best.front().var, old_color = int(x_[a] - model_.variables_[a].lower);
            color_tabu_[std::size_t(a) * 64 + old_color] =
                stats_.iterations + 5 + candidates.size() / 3 + rng_() % 7;
            evaluate_move(best);
            pending_ = true;
            commit_cursor_ = 0;
            finish_commit();
            return;
        }
        for (auto a : model_.constraints_[id].expressions) {
            tick();
            int u = a.var;
            if (u == v) continue;
            Move move;
            if (!allowed(v, x_[u])) continue;
            if (allowed(u, x_[v]) && (!domain_masks_[v] || rng_() % 8)) {
                add_change(move, v, x_[u]);
                add_change(move, u, x_[v]);
            } else
                move = permutation_cycle(v, u, allowed(u, x_[v]));
            if (move.empty()) continue;
            evaluate_move(move);
            long double score = move_score();
            Wide raw = raw_total_;
            for (int c : touched_) raw += next_value_[c] - cache_[c].value;
            bool aspiration = raw < result_.violation;
            if (score < score_fallback || (score == score_fallback && rng_() % 2)) {
                fallback = move;
                score_fallback = score;
            }
            bool tabu = false;
            for (auto change : move)
                if (tabu_[change.var] > stats_.iterations) tabu = true;
            if (!tabu || aspiration) {
                if (score < score_best || (score == score_best && rng_() % 2)) {
                    best = move;
                    score_best = score;
                }
            }
            rollback();
        }
        ++stats_.iterations;
        if (best.empty()) {
            best = std::move(fallback);
            score_best = score_fallback;
        }
        if (best.empty()) {
            ++stall_;
            return;
        }
        if (score_best >= 0 && stall_ % 16 == 0) bump_weights();
        evaluate_move(best);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
    }
    void step() {
        check();
        if (route_constraint_ >= 0 && stall_ > 100 && stats_.iterations % 512 == 0 && route_kick())
            return;
        if (route_constraint_ < 0 && stall_ > 500 && stats_.iterations % 128 == 0 &&
            !decisions_.empty() && options_.compound_moves) {
            Move kick;
            for (int k = 0; k < 6; ++k) {
                int v = decisions_[random_index(decisions_.size())];
                add_change(kick, v, random_value(v));
            }
            preserve_permutations(kick);
            if (!kick.empty()) {
                evaluate_move(kick);
                ++stats_.iterations;
                ++stats_.restarts;
                pending_ = true;
                commit_cursor_ = 0;
                finish_commit();
                return;
            }
        }
        Move best;
        long double best_score = std::numeric_limits<long double>::infinity();
        Int best_obj = 0;
        const bool shake = stall_ > 200 && rng_() % 32 == 0;
        for (int k = 0; k < options_.candidates; ++k) {
            tick();
            int id = -1;
            if (route_constraint_ >= 0 && violated_.empty() && options_.compound_moves)
                id = route_constraint_;
            else if (!violated_.empty() && rng_() % 5)
                id = violated_[random_index(violated_.size())];
            else if (model_.num_constraints() && rng_() % 3)
                id = random_index(model_.constraints_.size());
            Move move = propose(id);
            preserve_permutations(move);
            if (move.empty()) continue;
            evaluate_move(move);
            bool tabu = false, feasible = true;
            int bad = int(violated_.size());
            for (int c : touched_) bad += (next_value_[c] != 0) - (cache_[c].value != 0);
            feasible = bad == 0;
            for (auto a : move)
                if (a.value == undo_[a.var] && tabu_[a.var] > stats_.iterations) tabu = true;
            const bool aspiration =
                feasible && (!result_.feasible || Wide(model_.direction_) * next_objective_ <
                                                      Wide(model_.direction_) * result_.objective);
            long double score = move_score();
            bool select =
                (!tabu || aspiration || shake) &&
                (score < best_score - 1e-12L ||
                 (std::abs(score - best_score) <= 1e-12L &&
                  (Wide(model_.direction_) * next_objective_ < Wide(model_.direction_) * best_obj ||
                   rng_() % 3 == 0)));
            if (select) {
                best = std::move(move);
                best_score = score;
                best_obj = next_objective_;
            }
            rollback();
            if (shake && !best.empty()) break;
        }
        ++stats_.iterations;
        if (best_score >= -1e-12L) bump_weights();
        if (best.empty()) {
            ++stall_;
            return;
        }
        if (route_constraint_ >= 0 && violated_.empty() && best_score > 0) {
            long double delta =
                static_cast<long double>(Wide(model_.direction_) * (Wide(best_obj) - objective_));
            const long double temperature =
                std::max(0.1L, objective_scale_ * 0.15L) *
                std::pow(0.01L, static_cast<long double>(stats_.iterations % 1000) / 1000);
            if (std::generate_canonical<double, 53>(rng_) >= std::exp(-delta / temperature)) {
                ++stall_;
                return;
            }
        }
        evaluate_move(best);
        pending_ = true;
        commit_cursor_ = 0;
        finish_commit();
    }

public:
    explicit Solver(const Model &model, Options options = {})
        : model_(model), options_(std::move(options)), revision_(model.revision_),
          rng_(options_.seed) {
        if (options_.initial_solution.size() > std::size_t(model.num_variables()))
            throw std::invalid_argument("warm start too long");
        if (options_.candidates < 1 || options_.candidates > 1024)
            throw std::invalid_argument("candidates must be in [1,1024]");
        if (options_.clause_pair_samples < 1 || options_.clause_pair_samples > 64)
            throw std::invalid_argument("clause pair samples must be in [1,64]");
    }
    // The same object retains the assignment, caches, weights, RNG and incumbent.
    // Returned reference remains valid until the next call or object destruction.
    const Result &run(double seconds) {
        if (!std::isfinite(seconds) || seconds < 0)
            throw std::invalid_argument("time limit must be finite and nonnegative");
        if (model_.revision_ != revision_)
            throw std::logic_error("model changed while Solver is alive");
        started_ = Clock::now();
        const double safe_seconds = std::min(seconds, 1.0e8);
        deadline_ = std::min(
            options_.deadline,
            started_ + std::chrono::duration_cast<Clock::duration>(
                           std::chrono::duration<double>(std::max(0.0, safe_seconds - 0.00015))));
        const std::uint64_t before = stats_.iterations;
        result_.reason = StopReason::TimeLimit;
        try {
            check();
            prepare();
            if (pending_) finish_commit();
            try_schedule();
            while (stats_.iterations - before < options_.iteration_limit) {
                check();
                if (result_.feasible && !has_objective_) {
                    result_.reason = StopReason::Satisfied;
                    break;
                }
                if (result_.feasible && has_objective_ && result_.objective == objective_bound_) {
                    result_.reason = StopReason::ObjectiveBound;
                    break;
                }
                if (decisions_.empty()) {
                    result_.reason = StopReason::Fixed;
                    break;
                }
                if (schedule_mode_)
                    schedule_step();
                else if (cover_mode_ && result_.feasible)
                    cover_step();
                else if (clause_mode_)
                    clause_step();
                else if (packing_mode_)
                    packing_step();
                else if (permutation_mode_)
                    permutation_step();
                else
                    step();
            }
            if (stats_.iterations - before >= options_.iteration_limit)
                result_.reason = StopReason::IterationLimit;
        } catch (const Timeout &) {
            if (!pending_) rollback();
        }
        stats_.elapsed = std::chrono::duration<double>(Clock::now() - started_).count();
        stats_.total_elapsed += stats_.elapsed;
        result_.stats = stats_;
        return result_;
    }
    const Result &run() { return run(options_.time_limit); }
    const Result &result() const { return result_; }
    // Expensive diagnostic, intentionally outside run()'s time budget.
    bool check_consistency() const {
        if (phase_ != 4 || pending_) return false;
        auto e = model_.evaluate(x_);
        if (e.objective != objective_ || e.violation != static_cast<long double>(raw_total_) ||
            e.feasible != violated_.empty())
            return false;
        if (!schedule_mode_)
            for (int id = 0; id < model_.num_constraints(); ++id) {
                const auto &c = model_.constraints_[id];
                auto get = [&](int v) { return x_[v]; };
                Wide raw = detail::raw_violation(c, get, detail::NoTick{});
                if (raw != cache_[id].raw || (detail::active(c, get) ? raw : 0) != cache_[id].value)
                    return false;
                if (!recolor_mode_ && permutation_constraint_[id] && raw == 0) {
                    for (auto a : c.expressions) {
                        auto found = cache_[id].owner.find(detail::value(a, get));
                        if (found == cache_[id].owner.end() || found->second != a.var)
                            return false;
                    }
                }
            }
        if (clause_mode_ && clauses_.ready) {
            std::vector<Wide> gain(model_.num_variables(), 0);
            for (int c = 0; c < model_.num_constraints(); ++c) {
                const auto &row = model_.constraints_[c];
                int count = 0;
                for (auto t : row.terms)
                    if (t.first != clause_slack_[c]) count += (x_[t.first] != 0) == (t.second > 0);
                if ((clauses_.bad_pos[c] >= 0) != (count == 0)) return false;
                if (clauses_.bad_pos[c] >= 0 && clauses_.bad[clauses_.bad_pos[c]] != c)
                    return false;
                for (auto t : row.terms)
                    if (t.first != clause_slack_[c] &&
                        domains_[t.first].first < domains_[t.first].second) {
                        if (count == 0)
                            gain[t.first] += clause_weight(c);
                        else if (count == 1 && ((x_[t.first] != 0) == (t.second > 0)))
                            gain[t.first] -= clause_weight(c);
                    }
            }
            for (int v : decisions_) {
                if (gain[v] != clauses_.gain[v]) return false;
                if ((clauses_.good_pos[v] >= 0) != (clauses_.gain[v] > 0)) return false;
                if (clauses_.good_pos[v] >= 0 && clauses_.good[clauses_.good_pos[v]] != v)
                    return false;
            }
        }
        return true;
    }
    Result take_result() { return std::move(result_); }
};
inline Result solve(const Model &model, Options options = {}) {
    const auto started = Clock::now();
    if (!std::isfinite(options.time_limit) || options.time_limit < 0)
        throw std::invalid_argument("time limit must be finite and nonnegative");
    // Unlike persistent run(), one-shot solve() also has to release its caches.
    // Reserve a size-dependent cleanup margin; no uninterruptible final scan.
    double cleanup = Solver::cleanup_margin(model, options);
    auto end = started + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(
                             std::max(0.0, std::min(1e8, options.time_limit) - cleanup)));
    if (options.deadline > started) {
        const auto reserve =
            std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(cleanup));
        options.deadline -= std::min(reserve, options.deadline - started);
    }
    options.deadline = std::min(options.deadline, end);
    Result r;
    {
        Solver solver(model, std::move(options));
        solver.run();
        r = solver.take_result();
    }
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - started).count();
    r.stats.total_elapsed = r.stats.elapsed;
    return r;
}
} // namespace constraint_search
#endif
