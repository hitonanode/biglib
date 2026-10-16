#ifndef ASYMMETRIC_TSP_HPP
#define ASYMMETRIC_TSP_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace asymmetric_tsp {
using Clock = std::chrono::steady_clock;
using Cost = long double;

enum class StopReason { time_limit, trial_limit, trivial };

struct Options {
    double time_limit = 1.0; // Wall seconds, including preprocessing.
    Clock::time_point deadline = Clock::time_point::max(); // Optional outer deadline.
    double safety_margin = 0.002;                          // Reserved for returning the incumbent.
    std::uint64_t seed = 1;
    std::uint64_t max_trials = std::numeric_limits<std::uint64_t>::max();
    int candidates = 24;
    int max_exchange = 31;        // Sequential directed exchanges: 3..31.
    int search_breadth = 6;       // Branches beyond the first added arc.
    int search_nodes = 256;       // Per starting vertex, bounds deep search.
    int deep_search_nodes = 1024; // Node budget in the intensification phase.
    bool adaptive_depth = true;   // Deep search for sparsely split assignments.
    bool merge_tours = true;      // Exchange agreeing vertex-set subpaths.
    bool staged_search = true;    // Intensify only after the shallow descent.
    bool lookahead = true;        // Rank additions by the next removable arc.
    bool assignment = true;
    double assignment_fraction = 0.25;
    bool reversal = true;               // Includes all internal reverse-arc costs.
    bool ascent_candidates = true;      // Degree-penalized fixed-edge 1-tree candidates.
    bool structural_candidates = false; // Optional MST sensitivity ranking.
    bool active_nodes = true;           // Skip unchanged, exhausted search origins.
    bool patching = true;               // Reconnect improving two-cycle splits.
    bool guided_kicks = true;           // Also perturb expensive reduced-cost arcs.
    bool restarts = true;               // Diversify with candidate-based greedy tours.
    Cost relative_tolerance = 1e-12L;
    std::vector<int> initial_tour; // Permutation; omit repeated start vertex.
};

struct Statistics {
    double elapsed = 0;
    std::uint64_t trials = 0, moves = 0, improvements = 0, search_nodes = 0;
    std::uint64_t patching_moves = 0, merging_moves = 0;
    bool assignment_complete = false, assignment_repaired = false, candidates_complete = false;
    bool tree_candidates = false, one_tree_candidates = false;
    int assignment_cycles = 0, ascent_iterations = 0;
};

struct Result {
    std::vector<int> tour; // Closure tour.back() -> tour.front() implied.
    Cost cost = 0;
    Cost lower_bound =
        -std::numeric_limits<Cost>::infinity(); // Maximum of completed assignment/1-tree bounds.
    StopReason reason = StopReason::time_limit;
    Statistics stats;
};

namespace detail {
struct Budget {
    Clock::time_point start, stop;
    unsigned checks = 0;
    Budget(Clock::time_point started, const Options &o) : start(started) {
        if (std::isnan(o.time_limit) || o.time_limit < 0 || !std::isfinite(o.safety_margin) ||
            o.safety_margin < 0)
            throw std::invalid_argument("asymmetric_tsp: invalid time budget");
        auto end = o.deadline;
        const double representable =
            std::chrono::duration<double>(Clock::time_point::max() - start).count();
        if (o.time_limit < representable)
            end = std::min(end, start + std::chrono::duration_cast<Clock::duration>(
                                            std::chrono::duration<double>(o.time_limit)));
        if (end <= start) {
            stop = start;
            return;
        }
        if (end == Clock::time_point::max()) {
            stop = end;
            return;
        }
        const double available = std::chrono::duration<double>(end - start).count();
        const double reserve = std::min(o.safety_margin, available * 0.1);
        stop = end -
               std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(reserve));
    }
    bool expired() const { return stop != Clock::time_point::max() && Clock::now() >= stop; }
    bool poll() { return (++checks & 255U) == 0 && expired(); }
    bool phase_expired(Clock::time_point end) {
        return (++checks & 255U) == 0 && Clock::now() >= end;
    }
    Clock::time_point phase(double fraction) const {
        if (stop == Clock::time_point::max()) return stop;
        return start + std::chrono::duration_cast<Clock::duration>((stop - start) * fraction);
    }
};

template <class Distance> class Engine {
    friend struct EngineTestAccess;
    int n;
    const Cost edge_limit;
    Distance &distance;
    const Options &opt;
    Budget budget;
    Result result;
    std::mt19937_64 rng;
    std::vector<int> tour, next, prev, position, order, work;
    struct Candidate {
        int vertex;
        Cost cost, reduced;
    };
    struct Branch {
        Cost rank, added;
        int head;
        bool operator<(const Branch &b) const {
            return rank < b.rank || (!(b.rank < rank) && head < b.head);
        }
    };
    std::vector<std::vector<Candidate>> candidate;
    std::vector<unsigned char> active, cut_used;
    std::vector<Cost> row, col, reverse_prefix, candidate_out, candidate_in;
    std::vector<Cost> arc_cost, arc_reduced;
    std::vector<int> best_position, merge_difference;
    std::vector<Cost> best_arc;
    std::uint64_t merge_version = std::numeric_limits<std::uint64_t>::max();
    Cost root_candidate = 0;
    bool one_tree = false;
    Cost current = 0;
    std::array<int, 31> tails{}, heads{};
    std::array<std::vector<Branch>, 31> branches;
    int nodes_left = 0, exchange_limit = 5;

    Cost c(int a, int b) const {
        const Cost value = static_cast<Cost>(distance(a, b));
        using Input = typename std::decay<decltype(distance(a, b))>::type;
        // The entire integral input range fits below the limit for every int n.
        if constexpr (std::is_integral<Input>::value &&
                      std::numeric_limits<Cost>::max_exponent >
                          std::numeric_limits<Input>::digits +
                              2 * std::numeric_limits<int>::digits + 8)
            return value;
        // Checked on access: no uninterruptible O(n^2) validation pass.
        if (!std::isfinite(value) || std::abs(value) > edge_limit)
            throw std::invalid_argument("asymmetric_tsp: non-finite or excessive edge cost");
        return value;
    }
    Cost reduced(int a, int b) const { return c(a, b) - row[a] - col[b]; }
    bool improvement(Cost gain, Cost scale) const { return gain > opt.relative_tolerance * scale; }
    int advance(int i) const { return i + 1 == n ? 0 : i + 1; }
    Cost value(const std::vector<int> &t) const {
        if (n < 2) return 0;
        // Compensated summation also handles mixed-sign costs.
        Cost sum = 0, correction = 0;
        for (int i = 0; i < n; ++i) {
            const Cost x = c(t[i], t[advance(i)]);
            const Cost y = x - correction, s = sum + y;
            correction = (s - sum) - y;
            sum = s;
        }
        return sum;
    }
    void remember(const std::vector<int> &t, Cost v) {
        if (v < result.cost) {
            result.cost = v;
            result.tour = t;
            ++result.stats.improvements;
        }
    }
    void set_tour(const std::vector<int> &t) {
        for (int i = 0; i < n; ++i) next[t[i]] = t[advance(i)];
        std::fill(active.begin(), active.end(), 1);
        rebuild();
    }
    void rebuild(Cost known_cost = std::numeric_limits<Cost>::quiet_NaN()) {
#ifdef ASYMMETRIC_TSP_CHECK_INVARIANTS
        std::vector<unsigned char> seen(n, 0);
        int check = 0;
        for (int i = 0; i < n; ++i) {
            if (check < 0 || check >= n || seen[check])
                throw std::logic_error("asymmetric_tsp: broken tour");
            seen[check] = 1;
            check = next[check];
        }
        if (check != 0) throw std::logic_error("asymmetric_tsp: tour does not close");
#endif
        int v = 0;
        reverse_prefix[0] = 0;
        for (int i = 0; i < n; ++i) {
            tour[i] = v;
            position[v] = i;
            prev[next[v]] = v;
            arc_cost[v] = c(v, next[v]);
            arc_reduced[v] = arc_cost[v] - row[v] - col[next[v]];
            if (opt.reversal)
                reverse_prefix[i + 1] = reverse_prefix[i] + c(next[v], v) - arc_cost[v];
            v = next[v];
        }
        current = std::isnan(known_cost) ? value(tour) : known_cost;
    }
    bool accepted() {
        int v = 0;
        for (int &x : work) {
            x = v;
            v = next[v];
        }
        const Cost cost = value(work);
        // Recheck the actual objective: rounded incremental gains must not cause
        // cycling, including when time_limit=infinity and max_trials is finite.
        if (!(cost < current)) {
            for (int i = 0; i < n; ++i) next[tour[i]] = tour[advance(i)];
            return false;
        }
        ++result.stats.moves;
        if (opt.active_nodes) {
            // All changed directed arcs count, including reversed internal arcs.
            for (int i = 0; i < n; ++i) {
                const int a = tour[i], b = tour[advance(i)];
                if (next[a] != b) {
                    active[a] = active[b] = active[next[a]] = 1;
                    for (const auto &e : candidate[a]) active[e.vertex] = 1;
                    for (const auto &e : candidate[b]) active[e.vertex] = 1;
                }
            }
        }
        rebuild(cost);
        remember(tour, current);
        return true;
    }
    Result finish(StopReason why) {
        result.reason = why;
        result.stats.elapsed = std::chrono::duration<double>(Clock::now() - budget.start).count();
        return std::move(result);
    }

    void nearest() {
        std::vector<int> t;
        t.reserve(n);
        std::vector<unsigned char> used(n, 0);
        int a = int(rng() % std::uint64_t(n));
        t.push_back(a);
        used[a] = 1;
        bool stopped = false;
        while (int(t.size()) < n && !stopped) {
            int b = -1;
            Cost best = std::numeric_limits<Cost>::infinity();
            for (int j = 0; j < n; ++j) {
                if (!used[j]) {
                    const Cost w = c(a, j);
                    if (w < best) {
                        best = w;
                        b = j;
                    }
                }
                if (budget.poll()) {
                    stopped = true;
                    break;
                }
            }
            if (!stopped && b >= 0) {
                t.push_back(b);
                used[b] = 1;
                a = b;
            }
        }
        // A partially constructed path can always be completed on a complete graph.
        for (int j = 0; j < n; ++j)
            if (!used[j]) t.push_back(j);
        remember(t, value(t));
    }

    bool potentials() {
        std::vector<Cost> u(n), v(n, std::numeric_limits<Cost>::infinity());
        for (int i = 0; i < n; ++i) {
            Cost minimum = std::numeric_limits<Cost>::infinity();
            for (int j = 0; j < n; ++j) {
                if (i != j) minimum = std::min(minimum, c(i, j));
                if (budget.poll()) return false;
            }
            u[i] = minimum;
            for (int j = 0; j < n; ++j) {
                if (i != j) v[j] = std::min(v[j], c(i, j) - minimum);
                if (budget.poll()) return false;
            }
        }
        row.swap(u);
        col.swap(v);
        bound();
        return true;
    }
    void bound() {
        result.lower_bound = std::accumulate(row.begin(), row.end(), Cost(0)) +
                             std::accumulate(col.begin(), col.end(), Cost(0));
    }

    // Keep useful matched arcs if the assignment phase runs out of its share.
    // Complete free rows/columns greedily. A final forbidden diagonal can always
    // be repaired by swapping its destination with one already matched row.
    std::vector<int> complete_assignment(const std::vector<int> &owner) {
        int matched = 0;
        for (int j = 0; j < n; ++j) matched += owner[j] >= 0;
        // A mostly unfinished matching is too costly to repair on short budgets.
        if (5LL * matched < 4LL * n) return {};
        const auto end = budget.phase((1 + opt.assignment_fraction) * 0.5);
        std::vector<int> successor(n, -1);
        std::vector<unsigned char> used(n, 0);
        for (int j = 0; j < n; ++j)
            if (owner[j] >= 0) {
                successor[owner[j]] = j;
                used[j] = 1;
            }
        for (int i = 0; i < n; ++i)
            if (successor[i] < 0) {
                if (Clock::now() >= end) return {};
                int best = -1;
                Cost key = std::numeric_limits<Cost>::infinity();
                for (int j = 0; j < n; ++j) {
                    if (!used[j] && i != j) {
                        const Cost r = reduced(i, j);
                        if (r < key) {
                            key = r;
                            best = j;
                        }
                    }
                    if (budget.phase_expired(end)) return {};
                }
                if (best >= 0) {
                    successor[i] = best;
                    used[best] = 1;
                    continue;
                }
                int repair = -1;
                for (int j = 0; j < n; ++j) {
                    if (successor[j] >= 0) {
                        const Cost delta = c(j, i) + c(i, successor[j]) - c(j, successor[j]);
                        if (delta < key) {
                            key = delta;
                            repair = j;
                        }
                    }
                    if (budget.phase_expired(end)) return {};
                }
                successor[i] = successor[repair];
                successor[repair] = i;
                used[i] = 1;
            }
        result.stats.assignment_repaired = true;
        return successor;
    }

    // Primal-dual shortest augmenting paths for the assignment relaxation.
    // Dual updates and augmentation are atomic O(n); scans are interruptible.
    std::vector<int> assignment() {
        const auto end = budget.phase(opt.assignment_fraction);
        if (Clock::now() >= end) return {};
        std::vector<int> owner(n + 1, -1), parent(n);
        std::vector<Cost> slack(n);
        std::vector<unsigned char> seen(n + 1);
        for (int root = 0; root < n; ++root) {
            owner[n] = root;
            std::fill(slack.begin(), slack.end(), std::numeric_limits<Cost>::infinity());
            std::fill(seen.begin(), seen.end(), 0);
            int j0 = n;
            do {
                seen[j0] = 1;
                const int i = owner[j0];
                Cost delta = std::numeric_limits<Cost>::infinity();
                int j1 = -1;
                for (int j = 0; j < n; ++j) {
                    if (!seen[j]) {
                        if (j != i) {
                            const Cost r = reduced(i, j);
                            if (r < slack[j]) {
                                slack[j] = r;
                                parent[j] = j0;
                            }
                        }
                        if (slack[j] < delta) {
                            delta = slack[j];
                            j1 = j;
                        }
                    }
                    if (budget.phase_expired(end)) {
                        bound();
                        return complete_assignment(owner);
                    }
                }
                if (j1 < 0) {
                    bound();
                    return {};
                }
                for (int j = 0; j <= n; ++j) {
                    if (seen[j]) {
                        row[owner[j]] += delta;
                        if (j < n) col[j] -= delta;
                    } else if (j < n)
                        slack[j] -= delta;
                }
                j0 = j1;
            } while (owner[j0] >= 0);
            do {
                const int p = parent[j0];
                owner[j0] = owner[p];
                j0 = p;
            } while (j0 != n);
        }
        std::vector<int> successor(n);
        for (int j = 0; j < n; ++j) successor[owner[j]] = j;
        bound();
        result.stats.assignment_complete = true;
        return successor;
    }

    void patch(std::vector<int> successor) {
        if (successor.empty()) return;
        const auto end = budget.phase((1 + opt.assignment_fraction) * 0.5);
        std::vector<int> label(n, -1);
        std::vector<std::vector<int>> parts;
        for (int i = 0; i < n; ++i)
            if (label[i] < 0) {
                const int id = int(parts.size());
                parts.emplace_back();
                int j = i;
                do {
                    parts.back().push_back(j);
                    label[j] = id;
                    j = successor[j];
                } while (j != i);
            }
        using Entry = std::pair<int, int>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> heap;
        for (int i = 0; i < int(parts.size()); ++i) heap.emplace(int(parts[i].size()), i);
        int count = int(parts.size());
        result.stats.assignment_cycles = count;
        while (count > 1) {
            if (Clock::now() >= end) return;
            const auto e = heap.top();
            heap.pop();
            const int id = e.second;
            if (int(parts[id].size()) != e.first) continue;
            int a = -1, b = -1;
            Cost delta = std::numeric_limits<Cost>::infinity();
            for (int x : parts[id])
                for (int y = 0; y < n; ++y) {
                    if (label[y] != id) {
                        const Cost d = c(x, successor[y]) + c(y, successor[x]) -
                                       c(x, successor[x]) - c(y, successor[y]);
                        if (d < delta) {
                            delta = d;
                            a = x;
                            b = y;
                        }
                    }
                    if (budget.phase_expired(end)) return;
                }
            std::swap(successor[a], successor[b]);
            const int other = label[b];
            for (int x : parts[id]) {
                label[x] = other;
                parts[other].push_back(x);
            }
            std::vector<int>().swap(parts[id]); // Release merged storage; O(n) live members.
            heap.emplace(int(parts[other].size()), other);
            --count;
        }
        std::vector<int> t(n);
        int v = 0;
        for (int &x : t) {
            x = v;
            v = successor[v];
        }
        remember(t, value(t));
    }

    // A symmetric relaxation of directed reduced costs. Adding an edge to an
    // MST replaces the maximum-weight edge on its path. This sensitivity keeps
    // bridges between dense groups visible even when their absolute cost is high.
    bool candidate_tree(std::vector<std::vector<std::pair<int, Cost>>> &tree) {
        if (!opt.structural_candidates) return false;
        const auto end = budget.phase(0.75);
        if (Clock::now() >= end) return false;
        std::vector<Cost> key(n, std::numeric_limits<Cost>::infinity());
        std::vector<int> parent(n, -1);
        std::vector<int> pending(n - 1);
        std::iota(pending.begin(), pending.end(), 1);
        int count = n - 1;
        key[0] = 0;
        int a = 0;
        for (int step = 0; step < n; ++step) {
            if (parent[a] >= 0) {
                tree[a].emplace_back(parent[a], key[a]);
                tree[parent[a]].emplace_back(a, key[a]);
            }
            int following = -1;
            for (int j = 0; j < count; ++j) {
                const int v = pending[j];
                const Cost w = std::min(reduced(a, v), reduced(v, a));
                if (w < key[v]) {
                    key[v] = w;
                    parent[v] = a;
                }
                if (following < 0 || key[v] < key[pending[following]]) following = j;
                if (budget.phase_expired(end)) return false;
            }
            if (following >= 0) {
                a = pending[following];
                std::move(pending.begin() + following + 1, pending.begin() + count,
                          pending.begin() + following);
                --count;
            }
        }
        return true;
    }

    // Contract the mandatory in/out pair edges of an ATSP transformation.
    // Removing out(0) leaves in(0) and n-1 pairs; its tree has incoming arcs
    // only at vertex 0. Add one outgoing root arc to complete the fixed 1-tree.
    bool ascent_tree(std::vector<std::vector<std::pair<int, Cost>>> &tree) {
        if (!opt.ascent_candidates) return false;
        const Cost base = std::accumulate(row.begin(), row.end(), Cost(0)) +
                          std::accumulate(col.begin(), col.end(), Cost(0));
        const auto end = budget.phase(0.45);
        if (Clock::now() >= end) return false;
        std::vector<Cost> po(n, 0), pi(n, 0), best_out(n), best_in(n), key(n);
        std::vector<int> parent(n), edge_from(n), edge_to(n), dout(n), din(n);
        std::vector<int> pending(n - 1);
        std::vector<std::vector<std::pair<int, Cost>>> best_tree(n);
        Cost best = -std::numeric_limits<Cost>::infinity(), best_root = 0;
        Cost multiplier = 1;
        int stalled = 0;
        for (int iteration = 0; iteration < 80; ++iteration) {
            if (Clock::now() >= end) break;
            std::fill(key.begin(), key.end(), std::numeric_limits<Cost>::infinity());
            std::fill(parent.begin(), parent.end(), -1);
            std::iota(pending.begin(), pending.end(), 1);
            int count = n - 1;
            std::fill(dout.begin(), dout.end(), 0);
            std::fill(din.begin(), din.end(), 0);
            key[0] = 0;
            Cost total = 0, root = std::numeric_limits<Cost>::infinity();
            int root_head = -1;
            bool complete = true;
            int a = 0;
            for (int step = 0; step < n && complete; ++step) {
                total += key[a];
                if (parent[a] >= 0) {
                    ++dout[edge_from[a]];
                    ++din[edge_to[a]];
                }
                int following = -1;
                for (int j = 0; j < count; ++j) {
                    const int v = pending[j];
                    Cost w = reduced(v, a) + po[v] + pi[a];
                    int from = v, to = a;
                    if (a != 0) {
                        const Cost other = reduced(a, v) + po[a] + pi[v];
                        if (other < w) {
                            w = other;
                            from = a;
                            to = v;
                        }
                    }
                    if (w < key[v]) {
                        key[v] = w;
                        parent[v] = a;
                        edge_from[v] = from;
                        edge_to[v] = to;
                    }
                    if (following < 0 || key[v] < key[pending[following]]) following = j;
                    if (budget.phase_expired(end)) {
                        complete = false;
                        break;
                    }
                }
                if (following >= 0) {
                    a = pending[following];
                    // Ascending vertices preserve the original tie order and locality.
                    std::move(pending.begin() + following + 1, pending.begin() + count,
                              pending.begin() + following);
                    --count;
                }
            }
            if (!complete) break;
            for (int v = 1; v < n; ++v) {
                const Cost w = reduced(0, v) + po[0] + pi[v];
                if (w < root) {
                    root = w;
                    root_head = v;
                }
                if (budget.phase_expired(end)) {
                    complete = false;
                    break;
                }
            }
            if (!complete) break;
            ++result.stats.ascent_iterations;
            ++dout[0];
            ++din[root_head];
            const Cost lower = total + root - std::accumulate(po.begin(), po.end(), Cost(0)) -
                               std::accumulate(pi.begin(), pi.end(), Cost(0));
            if (lower > best) {
                best = lower;
                best_out = po;
                best_in = pi;
                best_root = root;
                for (auto &edges : best_tree) edges.clear();
                for (int v = 1; v < n; ++v) {
                    best_tree[v].emplace_back(parent[v], key[v]);
                    best_tree[parent[v]].emplace_back(v, key[v]);
                }
                stalled = 0;
            } else if (++stalled == 8) {
                multiplier *= Cost(0.5);
                stalled = 0;
            }
            Cost norm = 0;
            for (int v = 0; v < n; ++v)
                norm += Cost(dout[v] - 1) * (dout[v] - 1) + Cost(din[v] - 1) * (din[v] - 1);
            if (norm == 0 || !(result.cost - base > lower)) break;
            const Cost step = multiplier * (result.cost - base - lower) / norm;
            for (int v = 0; v < n; ++v) {
                po[v] += step * (dout[v] - 1);
                pi[v] += step * (din[v] - 1);
            }
        }
        if (!std::isfinite(best)) return false;
        candidate_out = std::move(best_out);
        candidate_in = std::move(best_in);
        root_candidate = best_root;
        tree = std::move(best_tree);
        one_tree = true;
        result.stats.one_tree_candidates = true;
        result.lower_bound = std::max(result.lower_bound, base + best);
        return true;
    }

    bool make_candidates() {
        const int k = std::min(n - 1, opt.candidates);
        std::vector<std::vector<std::pair<int, Cost>>> tree(n);
        const bool structural = ascent_tree(tree) || candidate_tree(tree);
        result.stats.tree_candidates = structural;
        std::vector<Cost> bottleneck(n);
        std::vector<int> stack, parent(n);
        stack.reserve(n);
        using Edge = std::pair<Cost, int>;
        std::vector<Edge> heap;
        heap.reserve(k);
        for (int i = 0; i < n; ++i) {
            if (structural) {
                stack.clear();
                stack.push_back(i);
                parent[i] = -1;
                bottleneck[i] = -std::numeric_limits<Cost>::infinity();
                for (std::size_t t = 0; t < stack.size(); ++t) {
                    const int a = stack[t];
                    for (const auto &e : tree[a])
                        if (e.first != parent[a]) {
                            parent[e.first] = a;
                            bottleneck[e.first] = std::max(bottleneck[a], e.second);
                            stack.push_back(e.first);
                        }
                    if (budget.poll()) return false;
                }
            }
            heap.clear();
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    const Cost rank = one_tree ? reduced(i, j) + candidate_out[i] + candidate_in[j] -
                                                     (i == 0 ? root_candidate : bottleneck[j])
                                               : reduced(i, j) - (structural ? bottleneck[j] : 0);
                    const Edge e{rank, j};
                    if (int(heap.size()) < k) {
                        heap.push_back(e);
                        std::push_heap(heap.begin(), heap.end());
                    } else if (e < heap.front()) {
                        std::pop_heap(heap.begin(), heap.end());
                        heap.back() = e;
                        std::push_heap(heap.begin(), heap.end());
                    }
                }
                if (budget.poll()) return false;
            }
            std::sort_heap(heap.begin(), heap.end());
            candidate[i].reserve(k);
            for (const auto &e : heap) {
                const Cost cost = c(i, e.second);
                candidate[i].push_back({e.second, cost, cost - row[i] - col[e.second]});
            }
        }
        result.stats.candidates_complete = true;
        return true;
    }

    // Move a forward block of 1..3 vertices after a candidate predecessor.
    bool relocate(int a) {
        int b = a;
        const int before = prev[a];
        for (int length = 1; length <= 3 && length < n - 1; ++length, b = next[b]) {
            const int after = next[b];
            // Choose the new b -> y arc, inserting the block just before y.
            for (const auto &e : candidate[b]) {
                const int y = e.vertex;
                const int x = prev[y];
                bool overlap = x == before;
                int z = a;
                for (int t = 0; t < length; ++t, z = next[z]) overlap |= x == z;
                if (overlap) continue;
                const Cost old = arc_cost[before] + arc_cost[b] + arc_cost[x];
                const Cost fresh = c(before, after) + c(x, a) + e.cost;
                if (improvement(old - fresh, std::abs(old) + std::abs(fresh))) {
                    next[before] = after;
                    next[x] = a;
                    next[b] = y;
                    if (accepted()) return true;
                }
                if (budget.poll()) return false;
            }
        }
        return false;
    }

    bool reverse(int a) {
        const int i = position[a], b = next[a];
        for (const auto &e : candidate[a]) {
            const int x = e.vertex;
            const int j = position[x];
            if (x == a || x == b) continue;
            const int y = next[x];
            const Cost old = arc_cost[a] + arc_cost[x];
            const Cost fresh = e.cost + c(b, y);
            const Cost internal =
                reverse_prefix[j] - reverse_prefix[i + 1] + (j <= i ? reverse_prefix[n] : 0);
            if (improvement(old - fresh - internal,
                            std::abs(old) + std::abs(fresh) + std::abs(internal))) {
                int v = b, after = y;
                while (v != y) {
                    const int old_next = next[v];
                    next[v] = after;
                    after = v;
                    v = old_next;
                }
                next[a] = x;
                if (accepted()) return true;
            }
            if (budget.poll()) return false;
        }
        return false;
    }

    // Compress each unchanged path between cuts to one edge. After assigning
    // tail[i] -> head[(i+1)%k], follow these k edges to reject subtours in O(k^2).
    bool single_cycle(int k) const {
        std::array<int, 31> following{};
        for (int i = 0; i < k; ++i) {
            int best = n + 1, which = -1;
            for (int j = 0; j < k; ++j) {
                int d = position[tails[j]] - position[tails[i]];
                if (d <= 0) d += n;
                if (d < best) {
                    best = d;
                    which = j;
                }
            }
            following[i] = which;
        }
        int v = 0;
        for (int step = 1; step <= k; ++step) {
            v = following[v + 1 == k ? 0 : v + 1];
            if (v == 0) return step == k;
        }
        return false;
    }

    bool exchange(int depth, Cost gain, Cost scale) {
        if (nodes_left-- <= 0 || budget.poll()) return false;
        ++result.stats.search_nodes;
        const int a = tails[depth - 1];
        struct CutMark {
            unsigned char &used;
            explicit CutMark(unsigned char &flag) : used(flag) { used = 1; }
            ~CutMark() { used = 0; }
        } mark(cut_used[a]);
        const Cost removed = arc_reduced[a];
        if (depth >= 3 && (depth & 1) && a != heads[0]) {
            const Cost closing = reduced(a, heads[0]);
            const Cost total = gain + removed - closing;
            const Cost s = scale + std::abs(removed) + std::abs(closing);
            if (improvement(total, s) && single_cycle(depth)) {
                for (int i = 0; i < depth; ++i) next[tails[i]] = heads[i + 1 == depth ? 0 : i + 1];
                if (accepted()) return true;
            }
        }
        if (depth == exchange_limit || nodes_left <= 0) return false;
        auto &choices = branches[depth - 1];
        choices.clear();
        for (const auto &e : candidate[a]) {
            const int h = e.vertex;
            if (budget.poll()) return false;
            const int t = prev[h];
            if (cut_used[t]) continue;
            const Cost added = e.reduced;
            if (gain + removed - added <
                -opt.relative_tolerance * (scale + std::abs(removed) + std::abs(added)))
                continue;
            choices.push_back({opt.lookahead ? added - arc_reduced[t] : added, added, h});
        }
        const int limit =
            std::min(int(choices.size()), depth == 1 ? opt.candidates : opt.search_breadth);
        std::partial_sort(choices.begin(), choices.begin() + limit, choices.end());
        for (int i = 0; i < limit; ++i) {
            const int h = choices[i].head;
            const Cost added = choices[i].added;
            tails[depth] = prev[h];
            heads[depth] = h;
            if (exchange(depth + 1, gain + removed - added,
                         scale + std::abs(removed) + std::abs(added)))
                return true;
            if (nodes_left <= 0) break;
        }
        return false;
    }

    void local_search() {
        exchange_limit = opt.staged_search ? std::min(5, opt.max_exchange) : opt.max_exchange;
        bool changed = true;
        while (changed && !budget.expired()) {
            changed = false;
            std::shuffle(order.begin(), order.end(), rng);
            for (int a : order) {
                if (budget.expired()) return;
                if (opt.active_nodes && !active[a]) continue;
                if (relocate(a) || (opt.reversal && reverse(a))) {
                    changed = true;
                    continue;
                }
                nodes_left = opt.search_nodes;
                tails[0] = a;
                heads[0] = next[a];
                if (exchange(1, 0, 0))
                    changed = true;
                else
                    active[a] = 0;
            }
        }
    }

    // Equal directed endpoints and equal vertex sets permit exchanging paths
    // without changing their boundary arcs or duplicating any vertex.
    bool merge() {
        if (!opt.merge_tours || !(result.cost < current)) return false;
        merge_difference.assign(n, 0);
        if (merge_version != result.stats.improvements) {
            best_position.resize(n);
            best_arc.resize(n);
            for (int i = 0; i < n; ++i) {
                best_position[result.tour[i]] = i;
                best_arc[i] = c(result.tour[i], result.tour[advance(i)]);
                if (budget.poll()) return false;
            }
            merge_version = result.stats.improvements;
        }
        const int limit = std::min(n - 1, 256);
        for (int a : order) {
            if (budget.expired()) return false;
            const int i = position[a], j = best_position[a];
            int mismatch = 0, bi = j, ci = i, new_head = a, length = 0;
            Cost old_cost = 0, new_cost = 0;
            auto mark = [&](int v, int d) {
                mismatch -= merge_difference[v] != 0;
                merge_difference[v] += d;
                mismatch += merge_difference[v] != 0;
            };
            for (length = 1; length <= limit; ++length) {
                old_cost += best_arc[bi];
                new_cost += arc_cost[new_head];
                bi = advance(bi);
                ci = advance(ci);
                const int b = result.tour[bi], d = tour[ci];
                mark(b, 1);
                mark(d, -1);
                new_head = d;
                if (!mismatch && b == d &&
                    improvement(old_cost - new_cost, std::abs(old_cost) + std::abs(new_cost))) {
                    auto combined = result.tour;
                    int bp = j, cp = i;
                    for (int t = 1; t < length; ++t) {
                        bp = advance(bp);
                        cp = advance(cp);
                        combined[bp] = tour[cp];
                    }
                    // Match rebuild()'s summation order before testing progress.
                    std::rotate(combined.begin(), std::find(combined.begin(), combined.end(), 0),
                                combined.end());
                    const Cost cost = value(combined);
                    if (cost < result.cost) {
                        set_tour(combined);
                        remember(tour, current);
                        ++result.stats.merging_moves;
                        ++result.stats.moves;
                        return true;
                    }
                }
                if (budget.poll()) return false;
            }
            bi = j;
            ci = i;
            for (int t = 1; t <= limit; ++t) {
                bi = advance(bi);
                ci = advance(ci);
                merge_difference[result.tour[bi]] = 0;
                merge_difference[tour[ci]] = 0;
            }
        }
        return false;
    }

    bool intensify() {
        if (!opt.staged_search || opt.max_exchange <= 5) return false;
        if (opt.adaptive_depth &&
            (result.stats.assignment_cycles == 0 || 16LL * result.stats.assignment_cycles > n))
            return false;
        exchange_limit = opt.max_exchange;
        for (int a : order) {
            if (budget.expired()) return false;
            nodes_left = opt.deep_search_nodes;
            tails[0] = a;
            heads[0] = next[a];
            if (exchange(1, 0, 0)) return true;
        }
        return false;
    }

    // Two successor swaps: the first splits the tour, the second reconnects
    // different components. Four distinct tails give a nonsequential 4-opt.
    bool patch_exchange() {
        if (!opt.patching) return false;
        for (int a : order)
            for (const auto &e : candidate[a]) {
                const int h = e.vertex;
                if (budget.poll()) return false;
                const int b = prev[h];
                if (a == b || next[a] == b || next[b] == a) continue;
                const Cost old = arc_reduced[a] + arc_reduced[b];
                const Cost added = e.reduced + reduced(b, next[a]);
                const Cost gain = old - added;
                if (!improvement(gain, std::abs(old) + std::abs(added))) continue;
                int lo = position[a], hi = position[b];
                if (lo > hi) std::swap(lo, hi);
                const bool inner = hi - lo <= n / 2;
                const int begin = inner ? lo + 1 : hi + 1;
                const int count = inner ? hi - lo : n - (hi - lo);
                for (int step = 0; step < count; ++step) {
                    const int c0 = tour[(begin + step) % n];
                    if (c0 == a || c0 == b) continue;
                    for (const auto &edge : candidate[c0]) {
                        const int head = edge.vertex;
                        if (budget.poll()) return false;
                        const int d = prev[head], p = position[d];
                        if (d == a || d == b || ((p > lo && p <= hi) == inner)) continue;
                        const Cost removed = arc_reduced[c0] + arc_reduced[d];
                        const Cost fresh = edge.reduced + reduced(d, next[c0]);
                        if (improvement(gain + removed - fresh, std::abs(old) + std::abs(added) +
                                                                    std::abs(removed) +
                                                                    std::abs(fresh))) {
                            std::swap(next[a], next[b]);
                            std::swap(next[c0], next[d]);
                            if (accepted()) {
                                ++result.stats.patching_moves;
                                return true;
                            }
                        }
                        if (budget.poll()) return false;
                    }
                }
            }
        return false;
    }

    void kick(int strength) {
        for (int r = 0; r < strength; ++r) {
            // Four nonempty paths, reordered A,D,C,B, without any reversal.
            std::array<int, 4> cuts{};
            for (int j = 0; j < 4; ++j) {
                bool duplicate;
                do {
                    cuts[j] = int(rng() % std::uint64_t(n));
                    duplicate = false;
                    for (int i = 0; i < j; ++i) duplicate |= cuts[i] == cuts[j];
                } while (duplicate);
                if (opt.guided_kicks && result.stats.trials % 2 == 0 && r == 0) {
                    Cost score = arc_reduced[tour[cuts[j]]];
                    for (int sample = 0; sample < 64; ++sample) {
                        const int p = int(rng() % std::uint64_t(n));
                        bool used = false;
                        for (int i = 0; i < j; ++i) used |= cuts[i] == p;
                        if (used) continue;
                        const Cost w = arc_reduced[tour[p]];
                        if (w > score) {
                            score = w;
                            cuts[j] = p;
                        }
                    }
                }
            }
            std::sort(cuts.begin(), cuts.end());
            const int a = tour[cuts[0]], b = tour[cuts[1]], d = tour[cuts[2]], e = tour[cuts[3]];
            std::swap(next[a], next[d]);
            std::swap(next[b], next[e]);
            rebuild();
        }
    }

    void restart() {
        std::vector<unsigned char> used(n, 0);
        std::vector<int> t;
        t.reserve(n);
        int a = int(rng() % std::uint64_t(n));
        while (int(t.size()) < n && !budget.expired()) {
            t.push_back(a);
            used[a] = 1;
            int best = -1, second = -1;
            Cost first_cost = std::numeric_limits<Cost>::infinity(), second_cost = first_cost;
            for (const auto &e : candidate[a])
                if (!used[e.vertex]) {
                    const int b = e.vertex;
                    const Cost w = e.reduced;
                    if (w < first_cost) {
                        second = best;
                        second_cost = first_cost;
                        best = b;
                        first_cost = w;
                    } else if (w < second_cost) {
                        second = b;
                        second_cost = w;
                    }
                }
            if (best < 0) {
                for (int b = 0; b < n; ++b) {
                    if (!used[b]) {
                        const Cost w = reduced(a, b);
                        if (w < first_cost) {
                            best = b;
                            first_cost = w;
                        }
                    }
                    if (budget.poll()) break;
                }
            }
            if (best < 0) break;
            a = second >= 0 && rng() % 5 == 0 ? second : best;
        }
        for (int b = 0; b < n; ++b)
            if (!used[b]) t.push_back(b);
        set_tour(t);
        remember(tour, current);
    }

public:
    Engine(int vertices, Distance &d, const Options &o, Clock::time_point start)
        : n(vertices), edge_limit(std::numeric_limits<Cost>::max() /
                                  (64 * Cost(std::max(vertices, 1)) * Cost(std::max(vertices, 1)))),
          distance(d), opt(o), budget(start, o), rng(o.seed) {}

    Result run() {
        if (n < 0 || opt.candidates < 1 || opt.max_exchange < 3 || opt.max_exchange > 31 ||
            opt.search_breadth < 1 || opt.search_nodes < 1 || opt.deep_search_nodes < 1 ||
            !std::isfinite(opt.assignment_fraction) || opt.assignment_fraction < 0 ||
            opt.assignment_fraction > 1 || !std::isfinite(opt.relative_tolerance) ||
            opt.relative_tolerance < 0)
            throw std::invalid_argument("asymmetric_tsp: invalid option or size");
        result.tour.resize(n);
        std::iota(result.tour.begin(), result.tour.end(), 0);
        if (!opt.initial_tour.empty()) {
            if (opt.initial_tour.size() != std::size_t(n))
                throw std::invalid_argument("asymmetric_tsp: initial_tour size");
            std::vector<unsigned char> seen(n, 0);
            for (int v : opt.initial_tour) {
                if (v < 0 || v >= n || seen[v])
                    throw std::invalid_argument(
                        "asymmetric_tsp: initial_tour is not a permutation");
                seen[v] = 1;
            }
            result.tour = opt.initial_tour;
        }
        result.cost = value(result.tour);
        if (n <= 2) {
            result.lower_bound = result.cost;
            return finish(StopReason::trivial);
        }
        if (budget.expired()) return finish(StopReason::time_limit);
        if (opt.max_trials == 0) return finish(StopReason::trial_limit);
        if (n == 3) {
            auto alternative = result.tour;
            std::swap(alternative[1], alternative[2]);
            remember(alternative, value(alternative));
            result.lower_bound = result.cost;
            return finish(StopReason::trivial);
        }
        nearest();
        if (budget.expired()) return finish(StopReason::time_limit);
        row.resize(n);
        col.resize(n);
        if (!potentials()) return finish(StopReason::time_limit);
        if (opt.assignment) patch(assignment());
        if (budget.expired()) return finish(StopReason::time_limit);
        candidate.resize(n);
        if (!make_candidates()) return finish(StopReason::time_limit);
        tour.resize(n);
        next.resize(n);
        prev.resize(n);
        position.resize(n);
        work.resize(n);
        reverse_prefix.resize(std::size_t(n) + 1);
        order.resize(n);
        std::iota(order.begin(), order.end(), 0);
        active.assign(n, 1);
        cut_used.assign(n, 0);
        arc_cost.resize(n);
        arc_reduced.resize(n);
        for (auto &b : branches) b.reserve(std::min(n - 1, opt.candidates));
        set_tour(result.tour);
        std::vector<int> home = tour;
        Cost home_cost = current;
        std::uint64_t failures = 0;
        while (!budget.expired() && result.stats.trials < opt.max_trials) {
            const Cost previous_best = result.cost;
            do {
                local_search();
            } while (!budget.expired() && (patch_exchange() || merge() || intensify()));
            ++result.stats.trials;
            if (current < home_cost || (result.stats.trials % 8 == 0)) {
                home = tour;
                home_cost = current;
            }
            failures = result.cost < previous_best ? 0 : failures + 1;
            if (budget.expired() || result.stats.trials >= opt.max_trials) break;
            if (opt.restarts && result.stats.trials % 12 == 0) {
                restart();
                home = tour;
                home_cost = current;
                continue;
            }
            if (failures % 7 == 0)
                set_tour(result.tour);
            else
                set_tour(home);
            kick(1 + int(std::min<std::uint64_t>(failures / 12, 5)));
        }
        return finish(budget.expired() ? StopReason::time_limit : StopReason::trial_limit);
    }
};
} // namespace detail

// distance(i,j) must be a deterministic, fast callable with finite costs for i!=j.
// n=0/1 has cost zero and does not read any self-loop. No matrix is copied.
template <class Distance> Result solve(int n, Distance &&distance, const Options &options = {}) {
    const auto start = Clock::now();
    detail::Engine<typename std::remove_reference<Distance>::type> engine(
        n, distance, options, start);
    return engine.run();
}

template <class T>
Result solve(const std::vector<std::vector<T>> &matrix, const Options &options = {}) {
    const auto start = Clock::now();
    if (matrix.size() > std::size_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("asymmetric_tsp: matrix too large");
    const int n = int(matrix.size());
    for (const auto &row : matrix)
        if (row.size() != matrix.size())
            throw std::invalid_argument("asymmetric_tsp: matrix must be square");
    auto distance = [&matrix](int i, int j) { return matrix[i][j]; };
    detail::Engine<decltype(distance)> engine(n, distance, options, start);
    return engine.run();
}
} // namespace asymmetric_tsp
#endif
