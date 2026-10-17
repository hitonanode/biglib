#ifndef SYMMETRIC_TSP_HPP
#define SYMMETRIC_TSP_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace symmetric_tsp {
using Clock = std::chrono::steady_clock;
using Cost = long double;
enum class StopReason { time_limit, trial_limit, trivial };

struct Options {
    double time_limit = 1.0;
    Clock::time_point deadline = Clock::time_point::max();
    double safety_margin = 0.002;
    std::uint64_t seed = 1;
    std::uint64_t max_trials = std::numeric_limits<std::uint64_t>::max(); // Per run.
    int candidates = 24;
    int search_candidates = 8; // Leading candidates considered by deep exchanges.
    int max_depth = 10;        // 2..20 deleted edges, distinct endpoints.
    int search_breadth = 5;
    int search_nodes = 160;     // Per origin in the deep phase.
    int ascent_iterations = 80; // 0 disables 1-tree/alpha preprocessing.
    bool sparse_ascent = true;  // Verify the final penalties on the complete graph.
    bool relocate = true;
    bool kicks = true;
    bool guided_kicks = false;     // Alternate local candidate-guided and global bridges.
    bool merge_tours = true;       // Exact equal-vertex-set path transfers.
    bool chained_rotations = true; // Gain-positive path rotations, including shared endpoints.
    bool patching = true;          // Nonsequential 4-opt across two cycles.
    bool partition_seed = true;    // Distance-only balanced partition, no metric assumption.
    bool active_search = true;     // Revisit deep origins near changed tour/candidate edges.
    bool staged_search = true;     // Search through depth 5 before deeper exchanges.
    Cost relative_tolerance = 1e-12L;
    std::vector<int> initial_tour;               // A partial ordering is allowed; no duplicates.
    std::vector<std::vector<int>> initial_paths; // Disjoint paths; alternative to initial_tour.
};

struct Statistics {
    double elapsed = 0; // This call; destruction after return excluded.
    std::uint64_t trials = 0, moves = 0, improvements = 0, search_nodes = 0, kicks = 0;
    std::uint64_t two_opt_moves = 0, relocate_moves = 0, exchange_moves = 0;
    std::uint64_t merging_moves = 0, patching_moves = 0, chained_moves = 0;
    int nearest_rows = 0, alpha_rows = 0, ascent_iterations = 0;
};
struct Result {
    std::vector<int> tour;
    Cost cost = 0;
    Cost lower_bound = -std::numeric_limits<Cost>::infinity();
    StopReason reason = StopReason::time_limit;
    Statistics stats; // Counters are cumulative across run() calls.
};

namespace detail {
inline Clock::time_point deadline(Clock::time_point start, double seconds, Clock::time_point outer) {
    if (std::isnan(seconds) || seconds < 0)
        throw std::invalid_argument("symmetric_tsp: invalid time_limit");
    const double range = std::chrono::duration<double>(Clock::time_point::max() - start).count();
    if (seconds < range)
        outer = std::min(outer, start + std::chrono::duration_cast<Clock::duration>(
                                            std::chrono::duration<double>(seconds)));
    return outer;
}
struct Budget {
    Clock::time_point start, stop;
    unsigned ticks = 0;
    Budget(Clock::time_point s, double seconds, Clock::time_point outer, double margin)
        : start(s) {
        stop = deadline(s, seconds, outer);
        if (stop <= s) {
            stop = s;
            return;
        }
        if (stop != Clock::time_point::max()) {
            const double remaining = std::chrono::duration<double>(stop - s).count();
            stop -= std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(std::min(margin, remaining * 0.1)));
        }
    }
    bool expired() const { return stop != Clock::time_point::max() && Clock::now() >= stop; }
    bool poll() { return (++ticks & 127U) == 0 && expired(); }
    bool poll(Clock::time_point until) { return (++ticks & 127U) == 0 && Clock::now() >= until; }
    Clock::time_point phase(double fraction) const {
        return stop == Clock::time_point::max()
                   ? stop
                   : start + std::chrono::duration_cast<Clock::duration>((stop - start) * fraction);
    }
};
} // namespace detail

// Distance is stored by value. A lambda referring to a matrix does not copy the matrix.
// The distance function and referenced input must remain unchanged for the session.
template <class Distance> class Solver {
    static constexpr int depth_limit = 20;
    static constexpr Cost infinity = std::numeric_limits<Cost>::infinity();
    struct Edge {
        int vertex;
        Cost cost, rank;
        bool operator<(const Edge &b) const {
            if (rank != b.rank) return rank < b.rank;
            if (cost != b.cost) return cost < b.cost;
            return vertex < b.vertex;
        }
    };
    struct Branch {
        int c, d;
        Cost gain, rank;
        bool operator<(const Branch &b) const {
            return rank != b.rank ? rank < b.rank : (c != b.c ? c < b.c : d < b.d);
        }
    };
    int n_;
    mutable Distance distance_;
    Options options_;
    std::uint64_t random_;
    bool initialized_ = false, seeded_ = false, alpha_ready_ = false, ascent_done_ = false;
    bool best_is_current_ = true, exact_integral_ = false;
    std::vector<int> tour_, best_, position_, next_, prev_, work_;
    std::vector<unsigned char> active_, deep_active_, used_;
    std::vector<std::vector<Edge>> candidates_;
    std::vector<Edge> row_heap_;
    int nearest_row_ = 0, nearest_column_ = 0, alpha_row_ = 0, alpha_column_ = 0;
    std::vector<Cost> pi_, best_pi_;
    Cost current_ = 0, best_cost_ = 0, bound_ = -infinity;
    Statistics stats_;
    // An interrupted Prim scan is retained, including its partially updated keys.
    struct Tree {
        std::vector<int> parent, degree, remaining;
        std::vector<Cost> key;
        int vertex = 1, scan = 0, next_index = -1, root_scan = 1, root_a = -1, root_b = -1;
        Cost sum = 0, first = infinity, second = infinity;
        bool started = false, mst_done = false;
    } tree_;
    std::vector<std::vector<Edge>> ascent_graph_;
    std::vector<Edge> ascent_heap_;
    std::vector<unsigned char> tree_seen_;
    std::vector<Cost> sparse_pi_, previous_gradient_;
    std::vector<int> ascent_next_;
    int graph_row_ = 0, sparse_count_ = 0, sparse_stall_ = 0;
    bool sparse_tree_started_ = false, sparse_finished_ = false;
    Cost sparse_bound_ = -infinity, sparse_step_scale_ = 1;
    std::vector<int> best_parent_;
    int best_root_a_ = -1, best_root_b_ = -1, stalled_ascent_ = 0;
    Cost best_root_second_ = 0, lambda_ = 1;
    std::vector<std::vector<std::pair<int, Cost>>> adjacency_;
    struct Walk {
        int vertex, parent;
        Cost maximum;
    };
    std::vector<Walk> stack_;
    std::vector<Cost> beta_;
    bool alpha_walk_started_ = false, alpha_walk_done_ = false;
    int origin_ = 0, stage_ = 0;
    bool pass_changed_ = false, pending_kick_ = false, locally_refined_ = false;
    std::uint64_t stagnation_ = 0;
    int merge_origin_ = 0;
    bool merge_changed_current_ = false;
    std::array<int, depth_limit> x_{}, y_{};
    std::array<std::vector<Branch>, depth_limit> branches_;
    int nodes_left_ = 0, exchange_depth_ = 0, patches_left_ = 0;

    std::uint64_t random() {
        // SplitMix64 integer mixer; no dependency on a random-number library.
        std::uint64_t z = (random_ += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    Cost cost(int a, int b) const {
        const Cost v = static_cast<Cost>(distance_(a, b));
        using Input = typename std::decay<decltype(distance_(a, b))>::type;
        if constexpr (std::is_integral<Input>::value && std::numeric_limits<Cost>::max_exponent >
                                                            std::numeric_limits<Input>::digits + 80)
            return v;
        const Cost limit = std::numeric_limits<Cost>::max() /
                           (1024 * static_cast<Cost>(std::max(1, n_)) * std::max(1, n_));
        if (!std::isfinite(v) || std::abs(v) > limit)
            throw std::invalid_argument("symmetric_tsp: non-finite or excessive cost");
        return v;
    }
    Cost transformed(int a, int b) const { return cost(a, b) + best_pi_[a] + best_pi_[b]; }
    int after(int i) const { return i + 1 == n_ ? 0 : i + 1; }
    Cost value(const std::vector<int> &t) const {
        if (n_ < 2) return 0;
        Cost sum = 0, correction = 0;
        for (int i = 0; i < n_; ++i) {
            const Cost z = cost(t[i], t[after(i)]) - correction;
            const Cost s = sum + z;
            correction = (s - sum) - z;
            sum = s;
        }
        return sum;
    }
    bool improves(Cost gain, Cost scale) const {
        return gain > options_.relative_tolerance * scale;
    }
    void index() {
        for (int i = 0; i < n_; ++i) {
            position_[tour_[i]] = i;
            next_[tour_[i]] = tour_[after(i)];
            prev_[tour_[after(i)]] = tour_[i];
        }
    }
    void remember() {
        merge_origin_ = 0;
        if (current_ < best_cost_) {
            best_is_current_ = true;
            best_cost_ = current_;
            ++stats_.improvements;
            stagnation_ = 0;
        }
    }
    void save_best() {
        if (best_is_current_) {
            best_ = tour_;
            best_is_current_ = false;
        }
    }
    void activate(int v) {
        active_[v] = active_[next_[v]] = active_[prev_[v]] = 1;
        deep_active_[v] = deep_active_[next_[v]] = deep_active_[prev_[v]] = 3;
        for (const auto &e : candidates_[v]) {
            active_[e.vertex] = 1;
            deep_active_[e.vertex] = 3;
        }
    }
    bool accept(Cost exact_cost = infinity) {
#ifdef SYMMETRIC_TSP_CHECK_INVARIANTS
        std::vector<unsigned char> seen(n_, 0);
        if (work_.size() != static_cast<std::size_t>(n_))
            throw std::logic_error("symmetric_tsp: incomplete exchange");
        for (int v : work_) {
            if (v < 0 || v >= n_ || seen[v]++)
                throw std::logic_error("symmetric_tsp: invalid exchange");
        }
#endif
        const Cost v = exact_integral_ && exact_cost != infinity ? exact_cost : value(work_);
#ifdef SYMMETRIC_TSP_CHECK_INVARIANTS
        if (exact_integral_ && v != value(work_))
            throw std::logic_error("symmetric_tsp: incorrect incremental cost");
#endif
        if (!(v < current_)) return false; // Re-evaluate to reject rounded phantom gains.
        for (int i = 0; i < n_; ++i) {
            const int a = work_[i], b = work_[after(i)];
            if (next_[a] != b && prev_[a] != b) {
                activate(a);
                activate(b);
            }
        }
        tour_.swap(work_);
        current_ = v;
        index();
        remember();
        ++stats_.moves;
        return true;
    }
    void initialize() {
        using Input = typename std::decay<decltype(distance_(0, 0))>::type;
        if constexpr (std::is_integral<Input>::value) {
            int bits = 0;
            for (unsigned v = static_cast<unsigned>(std::max(n_, 2 * depth_limit)); v; v >>= 1)
                ++bits;
            exact_integral_ =
                std::numeric_limits<Input>::digits + bits + 2 <= std::numeric_limits<Cost>::digits;
        }
        std::vector<unsigned char> seen(n_, 0);
        auto append = [&](int v) {
            if (v < 0 || v >= n_ || seen[v]++)
                throw std::invalid_argument("symmetric_tsp: invalid or duplicate initial vertex");
            tour_.push_back(v);
        };
        tour_.reserve(n_);
        for (int v : options_.initial_tour) append(v);
        for (const auto &p : options_.initial_paths)
            for (int v : p) append(v);
        for (int v = 0; v < n_; ++v)
            if (!seen[v]) tour_.push_back(v);
        best_ = tour_;
        current_ = best_cost_ = value(tour_);
        position_.resize(n_);
        next_.resize(n_);
        prev_.resize(n_);
        work_.reserve(n_);
        active_.assign(n_, 1);
        deep_active_.assign(n_, 3);
        used_.assign(n_, 0);
        candidates_.resize(n_);
        pi_.assign(n_, 0);
        best_pi_.assign(n_, 0);
        beta_.resize(n_);
        index();
        initialized_ = true;
    }
    void
    seed_candidates(const std::vector<int> &order, detail::Budget &budget, double fraction = 0.32) {
        const auto until = budget.phase(fraction);
        const int span = std::min(n_ - 1, std::max(16, options_.candidates * 2));
        for (int i = 0; i < n_; ++i) {
            std::vector<Edge> heap;
            for (int offset = 1; offset <= span; ++offset) {
                if (budget.poll(until)) return;
                const int shift = offset & 1 ? (offset + 1) / 2 : -offset / 2;
                const int j = (i + shift + n_) % n_;
                const Cost c = cost(order[i], order[j]);
                heap_edge(heap, {order[j], c, c});
            }
            std::sort(heap.begin(), heap.end());
            candidates_[order[i]].swap(heap);
        }
    }
    void partition_seed(detail::Budget &budget) {
        // Balanced, distance-only bisections give a cheap locality ordering. This
        // is merely a construction heuristic, including on nonmetric inputs.
        std::vector<int> order = tour_;
        std::vector<Cost> score(n_);
        std::vector<std::pair<int, int>> tasks{{0, n_}};
        std::vector<std::pair<int, int>> merges;
        const auto until = budget.phase(0.45);
        while (!tasks.empty() && Clock::now() < until) {
            const auto range = tasks.back();
            tasks.pop_back();
            const int l = range.first, r = range.second;
            if (r - l <= 24) {
                for (int i = l + 1; i < r; ++i) {
                    int best = i;
                    Cost d = infinity;
                    for (int j = i; j < r; ++j) {
                        const Cost v = cost(order[i - 1], order[j]);
                        if (v < d) {
                            d = v;
                            best = j;
                        }
                    }
                    std::swap(order[i], order[best]);
                }
                continue;
            }
            int a = order[l], b = a;
            for (int pass = 0; pass < 2; ++pass) {
                Cost farthest = -infinity;
                for (int j = l; j < r; j += std::max(1, (r - l) / 32)) {
                    const Cost d = order[j] == a ? 0 : cost(a, order[j]);
                    if (d > farthest) {
                        farthest = d;
                        b = order[j];
                    }
                }
                if (!pass) a = b;
            }
            bool stopped = false;
            for (int j = l; j < r; ++j) {
                if (budget.poll(until)) {
                    stopped = true;
                    break;
                }
                const int v = order[j];
                score[v] = (v == a ? 0 : cost(a, v)) - (v == b ? 0 : cost(b, v));
            }
            if (stopped) break;
            const int mid = (l + r) / 2;
            std::nth_element(
                order.begin() + l, order.begin() + mid, order.begin() + r,
                [&](int u, int v) { return score[u] != score[v] ? score[u] < score[v] : u < v; });
            merges.push_back({l, r});
            tasks.push_back({mid, r});
            tasks.push_back({l, mid});
        }
        // Orient the two child paths to minimize their connecting edge.
        for (auto it = merges.rbegin(); it != merges.rend(); ++it) {
            if (budget.poll()) break;
            const int l = it->first, r = it->second, mid = (l + r) / 2;
            Cost best = infinity;
            int orientation = 0;
            for (int mask = 0; mask < 4; ++mask) {
                const Cost c = cost(order[mask & 1 ? l : mid - 1], order[mask & 2 ? r - 1 : mid]);
                if (c < best) {
                    best = c;
                    orientation = mask;
                }
            }
            if (orientation & 1) std::reverse(order.begin() + l, order.begin() + mid);
            if (orientation & 2) std::reverse(order.begin() + mid, order.begin() + r);
        }
        work_ = order;
        accept();
        seed_candidates(order, budget, 0.60);
    }
    // Stitch supplied fragments without changing their internal order. Missing
    // vertices are singleton fragments. Every interruption still yields a full tour.
    void seed(detail::Budget &budget) {
        seeded_ = true;
        if (options_.initial_tour.size() == static_cast<std::size_t>(n_)) {
            seed_candidates(tour_, budget);
            return;
        }
        if (n_ < 4) return;
        if (n_ >= 512 && options_.partition_seed && options_.initial_tour.empty() &&
            options_.initial_paths.empty()) {
            partition_seed(budget);
            return;
        }
        std::vector<std::vector<int>> paths = options_.initial_paths;
        if (!options_.initial_tour.empty()) paths.push_back(options_.initial_tour);
        std::vector<unsigned char> seen(n_, 0);
        for (const auto &p : paths)
            for (int v : p) seen[v] = 1;
        paths.erase(
            std::remove_if(paths.begin(), paths.end(), [](const auto &p) { return p.empty(); }),
            paths.end());
        for (int v = 0; v < n_; ++v)
            if (!seen[v]) paths.push_back({v});
        std::vector<int> remaining(paths.size());
        std::iota(remaining.begin(), remaining.end(), 0);
        work_.clear();
        const auto until = budget.phase(0.18);
        int chosen = 0;
        while (!remaining.empty()) {
            const int p = remaining[chosen];
            work_.insert(work_.end(), paths[p].begin(), paths[p].end());
            remaining[chosen] = remaining.back();
            remaining.pop_back();
            if (remaining.empty()) break;
            if (Clock::now() >= until) {
                for (int q : remaining)
                    work_.insert(work_.end(), paths[q].begin(), paths[q].end());
                break;
            }
            Cost best = infinity;
            chosen = 0;
            const bool full = remaining.size() <= 512;
            const int scans = full ? static_cast<int>(remaining.size()) : 96;
            for (int k = 0; k < scans; ++k) {
                const int j = full ? k : static_cast<int>(random() % remaining.size());
                const Cost d = cost(work_.back(), paths[remaining[j]].front());
                if (d < best) {
                    best = d;
                    chosen = j;
                }
            }
        }
        accept();
        seed_candidates(tour_, budget);
    }
    void heap_edge(std::vector<Edge> &heap, Edge e) {
        const auto size = static_cast<std::size_t>(std::min(options_.candidates, n_ - 1));
        if (heap.size() < size) {
            heap.push_back(e);
            std::push_heap(heap.begin(), heap.end());
        } else if (e < heap.front()) {
            std::pop_heap(heap.begin(), heap.end());
            heap.back() = e;
            std::push_heap(heap.begin(), heap.end());
        }
    }
    void nearest_candidates(detail::Budget &budget, Clock::time_point until) {
        while (nearest_row_ < n_) {
            while (nearest_column_ < n_) {
                if (budget.poll(until)) return;
                const int j = nearest_column_++;
                if (j == nearest_row_) continue;
                const Cost v = cost(nearest_row_, j);
                heap_edge(row_heap_, {j, v, v});
            }
            std::sort(row_heap_.begin(), row_heap_.end());
            candidates_[nearest_row_].swap(row_heap_);
            row_heap_.clear();
            ++nearest_row_;
            nearest_column_ = 0;
            stats_.nearest_rows = nearest_row_;
        }
    }
    void start_tree() {
        tree_ = Tree{};
        tree_.parent.assign(n_, -1);
        tree_.degree.assign(n_, 0);
        tree_.key.assign(n_, infinity);
        for (int v = 2; v < n_; ++v) tree_.remaining.push_back(v);
        tree_.started = true;
    }
    bool build_tree(detail::Budget &budget, Clock::time_point until) {
        if (!tree_.started) start_tree();
        while (!tree_.mst_done) {
            while (tree_.scan < static_cast<int>(tree_.remaining.size())) {
                if (budget.poll(until)) return false;
                const int j = tree_.scan++, v = tree_.remaining[j], u = tree_.vertex;
                const Cost w = cost(u, v) + pi_[u] + pi_[v];
                if (w < tree_.key[v]) {
                    tree_.key[v] = w;
                    tree_.parent[v] = u;
                }
                if (tree_.next_index < 0 ||
                    tree_.key[v] < tree_.key[tree_.remaining[tree_.next_index]])
                    tree_.next_index = j;
            }
            if (tree_.remaining.empty()) {
                tree_.mst_done = true;
                break;
            }
            const int j = tree_.next_index, v = tree_.remaining[j];
            tree_.sum += tree_.key[v];
            ++tree_.degree[v];
            ++tree_.degree[tree_.parent[v]];
            tree_.remaining[j] = tree_.remaining.back();
            tree_.remaining.pop_back();
            tree_.vertex = v;
            tree_.scan = 0;
            tree_.next_index = -1;
        }
        return finish_tree_root(budget, until);
    }
    bool finish_tree_root(detail::Budget &budget, Clock::time_point until) {
        while (tree_.root_scan < n_) {
            if (budget.poll(until)) return false;
            const int v = tree_.root_scan++;
            const Cost w = cost(0, v) + pi_[0] + pi_[v];
            if (w < tree_.first) {
                tree_.second = tree_.first;
                tree_.root_b = tree_.root_a;
                tree_.first = w;
                tree_.root_a = v;
            } else if (w < tree_.second) {
                tree_.second = w;
                tree_.root_b = v;
            }
        }
        tree_.degree[0] = 2;
        ++tree_.degree[tree_.root_a];
        ++tree_.degree[tree_.root_b];
        return true;
    }
    void dense_ascent(detail::Budget &budget, Clock::time_point until) {
        while (stats_.ascent_iterations < options_.ascent_iterations) {
            if (Clock::now() >= until || !build_tree(budget, until)) return;
            Cost lb = tree_.sum + tree_.first + tree_.second;
            for (Cost p : pi_) lb -= 2 * p;
            if (lb > bound_) {
                bound_ = lb;
                best_pi_ = pi_;
                best_parent_ = tree_.parent;
                best_root_a_ = tree_.root_a;
                best_root_b_ = tree_.root_b;
                best_root_second_ = tree_.second;
                stalled_ascent_ = 0;
            } else if (++stalled_ascent_ == 4) {
                lambda_ *= 0.5L;
                stalled_ascent_ = 0;
            }
            ++stats_.ascent_iterations;
            Cost norm = 0;
            for (int d : tree_.degree) norm += static_cast<Cost>(d - 2) * (d - 2);
            if (norm == 0) {
                // A connected 1-tree of degree two is itself a Hamiltonian tour.
                std::vector<std::vector<int>> adj(n_);
                auto link = [&](int a, int b) {
                    adj[a].push_back(b);
                    adj[b].push_back(a);
                };
                for (int v = 2; v < n_; ++v) link(v, tree_.parent[v]);
                link(0, tree_.root_a);
                link(0, tree_.root_b);
                work_.clear();
                int p = -1, v = 0;
                for (int k = 0; k < n_; ++k) {
                    work_.push_back(v);
                    const int q = adj[v][adj[v][0] == p ? 1 : 0];
                    p = v;
                    v = q;
                }
                accept();
                ascent_done_ = true;
                break;
            }
            const Cost step = lambda_ * std::max(Cost(0), best_cost_ - lb) / norm;
            for (int v = 0; v < n_; ++v) pi_[v] += step * (tree_.degree[v] - 2);
            tree_.started = false;
        }
        ascent_done_ = true;
    }
    bool sparse_tree(detail::Budget &budget, Clock::time_point until) {
        if (ascent_graph_.empty()) {
            ascent_graph_.resize(n_);
            ascent_next_ = next_;
        }
        while (graph_row_ < n_) {
            if (budget.poll(until)) return false;
            const int a = graph_row_++;
            auto link = [&](int b, Cost w) {
                ascent_graph_[a].push_back({b, w, w});
                ascent_graph_[b].push_back({a, w, w});
            };
            for (const auto &e : candidates_[a]) link(e.vertex, e.cost);
            link(ascent_next_[a],
                 cost(a, ascent_next_[a])); // A fixed tour ensures connectivity after removing root 0.
        }
        const auto later = [](const Edge &a, const Edge &b) { return b < a; };
        if (!sparse_tree_started_) {
            tree_ = Tree{};
            tree_.parent.assign(n_, -1);
            tree_.degree.assign(n_, 0);
            tree_.key.assign(n_, infinity);
            tree_seen_.assign(n_, 0);
            tree_seen_[1] = 1;
            sparse_count_ = 1;
            ascent_heap_.clear();
            sparse_tree_started_ = true;
        }
        while (!tree_.mst_done) {
            const int a = tree_.vertex;
            while (tree_.scan < static_cast<int>(ascent_graph_[a].size())) {
                if (budget.poll(until)) return false;
                const Edge e = ascent_graph_[a][tree_.scan++];
                const int b = e.vertex;
                if (b == 0 || tree_seen_[b]) continue;
                const Cost w = e.cost + pi_[a] + pi_[b];
                if (w < tree_.key[b]) {
                    tree_.key[b] = w;
                    tree_.parent[b] = a;
                    ascent_heap_.push_back({b, w, w});
                    std::push_heap(ascent_heap_.begin(), ascent_heap_.end(), later);
                }
            }
            int v = -1;
            while (!ascent_heap_.empty()) {
                if (budget.poll(until)) return false;
                std::pop_heap(ascent_heap_.begin(), ascent_heap_.end(), later);
                const Edge e = ascent_heap_.back();
                ascent_heap_.pop_back();
                if (!tree_seen_[e.vertex] && e.cost == tree_.key[e.vertex]) {
                    v = e.vertex;
                    break;
                }
            }
            if (v < 0) {
                tree_.mst_done = true;
                break;
            }
            tree_seen_[v] = 1;
            tree_.sum += tree_.key[v];
            ++tree_.degree[v];
            ++tree_.degree[tree_.parent[v]];
            tree_.vertex = v;
            tree_.scan = 0;
            if (++sparse_count_ == n_ - 1) tree_.mst_done = true;
        }
        return finish_tree_root(budget, until);
    }
    void ascent(detail::Budget &budget, Clock::time_point until) {
        if (!options_.sparse_ascent) {
            dense_ascent(budget, until);
            return;
        }
        if (previous_gradient_.empty()) previous_gradient_.assign(n_, 0);
        const auto sparse_until = std::min(until, budget.phase(0.25));
        while (!sparse_finished_) {
            if (stats_.ascent_iterations >= options_.ascent_iterations ||
                (stats_.ascent_iterations >= 16 && Clock::now() >= sparse_until)) {
                sparse_finished_ = true;
                pi_ = sparse_pi_;
                tree_.started = false;
                break;
            }
            if (Clock::now() >= until || !sparse_tree(budget, until)) return;
            Cost lb = tree_.sum + tree_.first + tree_.second;
            for (Cost p : pi_) lb -= 2 * p;
            // Sparse trees are only a surrogate objective, NEVER a public lower bound.
            if (lb > sparse_bound_) {
                sparse_bound_ = lb;
                sparse_pi_ = pi_;
                sparse_stall_ = 0;
            } else if (++sparse_stall_ >= 12) {
                sparse_step_scale_ *= 0.5L;
                sparse_stall_ = 0;
            }
            ++stats_.ascent_iterations;
            Cost norm = 0;
            for (int d : tree_.degree) norm += static_cast<Cost>(d - 2) * (d - 2);
            if (norm == 0) {
                sparse_finished_ = true;
                pi_ = sparse_pi_;
                tree_.started = false;
                break;
            }
            const Cost step = sparse_step_scale_ * std::max(Cost(0), best_cost_ - lb) / norm;
            for (int v = 0; v < n_; ++v) {
                const Cost gradient = tree_.degree[v] - 2;
                pi_[v] += step * (0.7L * gradient + 0.3L * previous_gradient_[v]);
                previous_gradient_[v] = gradient;
            }
            sparse_tree_started_ = false;
        }
        if (Clock::now() >= until || !build_tree(budget, until)) return;
        Cost lb = tree_.sum + tree_.first + tree_.second;
        for (Cost p : pi_) lb -= 2 * p;
        bound_ = lb;
        best_pi_ = pi_;
        best_parent_ = tree_.parent;
        best_root_a_ = tree_.root_a;
        best_root_b_ = tree_.root_b;
        best_root_second_ = tree_.second;
        ascent_done_ = true;
    }
    void prepare_alpha() {
        adjacency_.resize(n_);
        for (int v = 2; v < n_; ++v) {
            const int p = best_parent_[v];
            const Cost w = transformed(v, p);
            adjacency_[v].push_back({p, w});
            adjacency_[p].push_back({v, w});
        }
        alpha_ready_ = true;
    }
    void alpha_candidates(detail::Budget &budget, Clock::time_point until) {
        if (best_parent_.empty()) return;
        if (!alpha_ready_) prepare_alpha();
        while (alpha_row_ < n_) {
            if (!alpha_walk_started_) {
                stack_.clear();
                stack_.push_back({alpha_row_, -1, -infinity});
                alpha_walk_started_ = true;
                alpha_walk_done_ = alpha_row_ == 0;
            }
            while (!alpha_walk_done_ && !stack_.empty()) {
                if (budget.poll(until)) return;
                const Walk w = stack_.back();
                stack_.pop_back();
                beta_[w.vertex] = w.maximum;
                for (const auto &e : adjacency_[w.vertex])
                    if (e.first != w.parent)
                        stack_.push_back({e.first, w.vertex, std::max(w.maximum, e.second)});
            }
            alpha_walk_done_ = true;
            while (alpha_column_ < n_) {
                if (budget.poll(until)) return;
                const int v = alpha_column_++;
                if (v == alpha_row_) continue;
                const Cost raw = cost(alpha_row_, v), w = raw + best_pi_[alpha_row_] + best_pi_[v];
                Cost alpha;
                if (alpha_row_ == 0 || v == 0) {
                    const int other = alpha_row_ == 0 ? v : alpha_row_;
                    alpha =
                        other == best_root_a_ || other == best_root_b_ ? 0 : w - best_root_second_;
                } else
                    alpha = w - beta_[v];
                heap_edge(row_heap_, {v, raw, std::max(Cost(0), alpha)});
            }
            std::sort(row_heap_.begin(), row_heap_.end());
            // Retain a few raw nearest neighbours as a second source of candidates.
            const auto &old = candidates_[alpha_row_];
            const int keep =
                std::min(static_cast<int>(old.size()), std::max(2, options_.candidates / 4));
            for (int i = 0; i < keep; ++i) {
                bool found = false;
                for (const auto &e : row_heap_)
                    if (e.vertex == old[i].vertex) {
                        found = true;
                        break;
                    }
                if (!found) row_heap_.push_back(old[i]);
            }
            candidates_[alpha_row_].swap(row_heap_);
            row_heap_.clear();
            ++alpha_row_;
            alpha_column_ = 0;
            alpha_walk_started_ = false;
            stats_.alpha_rows = alpha_row_;
        }
    }
    void preprocess(detail::Budget &budget) {
        const auto until = budget.phase(0.38);
        if (Clock::now() >= until) return;
        const int old_rows = nearest_row_ + alpha_row_;
        const Cost old_bound = bound_;
        nearest_candidates(budget, until);
        if (nearest_row_ == n_ && locally_refined_ && options_.ascent_iterations > 0 &&
            Clock::now() < until) {
            if (!ascent_done_) ascent(budget, until);
            if (ascent_done_ && Clock::now() < until) alpha_candidates(budget, until);
        }
        if (old_rows != nearest_row_ + alpha_row_ || old_bound != bound_) {
            // Keep the cursor: repeatedly restarting at vertex zero starves the
            // other origins when preprocessing is resumed in short slices.
            std::fill(active_.begin(), active_.end(), 1);
            std::fill(deep_active_.begin(), deep_active_.end(), 3);
            stage_ = 0;
            pass_changed_ = true;
        }
    }
    bool two_opt(int a, detail::Budget &budget) {
        for (int dir = 0; dir < 2; ++dir) {
            const int b = dir ? prev_[a] : next_[a];
            const Cost ab = cost(a, b);
            for (const auto &e : candidates_[a]) {
                if (budget.poll()) return false;
                const int c = e.vertex, d = dir ? prev_[c] : next_[c];
                if (c == b || c == a || d == a || d == b) continue;
                const Cost cd = cost(c, d), bd = cost(b, d);
                if (!improves(ab + cd - e.cost - bd,
                              std::abs(ab) + std::abs(cd) + std::abs(e.cost) + std::abs(bd)))
                    continue;
                int i = position_[dir ? b : a], j = position_[dir ? d : c];
                if (i > j) std::swap(i, j);
                if (exact_integral_) {
                    activate(a);
                    activate(b);
                    activate(c);
                    activate(d);
                    std::reverse(tour_.begin() + i + 1, tour_.begin() + j + 1);
                    for (int k = i; k <= j + 1 && k < n_; ++k) {
                        position_[tour_[k]] = k;
                        next_[tour_[k]] = tour_[after(k)];
                        prev_[tour_[after(k)]] = tour_[k];
                    }
                    current_ -= ab + cd - e.cost - bd;
#ifdef SYMMETRIC_TSP_CHECK_INVARIANTS
                    if (current_ != value(tour_))
                        throw std::logic_error("symmetric_tsp: incorrect 2-opt cost");
#endif
                    remember();
                    ++stats_.moves;
                    ++stats_.two_opt_moves;
                    return true;
                }
                work_ = tour_;
                std::reverse(work_.begin() + i + 1, work_.begin() + j + 1);
                if (accept()) {
                    ++stats_.two_opt_moves;
                    return true;
                }
            }
        }
        return false;
    }
    bool relocate(int a, detail::Budget &budget) {
        const int p = prev_[a];
        int last = a;
        for (int length = 1; length <= 3 && length < n_ - 1; ++length) {
            const int q = next_[last];
            for (const auto &e : candidates_[a])
                for (int reverse = 0; reverse < 2; ++reverse) {
                    if (budget.poll()) return false;
                    const int b = reverse ? prev_[e.vertex] : e.vertex, c = next_[b];
                    bool overlaps = b == p;
                    int v = a;
                    for (int k = 0; k < length; ++k) {
                        overlaps = overlaps || v == b;
                        v = next_[v];
                    }
                    if (overlaps) continue;
                    const Cost old = cost(p, a) + cost(last, q) + cost(b, c);
                    const Cost now = cost(p, q) + e.cost + (reverse ? cost(b, last) : cost(last, c));
                    if (!improves(old - now, std::abs(old) + std::abs(now))) continue;
                    work_.clear();
                    for (int u : tour_) {
                        int offset = position_[u] - position_[a];
                        if (offset < 0) offset += n_;
                        if (offset < length) continue;
                        work_.push_back(u);
                        if (u == b) {
                            v = reverse ? last : a;
                            for (int k = 0; k < length; ++k) {
                                work_.push_back(v);
                                v = reverse ? prev_[v] : next_[v];
                            }
                        }
                    }
                    if (accept(current_ - (old - now))) {
                        ++stats_.relocate_moves;
                        return true;
                    }
                }
            last = q;
        }
        return false;
    }
    // The original tour minus k disjoint edges consists of k unchanged paths.
    // Pair their endpoints in tour order, then alternate path/added-edge links.
    bool close_exchange(int k) {
        std::array<int, depth_limit> order{}, tail{};
        std::array<int, 2 * depth_limit> mate{};
        auto vertex = [&](int endpoint) {
            return endpoint & 1 ? y_[endpoint / 2] : x_[endpoint / 2];
        };
        for (int i = 0; i < k; ++i) {
            order[i] = i;
            tail[i] = next_[x_[i]] == y_[i] ? 2 * i : 2 * i + 1;
        }
        // Insertion sort is cheaper than generic sorting for these tiny arrays.
        for (int i = 1; i < k; ++i) {
            const int v = order[i];
            int j = i;
            while (j && position_[vertex(tail[order[j - 1]])] > position_[vertex(tail[v])]) {
                order[j] = order[j - 1];
                --j;
            }
            order[j] = v;
        }
        for (int i = 0; i < k; ++i) {
            const int a = tail[order[i]] ^ 1, b = tail[order[(i + 1) % k]];
            mate[a] = b;
            mate[b] = a;
        }
        auto added = [k](int e) { return e & 1 ? (e + 1) % (2 * k) : (e + 2 * k - 1) % (2 * k); };
        int endpoint = 0;
        for (int count = 1; count <= k; ++count) {
            endpoint = added(mate[endpoint]);
            if (endpoint == 0 && count != k) return false;
        }
        if (endpoint != 0) return false;
        work_.clear();
        for (int count = 0; count < k; ++count) {
            const int stop = vertex(mate[endpoint]);
            int v = vertex(endpoint);
            const bool forward = endpoint != tail[endpoint / 2];
            for (;;) {
                work_.push_back(v);
                if (v == stop) break;
                v = forward ? next_[v] : prev_[v];
            }
            endpoint = added(mate[endpoint]);
        }
        Cost gain = 0;
        if (exact_integral_)
            for (int i = 0; i < k; ++i) gain += cost(x_[i], y_[i]) - cost(y_[i], x_[(i + 1) % k]);
        if (accept(exact_integral_ ? current_ - gain : infinity)) {
            ++stats_.exchange_moves;
            return true;
        }
        return false;
    }
    // A positive, infeasible 2-opt can split the tour into two cycles. Join
    // them by replacing one unchanged edge in each; the join may lose gain.
    bool patch_two_cycles(int a, int b, int c, int d, Cost split_gain, detail::Budget &budget) {
        const bool forward = next_[a] == b;
        if ((forward ? next_[c] : prev_[c]) != d) return false;
        auto offset = [&](int v) {
            const int q = forward ? position_[v] - position_[b] : position_[b] - position_[v];
            return q < 0 ? q + n_ : q;
        };
        const int length = offset(c) + 1;
        const bool first = length <= n_ - length;
        const int count = first ? length : n_ - length;
        const int begin = first ? position_[b] : position_[d];
        // Spread a bounded scan over the smaller cycle, including its boundaries.
        const int scans = std::min(count, 32);
        for (int i = 0; i < scans && patches_left_ > 0; ++i) {
            if (budget.poll()) return false;
            const int delta = i * count / scans;
            const int u = tour_[(begin + (forward ? delta : n_ - delta)) % n_];
            const int v = forward ? next_[u] : prev_[u];
            if (u == a || u == b || u == c || u == d || v == a || v == b || v == c || v == d)
                continue;
            const Cost uv = cost(u, v);
            for (const auto &e : candidates_[u]) {
                const int w = e.vertex;
                if ((offset(w) < length) == first || w == a || w == b || w == c || w == d)
                    continue;
                for (int dir = 0; dir < 2 && patches_left_ > 0; ++dir) {
                    const int z = dir ? prev_[w] : next_[w];
                    if (z == a || z == b || z == c || z == d) continue;
                    --patches_left_;
                    const Cost wz = cost(w, z), vz = cost(v, z);
                    const Cost gain = split_gain + uv + wz - e.cost - vz;
                    if (!improves(gain, std::abs(split_gain) + std::abs(uv) + std::abs(wz) +
                                            std::abs(e.cost) + std::abs(vz)))
                        continue;
                    const std::array<int, 8> endpoints{{a, b, c, d, u, v, w, z}};
                    const std::array<int, 8> added{{d, c, b, a, w, z, u, v}};
                    work_.clear();
                    int previous = -1, vertex = a;
                    for (int j = 0; j < n_; ++j) {
                        work_.push_back(vertex);
                        int left = prev_[vertex], right = next_[vertex];
                        for (int k = 0; k < 8; ++k)
                            if (endpoints[k] == vertex) {
                                if (left == endpoints[k ^ 1])
                                    left = added[k];
                                else
                                    right = added[k];
                                break;
                            }
                        const int next = left == previous ? right : left;
                        previous = vertex;
                        vertex = next;
                        if (vertex == a) break;
                    }
                    if (work_.size() == static_cast<std::size_t>(n_) && vertex == a &&
                        accept(current_ - gain)) {
                        ++stats_.patching_moves;
                        return true;
                    }
                }
            }
        }
        return false;
    }
    bool exchange_dfs(int k, Cost gain, detail::Budget &budget) {
        if (k >= exchange_depth_ || nodes_left_ <= 0 || budget.poll()) return false;
        const int a = x_[0], b = y_[k - 1];
        auto &branches = branches_[k];
        branches.clear();
        int candidate_index = 0;
        for (const auto &e : candidates_[b]) {
            if (candidate_index++ >= options_.search_candidates) break;
            const int c = e.vertex;
            if (used_[c] || c == next_[b] || c == prev_[b]) continue;
            const Cost remaining = gain - e.cost - best_pi_[b] - best_pi_[c];
            if (!improves(remaining, std::abs(gain) + std::abs(e.cost + best_pi_[b] + best_pi_[c])))
                continue;
            for (int dir = 0; dir < 2; ++dir) {
                const int d = dir ? prev_[c] : next_[c];
                if (used_[d]) continue;
                const Cost g = remaining + transformed(c, d);
                branches.push_back({c, d, g, -g});
            }
        }
        const int width =
            std::min(static_cast<int>(branches.size()),
                     k == 1 ? std::max(2 * options_.search_breadth,
                                       std::min(options_.search_candidates, options_.candidates))
                            : options_.search_breadth);
        std::partial_sort(branches.begin(), branches.begin() + width, branches.end());
        for (int i = 0; i < width && nodes_left_ > 0; ++i) {
            if (budget.poll()) return false;
            const Branch bch = branches[i];
            --nodes_left_;
            ++stats_.search_nodes;
            x_[k] = bch.c;
            y_[k] = bch.d;
            used_[bch.c] = used_[bch.d] = 1;
            const Cost closing = transformed(bch.d, a);
            bool changed = false;
            if (bch.d != next_[a] && bch.d != prev_[a] &&
                improves(bch.gain - closing, std::abs(bch.gain) + std::abs(closing)))
                changed = close_exchange(k + 1);
            if (!changed && k == 1 && options_.patching && patches_left_ > 0) {
                const Cost split =
                    cost(a, y_[0]) + cost(bch.c, bch.d) - cost(y_[0], bch.c) - cost(bch.d, a);
                if (split > 0) changed = patch_two_cycles(a, y_[0], bch.c, bch.d, split, budget);
            }
            if (!changed) changed = exchange_dfs(k + 1, bch.gain, budget);
            used_[bch.c] = used_[bch.d] = 0;
            if (changed) return true;
        }
        return false;
    }
    // Open (a,b), then repeatedly connect the free endpoint to an interior
    // vertex and reverse the remaining suffix. Every state is a Hamiltonian
    // path. Only a closing edge that improves the ORIGINAL tour is committed.
    // Segments reference the unchanged tour, so rotations cost O(depth), not O(n).
    bool chained_rotation(int a, detail::Budget &budget) {
        struct Segment {
            int lo, hi;
            bool reverse;
        };
        constexpr int limit = 48;
        std::array<Segment, limit + 1> segments{};
        std::array<std::pair<int, int>, limit> additions{};
        std::array<int, limit + 1> endpoints{};
        for (int orientation : {-1, 1}) {
            auto base = [&](int v) {
                const int p = orientation * (position_[v] - position_[a]);
                return p < 0 ? p + n_ : p;
            };
            auto vertex = [&](int p) { return tour_[(position_[a] + orientation * p + n_) % n_]; };
            segments[0] = {0, n_ - 1, false};
            int count = 1;
            int b = vertex(n_ - 1);
            endpoints[0] = b;
            Cost raw_gain = cost(a, b), gain = transformed(a, b);
            for (int depth = 0; depth < limit; ++depth) {
                if (budget.poll()) return false;
                int chosen_c = -1, chosen_d = -1, chosen_segment = -1;
                Cost chosen_gain = -infinity, chosen_raw = 0;
                int candidate_index = 0;
                for (const auto &e : candidates_[b]) {
                    if (candidate_index++ >= options_.search_candidates) break;
                    if (budget.poll()) return false;
                    const int c = e.vertex;
                    if (c == a || c == b) continue;
                    const Cost remaining = gain - e.cost - best_pi_[b] - best_pi_[c];
                    if (!improves(remaining,
                                  std::abs(gain) + std::abs(e.cost + best_pi_[b] + best_pi_[c])))
                        continue;
                    const int p = base(c);
                    int segment = 0;
                    while (segment < count && (p < segments[segment].lo || p > segments[segment].hi))
                        ++segment;
                    if (segment == count) continue;
                    const auto s = segments[segment];
                    int successor;
                    if (p == (s.reverse ? s.lo : s.hi)) {
                        if (segment + 1 == count) continue;
                        const auto t = segments[segment + 1];
                        successor = t.reverse ? t.hi : t.lo;
                    } else
                        successor = p + (s.reverse ? -1 : 1);
                    const int d = vertex(successor);
                    bool forbidden = false;
                    for (int k = 0; k <= depth; ++k) forbidden = forbidden || endpoints[k] == d;
                    for (int k = 0; k < depth; ++k)
                        forbidden = forbidden ||
                                    (additions[k].first == c && additions[k].second == d) ||
                                    (additions[k].first == d && additions[k].second == c);
                    if (forbidden) continue;
                    const Cost cd = cost(c, d), raw = raw_gain - e.cost + cd;
                    const Cost g = remaining + cd + best_pi_[c] + best_pi_[d];
                    const Cost closing = cost(d, a);
                    if (improves(raw - closing, std::abs(raw) + std::abs(closing))) {
                        chosen_c = c;
                        chosen_d = d;
                        chosen_segment = segment;
                        chosen_gain = g;
                        chosen_raw = raw;
                        break;
                    }
                    if (g > chosen_gain) {
                        chosen_c = c;
                        chosen_d = d;
                        chosen_segment = segment;
                        chosen_gain = g;
                        chosen_raw = raw;
                    }
                }
                if (chosen_c < 0) break;
                const int p = base(chosen_c), j = chosen_segment;
                auto &s = segments[j];
                if (p != (s.reverse ? s.lo : s.hi)) {
                    for (int k = count; k > j + 1; --k) segments[k] = segments[k - 1];
                    segments[j + 1] = s;
                    if (s.reverse) {
                        s.lo = p;
                        segments[j + 1].hi = p - 1;
                    } else {
                        s.hi = p;
                        segments[j + 1].lo = p + 1;
                    }
                    ++count;
                }
                std::reverse(segments.begin() + j + 1, segments.begin() + count);
                for (int k = j + 1; k < count; ++k) segments[k].reverse = !segments[k].reverse;
                additions[depth] = {b, chosen_c};
                endpoints[depth + 1] = chosen_d;
                b = chosen_d;
                gain = chosen_gain;
                raw_gain = chosen_raw;
                const Cost closing = cost(b, a);
                if (!improves(raw_gain - closing, std::abs(raw_gain) + std::abs(closing)))
                    continue;
                work_.clear();
                for (int k = 0; k < count; ++k) {
                    const auto t = segments[k];
                    if (t.reverse)
                        for (int q = t.hi; q >= t.lo; --q) work_.push_back(vertex(q));
                    else
                        for (int q = t.lo; q <= t.hi; ++q) work_.push_back(vertex(q));
                }
                if (accept(current_ - (raw_gain - closing))) {
                    ++stats_.chained_moves;
                    return true;
                }
            }
        }
        return false;
    }
    bool exchange(int a, detail::Budget &budget) {
        nodes_left_ = options_.search_nodes;
        patches_left_ = 96;
        exchange_depth_ = options_.staged_search && stage_ == 1 ? std::min(5, options_.max_depth)
                                                                : options_.max_depth;
        for (int dir = 0; dir < 2; ++dir) {
            x_[0] = a;
            y_[0] = dir ? prev_[a] : next_[a];
            used_[a] = used_[y_[0]] = 1;
            const bool changed = exchange_dfs(1, transformed(a, y_[0]), budget);
            used_[a] = used_[y_[0]] = 0;
            if (changed) return true;
        }
        return false;
    }
    // Iterative partial transcription. For each endpoint pair, an interval
    // in each parent has exactly the same vertices iff the largest relative
    // position in the other parent is length - 1 (all positions are distinct).
    // This deterministic test needs no probabilistic hash or bitset comparison.
    bool merge(detail::Budget &budget) {
        merge_changed_current_ = false;
        if (!options_.merge_tours || best_is_current_) return false;
        std::vector<int> other_position(n_);
        for (int i = 0; i < n_; ++i) other_position[best_[i]] = i;
        bool changed = false;
        const int limit = std::min(n_ - 1, 512);
        for (int start = merge_origin_; start < n_; ++start) {
            merge_origin_ = start;
            if (budget.poll()) return changed;
            const int a = tour_[start], base = other_position[a];
            for (int direction : {1, -1}) {
                int maximum = 0;
                Cost cost_a = 0, cost_b = 0;
                for (int length = 2; length <= limit; ++length) {
                    if (budget.poll()) return changed;
                    const int end = (start + length - 1) % n_;
                    const int b = tour_[end];
                    int relative = direction * (other_position[b] - base);
                    if (relative < 0) relative += n_;
                    maximum = std::max(maximum, relative);
                    if (maximum >= limit) break;
                    cost_a += cost(tour_[(end + n_ - 1) % n_], b);
                    const int j = (base + direction * (length - 1) + n_) % n_;
                    const int previous = (j - direction + n_) % n_;
                    cost_b += cost(best_[previous], best_[j]);
                    if (maximum != length - 1 || relative != length - 1 || length < 4) continue;
                    if (!improves(std::abs(cost_a - cost_b), std::abs(cost_a) + std::abs(cost_b)))
                        continue;
                    if (cost_b < cost_a) {
                        work_ = tour_;
                        for (int k = 1; k < length - 1; ++k)
                            work_[(start + k) % n_] = best_[(base + direction * k + n_) % n_];
                        if (accept()) {
                            ++stats_.merging_moves;
                            merge_changed_current_ = true;
                            return true;
                        }
                    } else {
                        work_ = best_;
                        for (int k = 1; k < length - 1; ++k)
                            work_[(base + direction * k + n_) % n_] = tour_[(start + k) % n_];
                        const Cost v = value(work_);
                        if (v < best_cost_) {
                            best_.swap(work_);
                            best_cost_ = v;
                            ++stats_.improvements;
                            ++stats_.merging_moves;
                            stagnation_ = 0;
                            changed = true;
                            merge_origin_ = 0;
                            for (int k = 0; k < n_; ++k) other_position[best_[k]] = k;
                            // Restart after changing a parent; relative positions changed.
                            return true;
                        }
                    }
                }
            }
        }
        merge_origin_ = 0;
        return changed;
    }
    void kick() {
        save_best();
        ++stats_.kicks;
        ++stagnation_;
        // Keep an independent incumbent. Periodically restart near it.
        if ((stagnation_ & 3U) == 0) {
            tour_ = best_;
            current_ = best_cost_;
            index();
        }
        const int repeats = 1 + static_cast<int>(std::min<std::uint64_t>(3, stagnation_ / 12));
        for (int repeat = 0; repeat < repeats; ++repeat) {
            if (n_ >= 8) {
                std::array<int, 4> cut{};
                for (int i = 0; i < 4; ++i) {
                    bool duplicate;
                    do {
                        cut[i] = static_cast<int>(random() % static_cast<std::uint64_t>(n_));
                        duplicate = false;
                        for (int j = 0; j < i; ++j)
                            if (cut[j] == cut[i]) duplicate = true;
                    } while (duplicate);
                }
                if (options_.guided_kicks && repeat == 0 && (stats_.kicks & 1U)) {
                    // Focus some perturbations on nearby edges, with an expensive
                    // edge as the first cut. Uniform bridges still supply diversity.
                    auto excess = [&](int p) {
                        const int v = tour_[p];
                        Cost shortest = infinity;
                        for (const auto &e : candidates_[v]) shortest = std::min(shortest, e.cost);
                        return cost(prev_[v], v) - (shortest == infinity ? 0 : shortest);
                    };
                    Cost score = excess(cut[0]);
                    for (int sample = 0; sample < 12; ++sample) {
                        const int p = static_cast<int>(random() % static_cast<std::uint64_t>(n_));
                        const Cost w = excess(p);
                        if (w > score) {
                            score = w;
                            cut[0] = p;
                        }
                    }
                    for (int i = 1; i < 4; ++i) {
                        for (int attempt = 0;; ++attempt) {
                            int p;
                            const auto &edges = candidates_[tour_[cut[i - 1]]];
                            if (attempt < 16 && !edges.empty()) {
                                const int v =
                                    edges[random() % std::min<std::size_t>(16, edges.size())].vertex;
                                p = position_[v];
                                if (random() & 1U) p = after(p);
                            } else
                                p = static_cast<int>(random() % static_cast<std::uint64_t>(n_));
                            bool duplicate = false;
                            for (int j = 0; j < i; ++j) {
                                const int d = std::abs(cut[j] - p);
                                if (d < 2 || d > n_ - 2) duplicate = true;
                            }
                            // For n=8 an unlucky third cut can leave no fourth
                            // nonadjacent cut, so eventually require only distinctness.
                            if (attempt >= 32) {
                                duplicate = false;
                                for (int j = 0; j < i; ++j) duplicate = duplicate || cut[j] == p;
                            }
                            if (!duplicate) {
                                cut[i] = p;
                                break;
                            }
                        }
                    }
                }
                std::sort(cut.begin(), cut.end());
                work_.clear();
                auto append = [&](int l, int r) {
                    work_.insert(work_.end(), tour_.begin() + l, tour_.begin() + r);
                };
                append(0, cut[0]);
                append(cut[2], cut[3]);
                append(cut[1], cut[2]);
                append(cut[0], cut[1]);
                append(cut[3], n_);
                tour_.swap(work_);
            } else {
                const int a = static_cast<int>(random() % static_cast<std::uint64_t>(n_));
                const int b = static_cast<int>(random() % static_cast<std::uint64_t>(n_));
                std::swap(tour_[a], tour_[b]);
            }
        }
        current_ = value(tour_);
        index();
        remember();
        std::fill(active_.begin(), active_.end(), 1);
        std::fill(deep_active_.begin(), deep_active_.end(), 3);
    }
    Result finish(StopReason reason, Clock::time_point start) const {
        Result r;
        r.tour = best_is_current_ ? tour_ : best_;
        r.cost = best_cost_;
        r.lower_bound = bound_;
        r.reason = reason;
        r.stats = stats_;
        r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        return r;
    }

public:
    Solver(int vertices, Distance distance, const Options &options = {})
        : n_(vertices), distance_(std::move(distance)), options_(options), random_(options.seed) {
        if (n_ < 0 || options.search_candidates < 2 || options.search_candidates > 256 ||
            options.candidates < 2 || options.candidates > 256 || options.max_depth < 2 ||
            options.max_depth > depth_limit || options.search_breadth < 1 ||
            options.search_breadth > 256 || options.search_nodes < 1 ||
            options.ascent_iterations < 0 || options.ascent_iterations > 1000 ||
            !std::isfinite(options.safety_margin) || options.safety_margin < 0 ||
            !std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0 ||
            (!options.initial_tour.empty() && !options.initial_paths.empty()))
            throw std::invalid_argument("symmetric_tsp: invalid options");
    }
    Result run() { return run(options_.time_limit); }
    Result run(double seconds, Clock::time_point outer_deadline = Clock::time_point::max()) {
        const auto start = Clock::now();
        detail::Budget budget(
            start, seconds, std::min(options_.deadline, outer_deadline), options_.safety_margin);
        if (!initialized_) initialize();
        if (n_ <= 3) return finish(StopReason::trivial, start);
        if (!options_.max_trials) return finish(StopReason::trial_limit, start);
        if (budget.expired()) return finish(StopReason::time_limit, start);
        if (!seeded_) seed(budget);
        preprocess(budget);
        if (pending_kick_ && !budget.expired()) {
            if (options_.kicks) kick();
            pending_kick_ = false;
            stage_ = 0;
            origin_ = 0;
            pass_changed_ = false;
        }
        std::uint64_t trials = 0;
        while (!budget.expired() && trials < options_.max_trials) {
            if (origin_ == n_) {
                origin_ = 0;
                if (pass_changed_) {
                    pass_changed_ = false;
                    if (stage_ > 0) {
                        stage_ = 0;
                        std::fill(active_.begin(), active_.end(), 1);
                    }
                    continue;
                }
                if (stage_ == 0 && !locally_refined_) {
                    locally_refined_ = true;
                    preprocess(budget);
                    continue;
                }
                if (stage_ == 0 && options_.max_depth > 2) {
                    stage_ = 1;
                    continue;
                }
                if (stage_ == 1 && options_.staged_search && options_.max_depth > 5) {
                    stage_ = 2;
                    continue;
                }
                bool merged = false;
                while (!budget.expired() && merge(budget)) {
                    merged = true;
                    if (merge_changed_current_) break;
                }
                if (merged) {
                    if (!merge_changed_current_) {
                        tour_ = best_;
                        current_ = best_cost_;
                        best_is_current_ = true;
                        index();
                    }
                    stage_ = 0;
                    pass_changed_ = false;
                    std::fill(active_.begin(), active_.end(), 1);
                    std::fill(deep_active_.begin(), deep_active_.end(), 3);
                    continue;
                }
                if (budget.expired()) {
                    origin_ = n_;
                    break;
                }
                ++stats_.trials;
                ++trials;
                pending_kick_ = true;
                if (trials == options_.max_trials || !options_.kicks) break;
                kick();
                pending_kick_ = false;
                stage_ = 0;
                pass_changed_ = false;
            }
            const int a = origin_++;
            bool changed = false;
            if (stage_ == 0) {
                if (!active_[a]) continue;
                active_[a] = 0;
                changed = two_opt(a, budget);
                if (!changed && options_.relocate && !budget.expired())
                    changed = relocate(a, budget);
            } else {
                const unsigned char bit = stage_ == 1 ? 1 : 2;
                if (options_.active_search && !(deep_active_[a] & bit)) continue;
                deep_active_[a] &= static_cast<unsigned char>(~bit);
                changed = exchange(a, budget);
                if (!changed && options_.chained_rotations && !budget.expired())
                    changed = chained_rotation(a, budget);
            }
            if (changed)
                pass_changed_ = true;
            else if (budget.expired()) {
                // An interrupted origin has not been exhausted. Retry it on resume.
                --origin_;
                active_[a] = 1;
                deep_active_[a] = 3;
                break;
            }
        }
        return finish(budget.expired() ? StopReason::time_limit : StopReason::trial_limit, start);
    }
};

template <class Distance>
auto make_solver(int n, Distance &&distance, const Options &options = {}) {
    return Solver<typename std::decay<Distance>::type>(n, std::forward<Distance>(distance), options);
}
template <class T>
auto make_solver(const std::vector<std::vector<T>> &matrix, const Options &options = {}) {
    if (matrix.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("symmetric_tsp: matrix too large");
    for (const auto &row : matrix)
        if (row.size() != matrix.size())
            throw std::invalid_argument("symmetric_tsp: matrix is not square");
    return make_solver(
        static_cast<int>(matrix.size()), [&matrix](int a, int b) { return matrix[a][b]; }, options);
}
// Reject temporary matrices: a persistent session keeps a reference to its input.
template <class T> auto make_solver(std::vector<std::vector<T>> &&, const Options & = {}) = delete;
template <class T>
auto make_solver(const std::vector<std::vector<T>> &&, const Options & = {}) = delete;
template <class Distance> Result solve(int n, Distance &&distance, const Options &options = {}) {
    const auto start = Clock::now();
    const auto stop = detail::deadline(start, options.time_limit, options.deadline);
    auto solver = make_solver(n, std::forward<Distance>(distance), options);
    auto r = solver.run(options.time_limit, stop);
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    return r;
}
template <class T>
Result solve(const std::vector<std::vector<T>> &matrix, const Options &options = {}) {
    const auto start = Clock::now();
    const auto stop = detail::deadline(start, options.time_limit, options.deadline);
    auto solver = make_solver(matrix, options);
    auto r = solver.run(options.time_limit, stop);
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    return r;
}
} // namespace symmetric_tsp
#endif
