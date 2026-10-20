#ifndef STEINER_TREE_HPP
#define STEINER_TREE_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace steiner_tree {
using Clock = std::chrono::steady_clock;
enum class StopReason { time_limit, work_limit, optimal, infeasible };

struct Options {
    double time_limit = 1.0; // Wall seconds for this call, including result construction.
    Clock::time_point deadline = Clock::time_point::max();
    double safety_margin = 0.0002;
    std::uint64_t max_work = std::numeric_limits<std::uint64_t>::max();
    int exact_terminal_limit = 11; // 0 disables subset DP; trivial exact cases remain enabled.
    std::size_t exact_memory_limit = 64 * 1024 * 1024;
    bool local_search = true;
    bool recombination = true;
    bool vertex_moves = false; // Optional original-graph vertex insertion/removal neighborhoods.
    bool cascaded_recombination = true;
    bool vertex_perturbation = true;
    bool fast_key_paths = true;
    bool steiner_reconnection = true;
};

struct Statistics {
    double elapsed = 0; // This call. The following counters are cumulative.
    std::uint64_t work = 0, constructions = 0, candidates = 0;
    std::uint64_t improvements = 0, local_moves = 0, recombinations = 0;
    std::uint64_t vertex_insertions = 0, vertex_removals = 0, induced_improvements = 0;
    int shortest_path_rows = 0, floyd_pivots = 0, exact_masks = 0;
};

template <class Weight = std::int64_t> class Solver {
    static_assert(std::is_arithmetic<Weight>::value && !std::is_same<Weight, bool>::value,
                  "Weight must be an integer or floating point type");
    static_assert(!std::is_integral<Weight>::value || sizeof(Weight) <= 8,
                  "Integer edge weights must fit in 64 bits");

public:
    using cost_type =
        typename std::conditional<std::is_integral<Weight>::value, __int128_t, long double>::type;
    struct Edge {
        int u, v;
        Weight weight;
    };
    struct Result {
        std::vector<int> edges; // Original input edge IDs; sorted, no duplicates.
        cost_type cost = 0, lower_bound = 0;
        bool feasible = false;
        bool optimal = false; // Floating point: subject to arithmetic rounding.
        StopReason reason = StopReason::time_limit;
        Statistics stats;
    };

private:
    static constexpr cost_type inf() {
        if constexpr (std::is_integral<Weight>::value)
            return std::numeric_limits<cost_type>::max() / 8;
        else
            return std::numeric_limits<cost_type>::infinity();
    }
    struct DSU {
        std::vector<int> p, size;
        void reset(int n) {
            p.resize(n);
            std::iota(p.begin(), p.end(), 0);
            size.assign(n, 1);
        }
        int find(int x) {
            while (x != p[x]) {
                p[x] = p[p[x]];
                x = p[x];
            }
            return x;
        }
        bool join(int a, int b) {
            a = find(a);
            b = find(b);
            if (a == b) return false;
            if (size[a] < size[b]) std::swap(a, b);
            p[b] = a;
            size[a] += size[b];
            return true;
        }
    };
    struct Item {
        cost_type d;
        int v;
        bool operator<(const Item &b) const { return d != b.d ? d > b.d : v > b.v; }
    };
    template <class T> struct ReusableHeap : std::priority_queue<T> {
        void clear() { this->c.clear(); } // Keep storage across shortest paths/candidates.
    };
    using Heap = ReusableHeap<Item>;
    struct EdgeItem {
        Weight d;
        int v; // Normalization only orders original edge weights, not path sums.
        bool operator<(const EdgeItem &b) const { return d != b.d ? d > b.d : v > b.v; }
    };
    struct Shortest {
        std::vector<cost_type> d;
        std::vector<int> pred;
        Heap heap;
        int vertex = -1;
        std::size_t arc = 0, end = 0;
        cost_type vertex_distance = 0;
        bool done = false;
        void reset(int n) {
            d.assign(n, inf());
            pred.assign(n, -1);
            heap.clear();
            vertex = -1;
            arc = 0;
            done = false;
        }
        void source(int v, cost_type value = 0) {
            d[v] = value;
            pred[v] = -1;
            heap.push({value, v});
        }
    };
    struct Tree {
        std::vector<int> edges;
        cost_type cost = 0;
    };
    struct Normalize {
        bool active = false, done = false;
        std::vector<int> input;
        std::vector<int> base; // Already sorted by cost; merged with the heap.
        std::size_t cursor = 0;
        std::size_t base_cursor = 0;
        ReusableHeap<EdgeItem> heap;
        DSU dsu;
        Tree tree;
        void reset() {
            active = done = false;
            cursor = base_cursor = 0;
            input.clear();
            base.clear();
            heap.clear();
            tree.edges.clear();
            tree.cost = 0;
        }
    } norm_;
    struct Metric {
        bool active = false, done = false;
        std::vector<int> nodes, parent, remaining;
        std::vector<cost_type> key;
        std::vector<std::pair<int, int>> links;
        int count = 0, next = -1;
        std::size_t link = 0;
        int path_vertex = -1;
        cost_type total = 0;
        std::vector<int> edges;
        void reset() {
            active = done = false;
            count = 0;
            next = -1;
            link = 0;
            path_vertex = -1;
            total = 0;
            nodes.clear();
            parent.clear();
            remaining.clear();
            key.clear();
            links.clear();
            edges.clear();
        }
    } metric_;
    struct Cut {
        bool active = false, done = false;
        std::vector<int> retained, vertices, component, a, b;
        DSU dsu;
        std::vector<int> removed, index;
        std::vector<unsigned char> present;
        std::vector<cost_type> distance;
        std::vector<cost_type> to_vertex;
        std::vector<int> nearest;
        std::size_t row = 0;
        int components = 0;
        bool star_checked = false;
        bool star_enabled = false;
        int star_vertex = 0, star_center = -1;
        cost_type star_best = inf();
        void reset() {
            active = done = star_checked = star_enabled = false;
            row = 0;
            components = star_vertex = 0;
            star_center = -1;
            star_best = inf();
            retained.clear();
            vertices.clear();
            component.clear();
            a.clear();
            b.clear();
            distance.clear();
            to_vertex.clear();
            nearest.clear();
            removed.clear();
            index.clear();
            present.clear();
        }
    } cut_;
    struct Induced {
        int vertex = -1; // -1: optimize induced graph; otherwise insert/remove this vertex.
        bool removal = false, built = false;
        std::size_t cursor = 0;
        std::vector<int> edges;
        void reset() {
            vertex = -1;
            removal = built = false;
            cursor = 0;
            edges.clear();
        }
    } induced_;
    struct Exchange {
        std::vector<int> vertices, parent, depth, parent_path, interior_path;
        std::vector<cost_type> path_cost;
        std::size_t i = 0, j = 1;
        int path = -1, a = -1, b = -1;
        cost_type gain = 0;
        std::vector<cost_type> bound;
        bool done = false;
        void reset() {
            i = 0;
            j = 1;
            path = a = b = -1;
            gain = 0;
            done = false;
            vertices.clear();
            parent.clear();
            depth.clear();
            parent_path.clear();
            interior_path.clear();
            path_cost.clear();
        }
    } exchange_;
    struct Cell {
        cost_type cost = inf();
        int split = 0, edge = -1;
    };
    struct Exact {
        bool active = false, finished = false;
        int count = 0, full = 0, mask = 1, sub = 0, vertex = 0, stage = 0;
        std::vector<std::vector<Cell>> table;
        Shortest shortest;
        std::vector<std::pair<int, int>> trace;
        std::vector<int> edges;
    } exact_;

    int n_ = 0;
    std::vector<Edge> edges_;
    struct Arc {
        int edge, vertex;
        Weight weight;
    };
    std::vector<std::vector<Arc>> adj_;
    std::vector<int> terminals_;
    std::vector<unsigned char> terminal_;
    bool feasible_ = false, optimal_ = false, all_terminals_ = false;
    Tree best_, current_;
    Tree cascade_best_, fresh_;
    bool cascading_ = false;
    int cascade_failures_ = 0;
    cost_type lower_bound_ = 0;
    Statistics stats_;
    std::uint64_t rng_;
    int phase_ = 0; // 0 construction, 1 all-pairs, 2 terminal MST, 3 local search, 4 MST exact.
    int construction_mode_ = 0;
    std::size_t weight_cursor_ = 0;
    std::vector<cost_type> perturbed_;
    std::vector<unsigned char> bias_, in_tree_;
    std::vector<int> vertex_factor_;
    int missing_ = 0;
    Shortest shortest_;
    Tree constructed_;
    bool construction_started_ = false, construction_normalizing_ = false;
    std::vector<std::vector<cost_type>> distance_;
    std::vector<std::vector<int>> predecessor_;
    int source_ = 0;
    bool dense_apsp_ = false;
    bool narrow_distances_ = false;
    std::size_t apsp_edge_ = 0;
    int floyd_pivot_ = 0, floyd_row_ = 0;
    bool row_started_ = false, terminal_mst_started_ = false;
    bool local_ready_ = false;
    int local_stage_ = 0;
    std::size_t candidate_ = 0;
    std::vector<int> keys_, additions_, removals_, vertex_additions_, vertex_removals_;
    std::vector<std::vector<int>> key_paths_, key_incident_;
    std::vector<int> degree_;
    std::vector<int> tree_sorted_;
    std::vector<int> edge_path_;
    std::vector<int> prune_degree_, prune_neighbor_, prune_queue_;
    std::vector<std::vector<int>> local_tree_;
    std::vector<int> path_buffer_;
    std::vector<unsigned char> key_flags_;
    std::vector<Tree> elite_;
    int candidate_kind_ = 0; // 1 metric, 2 cut, 3 induced graph; heuristic normalizer is shared.
    std::uint64_t generation_ = 0;
    bool exact_decided_ = false;

    int other(int e, int v) const { return edges_[e].u ^ edges_[e].v ^ v; }
    cost_type weight(int e) const { return static_cast<cost_type>(edges_[e].weight); }
    std::uint64_t random() {
        std::uint64_t z = (rng_ += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    bool better(cost_type a, cost_type b) const {
        if constexpr (std::is_integral<Weight>::value)
            return a < b;
        else
            return a < b && b - a > std::max(std::fabs(a), std::fabs(b)) *
                                        (32 * std::numeric_limits<cost_type>::epsilon());
    }
    void shuffle(std::vector<int> &a) {
        for (std::size_t i = a.size(); i > 1; --i) std::swap(a[i - 1], a[random() % i]);
    }
    // Input is a forest. Keep the terminal component and remove nonterminal leaves.
    void prune(Tree &tree) {
        tree.cost = 0;
        if (terminals_.size() <= 1) {
            tree.edges.clear();
            return;
        }
        auto &degree = prune_degree_;
        auto &neighbor = prune_neighbor_;
        auto &queue = prune_queue_;
        degree.assign(n_, 0);
        neighbor.assign(n_, 0);
        queue.clear();
        // In a forest, a leaf's XOR of remaining neighbors is its sole neighbor.
        // This removes per-vertex adjacency allocations without changing leaf order.
        for (int e : tree.edges) {
            int u = edges_[e].u, v = edges_[e].v;
            ++degree[u];
            ++degree[v];
            neighbor[u] ^= v;
            neighbor[v] ^= u;
        }
        for (int v = 0; v < n_; ++v)
            if (!terminal_[v] && degree[v] <= 1) queue.push_back(v);
        for (std::size_t i = 0; i < queue.size(); ++i) {
            int v = queue[i];
            if (degree[v] == 0) continue;
            int u = neighbor[v];
            degree[v] = 0;
            neighbor[u] ^= v;
            if (--degree[u] == 1 && !terminal_[u]) queue.push_back(u);
        }
        std::size_t kept = 0;
        for (int e : tree.edges)
            if (degree[edges_[e].u] && degree[edges_[e].v]) {
                tree.edges[kept++] = e;
                tree.cost += weight(e);
            }
        tree.edges.resize(kept);
        std::sort(tree.edges.begin(), tree.edges.end());
    }
    void consider(const Tree &t) {
        if (better(t.cost, best_.cost)) {
            best_ = t;
            ++stats_.improvements;
        }
        if (best_.cost == 0) {
            optimal_ = true;
            lower_bound_ = 0;
        }
    }
    void begin_normalize(std::vector<int> &&ids) {
        norm_.reset();
        norm_.active = true;
        norm_.input.swap(ids);
        norm_.dsu.reset(n_);
    }
    void normalize_step() {
        if (norm_.cursor < norm_.input.size()) {
            int e = norm_.input[norm_.cursor++];
            norm_.heap.push({edges_[e].weight, e});
            return;
        }
        if (!norm_.heap.empty() || norm_.base_cursor < norm_.base.size()) {
            int e;
            if (norm_.base_cursor < norm_.base.size() &&
                (norm_.heap.empty() || weight(norm_.base[norm_.base_cursor]) <= norm_.heap.top().d))
                e = norm_.base[norm_.base_cursor++];
            else {
                e = norm_.heap.top().v;
                norm_.heap.pop();
            }
            if (norm_.dsu.join(edges_[e].u, edges_[e].v)) norm_.tree.edges.push_back(e);
            return;
        }
        prune(norm_.tree);
        norm_.done = true;
        norm_.active = false;
    }
    // One heap pop or one edge relaxation. The full scan position survives interruptions.
    template <bool Perturb = false> void shortest_step(Shortest &s) {
        if (s.vertex < 0) {
            if (s.heap.empty()) {
                s.done = true;
                return;
            }
            Item item = s.heap.top();
            s.heap.pop();
            if (item.d != s.d[item.v]) return;
            s.vertex = item.v;
            s.arc = 0;
            s.end = adj_[item.v].size();
            s.vertex_distance = item.d;
            return;
        }
        if (s.arc == s.end) {
            s.vertex = -1;
            return;
        }
        const auto &arc = adj_[s.vertex][s.arc++];
        int e = arc.edge, v = arc.vertex;
        cost_type nd =
            s.vertex_distance + (Perturb ? perturbed_[e] : static_cast<cost_type>(arc.weight));
        if (nd < s.d[v]) {
            s.d[v] = nd;
            s.pred[v] = e;
            s.heap.push({nd, v});
        }
    }
    void begin_construction(bool perturb, bool recombine, bool vertex_perturbation = false) {
        construction_started_ = true;
        construction_normalizing_ = false;
        construction_mode_ = perturb ? (cascading_ ? 2 : 1) : 0;
        weight_cursor_ = 0;
        constructed_.edges.clear();
        constructed_.cost = 0;
        in_tree_.assign(n_, 0);
        missing_ = static_cast<int>(terminals_.size());
        shortest_.reset(n_);
        int root = terminals_[perturb ? random() % terminals_.size() : 0];
        in_tree_[root] = 1;
        --missing_;
        shortest_.source(root);
        if (perturb) {
            perturbed_.resize(edges_.size());
            bias_.assign(edges_.size(), 0);
            vertex_factor_.clear();
            if (!cascading_ && vertex_perturbation && random() % 2) {
                vertex_factor_.resize(n_);
                for (int &f : vertex_factor_) f = 512 + static_cast<int>(random() % 1025);
            }
            if (cascading_) {
                std::vector<int> choices;
                for (std::size_t i = 0; i < elite_.size(); ++i)
                    if (elite_[i].edges != cascade_best_.edges)
                        choices.push_back(static_cast<int>(i));
                const Tree &b = elite_[choices[random() % choices.size()]];
                for (int e : cascade_best_.edges) ++bias_[e];
                for (int e : b.edges) ++bias_[e];
                ++stats_.recombinations;
            } else if (recombine && elite_.size() >= 2 && generation_ % 3 != 0) {
                const Tree &a = elite_[random() % elite_.size()];
                const Tree &b = elite_[random() % elite_.size()];
                for (int e : a.edges) ++bias_[e];
                for (int e : b.edges) ++bias_[e];
                ++stats_.recombinations;
            } else if (generation_ % 3 == 1) {
                for (int e : best_.edges) bias_[e] = 1;
            }
        }
    }
    void construction_step() {
        if (construction_mode_ && weight_cursor_ < edges_.size()) {
            std::size_t e = weight_cursor_++;
            int multiplier = 512 + static_cast<int>(random() % 1025);
            if (!vertex_factor_.empty())
                multiplier = (vertex_factor_[edges_[e].u] + vertex_factor_[edges_[e].v]) / 2;
            if (construction_mode_ == 2)
                multiplier = bias_[e] == 2   ? 1
                             : bias_[e] == 1 ? 100 + static_cast<int>(random() % 401)
                                             : 1000;
            else if (bias_[e] == 1)
                multiplier /= 2;
            else if (bias_[e] == 2)
                multiplier /= 4;
            perturbed_[e] = weight(static_cast<int>(e)) * multiplier;
            return;
        }
        if (construction_normalizing_) {
            normalize_step();
            if (norm_.done) {
                if (construction_mode_ || better(norm_.tree.cost, current_.cost))
                    current_ = norm_.tree;
                consider(current_);
                ++stats_.constructions;
                construction_started_ = false;
                phase_ = source_ == n_ ? 3 : 1;
                local_ready_ = false;
            }
            return;
        }
        if (!missing_) {
            begin_normalize(std::move(constructed_.edges));
            construction_normalizing_ = true;
            return;
        }
        if (shortest_.vertex >= 0 && shortest_.arc == 0 && terminal_[shortest_.vertex] &&
            !in_tree_[shortest_.vertex]) {
            int v = shortest_.vertex;
            // Strict relaxations ensure acyclic predecessor paths, including zero-weight edges.
            while (!in_tree_[v]) {
                int e = shortest_.pred[v];
                if (e < 0) throw std::logic_error("steiner_tree: broken construction path");
                constructed_.edges.push_back(e);
                in_tree_[v] = 1;
                if (terminal_[v]) --missing_;
                shortest_.source(v);
                v = other(e, v);
            }
            shortest_.vertex = -1;
            return;
        }
        if (construction_mode_)
            shortest_step<true>(shortest_);
        else
            shortest_step(shortest_);
        if (shortest_.done && missing_)
            throw std::logic_error("steiner_tree: unreachable terminal");
    }
    void apsp_step() {
        if (source_ == n_) {
            phase_ = 2;
            return;
        }
        if (dense_apsp_) {
            if (distance_.size() < static_cast<std::size_t>(n_)) {
                int row = static_cast<int>(distance_.size());
                distance_.emplace_back(n_, inf());
                distance_.back()[row] = 0;
                predecessor_.emplace_back(n_, -1);
                return;
            }
            if (apsp_edge_ < edges_.size()) {
                int id = static_cast<int>(apsp_edge_++);
                const auto &e = edges_[id];
                if (weight(id) < distance_[e.u][e.v]) {
                    distance_[e.u][e.v] = distance_[e.v][e.u] = weight(id);
                    predecessor_[e.u][e.v] = predecessor_[e.v][e.u] = id;
                }
                return;
            }
            // A Floyd row is the interruption unit. No row is exposed until all pivots finish.
            int i = floyd_row_, k = floyd_pivot_;
            cost_type via = distance_[i][k];
            auto *row = distance_[i].data();
            const auto *pivot = distance_[k].data();
            auto *pred = predecessor_[i].data();
            const auto *pivot_pred = predecessor_[k].data();
            if (i != k && via < inf())
                for (int j = 0; j < n_; ++j) {
                    cost_type value = via + pivot[j];
                    if (value < row[j]) {
                        row[j] = value;
                        pred[j] = pivot_pred[j];
                    }
                }
            if (++floyd_row_ == n_) {
                floyd_row_ = 0;
                ++floyd_pivot_;
                ++stats_.floyd_pivots;
                if (floyd_pivot_ == n_) {
                    source_ = n_;
                    stats_.shortest_path_rows = n_;
                }
            }
            return;
        }
        if (!row_started_) {
            shortest_.reset(n_);
            shortest_.source(source_);
            row_started_ = true;
        }
        shortest_step(shortest_);
        if (shortest_.done) {
            distance_.push_back(std::move(shortest_.d));
            predecessor_.push_back(std::move(shortest_.pred));
            ++source_;
            ++stats_.shortest_path_rows;
            row_started_ = false;
        }
    }
    void begin_metric(std::vector<int> nodes) {
        metric_.reset();
        metric_.active = true;
        metric_.nodes = std::move(nodes);
        int size = static_cast<int>(metric_.nodes.size());
        metric_.key.assign(size, inf());
        metric_.parent.assign(size, -1);
        metric_.remaining.resize(size);
        std::iota(metric_.remaining.begin(), metric_.remaining.end(), 0);
        if (size) {
            metric_.key[0] = 0;
            metric_.next = 0;
        }
    }
    template <bool ComponentMetric, class Distance> void metric_step_impl() {
        int size = static_cast<int>(metric_.nodes.size());
        if (metric_.count < size) {
            int position = metric_.next;
            if (position < 0) {
                metric_.done = true;
                return;
            }
            int p = metric_.remaining[position];
            metric_.remaining[position] = metric_.remaining.back();
            metric_.remaining.pop_back();
            ++metric_.count;
            metric_.total += metric_.key[p];
            int par = metric_.parent[p];
            if (par >= 0) {
                if constexpr (ComponentMetric) {
                    std::size_t pos = static_cast<std::size_t>(par) * size + p;
                    metric_.links.emplace_back(cut_.a[pos], cut_.b[pos]);
                } else
                    metric_.links.emplace_back(metric_.nodes[par], metric_.nodes[p]);
            }
            const auto *row = ComponentMetric
                                  ? cut_.distance.data() + static_cast<std::size_t>(p) * size
                                  : distance_[metric_.nodes[p]].data();
            // Select the next vertex while relaxing this row. Index tie breaking
            // matches the original ordered scan even though removal is now O(1).
            Distance minimum = static_cast<Distance>(inf());
            int next_vertex = -1;
            metric_.next = -1;
            for (std::size_t at = 0; at < metric_.remaining.size(); ++at) {
                int j = metric_.remaining[at];
                Distance value = static_cast<Distance>(row[ComponentMetric ? j : metric_.nodes[j]]);
                Distance key = static_cast<Distance>(metric_.key[j]);
                if (value < key) {
                    key = value;
                    metric_.key[j] = value;
                    metric_.parent[j] = p;
                }
                if (key < minimum || (key == minimum && j < next_vertex)) {
                    minimum = key;
                    next_vertex = j;
                    metric_.next = static_cast<int>(at);
                }
            }
            return;
        }
        if (metric_.link < metric_.links.size()) {
            auto link = metric_.links[metric_.link];
            if (metric_.path_vertex < 0) metric_.path_vertex = link.second;
            if (metric_.path_vertex == link.first) {
                ++metric_.link;
                metric_.path_vertex = -1;
                return;
            }
            int e = predecessor_[link.first][metric_.path_vertex];
            if (e < 0) throw std::logic_error("steiner_tree: broken metric path");
            metric_.edges.push_back(e);
            metric_.path_vertex = other(e, metric_.path_vertex);
            return;
        }
        metric_.done = true;
        metric_.active = false;
    }
    template <bool ComponentMetric = false> void metric_step() {
        if constexpr (std::is_integral<Weight>::value) {
            if (narrow_distances_) {
                metric_step_impl<ComponentMetric, std::uint64_t>();
                return;
            }
        }
        metric_step_impl<ComponentMetric, cost_type>();
    }
    void terminal_mst_step() {
        if (!terminal_mst_started_) {
            begin_metric(terminals_);
            terminal_mst_started_ = true;
        }
        if (!metric_.done) {
            metric_step();
            return;
        }
        if (!norm_.active && !norm_.done) begin_normalize(std::move(metric_.edges));
        if (!norm_.done) {
            normalize_step();
            return;
        }
        // MST(T) <= (2 - 2/|T|) OPT. Integer division rounds the lower bound upwards.
        const cost_type k = static_cast<cost_type>(terminals_.size());
        if constexpr (std::is_integral<Weight>::value)
            lower_bound_ = (metric_.total * k + 2 * (k - 1) - 1) / (2 * (k - 1));
        else
            lower_bound_ = metric_.total * k / (2 * (k - 1));
        for (int u : terminals_)
            for (int v : terminals_) lower_bound_ = std::max(lower_bound_, distance_[u][v]);
        consider(norm_.tree);
        if (better(norm_.tree.cost, current_.cost)) current_ = norm_.tree;
        if constexpr (std::is_integral<Weight>::value)
            if (lower_bound_ == best_.cost) optimal_ = true;
        norm_.reset();
        phase_ = 3;
        local_ready_ = false;
    }
    void prepare_local(const Options &options) {
        if (options.vertex_moves) {
            tree_sorted_ = current_.edges;
            std::sort(tree_sorted_.begin(), tree_sorted_.end(), [this](int a, int b) {
                return weight(a) != weight(b) ? weight(a) < weight(b) : a < b;
            });
        }
        degree_.assign(n_, 0);
        auto &tree = local_tree_;
        tree.resize(n_);
        for (auto &list : tree) list.clear();
        for (int e : current_.edges) {
            ++degree_[edges_[e].u];
            ++degree_[edges_[e].v];
            tree[edges_[e].u].push_back(e);
            tree[edges_[e].v].push_back(e);
        }
        keys_ = terminals_;
        additions_.clear();
        removals_.clear();
        auto &is_key = key_flags_;
        is_key = terminal_;
        for (int v = 0; v < n_; ++v)
            if (!terminal_[v] && degree_[v] >= 3) {
                keys_.push_back(v);
                removals_.push_back(v);
                is_key[v] = 1;
            }
        for (int v = 0; v < n_; ++v)
            if (!is_key[v] && distance_[terminals_[0]][v] != inf()) additions_.push_back(v);
        shuffle(additions_);
        shuffle(removals_);
        vertex_additions_.clear();
        vertex_removals_.clear();
        for (int v : additions_)
            if (options.vertex_moves && !degree_[v]) {
                int neighbors = 0;
                for (const auto &arc : adj_[v])
                    if (degree_[arc.vertex]) ++neighbors;
                if (neighbors >= 2) vertex_additions_.push_back(v);
            }
        for (int v = 0; v < n_; ++v)
            if (options.vertex_moves && degree_[v] && !terminal_[v]) vertex_removals_.push_back(v);
        shuffle(vertex_removals_);
        std::size_t path_count = 0;
        key_incident_.resize(n_);
        for (auto &list : key_incident_) list.clear();
        edge_path_.resize(edges_.size());
        for (int v : keys_)
            for (int first : tree[v]) {
                int e = first, u = other(e, v);
                auto &path = path_buffer_;
                path.clear();
                path.push_back(e);
                while (!is_key[u]) {
                    int next = tree[u][0] == e ? tree[u][1] : tree[u][0];
                    e = next;
                    u = other(e, u);
                    path.push_back(e);
                }
                if (v < u) {
                    int id = static_cast<int>(path_count++);
                    for (int edge : path) edge_path_[edge] = id;
                    if (key_paths_.size() < path_count) key_paths_.emplace_back();
                    key_paths_[id].swap(path);
                    key_incident_[v].push_back(id);
                    key_incident_[u].push_back(id);
                }
            }
        key_paths_.resize(path_count);
        exchange_.reset();
        auto &x = exchange_;
        x.parent.assign(n_, -1);
        x.depth.assign(n_, 0);
        x.parent_path.assign(n_, -1);
        x.interior_path.assign(n_, -1);
        x.path_cost.assign(key_paths_.size(), 0);
        for (std::size_t p = 0; p < key_paths_.size(); ++p)
            for (int e : key_paths_[p]) {
                x.path_cost[p] += weight(e);
                if (!is_key[edges_[e].u]) x.interior_path[edges_[e].u] = static_cast<int>(p);
                if (!is_key[edges_[e].v]) x.interior_path[edges_[e].v] = static_cast<int>(p);
            }
        if constexpr (std::is_integral<Weight>::value) x.bound.resize(n_);
        x.vertices.push_back(terminals_[0]);
        x.parent[terminals_[0]] = terminals_[0];
        if constexpr (std::is_integral<Weight>::value) x.bound[terminals_[0]] = 0;
        for (std::size_t i = 0; i < x.vertices.size(); ++i) {
            int v = x.vertices[i];
            for (int e : tree[v]) {
                int u = other(e, v);
                if (x.parent[u] >= 0) continue;
                x.parent[u] = v;
                x.depth[u] = x.depth[v] + 1;
                x.parent_path[u] = edge_path_[e];
                if constexpr (std::is_integral<Weight>::value)
                    x.bound[u] = std::max(x.bound[v], x.path_cost[edge_path_[e]]);
                x.vertices.push_back(u);
            }
        }
        local_stage_ = 0;
        candidate_ = 0;
        candidate_kind_ = 0;
        local_ready_ = true;
        metric_.reset();
        norm_.reset();
        cut_.reset();
    }
    void begin_induced(int vertex, bool removal = false) {
        induced_.reset();
        induced_.vertex = vertex;
        induced_.removal = removal;
        candidate_kind_ = 3;
    }
    void induced_step() {
        auto &job = induced_;
        bool insertion = job.vertex >= 0 && !job.removal;
        std::size_t count = insertion ? adj_[job.vertex].size() : edges_.size();
        if (job.cursor < count) {
            int id =
                insertion ? adj_[job.vertex][job.cursor++].edge : static_cast<int>(job.cursor++);
            const auto &e = edges_[id];
            if (insertion ? degree_[other(id, job.vertex)] != 0
                          : degree_[e.u] && degree_[e.v] && e.u != job.vertex && e.v != job.vertex)
                job.edges.push_back(id);
            return;
        }
        begin_normalize(std::move(job.edges));
        if (insertion) norm_.base = tree_sorted_;
        job.built = true;
    }
    bool normalized_feasible() {
        int root = norm_.dsu.find(terminals_[0]);
        for (int t : terminals_)
            if (norm_.dsu.find(t) != root) return false;
        return true;
    }
    // Evaluate every key-path cut in one pass over tree vertex pairs. The tree path
    // tells which cuts separate the endpoints; endpoints inside a deleted path are excluded.
    void exchange_step() {
        auto &x = exchange_;
        if (x.i < x.vertices.size()) {
            if (x.j == x.vertices.size()) {
                ++x.i;
                x.j = x.i + 1;
                return;
            }
            int a = x.vertices[x.i], b = x.vertices[x.j++], u = a, v = b;
            cost_type replacement = distance_[a][b];
            // Every key path on a--b also lies on root--a or root--b. Their
            // prefix maxima give a cheap upper bound, including ineligible paths.
            // Keep the original walk (and its tie order) for potentially improving pairs.
            if constexpr (std::is_integral<Weight>::value)
                if (std::max(x.bound[a], x.bound[b]) - replacement <= x.gain) return;
            while (u != v) {
                if (x.depth[u] < x.depth[v]) std::swap(u, v);
                int path = x.parent_path[u];
                if (path != x.interior_path[a] && path != x.interior_path[b]) {
                    cost_type gain = x.path_cost[path] - replacement;
                    if (better(x.gain, gain)) {
                        x.gain = gain;
                        x.path = path;
                        x.a = a;
                        x.b = b;
                    }
                }
                u = x.parent[u];
            }
            return;
        }
        x.done = true;
        if (x.path < 0) {
            candidate_ = 1;
            return;
        }
        begin_metric({});
        metric_.links.emplace_back(x.a, x.b);
        for (int e : current_.edges)
            if (edge_path_[e] != x.path) metric_.edges.push_back(e);
        candidate_kind_ = 1;
    }
    void begin_cut(const std::vector<int> &paths, bool steiner_reconnection) {
        cut_.reset();
        cut_.active = true;
        cut_.star_enabled = steiner_reconnection;
        auto &removed = cut_.removed;
        for (int p : paths)
            removed.insert(removed.end(), key_paths_[p].begin(), key_paths_[p].end());
        std::sort(removed.begin(), removed.end());
        auto &d = cut_.dsu;
        d.reset(n_);
        auto &active = cut_.present;
        active = terminal_;
        std::size_t gone = 0;
        for (int e : current_.edges) {
            while (gone < removed.size() && removed[gone] < e) ++gone;
            if (gone < removed.size() && removed[gone] == e) continue;
            cut_.retained.push_back(e);
            d.join(edges_[e].u, edges_[e].v);
            active[edges_[e].u] = active[edges_[e].v] = 1;
        }
        auto &index = cut_.index;
        index.assign(n_, -1);
        cut_.component.assign(n_, -1);
        for (int v = 0; v < n_; ++v)
            if (active[v]) {
                int root = d.find(v);
                if (index[root] < 0) index[root] = cut_.components++;
                cut_.component[v] = index[root];
                cut_.vertices.push_back(v);
            }
        std::size_t size = static_cast<std::size_t>(cut_.components) * cut_.components;
        cut_.distance.assign(size, inf());
        cut_.a.assign(size, -1);
        cut_.b.assign(size, -1);
        if (cut_.star_enabled && cut_.components >= 3) {
            cut_.to_vertex.assign(static_cast<std::size_t>(cut_.components) * n_, inf());
            cut_.nearest.assign(static_cast<std::size_t>(cut_.components) * n_, -1);
        }
    }
    template <class Distance> void cut_step_impl() {
        if (cut_.row < cut_.vertices.size()) {
            int u = cut_.vertices[cut_.row++], a = cut_.component[u], k = cut_.components;
            const auto *source = distance_[u].data();
            auto *row = cut_.distance.data() + static_cast<std::size_t>(a) * k;
            auto *from = cut_.a.data() + static_cast<std::size_t>(a) * k;
            auto *to = cut_.b.data() + static_cast<std::size_t>(a) * k;
            for (int v : cut_.vertices) {
                int b = cut_.component[v];
                if (a == b) continue;
                if (static_cast<Distance>(source[v]) < static_cast<Distance>(row[b])) {
                    row[b] = source[v];
                    from[b] = u;
                    to[b] = v;
                }
            }
            if (cut_.star_enabled && k >= 3) {
                auto *star = cut_.to_vertex.data() + static_cast<std::size_t>(a) * n_;
                auto *nearest = cut_.nearest.data() + static_cast<std::size_t>(a) * n_;
                for (int v = 0; v < n_; ++v)
                    if (static_cast<Distance>(source[v]) < static_cast<Distance>(star[v])) {
                        star[v] = source[v];
                        nearest[v] = u;
                    }
            }
            return;
        }
        if (!metric_.active && !metric_.done) {
            std::vector<int> nodes(cut_.components);
            std::iota(nodes.begin(), nodes.end(), 0);
            begin_metric(std::move(nodes));
            metric_.edges = cut_.retained;
        }
        if (cut_.star_enabled && cut_.components >= 3 && !cut_.star_checked &&
            metric_.count == cut_.components) {
            if (cut_.star_vertex == 0) cut_.star_best = metric_.total;
            if (cut_.star_vertex < n_) {
                int v = cut_.star_vertex++;
                cost_type total = 0;
                for (int a = 0; a < cut_.components; ++a) {
                    cost_type d = cut_.to_vertex[static_cast<std::size_t>(a) * n_ + v];
                    if (d == inf()) {
                        total = inf();
                        break;
                    }
                    total += d;
                }
                if (better(total, cut_.star_best)) {
                    cut_.star_best = total;
                    cut_.star_center = v;
                }
                return;
            }
            int center = cut_.star_center;
            if (center >= 0) {
                metric_.links.clear();
                for (int a = 0; a < cut_.components; ++a)
                    metric_.links.emplace_back(
                        center, cut_.nearest[static_cast<std::size_t>(a) * n_ + center]);
            }
            cut_.star_checked = true;
            return;
        }
        if (!metric_.done) {
            metric_step<true>();
            return;
        }
        cut_.done = true;
    }
    void cut_step() {
        if constexpr (std::is_integral<Weight>::value) {
            if (narrow_distances_) {
                cut_step_impl<std::uint64_t>();
                return;
            }
        }
        cut_step_impl<cost_type>();
    }
    void add_elite(const Tree &tree) {
        for (const auto &e : elite_)
            if (e.edges == tree.edges) return;
        if (elite_.size() < 16)
            elite_.push_back(tree);
        else {
            std::size_t replace = elite_.size(), closest = std::numeric_limits<std::size_t>::max();
            for (std::size_t i = 0; i < elite_.size(); ++i)
                if (tree.cost <= elite_[i].cost) {
                    std::size_t common = 0, a = 0, b = 0;
                    while (a < tree.edges.size() && b < elite_[i].edges.size()) {
                        if (tree.edges[a] == elite_[i].edges[b]) {
                            ++common;
                            ++a;
                            ++b;
                        } else if (tree.edges[a] < elite_[i].edges[b])
                            ++a;
                        else
                            ++b;
                    }
                    std::size_t delta = tree.edges.size() + elite_[i].edges.size() - 2 * common;
                    if (delta < closest) {
                        closest = delta;
                        replace = i;
                    }
                }
            if (replace < elite_.size()) elite_[replace] = tree;
        }
    }
    void restart(const Options &options) {
        if (options.recombination && options.cascaded_recombination) {
            if (!cascading_) {
                fresh_ = cascade_best_ = current_;
                cascade_failures_ = 0;
            } else if (better(current_.cost, cascade_best_.cost))
                cascade_best_ = current_;
            else
                ++cascade_failures_;
            bool different = false;
            for (const auto &e : elite_)
                if (e.edges != cascade_best_.edges) different = true;
            if (different && cascade_failures_ < 3) {
                cascading_ = true;
                phase_ = 0;
                begin_construction(true, true);
                return;
            }
            add_elite(fresh_);
            add_elite(cascade_best_);
            cascading_ = false;
            ++generation_;
            phase_ = 0;
            begin_construction(true, false, options.vertex_perturbation);
            return;
        }
        cascading_ = false;
        add_elite(current_);
        add_elite(best_);
        ++generation_;
        phase_ = 0;
        begin_construction(true, options.recombination, options.vertex_perturbation);
    }
    void local_step(const Options &options) {
        if (!options.local_search) {
            restart(options);
            return;
        }
        if (!local_ready_) {
            prepare_local(options);
            return;
        }
        if (candidate_kind_) {
            if (candidate_kind_ == 1 && !metric_.done) {
                metric_step();
                return;
            }
            if (candidate_kind_ == 2 && !cut_.done) {
                cut_step();
                return;
            }
            if (candidate_kind_ == 3 && !induced_.built) {
                induced_step();
                return;
            }
            if (!norm_.active && !norm_.done) begin_normalize(std::move(metric_.edges));
            if (!norm_.done) {
                normalize_step();
                return;
            }
            ++stats_.candidates;
            if (better(norm_.tree.cost, current_.cost) &&
                (candidate_kind_ != 3 || normalized_feasible())) {
                if (candidate_kind_ == 3) {
                    if (induced_.vertex < 0)
                        ++stats_.induced_improvements;
                    else if (induced_.removal)
                        ++stats_.vertex_removals;
                    else
                        ++stats_.vertex_insertions;
                }
                current_ = norm_.tree;
                consider(current_);
                ++stats_.local_moves;
                local_ready_ = false;
            }
            ++candidate_;
            candidate_kind_ = 0;
            metric_.reset();
            norm_.reset();
            cut_.reset();
            return;
        }
        if (local_stage_ == 0) {
            if (options.vertex_moves) {
                if (!candidate_) {
                    begin_induced(-1);
                    return;
                }
                if (candidate_ <= vertex_additions_.size()) {
                    begin_induced(vertex_additions_[candidate_ - 1]);
                    return;
                }
            }
        } else if (local_stage_ == 1) {
            if (options.fast_key_paths) {
                if (!exchange_.done) exchange_step();
                if (!exchange_.done || candidate_kind_) return;
            } else if (candidate_ < key_paths_.size()) {
                begin_cut({static_cast<int>(candidate_)}, options.steiner_reconnection);
                candidate_kind_ = 2;
                return;
            }
        } else if (local_stage_ == 2) {
            if (candidate_ < removals_.size()) {
                begin_cut(key_incident_[removals_[candidate_]], options.steiner_reconnection);
                candidate_kind_ = 2;
                return;
            }
        } else if (local_stage_ == 3) {
            if (candidate_ < additions_.size()) {
                auto nodes = keys_;
                nodes.push_back(additions_[candidate_]);
                begin_metric(std::move(nodes));
                candidate_kind_ = 1;
                return;
            }
        } else if (local_stage_ == 4) {
            if (candidate_ < removals_.size()) {
                auto nodes = keys_;
                nodes.erase(std::find(nodes.begin(), nodes.end(), removals_[candidate_]));
                begin_metric(std::move(nodes));
                candidate_kind_ = 1;
                return;
            }
        } else if (local_stage_ == 5) {
            if (options.vertex_moves && candidate_ < vertex_removals_.size()) {
                begin_induced(vertex_removals_[candidate_], true);
                return;
            }
        } else {
            restart(options);
            return;
        }
        ++local_stage_;
        candidate_ = 0;
    }
    void decide_exact(const Options &o) {
        // Once started, the DP survives calls; disabling it merely pauses it.
        if (exact_decided_ || o.exact_terminal_limit == 0) return;
        exact_decided_ = true;
        int k = static_cast<int>(terminals_.size()) - 1;
        if (k < 1 || static_cast<int>(terminals_.size()) > o.exact_terminal_limit || k > 20)
            return;
        std::size_t masks = std::size_t(1) << k;
        if (static_cast<std::size_t>(n_) > o.exact_memory_limit / sizeof(Cell) / masks) return;
        long double work = n_;
        for (int i = 0; i < k; ++i) work *= 3;
        if (work > 40000000) return;
        exact_.active = true;
        exact_.count = k;
        exact_.full = static_cast<int>(masks) - 1;
        exact_.table.resize(masks);
    }
    void exact_step() {
        auto &x = exact_;
        if (x.stage == 0) {
            x.table[x.mask].resize(n_);
            x.vertex = 0;
            if ((x.mask & (x.mask - 1)) == 0) {
                int bit = 0;
                while ((1 << bit) != x.mask) ++bit;
                x.table[x.mask][terminals_[bit + 1]].cost = 0;
                x.stage = 2;
            } else {
                x.sub = (x.mask - 1) & x.mask;
                x.stage = 1;
            }
            return;
        }
        if (x.stage == 1) {
            if (!x.sub) {
                x.stage = 2;
                return;
            }
            int rest = x.mask ^ x.sub;
            if (x.sub > rest) {
                x.sub = (x.sub - 1) & x.mask;
                x.vertex = 0;
                return;
            }
            int v = x.vertex++;
            cost_type value = x.table[x.sub][v].cost + x.table[rest][v].cost;
            if (value < x.table[x.mask][v].cost) {
                auto &cell = x.table[x.mask][v];
                cell.cost = value;
                cell.split = x.sub;
            }
            if (x.vertex == n_) {
                x.vertex = 0;
                x.sub = (x.sub - 1) & x.mask;
            }
            return;
        }
        if (x.stage == 2) {
            x.shortest.reset(n_);
            for (int v = 0; v < n_; ++v)
                if (x.table[x.mask][v].cost < inf()) x.shortest.source(v, x.table[x.mask][v].cost);
            x.stage = 3;
            return;
        }
        if (x.stage == 3) {
            shortest_step(x.shortest);
            if (!x.shortest.done) return;
            for (int v = 0; v < n_; ++v)
                if (x.shortest.pred[v] >= 0) {
                    auto &cell = x.table[x.mask][v];
                    cell.cost = x.shortest.d[v];
                    cell.split = 0;
                    cell.edge = x.shortest.pred[v];
                }
            ++stats_.exact_masks;
            if (x.mask == x.full) {
                x.trace.emplace_back(x.full, terminals_[0]);
                x.stage = 4;
            } else {
                ++x.mask;
                x.stage = 0;
            }
            return;
        }
        if (x.stage == 4) {
            if (!x.trace.empty()) {
                auto pair = x.trace.back();
                x.trace.pop_back();
                const Cell &cell = x.table[pair.first][pair.second];
                if (cell.edge >= 0) {
                    x.edges.push_back(cell.edge);
                    x.trace.emplace_back(pair.first, other(cell.edge, pair.second));
                } else if (cell.split) {
                    x.trace.emplace_back(cell.split, pair.second);
                    x.trace.emplace_back(pair.first ^ cell.split, pair.second);
                }
                return;
            }
            // A DP decomposition is a union of paths. Remove cycles and duplicate edges.
            DSU d;
            d.reset(n_);
            std::vector<int> forest;
            std::sort(x.edges.begin(), x.edges.end(), [this](int a, int b) {
                return weight(a) != weight(b) ? weight(a) < weight(b) : a < b;
            });
            for (int e : x.edges)
                if (d.join(edges_[e].u, edges_[e].v)) forest.push_back(e);
            Tree tree;
            tree.edges = std::move(forest);
            prune(tree);
            consider(tree);
            lower_bound_ = x.table[x.full][terminals_[0]].cost;
            if constexpr (std::is_integral<Weight>::value) {
                if (tree.cost != lower_bound_)
                    throw std::logic_error("steiner_tree: inconsistent exact reconstruction");
            }
            optimal_ = true;
            x.finished = true;
            x.active = false;
        }
    }
    void heuristic_step(const Options &o) {
        if (phase_ == 0) {
            if (!construction_started_) begin_construction(false, false);
            construction_step();
        } else if (phase_ == 1)
            apsp_step();
        else if (phase_ == 2) {
            // The previous construction's normalizer must not masquerade as this one.
            if (!terminal_mst_started_) norm_.reset();
            terminal_mst_step();
        } else if (phase_ == 3)
            local_step(o);
        else {
            normalize_step();
            if (norm_.done) {
                consider(norm_.tree);
                lower_bound_ = best_.cost;
                optimal_ = true;
            }
        }
    }

public:
    // Construction validates and copies the graph, and repairs the supplied edge fragments.
    // Setup uses DSU plus O(n log n) output ordering, and is outside solve()'s time limit.
    Solver(int n, const std::vector<Edge> &edges, const std::vector<int> &terminals,
           const std::vector<int> &initial_edges = {}, std::uint64_t seed = 1)
        : n_(n), edges_(edges), rng_(seed) {
        if (n < 0 || edges.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("steiner_tree: invalid graph size");
        adj_.resize(n);
        terminal_.assign(n, 0);
        std::vector<int> adjacency_size(n, 0);
        dense_apsp_ = edges.size() > static_cast<std::size_t>(n) * static_cast<std::size_t>(n) / 8;
        DSU graph;
        graph.reset(n);
        cost_type max_edge_weight = 0;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto &e = edges[i];
            if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n ||
                !std::isfinite(static_cast<long double>(e.weight)) || e.weight < 0)
                throw std::invalid_argument(
                    "steiner_tree: invalid endpoint or nonnegative finite weight required");
            if constexpr (!std::is_integral<Weight>::value)
                if (e.weight > std::numeric_limits<cost_type>::max() / (65536.0L * std::max(1, n)))
                    throw std::overflow_error("steiner_tree: weight too large for safe sums");
            if (e.u != e.v) {
                if constexpr (std::is_integral<Weight>::value)
                    max_edge_weight = std::max(max_edge_weight, static_cast<cost_type>(e.weight));
                ++adjacency_size[e.u];
                ++adjacency_size[e.v];
                graph.join(e.u, e.v);
            }
        }
        // Finite shortest paths have at most n-1 edges. In this range their low
        // 64 bits preserve ordering; inf() maps to UINT64_MAX. Sums and storage
        // remain wide, and large integer / floating weights use the wide path.
        if constexpr (std::is_integral<Weight>::value)
            narrow_distances_ =
                max_edge_weight * std::max(0, n - 1) < std::numeric_limits<std::uint64_t>::max();
        for (int v = 0; v < n; ++v) adj_[v].reserve(adjacency_size[v]);
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto &e = edges[i];
            if (e.u != e.v) {
                adj_[e.u].push_back({static_cast<int>(i), e.v, e.weight});
                adj_[e.v].push_back({static_cast<int>(i), e.u, e.weight});
            }
        }
        for (int t : terminals) {
            if (t < 0 || t >= n) throw std::invalid_argument("steiner_tree: invalid terminal");
            if (!terminal_[t]) {
                terminal_[t] = 1;
                terminals_.push_back(t);
            }
        }
        for (int e : initial_edges)
            if (e < 0 || static_cast<std::size_t>(e) >= edges.size())
                throw std::invalid_argument("steiner_tree: invalid initial edge ID");
        feasible_ = true;
        if (terminals_.size() <= 1) {
            optimal_ = true;
            return;
        }
        int root = graph.find(terminals_[0]);
        for (int t : terminals_)
            if (graph.find(t) != root) {
                feasible_ = false;
                return;
            }
        DSU forest;
        forest.reset(n);
        std::vector<int> ids;
        for (int e : initial_edges)
            if (graph.find(edges_[e].u) == root && forest.join(edges_[e].u, edges_[e].v))
                ids.push_back(e);
        for (std::size_t e = 0; e < edges.size(); ++e)
            if (graph.find(edges[e].u) == root && forest.join(edges[e].u, edges[e].v))
                ids.push_back(static_cast<int>(e));
        current_.edges = std::move(ids);
        prune(current_);
        best_ = current_;
        int vertices = 0;
        std::size_t component_edges = 0;
        all_terminals_ = true;
        for (int v = 0; v < n; ++v)
            if (graph.find(v) == root) {
                ++vertices;
                if (!terminal_[v]) all_terminals_ = false;
                component_edges += adj_[v].size();
            }
        if (best_.cost == 0 || component_edges / 2 == static_cast<std::size_t>(vertices - 1)) {
            optimal_ = true;
            lower_bound_ = best_.cost;
        } else if (all_terminals_) {
            std::vector<int> all(edges.size());
            std::iota(all.begin(), all.end(), 0);
            begin_normalize(std::move(all));
            phase_ = 4;
        }
    }

    Result result() const {
        Result r;
        r.edges = best_.edges;
        r.cost = best_.cost;
        r.lower_bound = lower_bound_;
        r.feasible = feasible_;
        r.optimal = optimal_;
        r.stats = stats_;
        r.reason = !feasible_ ? StopReason::infeasible
                   : optimal_ ? StopReason::optimal
                              : StopReason::time_limit;
        return r;
    }

    Result solve(const Options &o = {}) {
        const auto start = Clock::now();
        if (std::isnan(o.time_limit) || o.time_limit < 0 || !std::isfinite(o.safety_margin) ||
            o.safety_margin < 0 || o.exact_terminal_limit < 0)
            throw std::invalid_argument("steiner_tree: invalid solve options");
        if (!exact_.active && !exact_.finished) exact_decided_ = false;
        auto stop = o.deadline;
        const double capacity =
            std::chrono::duration<double>(Clock::time_point::max() - start).count();
        if (o.time_limit < capacity)
            stop = std::min(stop, start + std::chrono::duration_cast<Clock::duration>(
                                              std::chrono::duration<double>(o.time_limit)));
        if (stop > start && stop != Clock::time_point::max()) {
            double seconds = std::chrono::duration<double>(stop - start).count();
            double reserve =
                std::min(seconds * 0.1, std::max(o.safety_margin, std::min(0.002, seconds * 0.01)));
            stop -=
                std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(reserve));
        }
        std::uint64_t done = 0;
        StopReason reason = StopReason::time_limit;
        while (feasible_ && !optimal_) {
            if (done == o.max_work) {
                reason = StopReason::work_limit;
                break;
            }
            if (Clock::now() >= stop) break;
            const auto end = done + std::min<std::uint64_t>(32, o.max_work - done);
            do {
                if (stats_.constructions && !all_terminals_) decide_exact(o);
                if (exact_.active && o.exact_terminal_limit && (stats_.work & 3U) == 0)
                    exact_step();
                else
                    heuristic_step(o);
                ++done;
                ++stats_.work;
            } while (done < end && !optimal_);
        }
        Result r = result();
        if (feasible_ && !optimal_) r.reason = reason;
        r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        stats_.elapsed = r.stats.elapsed;
        return r;
    }
    Result solve(double seconds) {
        Options o;
        o.time_limit = seconds;
        return solve(o);
    }
};

// One-shot helper: graph setup consumes the same absolute deadline as search. Setup itself
// is not preemptible; reuse Solver for tight repeated calls. See the timing measurements.
template <class Weight>
typename Solver<Weight>::Result
solve(int n, const std::vector<typename Solver<Weight>::Edge> &edges,
      const std::vector<int> &terminals, Options options = {},
      const std::vector<int> &initial_edges = {}, std::uint64_t seed = 1) {
    auto start = Clock::now();
    if (std::isnan(options.time_limit) || options.time_limit < 0)
        throw std::invalid_argument("steiner_tree: invalid time_limit");
    double cap = std::chrono::duration<double>(Clock::time_point::max() - start).count();
    if (options.time_limit < cap)
        options.deadline = std::min(
            options.deadline, start + std::chrono::duration_cast<Clock::duration>(
                                          std::chrono::duration<double>(options.time_limit)));
    Solver<Weight> solver(n, edges, terminals, initial_edges, seed);
    auto r = solver.solve(options);
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    return r;
}
} // namespace steiner_tree
#endif
