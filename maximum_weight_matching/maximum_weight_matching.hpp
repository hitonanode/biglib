#ifndef MAXIMUM_WEIGHT_MATCHING_HPP
#define MAXIMUM_WEIGHT_MATCHING_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace maximum_weight_matching {

using Integer = __int128_t;

inline std::string to_string(Integer value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    __uint128_t magnitude = negative ? __uint128_t(-(value + 1)) + 1 : __uint128_t(value);
    std::string result;
    while (magnitude) {
        result.push_back(char('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

enum class Status { Optimal, Limit, NumericalError };
enum class Algorithm { Automatic, Dense, Sparse };

struct Options {
    double time_limit = std::numeric_limits<double>::infinity(); // Per call, seconds.
    std::uint64_t step_limit = std::numeric_limits<std::uint64_t>::max();
};

struct Statistics {
    std::uint64_t steps = 0, phases = 0, augmentations = 0, relocations = 0;
    std::uint64_t contractions = 0, expansions = 0, edge_scans = 0, dual_updates = 0;
    int retained_seed_edges = 0, initial_search_edges = 0;
    Algorithm algorithm = Algorithm::Automatic;
    double elapsed = 0; // Cumulative time spent inside solve().
};

namespace detail {

template <class Weight> struct Types {
    static_assert((std::is_integral<Weight>::value && std::is_signed<Weight>::value &&
                   sizeof(Weight) <= 8) ||
                      std::is_floating_point<Weight>::value,
                  "Use a signed integer up to 64 bits, or a floating point type");

public:
    static constexpr bool exact = std::is_integral<Weight>::value;
    using Value = typename std::conditional<exact, Integer, long double>::type;
    struct Edge {
        int u, v;
        Weight weight;
    };
    struct Result {
        Status status = Status::Limit;
        std::vector<int> mate, mate_edge, edges;
        Value weight = 0;
        Value dual_bound_twice = 0; // Exact upper bound / 2 for integer inputs.
        Value tolerance = 0;        // Floating point arithmetic only.
        Statistics stats;
    };
    struct OddSet {
        std::vector<int> vertices;
        Value dual_twice;
    };
    struct Certificate {
        std::vector<Value> vertex_dual_twice;
        std::vector<OddSet> odd_sets;
    };
};

// A deterministic AVL sequence. Leaves have stable vertex handles; joins and
// splits change only O(log n) internal nodes. Both minima share the lazy offset.
template <class Number> class VertexSequence {
public:
    struct Node {
        int left = -1, right = -1, parent = -1, height = 1, size = 1, owner = -1;
        int dual_arg = -1, incoming_arg = -1;
        Number dual = 0, incoming = 0, add = 0;
    };
    std::vector<Node> data;
    std::vector<int> free_nodes;
    void reset(const std::vector<Number> &values) {
        const int n = int(values.size());
        data.assign(n ? 2 * n - 1 : 0, Node{});
        free_nodes.clear();
        for (int v = 0; v < n; ++v) {
            data[v].dual = values[v];
            data[v].dual_arg = v;
            data[v].owner = v;
        }
        for (int v = int(data.size()) - 1; v >= n; --v) free_nodes.push_back(v);
    }
    int root(int v) const {
        while (data[v].parent >= 0) v = data[v].parent;
        return v;
    }
    int rank(int v) const {
        int position = 0;
        while (data[v].parent >= 0) {
            const int p = data[v].parent;
            if (data[p].right == v) position += data[data[p].left].size;
            v = p;
        }
        return position;
    }
    Number vertex_dual(int v) const {
        Number value = data[v].dual;
        while (data[v].parent >= 0) {
            v = data[v].parent;
            value += data[v].add;
        }
        return value;
    }
    void add_to(int root, Number value) {
        if (root < 0 || value == 0) return;
        Node &p = data[root];
        p.dual += value;
        if (p.incoming_arg >= 0) p.incoming += value;
        p.add += value;
    }
    int set_incoming(int v, bool present, Number source_key = 0) {
        std::array<int, 96> path;
        int count = 0;
        for (int x = v; x >= 0; x = data[x].parent) path[count++] = x;
        for (int i = count - 1; i >= 0; --i) push(path[i]);
        data[v].incoming_arg = present ? v : -1;
        data[v].incoming = present ? data[v].dual + source_key : 0;
        for (int i = 1; i < count; ++i) pull(path[i]);
        return path[count - 1];
    }
    template <class Function> void for_each(int root, Function function) const {
        if (root < 0) return;
        std::array<int, 96> stack;
        int count = 0;
        stack[count++] = root;
        while (count) {
            const int x = stack[--count];
            if (data[x].left < 0)
                function(x);
            else {
                stack[count++] = data[x].right;
                stack[count++] = data[x].left;
            }
        }
    }
    int concatenate(int a, int b) {
        if (a < 0) return b;
        if (b < 0) return a;
        if (height(a) > height(b) + 1) {
            push(a);
            const int l = data[a].left, r = data[a].right;
            data[l].parent = data[r].parent = -1;
            const int joined = concatenate(r, b);
            children(a, l, joined);
            return balance(a);
        }
        if (height(b) > height(a) + 1) {
            push(b);
            const int l = data[b].left, r = data[b].right;
            data[l].parent = data[r].parent = -1;
            const int joined = concatenate(a, l);
            children(b, joined, r);
            return balance(b);
        }
        const int p = free_nodes.back();
        free_nodes.pop_back();
        data[p] = Node{};
        children(p, a, b);
        return p;
    }
    std::pair<int, int> split(int root, int prefix_size) {
        if (prefix_size == 0) return {-1, root};
        if (prefix_size == data[root].size) return {root, -1};
        push(root);
        const int l = data[root].left, r = data[root].right, left_size = data[l].size;
        data[l].parent = data[r].parent = -1;
        free_nodes.push_back(root);
        if (prefix_size == left_size) return {l, r};
        if (prefix_size < left_size) {
            const auto parts = split(l, prefix_size);
            return {parts.first, concatenate(parts.second, r)};
        }
        const auto parts = split(r, prefix_size - left_size);
        return {concatenate(l, parts.first), parts.second};
    }

private:
    int height(int x) const { return x < 0 ? 0 : data[x].height; }
    void push(int x) {
        if (data[x].left >= 0) {
            add_to(data[x].left, data[x].add);
            add_to(data[x].right, data[x].add);
        }
        data[x].add = 0;
    }
    void pull(int x) {
        Node &p = data[x];
        const Node &a = data[p.left];
        const Node &b = data[p.right];
        p.size = a.size + b.size;
        p.height = 1 + std::max(a.height, b.height);
        const Node &d = (a.dual < b.dual || (a.dual == b.dual && a.dual_arg < b.dual_arg)) ? a : b;
        p.dual = d.dual;
        p.dual_arg = d.dual_arg;
        const Node &c =
            b.incoming_arg < 0 || (a.incoming_arg >= 0 &&
                                   (a.incoming < b.incoming ||
                                    (a.incoming == b.incoming && a.incoming_arg < b.incoming_arg)))
                ? a
                : b;
        p.incoming = c.incoming;
        p.incoming_arg = c.incoming_arg;
    }
    void children(int x, int a, int b) {
        data[x].left = a;
        data[x].right = b;
        data[x].parent = -1;
        data[a].parent = data[b].parent = x;
        pull(x);
    }
    int rotate_left(int x) {
        push(x);
        const int y = data[x].right;
        push(y);
        const int a = data[x].left, b = data[y].left, c = data[y].right;
        children(x, a, b);
        children(y, x, c);
        return y;
    }
    int rotate_right(int x) {
        push(x);
        const int y = data[x].left;
        push(y);
        const int a = data[y].left, b = data[y].right, c = data[x].right;
        children(x, b, c);
        children(y, a, x);
        return y;
    }
    int balance(int x) {
        const int l = data[x].left, r = data[x].right;
        if (height(l) > height(r) + 1) {
            if (height(data[l].left) < height(data[l].right)) {
                const int rotated = rotate_left(l);
                children(x, rotated, r);
            }
            return rotate_right(x);
        }
        if (height(r) > height(l) + 1) {
            if (height(data[r].right) < height(data[r].left)) {
                const int rotated = rotate_right(r);
                children(x, l, rotated);
            }
            return rotate_left(x);
        }
        return x;
    }
};

template <class Weight, class Number> class Engine {
public:
    static constexpr bool exact = Types<Weight>::exact;
    using Value = Number;
    using Edge = typename Types<Weight>::Edge;
    using Result = typename Types<Weight>::Result;
    using Certificate = typename Types<Weight>::Certificate;

    // Edge IDs are indices in edges. Parallel edges are allowed; self-loops are rejected.
    // Construction/reset validate, greedily complete the seed and initialize prices.
    // Initialization: O(n*n+m) for fixed-width integers, O(n*n+m log(m+1))
    // for floating weights. It is outside solve's budget.
    Engine(int vertex_count, std::vector<Edge> edges, const std::vector<int> &initial_matching = {})
        : n_(checked_count(vertex_count)), edges_(std::move(edges)), adjacency_(n_) {
        if (edges_.size() > std::size_t(std::numeric_limits<int>::max()))
            throw std::length_error("too many edges");
        for (int e = 0; e < int(edges_.size()); ++e) {
            const auto &a = edges_[e];
            if (a.u < 0 || a.u >= n_ || a.v < 0 || a.v >= n_ || a.u == a.v)
                throw std::invalid_argument("invalid edge endpoints");
            if constexpr (!exact) {
                if (!std::isfinite(a.weight)) throw std::invalid_argument("non-finite weight");
                if (std::abs(Value(a.weight)) >
                    std::numeric_limits<Value>::max() / (Value(64) * std::max(1, n_)))
                    throw std::overflow_error("weight too large");
            }
            adjacency_[a.u].push_back(e);
            adjacency_[a.v].push_back(e);
        }
        reset(initial_matching);
    }

    // A fresh optimization with a feasible seed. Previously saved search state is discarded.
    // The entire seed is an incumbent; a dual-feasible tight subset also seeds the search.
    void reset(const std::vector<int> &initial_matching = {}) {
        std::vector<int> seed(n_, -1);
        Value seed_value = 0;
        for (int e : initial_matching) {
            if (e < 0 || e >= int(edges_.size())) throw std::invalid_argument("invalid seed edge");
            const auto &a = edges_[e];
            if (seed[a.u] >= 0 || seed[a.v] >= 0)
                throw std::invalid_argument("seed is not a matching");
            seed[a.u] = seed[a.v] = e;
            seed_value += Value(a.weight);
        }
        Value largest = 0;
        for (const auto &e : edges_) largest = std::max(largest, Value(e.weight));
        tolerance_ = 0;
        if constexpr (!exact)
            tolerance_ = std::max(
                largest * (Value(64) * std::numeric_limits<Value>::epsilon() * std::max(1, n_)),
                Value(8) * std::numeric_limits<Value>::denorm_min());
        stats_ = {};
        stats_.algorithm = Algorithm::Dense;
        optimal_ = failed_ = phase_ = false;
        nodes_.assign(2 * n_, Node{});
        labels_.assign(n_, Unseen);
        blossoms_.clear();
        pair_edges_.assign(std::size_t(n_) * n_, -1);
        free_slots_.clear();
        for (int slot = n_ - 1; slot >= 0; --slot) free_slots_.push_back(slot);
        owner_.resize(n_);
        dual_.assign(n_, largest);
        matching_.assign(n_, -1);
        incumbent_.assign(n_, -1);
        incumbent_value_ = 0;
        if (seed_value > 0) {
            incumbent_ = seed;
            incumbent_value_ = seed_value;
        }
        std::vector<int> greedy = seed;
        std::vector<int> candidates;
        for (int e : initial_matching)
            if (edges_[e].weight < 0) greedy[edges_[e].u] = greedy[edges_[e].v] = -1;
        for (int e = 0; e < int(edges_.size()); ++e)
            if (edges_[e].weight > 0 && greedy[edges_[e].u] < 0 && greedy[edges_[e].v] < 0)
                candidates.push_back(e);
        sort_candidates(candidates, largest);
        std::vector<int> starting_edges;
        for (int e : initial_matching)
            if (edges_[e].weight > 0) starting_edges.push_back(e);
        for (int e : candidates)
            if (greedy[edges_[e].u] < 0 && greedy[edges_[e].v] < 0) {
                greedy[edges_[e].u] = greedy[edges_[e].v] = e;
                starting_edges.push_back(e);
            }
        save_incumbent(greedy);
        available_.clear();
        for (int i = 0; i < n_; ++i) {
            owner_[i] = i;
            nodes_[i].alive = true;
            nodes_[i].base = i;
            nodes_[i].leaves = {i};
            available_.push_back(2 * n_ - 1 - i);
        }
        // Lower prices without violating any edge constraint. Only successful pairs
        // are matched. Retain lowered prices also for unsuccessful pairs; the
        // positive exposed integer prices must retain a common parity.
        for (int e : starting_edges) {
            const auto &a = edges_[e];
            if (a.weight <= 0) continue;
            for (int v : {a.u, a.v}) {
                Value lower = 0;
                for (int f : adjacency_[v])
                    lower = std::max(lower, Value(2) * Value(edges_[f].weight) - dual_[other(f, v)]);
                dual_[v] = lower;
            }
            if (slack(e) <= tolerance_) {
                matching_[a.u] = matching_[a.v] = e;
                ++stats_.initial_search_edges;
                if (seed[a.u] == e) ++stats_.retained_seed_edges;
            } else if constexpr (exact) {
                // Positive exposed prices share parity, so S-S slack is even.
                for (int v : {a.u, a.v})
                    if (dual_[v] > 0 && dual_[v] % 2 != largest % 2) ++dual_[v];
            }
        }
        incoming_.assign(n_, Arc{});
        incoming_key_.assign(n_, 0);
        cache_candidates_ = edges_.size() > std::size_t(n_) * 32;
        reserves_.assign(n_, {});
        frontier_.assign(n_, std::numeric_limits<Value>::max());
        refresh_.assign(n_, 0);
        versions_.assign(n_, 0);
        queue_.clear();
        scan_vertex_ = -1;
        scan_position_ = 0;
        marks_.assign(2 * n_, 0);
        generation_ = 0;
        dual_objective_ = bound_twice();
        if constexpr (exact) optimal_ = dual_objective_ == 2 * incumbent_value_;
    }

    // Repeated calls resume the exact forest, blossoms, dual prices and edge cursor.
    Result solve(const Options &options = {}) {
        if (std::isnan(options.time_limit) || options.time_limit < 0)
            throw std::invalid_argument("time_limit must be nonnegative");
        const auto start = Clock::now();
        std::uint64_t steps = 0;
        while (!optimal_ && !failed_ && steps < options.step_limit) {
            if ((steps & 255) == 0 && std::isfinite(options.time_limit) &&
                std::chrono::duration<double>(Clock::now() - start).count() >= options.time_limit)
                break;
            const auto taken =
                advance(std::min<std::uint64_t>(options.step_limit - steps, 256 - (steps & 255)));
            steps += taken;
            stats_.steps += taken;
        }
        stats_.elapsed += std::chrono::duration<double>(Clock::now() - start).count();
        return result();
    }

    Result result() const {
        Result r;
        r.status = failed_ ? Status::NumericalError : optimal_ ? Status::Optimal : Status::Limit;
        r.mate.assign(n_, -1);
        r.mate_edge = incumbent_;
        r.weight = incumbent_value_;
        for (int v = 0; v < n_; ++v)
            if (incumbent_[v] >= 0) {
                const int e = incumbent_[v];
                r.mate[v] = other(e, v);
                if (v == edges_[e].u) r.edges.push_back(e);
            }
        r.dual_bound_twice = bound_twice();
        r.tolerance = tolerance_;
        r.stats = stats_;
        return r;
    }

    Certificate certificate() const {
        Certificate c;
        c.vertex_dual_twice.assign(dual_.begin(), dual_.end());
        for (int b : blossoms_)
            if (nodes_[b].alive && nodes_[b].price > 0)
                c.odd_sets.push_back({nodes_[b].leaves, nodes_[b].price});
        return c;
    }

private:
    using Clock = std::chrono::steady_clock;
    struct Arc {
        int edge = -1, from = -1, to = -1;
        Arc reversed() const { return {edge, to, from}; }
    };
    enum Label { Unseen, Outer, Inner };
    struct Node {
        bool alive = false;
        int parent = -1, base = -1, tree = -1, slot = -1, best_pair = -1;
        Label label = Unseen;
        Arc entry;
        Value price = 0; // Twice the odd-set dual variable.
        std::vector<int> children, leaves;
        std::vector<Arc> ring; // children[i] -> children[(i+1)%size].
    };
    struct Queued {
        int vertex;
        std::uint64_t version;
    };
    int n_;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> adjacency_;
    std::vector<Node> nodes_;
    std::vector<Label> labels_; // Top-level labels, cached per original vertex.
    std::vector<int> pair_edges_, free_slots_;
    std::vector<int> blossoms_; // Only live odd sets, including nested ones.
    std::vector<int> owner_, matching_, incumbent_, available_;
    std::vector<Queued> queue_;
    std::vector<std::uint64_t> versions_;
    std::vector<Value> dual_;
    std::vector<Arc> incoming_;
    std::vector<Value> incoming_key_;
    struct Candidate {
        Arc arc;
        Value key = 0;
    };
    std::vector<std::array<Candidate, 3>> reserves_;
    std::vector<Value> frontier_;
    std::vector<unsigned char> refresh_;
    Value incumbent_value_ = 0, tolerance_ = 0, phase_time_ = 0, dual_objective_ = 0;
    int root_count_ = 0;
    Statistics stats_;
    bool optimal_ = false, failed_ = false, phase_ = false;
    bool cache_candidates_ = false;
    int scan_vertex_ = -1;
    std::uint64_t scan_version_ = 0;
    std::size_t scan_position_ = 0;
    std::vector<std::uint64_t> marks_;
    std::uint64_t generation_ = 0;

    static int checked_count(int n) {
        if (n < 0 || n > std::numeric_limits<int>::max() / 2)
            throw std::invalid_argument("invalid vertex count");
        return n;
    }
    void sort_candidates(std::vector<int> &ids, Value largest) const {
        if constexpr (exact) {
            if (ids.size() >= 2048) {
                // Stable radix passes preserve ascending edge ID for equal weights.
                // All candidates have positive weights, so the reversed key is unsigned.
                constexpr unsigned bits = 11, buckets = 1u << bits;
                const auto maximum = std::uint64_t(largest);
                std::vector<int> scratch(ids.size());
                std::array<std::size_t, buckets> offsets;
                for (unsigned shift = 0; shift < 63 && (maximum >> shift); shift += bits) {
                    offsets.fill(0);
                    for (int e : ids)
                        ++offsets[((maximum - std::uint64_t(edges_[e].weight)) >> shift) & (buckets - 1)];
                    std::size_t position = 0;
                    for (auto &count : offsets) {
                        const auto next = position + count;
                        count = position;
                        position = next;
                    }
                    for (int e : ids)
                        scratch[offsets[((maximum - std::uint64_t(edges_[e].weight)) >> shift) &
                                        (buckets - 1)]++] = e;
                    ids.swap(scratch);
                }
                return;
            }
        }
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            return edges_[a].weight != edges_[b].weight ? edges_[a].weight > edges_[b].weight
                                                        : a < b;
        });
    }
    int other(int edge, int v) const { return edges_[edge].u ^ edges_[edge].v ^ v; }
    Arc arc(int edge, int from) const { return {edge, from, other(edge, from)}; }
    Value slack(int edge) const {
        const auto &a = edges_[edge];
        return dual_[a.u] + dual_[a.v] - Value(2) * Value(a.weight);
    }
    Value incoming_slack(int v) const {
        if constexpr (exact)
            return incoming_key_[v] + dual_[v] - phase_time_;
        else
            return slack(incoming_[v].edge);
    }
    void clear_incoming(int v) {
        incoming_[v] = {};
        reserves_[v] = {};
        frontier_[v] = std::numeric_limits<Value>::max();
    }
    void remember_incoming(int v, Arc a) {
        if constexpr (!exact) {
            if (incoming_[v].edge < 0 || slack(a.edge) < slack(incoming_[v].edge))
                incoming_[v] = a;
        } else {
            Value key = dual_[a.from] + phase_time_ - Value(2) * Value(edges_[a.edge].weight);
            if (!cache_candidates_) {
                if (incoming_[v].edge < 0 || key < incoming_key_[v]) {
                    incoming_[v] = a;
                    incoming_key_[v] = key;
                }
                return;
            }
            // An unretained old candidate may be smaller than a new candidate beyond
            // this frontier. Keep the frontier when invalid candidates leave holes.
            if (key > frontier_[v] || a.edge == incoming_[v].edge) return;
            for (const auto &b : reserves_[v])
                if (a.edge == b.arc.edge) return;
            if (incoming_[v].edge < 0 || key < incoming_key_[v]) {
                std::swap(a, incoming_[v]);
                std::swap(key, incoming_key_[v]);
            }
            if (a.edge >= 0)
                for (auto &b : reserves_[v]) {
                    if (b.arc.edge < 0 || key < b.key) {
                        std::swap(a, b.arc);
                        std::swap(key, b.key);
                    }
                    if (a.edge < 0) break;
                }
            const auto &last = reserves_[v].back();
            if (last.arc.edge >= 0) frontier_[v] = std::min(frontier_[v], last.key);
        }
    }
    bool repair_incoming(int v) {
        if (refresh_[v]) {
            refresh_[v] = 0;
            clear_incoming(v);
            return true;
        }
        if (incoming_[v].edge < 0) return false;
        if constexpr (exact) {
            if (cache_candidates_) {
                auto &list = reserves_[v];
                std::size_t count = 0;
                for (auto b : list)
                    if (b.arc.edge >= 0 && labels_[b.arc.from] == Outer) list[count++] = b;
                while (count < list.size()) list[count++] = {};
                if (labels_[incoming_[v].from] == Outer) return false;
                if (list[0].arc.edge >= 0) {
                    incoming_[v] = list[0].arc;
                    incoming_key_[v] = list[0].key;
                    for (std::size_t i = 1; i < list.size(); ++i) list[i - 1] = list[i];
                    list.back() = {};
                    return false;
                }
            } else if (labels_[incoming_[v].from] == Outer)
                return false;
        } else if (labels_[incoming_[v].from] == Outer)
            return false;
        clear_incoming(v);
        return true;
    }
    Value bound_twice() const {
        Value bound = 0;
        for (Value y : dual_) bound += y;
        for (int b : blossoms_)
            if (nodes_[b].alive)
                bound += nodes_[b].price * Value((nodes_[b].leaves.size() - 1) / 2);
        return bound;
    }
    void save_incumbent(const std::vector<int> &matching) {
        Value weight = 0;
        for (int v = 0; v < n_; ++v)
            if (matching[v] >= 0 && edges_[matching[v]].u == v)
                weight += Value(edges_[matching[v]].weight);
        if (weight > incumbent_value_) {
            incumbent_value_ = weight;
            incumbent_ = matching;
        }
    }
    void enqueue(int v) { queue_.push_back({v, ++versions_[v]}); }
    int &pair_edge(int a, int b) { return pair_edges_[std::size_t(a) * n_ + b]; }
    void consider_pair(int b, int edge) {
        if (edge >= 0 && (nodes_[b].best_pair < 0 || slack(edge) < slack(nodes_[b].best_pair)))
            nodes_[b].best_pair = edge;
    }
    void allocate_slot(int b) {
        const int slot = free_slots_.back();
        free_slots_.pop_back();
        nodes_[b].slot = slot;
        nodes_[b].best_pair = -1;
        for (int i = 0; i < n_; ++i) pair_edge(slot, i) = pair_edge(i, slot) = -1;
    }
    void release_slot(int b) {
        free_slots_.push_back(nodes_[b].slot);
        nodes_[b].slot = -1;
        nodes_[b].best_pair = -1;
    }
    void remember_pair(int a, int b, int edge) {
        const int x = nodes_[a].slot, y = nodes_[b].slot;
        int &best = pair_edge(x, y);
        if (best < 0 || slack(edge) < slack(best)) best = pair_edge(y, x) = edge;
        consider_pair(a, edge);
        consider_pair(b, edge);
    }
    void merge_pairs(int blossom) {
        Node &node = nodes_[blossom];
        const int slot = nodes_[node.children.front()].slot;
        node.slot = slot;
        node.best_pair = -1;
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0 && b != blossom && nodes_[b].label == Outer) {
                const int other_slot = nodes_[b].slot;
                int best = -1;
                for (int child : node.children)
                    if (nodes_[child].label == Outer) {
                        const int e = pair_edge(nodes_[child].slot, other_slot);
                        if (e >= 0 && (best < 0 || slack(e) < slack(best))) best = e;
                    }
                pair_edge(slot, other_slot) = pair_edge(other_slot, slot) = best;
                consider_pair(blossom, best);
                consider_pair(b, best);
            }
        for (int child : node.children)
            if (nodes_[child].slot >= 0) {
                if (nodes_[child].slot != slot) free_slots_.push_back(nodes_[child].slot);
                nodes_[child].slot = -1;
                nodes_[child].best_pair = -1;
            }
    }
    void repair_pairs() {
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0 && nodes_[b].label == Outer) {
                const int e = nodes_[b].best_pair;
                if (e < 0 || (labels_[edges_[e].u] == Outer && labels_[edges_[e].v] == Outer))
                    continue;
                nodes_[b].best_pair = -1;
                for (int c = 0; c < 2 * n_; ++c)
                    if (c != b && nodes_[c].alive && nodes_[c].parent < 0 && nodes_[c].label == Outer)
                        consider_pair(b, pair_edge(nodes_[b].slot, nodes_[c].slot));
            }
    }
    void make_outer(int b, Arc entry, int tree) {
        allocate_slot(b);
        nodes_[b].label = Outer;
        nodes_[b].entry = entry;
        nodes_[b].tree = tree;
        for (int v : nodes_[b].leaves) {
            labels_[v] = Outer;
            enqueue(v);
        }
    }
    void begin_phase() {
        for (auto &node : nodes_) {
            node.label = Unseen;
            node.entry = {};
        }
        for (int v = 0; v < n_; ++v) clear_incoming(v);
        std::fill(labels_.begin(), labels_.end(), Unseen);
        queue_.clear();
        scan_vertex_ = -1;
        phase_time_ = 0;
        root_count_ = 0;
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0 && matching_[nodes_[b].base] < 0 &&
                dual_[nodes_[b].base] > tolerance_) {
                make_outer(b, {}, nodes_[b].base);
                ++root_count_;
            }
        if (queue_.empty()) {
            optimal_ = true;
            save_incumbent(matching_);
            return;
        }
        phase_ = true;
        ++stats_.phases;
    }
    void finish_phase(int first_tree, int second_tree = -1) {
        save_incumbent(matching_);
        root_count_ -= second_tree >= 0 ? 2 : 1;
        // Remove only trees touched by this augmentation/exposure move. Prices and
        // all other alternating trees remain valid; re-scan only invalidated minima.
        bool roots_remain = false;
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0) {
                Node &node = nodes_[b];
                if (node.label != Unseen && (node.tree == first_tree || node.tree == second_tree)) {
                    if (node.label == Outer) release_slot(b);
                    if (node.label == Outer)
                        for (int v : node.leaves) {
                            ++versions_[v];
                            incoming_[v] = {};
                            refresh_[v] = 1;
                        }
                    for (int v : node.leaves) labels_[v] = Unseen;
                    node.label = Unseen;
                    node.entry = {};
                    node.tree = -1;
                }
                if (node.label == Outer && matching_[node.base] < 0 && dual_[node.base] > tolerance_)
                    roots_remain = true;
            }
        if (!roots_remain) {
            optimal_ = true;
            return;
        }
        for (int v = 0; v < n_; ++v)
            if (labels_[v] != Outer) {
                if (!repair_incoming(v)) continue;
                for (int e : adjacency_[v]) {
                    ++stats_.edge_scans;
                    const int u = other(e, v);
                    if (labels_[u] == Outer) remember_incoming(v, arc(e, u));
                }
            }
        repair_pairs();
        ++stats_.phases;
    }
    int predecessor(int b) const {
        return nodes_[b].entry.edge < 0 ? -1 : owner_[nodes_[b].entry.from];
    }
    int common_ancestor(int a, int b) {
        if (++generation_ == 0) {
            std::fill(marks_.begin(), marks_.end(), 0);
            ++generation_;
        }
        for (int v = a; v >= 0; v = predecessor(v)) marks_[v] = generation_;
        for (int v = b; v >= 0; v = predecessor(v))
            if (marks_[v] == generation_) return v;
        return -1;
    }
    void contract(Arc bridge, int common) {
        std::vector<int> left, right;
        std::vector<Arc> left_links, right_links;
        for (int b = owner_[bridge.from]; b != common; b = predecessor(b)) {
            left.push_back(b);
            left_links.push_back(nodes_[b].entry);
        }
        for (int b = owner_[bridge.to]; b != common; b = predecessor(b)) {
            right.push_back(b);
            right_links.push_back(nodes_[b].entry.reversed());
        }
        const int id = available_.back();
        available_.pop_back();
        blossoms_.push_back(id);
        Node &node = nodes_[id];
        node = Node{};
        node.alive = true;
        node.label = Outer;
        node.base = nodes_[common].base;
        node.entry = nodes_[common].entry;
        node.tree = nodes_[common].tree;
        node.children.push_back(common);
        node.children.insert(node.children.end(), left.rbegin(), left.rend());
        node.children.insert(node.children.end(), right.begin(), right.end());
        node.ring.insert(node.ring.end(), left_links.rbegin(), left_links.rend());
        node.ring.push_back(bridge);
        node.ring.insert(node.ring.end(), right_links.begin(), right_links.end());
        for (int child : node.children) {
            nodes_[child].parent = id;
            for (int v : nodes_[child].leaves) {
                if (nodes_[child].label != Outer) enqueue(v);
                owner_[v] = id;
                labels_[v] = Outer;
                node.leaves.push_back(v);
            }
        }
        merge_pairs(id);
        ++stats_.contractions;
    }
    int child_position(int blossom, int vertex) const {
        int child = vertex;
        while (nodes_[child].parent != blossom) child = nodes_[child].parent;
        const auto &cycle = nodes_[blossom].children;
        return int(std::find(cycle.begin(), cycle.end(), child) - cycle.begin());
    }
    Arc ring_step(int b, int position, int direction) const {
        const auto &ring = nodes_[b].ring;
        return direction > 0 ? ring[position]
                             : ring[(position + int(ring.size()) - 1) % ring.size()].reversed();
    }
    // Rotate the canonical near-perfect matching inside a blossom so x is its base.
    void expose(int b, int x) {
        if (b < n_ || nodes_[b].base == x) return;
        Node &node = nodes_[b];
        const int start = child_position(b, x), k = int(node.children.size());
        const int direction = (start % 2 == 0) ? -1 : 1;
        expose(node.children[start], x);
        int position = start;
        while (position != 0) {
            const int middle = (position + direction + k) % k;
            const int next = (middle + direction + k) % k;
            const Arc link = ring_step(b, middle, direction);
            expose(node.children[middle], link.from);
            expose(node.children[next], link.to);
            matching_[link.from] = matching_[link.to] = link.edge;
            position = next;
        }
        std::rotate(node.children.begin(), node.children.begin() + start, node.children.end());
        std::rotate(node.ring.begin(), node.ring.begin() + start, node.ring.end());
        node.base = x;
    }
    // Flip an arm from x to its exposed root. edge==-1 relocates exposure to x.
    void flip_arm(int x, int edge) {
        for (;;) {
            const int b = owner_[x];
            const Arc back = nodes_[b].entry;
            expose(b, x);
            matching_[x] = edge;
            if (back.edge < 0) break;
            const int inner = owner_[back.from];
            const Arc next = nodes_[inner].entry;
            expose(inner, next.to);
            matching_[next.to] = next.edge;
            x = next.from;
            edge = next.edge;
        }
    }
    void use_edge(Arc a) {
        const int from = owner_[a.from], to = owner_[a.to];
        if (from == to) return;
        if (nodes_[to].label == Unseen) {
            const int base = nodes_[to].base;
            if (matching_[base] < 0) {
                flip_arm(a.from, a.edge);
                flip_arm(a.to, a.edge);
                ++stats_.augmentations;
                finish_phase(nodes_[from].tree);
            } else {
                nodes_[to].label = Inner;
                for (int v : nodes_[to].leaves) labels_[v] = Inner;
                nodes_[to].entry = a;
                nodes_[to].tree = nodes_[from].tree;
                const Arc matched = arc(matching_[base], base);
                make_outer(owner_[matched.to], matched, nodes_[from].tree);
            }
        } else if (nodes_[to].label == Outer) {
            const int common = common_ancestor(from, to);
            if (common >= 0)
                contract(a, common);
            else {
                flip_arm(a.from, a.edge);
                flip_arm(a.to, a.edge);
                ++stats_.augmentations;
                finish_phase(nodes_[from].tree, nodes_[to].tree);
            }
        }
    }
    void expand(int b) {
        const Arc entry = nodes_[b].entry;
        const int tree = nodes_[b].tree;
        const int start = child_position(b, entry.to);
        const int k = int(nodes_[b].children.size());
        const int direction = start % 2 == 0 ? -1 : 1;
        for (int child : nodes_[b].children) {
            nodes_[child].parent = -1;
            nodes_[child].label = Unseen;
            nodes_[child].entry = {};
            nodes_[child].tree = -1;
            for (int v : nodes_[child].leaves) {
                owner_[v] = child;
                labels_[v] = Unseen;
            }
        }
        int position = start;
        Arc incoming = entry;
        for (;;) {
            const int child = nodes_[b].children[position];
            nodes_[child].label = Inner;
            for (int v : nodes_[child].leaves) labels_[v] = Inner;
            nodes_[child].entry = incoming;
            nodes_[child].tree = tree;
            if (position == 0) break;
            const Arc matched = ring_step(b, position, direction);
            position = (position + direction + k) % k;
            make_outer(nodes_[b].children[position], matched, tree);
            incoming = ring_step(b, position, direction);
            position = (position + direction + k) % k;
        }
        nodes_[b] = Node{};
        available_.push_back(b);
        auto it = std::find(blossoms_.begin(), blossoms_.end(), b);
        *it = blossoms_.back();
        blossoms_.pop_back();
        ++stats_.expansions;
    }
    void scan_edge(int v, int e) {
        ++stats_.edge_scans;
        const int w = other(e, v), bv = owner_[v], bw = owner_[w];
        if (bv == bw) return;
        const Value gap = slack(e);
        if (gap < -tolerance_) {
            failed_ = true;
            return;
        }
        if (labels_[w] == Outer) {
            if (gap <= tolerance_)
                use_edge(arc(e, v));
            else
                remember_pair(bv, bw, e);
        } else {
            remember_incoming(w, arc(e, v));
            if (labels_[w] == Unseen && gap <= tolerance_) use_edge(arc(e, v));
        }
    }
    void dual_event() {
        // One minimum per outer top-level blossom. Contraction merges its rows
        // in O(n * number_of_children), eliminating per-edge heap operations.
        enum Kind { Expose, Grow, Join, Expand } kind = Expose;
        int index = -1;
        Value delta = 0;
        for (int v = 0; v < n_; ++v)
            if (labels_[v] == Outer && (index < 0 || dual_[v] < delta)) {
                delta = dual_[v];
                index = v;
            }
        if (index < 0) {
            failed_ = true;
            return;
        }
        for (int v = 0; v < n_; ++v)
            if (labels_[v] == Unseen && incoming_[v].edge >= 0) {
                const Value candidate = incoming_slack(v);
                if (candidate < delta) {
                    delta = candidate;
                    kind = Grow;
                    index = v;
                }
            }
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0 && nodes_[b].label == Outer &&
                nodes_[b].best_pair >= 0) {
                const int e = nodes_[b].best_pair;
                const Value candidate = slack(e) / 2;
                if (candidate < delta) {
                    delta = candidate;
                    kind = Join;
                    index = e;
                }
            }
        for (int b : blossoms_)
            if (nodes_[b].alive && nodes_[b].parent < 0 && nodes_[b].label == Inner &&
                nodes_[b].price / 2 < delta) {
                delta = nodes_[b].price / 2;
                kind = Expand;
                index = b;
            }
        if (delta < -tolerance_) {
            failed_ = true;
            return;
        }
        delta = std::max(Value(0), delta);
        for (int v = 0; v < n_; ++v) {
            if (labels_[v] == Outer)
                dual_[v] -= delta;
            else if (labels_[v] == Inner)
                dual_[v] += delta;
        }
        for (int b : blossoms_)
            if (nodes_[b].alive && nodes_[b].parent < 0) {
                if (nodes_[b].label == Outer)
                    nodes_[b].price += 2 * delta;
                else if (nodes_[b].label == Inner)
                    nodes_[b].price -= 2 * delta;
            }
        phase_time_ += delta;
        // Each alternating tree has one more outer than inner top-level node.
        // Odd-set coefficients cancel the additional leaf terms exactly.
        if constexpr (exact) dual_objective_ -= delta * root_count_;
        ++stats_.dual_updates;
        // A supplied incumbent may already attain the dual bound even while the
        // repaired internal matching is incomplete. Equality is an exact proof.
        if constexpr (exact)
            if (dual_objective_ == 2 * incumbent_value_) {
                optimal_ = true;
                return;
            }
        if (kind == Expose) {
            const int tree = nodes_[owner_[index]].tree;
            flip_arm(index, -1);
            ++stats_.relocations;
            finish_phase(tree);
        } else if (kind == Grow)
            use_edge(incoming_[index]);
        else if (kind == Join)
            use_edge(arc(index, edges_[index].u));
        else
            expand(index);
    }
    std::uint64_t advance(std::uint64_t limit) {
        if (!phase_) {
            begin_phase();
            return 1;
        }
        if (scan_vertex_ >= 0 &&
            (versions_[scan_vertex_] != scan_version_ || labels_[scan_vertex_] != Outer))
            scan_vertex_ = -1;
        if (scan_vertex_ < 0 || scan_position_ == adjacency_[scan_vertex_].size()) {
            scan_vertex_ = -1;
            while (!queue_.empty()) {
                const auto next = queue_.back();
                queue_.pop_back();
                if (versions_[next.vertex] != next.version || labels_[next.vertex] != Outer)
                    continue;
                scan_vertex_ = next.vertex;
                scan_version_ = next.version;
                scan_position_ = 0;
                break;
            }
        }
        if (scan_vertex_ < 0) {
            dual_event();
            return 1;
        }
        const int v = scan_vertex_;
        const auto &adjacent = adjacency_[v];
        const std::size_t first = scan_position_,
                          stop = std::min(adjacent.size(), first + std::size_t(limit));
        std::size_t position = first;
        while (position < stop) {
            scan_edge(v, adjacent[position++]);
            if (optimal_ || failed_ || versions_[v] != scan_version_ || labels_[v] != Outer) break;
        }
        scan_position_ = position;
        return std::max<std::uint64_t>(1, position - first);
    }
};

template <class Weight, class Number> class PriorityEngine {
public:
    static constexpr bool exact = Types<Weight>::exact;
    using Value = Number;
    using Edge = typename Types<Weight>::Edge;
    using Result = typename Types<Weight>::Result;
    using Certificate = typename Types<Weight>::Certificate;

    // Edge IDs are indices in edges. Parallel edges are allowed; self-loops are rejected.
    // Construction/reset validate, greedily complete the seed and initialize prices.
    // Initialization: O(n+m) for fixed-width integers, O(n+m log(m+1))
    // for floating weights. It is outside solve's budget.
    PriorityEngine(int vertex_count, std::vector<Edge> edges,
                   const std::vector<int> &initial_matching = {})
        : n_(checked_count(vertex_count)), edges_(std::move(edges)), adjacency_(n_) {
        if (edges_.size() > std::size_t(std::numeric_limits<int>::max()))
            throw std::length_error("too many edges");
        for (int e = 0; e < int(edges_.size()); ++e) {
            const auto &a = edges_[e];
            if (a.u < 0 || a.u >= n_ || a.v < 0 || a.v >= n_ || a.u == a.v)
                throw std::invalid_argument("invalid edge endpoints");
            if constexpr (!exact) {
                if (!std::isfinite(a.weight)) throw std::invalid_argument("non-finite weight");
                if (std::abs(Value(a.weight)) >
                    std::numeric_limits<Value>::max() / (Value(64) * std::max(1, n_)))
                    throw std::overflow_error("weight too large");
            }
            adjacency_[a.u].push_back(e);
            adjacency_[a.v].push_back(e);
        }
        reset(initial_matching);
    }

    // A fresh optimization with a feasible seed. Previously saved search state is discarded.
    // The entire seed is an incumbent; a dual-feasible tight subset also seeds the search.
    void reset(const std::vector<int> &initial_matching = {}) {
        std::vector<int> seed(n_, -1);
        Value seed_value = 0;
        for (int e : initial_matching) {
            if (e < 0 || e >= int(edges_.size())) throw std::invalid_argument("invalid seed edge");
            const auto &a = edges_[e];
            if (seed[a.u] >= 0 || seed[a.v] >= 0)
                throw std::invalid_argument("seed is not a matching");
            seed[a.u] = seed[a.v] = e;
            seed_value += Value(a.weight);
        }
        Value largest = 0;
        for (const auto &e : edges_) largest = std::max(largest, Value(e.weight));
        tolerance_ = 0;
        if constexpr (!exact)
            tolerance_ = std::max(
                largest * (Value(64) * std::numeric_limits<Value>::epsilon() * std::max(1, n_)),
                Value(8) * std::numeric_limits<Value>::denorm_min());
        stats_ = {};
        stats_.algorithm = Algorithm::Sparse;
        optimal_ = failed_ = phase_ = false;
        nodes_.assign(2 * n_, Node{});
        blossoms_.clear();
        dual_.assign(n_, largest);
        matching_.assign(n_, -1);
        incumbent_.assign(n_, -1);
        incumbent_value_ = 0;
        if (seed_value > 0) {
            incumbent_ = seed;
            incumbent_value_ = seed_value;
        }
        std::vector<int> greedy = seed;
        std::vector<int> candidates;
        for (int e : initial_matching)
            if (edges_[e].weight < 0) greedy[edges_[e].u] = greedy[edges_[e].v] = -1;
        for (int e = 0; e < int(edges_.size()); ++e)
            if (edges_[e].weight > 0 && greedy[edges_[e].u] < 0 && greedy[edges_[e].v] < 0)
                candidates.push_back(e);
        sort_candidates(candidates, largest);
        std::vector<int> starting_edges;
        for (int e : initial_matching)
            if (edges_[e].weight > 0) starting_edges.push_back(e);
        for (int e : candidates)
            if (greedy[edges_[e].u] < 0 && greedy[edges_[e].v] < 0) {
                greedy[edges_[e].u] = greedy[edges_[e].v] = e;
                starting_edges.push_back(e);
            }
        save_incumbent(greedy);
        available_.clear();
        for (int i = 0; i < n_; ++i) {
            nodes_[i].alive = true;
            nodes_[i].base = i;
            nodes_[i].sequence = nodes_[i].first = i;
            available_.push_back(2 * n_ - 1 - i);
        }
        // Lower prices without violating any edge constraint. Only successful pairs
        // are matched. Retain lowered prices also for unsuccessful pairs; the
        // positive exposed integer prices must retain a common parity.
        for (int e : starting_edges) {
            const auto &a = edges_[e];
            if (a.weight <= 0) continue;
            for (int v : {a.u, a.v}) {
                Value lower = 0;
                for (int f : adjacency_[v])
                    lower = std::max(lower, Value(2) * Value(edges_[f].weight) - dual_[other(f, v)]);
                dual_[v] = lower;
            }
            if (dual_[a.u] + dual_[a.v] - Value(2) * Value(a.weight) <= tolerance_) {
                matching_[a.u] = matching_[a.v] = e;
                ++stats_.initial_search_edges;
                if (seed[a.u] == e) ++stats_.retained_seed_edges;
            } else if constexpr (exact) {
                // Positive exposed prices share parity, so S-S slack is even.
                for (int v : {a.u, a.v})
                    if (dual_[v] > 0 && dual_[v] % 2 != largest % 2) ++dual_[v];
            }
        }
        incoming_.assign(n_, Arc{});
        incoming_key_.assign(n_, 0);
        refresh_.assign(n_, 0);
        versions_.assign(n_, 0);
        queue_.clear();
        events_ = EventQueue{};
        scan_vertex_ = -1;
        scan_position_ = 0;
        marks_.assign(2 * n_, 0);
        generation_ = 0;
        phase_time_ = 0;
        sequence_.reset(dual_);
        node_events_.reset(2 * n_);
        dual_objective_ = bound_twice();
        if constexpr (exact) optimal_ = dual_objective_ == 2 * incumbent_value_;
    }

    // Repeated calls resume the exact forest, blossoms, dual prices and edge cursor.
    Result solve(const Options &options = {}) {
        if (std::isnan(options.time_limit) || options.time_limit < 0)
            throw std::invalid_argument("time_limit must be nonnegative");
        const auto start = Clock::now();
        std::uint64_t steps = 0;
        while (!optimal_ && !failed_ && steps < options.step_limit) {
            if ((steps & 255) == 0 && std::isfinite(options.time_limit) &&
                std::chrono::duration<double>(Clock::now() - start).count() >= options.time_limit)
                break;
            const auto taken =
                advance(std::min<std::uint64_t>(options.step_limit - steps, 256 - (steps & 255)));
            steps += taken;
            stats_.steps += taken;
        }
        stats_.elapsed += std::chrono::duration<double>(Clock::now() - start).count();
        return result();
    }

    Result result() const {
        Result r;
        r.status = failed_ ? Status::NumericalError : optimal_ ? Status::Optimal : Status::Limit;
        r.mate.assign(n_, -1);
        r.mate_edge = incumbent_;
        r.weight = incumbent_value_;
        for (int v = 0; v < n_; ++v)
            if (incumbent_[v] >= 0) {
                const int e = incumbent_[v];
                r.mate[v] = other(e, v);
                if (v == edges_[e].u) r.edges.push_back(e);
            }
        r.dual_bound_twice = bound_twice();
        r.tolerance = tolerance_;
        r.stats = stats_;
        return r;
    }

    Certificate certificate() const {
        Certificate c;
        c.vertex_dual_twice.reserve(n_);
        for (int v = 0; v < n_; ++v) c.vertex_dual_twice.push_back(vertex_dual(v));
        for (int b : blossoms_)
            if (blossom_dual(b) > 0) {
                typename Types<Weight>::OddSet set;
                set.dual_twice = blossom_dual(b);
                std::vector<int> stack{b};
                while (!stack.empty()) {
                    const int x = stack.back();
                    stack.pop_back();
                    if (x < n_)
                        set.vertices.push_back(x);
                    else
                        stack.insert(
                            stack.end(), nodes_[x].children.begin(), nodes_[x].children.end());
                }
                c.odd_sets.push_back(std::move(set));
            }
        return c;
    }

private:
    using Clock = std::chrono::steady_clock;
    struct Arc {
        int edge = -1, from = -1, to = -1;
        Arc reversed() const { return {edge, to, from}; }
    };
    enum Label { Unseen, Outer, Inner };
    struct Node {
        bool alive = false;
        int parent = -1, base = -1, tree = -1;
        Label label = Unseen;
        Arc entry;
        Value price = 0; // Twice the odd-set dual variable.
        int sequence = -1, first = -1, size = 1, list_position = -1;
        std::vector<int> children;
        std::vector<Arc> ring; // children[i] -> children[(i+1)%size].
    };
    struct Event {
        Value time;
        int edge;
        std::uint64_t u_version, v_version;
    };
    struct Later {
        bool operator()(const Event &a, const Event &b) const {
            return a.time != b.time ? a.time > b.time : a.edge > b.edge;
        }
    };
    struct EventQueue {
        std::vector<Event> entries;
        bool empty() const { return entries.empty(); }
        const Event &top() const { return entries.front(); }
        void push(Event e) {
            entries.push_back(e);
            std::push_heap(entries.begin(), entries.end(), Later{});
        }
        void pop() {
            std::pop_heap(entries.begin(), entries.end(), Later{});
            entries.pop_back();
        }
    };
    struct NodeEvents {
        std::vector<int> heap, position;
        std::vector<Value> keys;
        void reset(int n) {
            heap.clear();
            position.assign(n, -1);
            keys.resize(n);
        }
        bool less(int a, int b) const { return keys[a] != keys[b] ? keys[a] < keys[b] : a < b; }
        void swap_positions(int a, int b) {
            std::swap(heap[a], heap[b]);
            position[heap[a]] = a;
            position[heap[b]] = b;
        }
        int up(int at) {
            while (at && less(heap[at], heap[(at - 1) / 2])) {
                swap_positions(at, (at - 1) / 2);
                at = (at - 1) / 2;
            }
            return at;
        }
        void down(int at) {
            for (;;) {
                int next = 2 * at + 1;
                if (next >= int(heap.size())) return;
                if (next + 1 < int(heap.size()) && less(heap[next + 1], heap[next])) ++next;
                if (!less(heap[next], heap[at])) return;
                swap_positions(at, next);
                at = next;
            }
        }
        void set(int id, Value value) {
            keys[id] = value;
            if (position[id] < 0) {
                position[id] = int(heap.size());
                heap.push_back(id);
            }
            down(up(position[id]));
        }
        void erase(int id) {
            const int at = position[id];
            if (at < 0) return;
            swap_positions(at, int(heap.size()) - 1);
            heap.pop_back();
            position[id] = -1;
            if (at < int(heap.size())) down(up(at));
        }
    } node_events_;
    struct Queued {
        int vertex;
        std::uint64_t version;
    };
    int n_;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> adjacency_;
    std::vector<Node> nodes_;
    std::vector<int> blossoms_; // Only live odd sets, including nested ones.
    std::vector<int> matching_, incumbent_, available_;
    std::vector<Queued> queue_;
    std::vector<std::uint64_t> versions_;
    std::vector<Value> dual_;
    std::vector<Arc> incoming_;
    std::vector<Value> incoming_key_;
    VertexSequence<Value> sequence_;
    std::vector<unsigned char> refresh_;
    EventQueue events_;
    Value incumbent_value_ = 0, tolerance_ = 0, phase_time_ = 0, dual_objective_ = 0;
    int root_count_ = 0;
    Statistics stats_;
    bool optimal_ = false, failed_ = false, phase_ = false;
    int scan_vertex_ = -1;
    std::uint64_t scan_version_ = 0;
    std::size_t scan_position_ = 0;
    std::vector<std::uint64_t> marks_;
    std::uint64_t generation_ = 0;

    static int checked_count(int n) {
        if (n < 0 || n > std::numeric_limits<int>::max() / 2)
            throw std::invalid_argument("invalid vertex count");
        return n;
    }
    void sort_candidates(std::vector<int> &ids, Value largest) const {
        if constexpr (exact) {
            if (ids.size() >= 2048) {
                // Stable radix passes preserve ascending edge ID for equal weights.
                // All candidates have positive weights, so the reversed key is unsigned.
                constexpr unsigned bits = 11, buckets = 1u << bits;
                const auto maximum = std::uint64_t(largest);
                std::vector<int> scratch(ids.size());
                std::array<std::size_t, buckets> offsets;
                for (unsigned shift = 0; shift < 63 && (maximum >> shift); shift += bits) {
                    offsets.fill(0);
                    for (int e : ids)
                        ++offsets[((maximum - std::uint64_t(edges_[e].weight)) >> shift) & (buckets - 1)];
                    std::size_t position = 0;
                    for (auto &count : offsets) {
                        const auto next = position + count;
                        count = position;
                        position = next;
                    }
                    for (int e : ids)
                        scratch[offsets[((maximum - std::uint64_t(edges_[e].weight)) >> shift) &
                                        (buckets - 1)]++] = e;
                    ids.swap(scratch);
                }
                return;
            }
        }
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            return edges_[a].weight != edges_[b].weight ? edges_[a].weight > edges_[b].weight
                                                        : a < b;
        });
    }
    int other(int edge, int v) const { return edges_[edge].u ^ edges_[edge].v ^ v; }
    Arc arc(int edge, int from) const { return {edge, from, other(edge, from)}; }
    int owner(int v) const { return sequence_.data[sequence_.root(v)].owner; }
    static int rate(Label label) { return label == Outer ? -1 : label == Inner ? 1 : 0; }
    Value vertex_dual(int v) const {
        return sequence_.vertex_dual(v) + Value(rate(nodes_[owner(v)].label)) * phase_time_;
    }
    Value blossom_dual(int b) const {
        const Node &node = nodes_[b];
        return node.price - (node.parent < 0 ? Value(2 * rate(node.label)) * phase_time_ : Value(0));
    }
    Value slack(int edge) const {
        const auto &a = edges_[edge];
        return vertex_dual(a.u) + vertex_dual(a.v) - Value(2) * Value(a.weight);
    }
    void update_event(int b) {
        const Node &node = nodes_[b];
        if (!node.alive || node.parent >= 0) {
            node_events_.erase(b);
            return;
        }
        const auto &root = sequence_.data[node.sequence];
        if (node.label == Outer)
            node_events_.set(b, root.dual);
        else if (node.label == Inner && b >= n_)
            node_events_.set(b, node.price / 2);
        else if (node.label == Unseen && root.incoming_arg >= 0)
            node_events_.set(b, root.incoming);
        else
            node_events_.erase(b);
    }
    void change_label(int b, Label next) {
        Node &node = nodes_[b];
        const Value price = blossom_dual(b);
        sequence_.add_to(node.sequence, Value(rate(node.label) - rate(next)) * phase_time_);
        node.label = next;
        node.price = price + Value(2 * rate(next)) * phase_time_;
        update_event(b);
    }
    Value incoming_slack(int v) const { return incoming_key_[v] + vertex_dual(v) - phase_time_; }
    void clear_incoming(int v) {
        incoming_[v] = {};
        const int root = sequence_.set_incoming(v, false);
        update_event(sequence_.data[root].owner);
    }
    void remember_incoming(int v, Arc a) {
        const Value key =
            vertex_dual(a.from) + phase_time_ - Value(2) * Value(edges_[a.edge].weight);
        if (incoming_[v].edge >= 0 && key >= incoming_key_[v]) return;
        incoming_[v] = a;
        incoming_key_[v] = key;
        const int root = sequence_.set_incoming(v, true, key);
        update_event(sequence_.data[root].owner);
    }
    bool repair_incoming(int v) {
        if (refresh_[v]) {
            refresh_[v] = 0;
            clear_incoming(v);
            return true;
        }
        if (incoming_[v].edge < 0 || nodes_[owner(incoming_[v].from)].label == Outer) return false;
        clear_incoming(v);
        return true;
    }
    Value bound_twice() const {
        Value bound = 0;
        for (int v = 0; v < n_; ++v) bound += vertex_dual(v);
        for (int b : blossoms_) bound += blossom_dual(b) * Value((nodes_[b].size - 1) / 2);
        return bound;
    }
    void save_incumbent(const std::vector<int> &matching) {
        Value weight = 0;
        for (int v = 0; v < n_; ++v)
            if (matching[v] >= 0 && edges_[matching[v]].u == v)
                weight += Value(edges_[matching[v]].weight);
        if (weight > incumbent_value_) {
            incumbent_value_ = weight;
            incumbent_ = matching;
        }
    }
    void enqueue(int v) { queue_.push_back({v, ++versions_[v]}); }
    void make_outer(int b, Arc entry, int tree) {
        change_label(b, Outer);
        nodes_[b].entry = entry;
        nodes_[b].tree = tree;
        sequence_.for_each(nodes_[b].sequence, [&](int v) { enqueue(v); });
    }
    void begin_phase() {
        // Only reset() starts a new forest. Unaffected trees survive augmentations.
        queue_.clear();
        events_ = EventQueue{};
        scan_vertex_ = -1;
        root_count_ = 0;
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0 && matching_[nodes_[b].base] < 0 &&
                vertex_dual(nodes_[b].base) > tolerance_) {
                make_outer(b, {}, nodes_[b].base);
                ++root_count_;
            }
        if (queue_.empty()) {
            optimal_ = true;
            save_incumbent(matching_);
            return;
        }
        phase_ = true;
        ++stats_.phases;
    }
    bool valid_event(const Event &e) const {
        const auto &a = edges_[e.edge];
        if (versions_[a.u] != e.u_version || versions_[a.v] != e.v_version) return false;
        const int u = owner(a.u), v = owner(a.v);
        return u != v && nodes_[u].label == Outer && nodes_[v].label == Outer;
    }
    void finish_phase(int first_tree, int second_tree = -1) {
        save_incumbent(matching_);
        root_count_ -= second_tree >= 0 ? 2 : 1;
        bool roots_remain = false;
        for (int b = 0; b < 2 * n_; ++b)
            if (nodes_[b].alive && nodes_[b].parent < 0) {
                Node &node = nodes_[b];
                if (node.label != Unseen && (node.tree == first_tree || node.tree == second_tree)) {
                    if (node.label == Outer)
                        sequence_.for_each(node.sequence, [&](int v) {
                            ++versions_[v];
                            refresh_[v] = 1;
                        });
                    change_label(b, Unseen);
                    node.entry = {};
                    node.tree = -1;
                }
                if (node.label == Outer && matching_[node.base] < 0 &&
                    vertex_dual(node.base) > tolerance_)
                    roots_remain = true;
            }
        if (!roots_remain) {
            optimal_ = true;
            return;
        }
        // Removing sources cannot create a new minimum for an empty cache. Repair
        // all affected vertices before any further edge scan can re-label a source.
        for (int v = 0; v < n_; ++v)
            if (nodes_[owner(v)].label != Outer && repair_incoming(v)) {
                Arc best;
                Value key = 0;
                for (int e : adjacency_[v]) {
                    ++stats_.edge_scans;
                    const int u = other(e, v);
                    if (nodes_[owner(u)].label != Outer) continue;
                    const Value candidate =
                        vertex_dual(u) + phase_time_ - Value(2) * Value(edges_[e].weight);
                    if (best.edge < 0 || candidate < key) {
                        best = arc(e, u);
                        key = candidate;
                    }
                }
                if (best.edge >= 0) remember_incoming(v, best);
            }
        if (queue_.size() > std::size_t(2) * n_) {
            queue_.erase(std::remove_if(queue_.begin(), queue_.end(),
                                        [&](const Queued &item) {
                                            return versions_[item.vertex] != item.version ||
                                                   nodes_[owner(item.vertex)].label != Outer;
                                        }),
                         queue_.end());
        }
        if (events_.entries.size() > 4 * edges_.size()) {
            auto &entries = events_.entries;
            entries.erase(std::remove_if(entries.begin(), entries.end(),
                                         [&](const Event &e) { return !valid_event(e); }),
                          entries.end());
            std::make_heap(entries.begin(), entries.end(), Later{});
        }
        ++stats_.phases;
    }
    int predecessor(int b) const {
        return nodes_[b].entry.edge < 0 ? -1 : owner(nodes_[b].entry.from);
    }
    int common_ancestor(int a, int b) {
        if (++generation_ == 0) {
            std::fill(marks_.begin(), marks_.end(), 0);
            ++generation_;
        }
        // Alternating the two upward walks charges a contraction to its cycle.
        while (a >= 0 || b >= 0) {
            if (a >= 0) {
                if (marks_[a] == generation_) return a;
                marks_[a] = generation_;
                a = predecessor(a);
            }
            std::swap(a, b);
        }
        return -1;
    }
    void contract(Arc bridge, int common) {
        std::vector<int> left, right;
        std::vector<Arc> left_links, right_links;
        for (int b = owner(bridge.from); b != common; b = predecessor(b)) {
            left.push_back(b);
            left_links.push_back(nodes_[b].entry);
        }
        for (int b = owner(bridge.to); b != common; b = predecessor(b)) {
            right.push_back(b);
            right_links.push_back(nodes_[b].entry.reversed());
        }
        const int id = available_.back();
        available_.pop_back();
        blossoms_.push_back(id);
        Node &node = nodes_[id];
        node = Node{};
        node.list_position = int(blossoms_.size()) - 1;
        node.alive = true;
        node.label = Outer;
        node.base = nodes_[common].base;
        node.entry = nodes_[common].entry;
        node.tree = nodes_[common].tree;
        node.price = -Value(2) * phase_time_;
        node.size = 0;
        node.first = nodes_[common].first;
        node.children.push_back(common);
        node.children.insert(node.children.end(), left.rbegin(), left.rend());
        node.children.insert(node.children.end(), right.begin(), right.end());
        node.ring.insert(node.ring.end(), left_links.rbegin(), left_links.rend());
        node.ring.push_back(bridge);
        node.ring.insert(node.ring.end(), right_links.begin(), right_links.end());
        for (int child : node.children) {
            Node &sub = nodes_[child];
            node_events_.erase(child);
            sub.price = blossom_dual(child);
            sequence_.add_to(sub.sequence, Value(rate(sub.label) + 1) * phase_time_);
            if (sub.label != Outer) sequence_.for_each(sub.sequence, [&](int v) { enqueue(v); });
            sub.parent = id;
            node.size += sub.size;
            node.sequence = sequence_.concatenate(node.sequence, sub.sequence);
        }
        sequence_.data[node.sequence].owner = id;
        update_event(id);
        ++stats_.contractions;
    }
    int child_position(int blossom, int vertex) const {
        const int target = sequence_.rank(vertex);
        const auto &children = nodes_[blossom].children;
        for (int i = 0; i < int(children.size()); ++i) {
            const Node &child = nodes_[children[i]];
            const int first = sequence_.rank(child.first);
            if (first <= target && target < first + child.size) return i;
        }
        throw std::logic_error("broken blossom interval");
    }
    Arc ring_step(int b, int position, int direction) const {
        const auto &ring = nodes_[b].ring;
        return direction > 0 ? ring[position]
                             : ring[(position + int(ring.size()) - 1) % ring.size()].reversed();
    }
    // Rotate the canonical near-perfect matching inside a blossom so x is its base.
    void expose(int b, int x) {
        if (b < n_ || nodes_[b].base == x) return;
        Node &node = nodes_[b];
        const int start = child_position(b, x), k = int(node.children.size());
        const int direction = (start % 2 == 0) ? -1 : 1;
        expose(node.children[start], x);
        int position = start;
        while (position != 0) {
            const int middle = (position + direction + k) % k;
            const int next = (middle + direction + k) % k;
            const Arc link = ring_step(b, middle, direction);
            expose(node.children[middle], link.from);
            expose(node.children[next], link.to);
            matching_[link.from] = matching_[link.to] = link.edge;
            position = next;
        }
        std::rotate(node.children.begin(), node.children.begin() + start, node.children.end());
        std::rotate(node.ring.begin(), node.ring.begin() + start, node.ring.end());
        node.base = x;
    }
    // Flip an arm from x to its exposed root. edge==-1 relocates exposure to x.
    void flip_arm(int x, int edge) {
        for (;;) {
            const int b = owner(x);
            const Arc back = nodes_[b].entry;
            expose(b, x);
            matching_[x] = edge;
            if (back.edge < 0) break;
            const int inner = owner(back.from);
            const Arc next = nodes_[inner].entry;
            expose(inner, next.to);
            matching_[next.to] = next.edge;
            x = next.from;
            edge = next.edge;
        }
    }
    void use_edge(Arc a) {
        const int from = owner(a.from), to = owner(a.to);
        if (from == to) return;
        if (nodes_[to].label == Unseen) {
            const int base = nodes_[to].base;
            if (matching_[base] < 0) {
                flip_arm(a.from, a.edge);
                flip_arm(a.to, a.edge);
                ++stats_.augmentations;
                finish_phase(nodes_[from].tree);
            } else {
                change_label(to, Inner);
                nodes_[to].entry = a;
                nodes_[to].tree = nodes_[from].tree;
                const Arc matched = arc(matching_[base], base);
                make_outer(owner(matched.to), matched, nodes_[from].tree);
            }
        } else if (nodes_[to].label == Outer) {
            const int common = common_ancestor(from, to);
            if (common >= 0)
                contract(a, common);
            else {
                flip_arm(a.from, a.edge);
                flip_arm(a.to, a.edge);
                ++stats_.augmentations;
                finish_phase(nodes_[from].tree, nodes_[to].tree);
            }
        }
    }
    void expand(int b) {
        const Arc entry = nodes_[b].entry;
        const int tree = nodes_[b].tree, start = child_position(b, entry.to);
        const int k = int(nodes_[b].children.size()), direction = start % 2 == 0 ? -1 : 1;
        std::vector<std::pair<int, int>> order;
        for (int child : nodes_[b].children)
            order.push_back({sequence_.rank(nodes_[child].first), child});
        std::sort(order.begin(), order.end());
        node_events_.erase(b);
        int remaining = nodes_[b].sequence;
        for (auto item : order) {
            const int child = item.second;
            auto parts = sequence_.split(remaining, nodes_[child].size);
            remaining = parts.second;
            Node &node = nodes_[child];
            node.sequence = parts.first;
            sequence_.data[node.sequence].owner = child;
            sequence_.add_to(node.sequence, phase_time_); // Parent was inner: Y=Y'+time.
            node.parent = -1;
            node.label = Unseen;
            node.entry = {};
            node.tree = -1;
            update_event(child);
        }
        int position = start;
        Arc incoming = entry;
        for (;;) {
            const int child = nodes_[b].children[position];
            change_label(child, Inner);
            nodes_[child].entry = incoming;
            nodes_[child].tree = tree;
            if (position == 0) break;
            const Arc matched = ring_step(b, position, direction);
            position = (position + direction + k) % k;
            make_outer(nodes_[b].children[position], matched, tree);
            incoming = ring_step(b, position, direction);
            position = (position + direction + k) % k;
        }
        const int slot = nodes_[b].list_position, last = blossoms_.back();
        blossoms_[slot] = last;
        nodes_[last].list_position = slot;
        blossoms_.pop_back();
        nodes_[b] = Node{};
        available_.push_back(b);
        ++stats_.expansions;
    }
    void scan_edge(int v, int e) {
        ++stats_.edge_scans;
        const int w = other(e, v), bv = owner(v), bw = owner(w);
        if (bv == bw) return;
        const Value gap = slack(e);
        if (gap < -tolerance_) {
            failed_ = true;
            return;
        }
        if (nodes_[bw].label == Outer) {
            if (gap <= tolerance_)
                use_edge(arc(e, v));
            else
                events_.push(
                    {phase_time_ + gap / 2, e, versions_[edges_[e].u], versions_[edges_[e].v]});
        } else {
            remember_incoming(w, arc(e, v));
            if (nodes_[bw].label == Unseen && gap <= tolerance_) use_edge(arc(e, v));
        }
    }
    void dual_event() {
        while (!events_.empty() && !valid_event(events_.top())) events_.pop();
        if (node_events_.heap.empty()) {
            failed_ = true;
            return;
        }
        const int b = node_events_.heap.front();
        Value deadline = node_events_.keys[b];
        bool join = false;
        if (!events_.empty() && events_.top().time < deadline) {
            deadline = events_.top().time;
            join = true;
        }
        Value delta = deadline - phase_time_;
        if (delta < -tolerance_) {
            failed_ = true;
            return;
        }
        delta = std::max(Value(0), delta);
        phase_time_ += delta;
        if constexpr (exact) dual_objective_ -= delta * root_count_;
        ++stats_.dual_updates;
        if constexpr (exact)
            if (dual_objective_ == 2 * incumbent_value_) {
                optimal_ = true;
                return;
            }
        if (join) {
            const int e = events_.top().edge;
            events_.pop();
            use_edge(arc(e, edges_[e].u));
        } else if (nodes_[b].label == Outer) {
            const int v = sequence_.data[nodes_[b].sequence].dual_arg, tree = nodes_[b].tree;
            flip_arm(v, -1);
            ++stats_.relocations;
            finish_phase(tree);
        } else if (nodes_[b].label == Unseen) {
            const int v = sequence_.data[nodes_[b].sequence].incoming_arg;
            use_edge(incoming_[v]);
        } else
            expand(b);
    }
    std::uint64_t advance(std::uint64_t limit) {
        if (!phase_) {
            begin_phase();
            return 1;
        }
        if (scan_vertex_ >= 0 && (versions_[scan_vertex_] != scan_version_ ||
                                  nodes_[owner(scan_vertex_)].label != Outer))
            scan_vertex_ = -1;
        if (scan_vertex_ < 0 || scan_position_ == adjacency_[scan_vertex_].size()) {
            scan_vertex_ = -1;
            while (!queue_.empty()) {
                const auto next = queue_.back();
                queue_.pop_back();
                if (versions_[next.vertex] != next.version ||
                    nodes_[owner(next.vertex)].label != Outer)
                    continue;
                scan_vertex_ = next.vertex;
                scan_version_ = next.version;
                scan_position_ = 0;
                break;
            }
        }
        if (scan_vertex_ < 0) {
            dual_event();
            return 1;
        }
        const int v = scan_vertex_;
        const auto &adjacent = adjacency_[v];
        const std::size_t first = scan_position_,
                          stop = std::min(adjacent.size(), first + std::size_t(limit));
        std::size_t position = first;
        while (position < stop) {
            scan_edge(v, adjacent[position++]);
            if (optimal_ || failed_ || versions_[v] != scan_version_ ||
                nodes_[owner(v)].label != Outer)
                break;
        }
        scan_position_ = position;
        return std::max<std::uint64_t>(1, position - first);
    }
};

} // namespace detail

template <class Weight = std::int64_t> class Solver : public detail::Types<Weight> {
    using Types = detail::Types<Weight>;
    using Narrow = detail::Engine<Weight, std::int64_t>;
    using Wide = detail::Engine<Weight, typename Types::Value>;
    using SparseNarrow = detail::PriorityEngine<Weight, std::int64_t>;
    using SparseWide = detail::PriorityEngine<Weight, typename Types::Value>;
    using Storage =
        std::conditional_t<Types::exact, std::variant<Narrow, Wide, SparseNarrow, SparseWide>,
                           std::variant<Wide, SparseWide>>;
    Storage engine_;

    static Storage initialize(int n, std::vector<typename Types::Edge> edges,
                              const std::vector<int> &seed, Algorithm algorithm) {
        if (n < 0 || n > std::numeric_limits<int>::max() / 2)
            throw std::invalid_argument("invalid vertex count");
        unsigned logarithm = 1;
        for (unsigned k = unsigned(n); k > 1; k >>= 1) ++logarithm;
        const bool sparse =
            algorithm == Algorithm::Sparse ||
            (algorithm == Algorithm::Automatic && edges.size() < std::size_t(n) * n / logarithm);
        if constexpr (Types::exact) {
            Integer magnitude = 0;
            for (const auto &e : edges) {
                const Integer w = e.weight;
                magnitude = std::max(magnitude, w < 0 ? -w : w);
            }
            // Feasible dual objective <= n*max|w|, all dual terms nonnegative.
            // Total dual time <= that bound (at least one root per update).
            // The margin also covers slacks, lazy offsets and heap deadlines.
            const bool narrow = Integer(32) * std::max(1, n) * magnitude <=
                                std::numeric_limits<std::int64_t>::max();
            if (sparse) {
                if (narrow) return Storage(std::in_place_index<2>, n, std::move(edges), seed);
                return Storage(std::in_place_index<3>, n, std::move(edges), seed);
            }
            if (narrow) return Storage(std::in_place_index<0>, n, std::move(edges), seed);
            return Storage(std::in_place_index<1>, n, std::move(edges), seed);
        } else {
            if (sparse) return Storage(std::in_place_index<1>, n, std::move(edges), seed);
            return Storage(std::in_place_index<0>, n, std::move(edges), seed);
        }
    }

public:
    using typename Types::Certificate;
    using typename Types::Edge;
    using typename Types::Result;

    Solver(int vertex_count, std::vector<Edge> edges, const std::vector<int> &initial_matching = {},
           Algorithm algorithm = Algorithm::Automatic)
        : engine_(initialize(vertex_count, std::move(edges), initial_matching, algorithm)) {}

    void reset(const std::vector<int> &initial_matching = {}) {
        std::visit([&](auto &engine) { engine.reset(initial_matching); }, engine_);
    }
    Result solve(const Options &options = {}) {
        return std::visit([&](auto &engine) { return engine.solve(options); }, engine_);
    }
    Result result() const {
        return std::visit([](const auto &engine) { return engine.result(); }, engine_);
    }
    Certificate certificate() const {
        return std::visit([](const auto &engine) { return engine.certificate(); }, engine_);
    }
};

} // namespace maximum_weight_matching
#endif
