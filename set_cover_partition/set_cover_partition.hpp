#ifndef SET_COVER_PARTITION_HPP
#define SET_COVER_PARTITION_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace set_cover_partition {
using Cost = long double;
using Clock = std::chrono::steady_clock;
enum class Mode { Cover, Partition };
enum class StopReason { TimeLimit, IterationLimit, NoCoveringSet, Exhausted };

struct Options {
    double time_limit = 1.0; // Seconds per call.
    Clock::time_point deadline = Clock::time_point::max();
    double safety_margin = 0.0005;
    std::uint64_t seed = 1;
    std::uint64_t iteration_limit = std::numeric_limits<std::uint64_t>::max();
    int swap_origins = 16;             // Selected columns per 2-flip scan.
    bool composite_moves = true;       // Combine promising swaps into 3-/4-flips.
    bool adaptive_penalties = true;    // Both increase and decrease row penalties.
    bool priced_restarts = true;       // Construct covers from dual reduced costs.
    bool fast_cover = true;            // Row-driven incremental search for Cover.
    bool cover_cost_bound = true;      // Search below the incumbent cover cost.
    bool incremental_penalties = true; // Update only the changed penalty terms.
    int ejection_moves = 256;          // Partition-only batch; 0 disables it.
    int repair_interval = 64;          // 0 disables large neighborhoods.
    int repair_nodes = 1500;           // Per neighborhood, across run calls.
    int repair_columns_per_row = 64;   // 0: all eligible columns.
    int repair_max_remove = 32;
    int dual_iterations = 128;     // 0 disables Lagrangian pricing.
    std::vector<int> initial_sets; // A hint, never fixed constraints.
};

struct Statistics {
    double elapsed = 0;
    std::uint64_t iterations = 0, flips = 0, swaps = 0, improvements = 0;
    std::uint64_t penalty_updates = 0, repairs = 0, repair_nodes = 0;
    std::uint64_t dual_iterations = 0;
    std::uint64_t composite_moves = 0, penalty_decreases = 0;
    std::uint64_t cover_cycles = 0;
    std::uint64_t ejection_steps = 0, ejection_batches = 0, ejection_improvements = 0;
};
struct Result {
    std::vector<int> selected;
    std::vector<int> coverage;
    Cost cost = 0; // Also defined for an infeasible result.
    Cost lower_bound = -std::numeric_limits<Cost>::infinity();
    std::int64_t violation = 0; // Missing + (Partition only) excess.
    int uncovered = 0;
    bool feasible = false;
    StopReason reason = StopReason::TimeLimit;
    Statistics stats; // Counts cumulative; elapsed for this call.
};

// Build once, then keep unchanged and alive while any Solver refers to it.
// Model construction/validation is intentionally outside the search time budget.
class Model {
public:
    struct Set {
        Cost cost;
        std::vector<int> elements;
        std::vector<std::uint64_t> bits;
    };

private:
    friend class Solver;
    int n_;
    Mode mode_;
    std::vector<Set> sets_;
    std::vector<std::vector<int>> rows_;
    std::unordered_map<std::uint64_t, std::vector<int>> hash_groups_;
    std::vector<int> group_, representative_;
    Cost magnitude_ = 0, scale_ = 0;
    std::uint64_t revision_ = 0;

public:
    explicit Model(int elements, Mode mode = Mode::Cover) : n_(elements), mode_(mode) {
        if (n_ < 0 || (mode != Mode::Cover && mode != Mode::Partition))
            throw std::invalid_argument("set_cover_partition: invalid model");
        rows_.resize(static_cast<std::size_t>(n_));
    }
    int add_set(Cost cost, std::vector<int> elements) {
        if (!std::isfinite(cost) ||
            std::abs(cost) > std::numeric_limits<Cost>::max() / 4096 - magnitude_)
            throw std::invalid_argument("set_cover_partition: non-finite or excessive cost");
        if (sets_.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("set_cover_partition: too many sets");
        std::sort(elements.begin(), elements.end());
        if (std::adjacent_find(elements.begin(), elements.end()) != elements.end())
            throw std::invalid_argument("set_cover_partition: repeated element");
        for (int i : elements)
            if (i < 0 || i >= n_)
                throw std::invalid_argument("set_cover_partition: element out of range");
        const int id = static_cast<int>(sets_.size());
        std::uint64_t hash = UINT64_C(1469598103934665603);
        for (int i : elements) {
            hash ^= static_cast<std::uint64_t>(i) + 1;
            hash *= UINT64_C(1099511628211);
        }
        auto &bucket = hash_groups_[hash];
        int group = -1;
        for (int g : bucket)
            if (sets_[representative_[g]].elements == elements) {
                group = g;
                break;
            }
        if (group_.size() == group_.capacity())
            group_.reserve(std::max<std::size_t>(4, 2 * group_.size()));
        if (group < 0) {
            if (bucket.size() == bucket.capacity())
                bucket.reserve(std::max<std::size_t>(4, 2 * bucket.size()));
            if (representative_.size() == representative_.capacity())
                representative_.reserve(std::max<std::size_t>(4, 2 * representative_.size()));
        }
        Set s{cost, std::move(elements),
              std::vector<std::uint64_t>((static_cast<std::size_t>(n_) + 63) / 64)};
        for (int i : s.elements) s.bits[i / 64] |= UINT64_C(1) << (i % 64);
        // Reserve first: allocation failure must not leave half-inserted incidences.
        for (int i : s.elements)
            if (rows_[i].size() == rows_[i].capacity())
                rows_[i].reserve(std::max<std::size_t>(4, 2 * rows_[i].size()));
        sets_.push_back(std::move(s));
        if (group < 0) {
            group = static_cast<int>(representative_.size());
            bucket.push_back(group);
            representative_.push_back(id);
        } else if (cost < sets_[representative_[group]].cost)
            representative_[group] = id;
        group_.push_back(group);
        for (int i : sets_.back().elements) rows_[i].push_back(id);
        magnitude_ += std::abs(cost);
        scale_ = std::max(scale_, std::abs(cost));
        ++revision_;
        return id;
    }
    int elements() const { return n_; }
    int sets() const { return static_cast<int>(sets_.size()); }
    Mode mode() const { return mode_; }
    const Set &set(int j) const { return sets_.at(static_cast<std::size_t>(j)); }
    const std::vector<int> &row(int i) const { return rows_.at(static_cast<std::size_t>(i)); }
    Cost scale() const { return scale_ > 0 ? scale_ : 1; }
    std::uint64_t revision() const { return revision_; }
    int cheapest_equivalent(int j) const { return representative_[group_[j]]; }

    Result evaluate(const std::vector<int> &selected) const {
        Result r;
        r.selected = selected;
        r.coverage.assign(n_, 0);
        std::vector<unsigned char> seen(sets_.size(), 0);
        Cost correction = 0;
        for (int j : selected) {
            if (j < 0 || j >= sets() || seen[j])
                throw std::invalid_argument("set_cover_partition: invalid selected set");
            seen[j] = 1;
            const Cost y = sets_[j].cost - correction, t = r.cost + y;
            correction = (t - r.cost) - y;
            r.cost = t;
            for (int i : sets_[j].elements) ++r.coverage[i];
        }
        for (int k : r.coverage) {
            if (!k) {
                ++r.uncovered;
                ++r.violation;
            } else if (mode_ == Mode::Partition)
                r.violation += k - 1;
        }
        r.feasible = r.violation == 0;
        return r;
    }
};

namespace detail {
struct Budget {
    Clock::time_point start, stop;
    unsigned ticks = 0;
    Budget(Clock::time_point begin, double seconds, Clock::time_point outer, double margin)
        : start(begin), stop(outer) {
        if (std::isnan(seconds) || seconds < 0)
            throw std::invalid_argument("set_cover_partition: invalid time limit");
        const double range = std::chrono::duration<double>(Clock::time_point::max() - begin).count();
        if (seconds < range)
            stop = std::min(stop, begin + std::chrono::duration_cast<Clock::duration>(
                                              std::chrono::duration<double>(seconds)));
        if (stop <= begin) {
            stop = begin;
            return;
        }
        if (stop != Clock::time_point::max()) {
            const double left = std::chrono::duration<double>(stop - begin).count();
            stop -= std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(std::min(left * 0.05, margin)));
        }
    }
    bool expired() const { return stop != Clock::time_point::max() && Clock::now() >= stop; }
    bool poll() { return (++ticks & 63U) == 0 && expired(); }
    // Tiny incidence updates need less frequent clock reads than whole columns.
    bool sparse_poll() { return (++ticks & 1023U) == 0 && expired(); }
    // Same clock-read positions as poll()/sparse_poll(), without a per-item
    // branch inside a tight loop. Commit its cursor after processing the batch.
    template <class Index> Index batch_end(Index pos, Index limit, unsigned mask = 63U) {
        if ((++ticks & mask) == 0 && expired()) return pos;
        const Index count = std::min<Index>(limit - pos, mask + 1U - (ticks & mask));
        ticks += static_cast<unsigned>(count) - 1U;
        return pos + count;
    }
};
inline int popcount(std::uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(v);
#else
    v -= (v >> 1) & UINT64_C(0x5555555555555555);
    v = (v & UINT64_C(0x3333333333333333)) + ((v >> 2) & UINT64_C(0x3333333333333333));
    return static_cast<int>(
        (((v + (v >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f)) * UINT64_C(0x0101010101010101)) >> 56);
#endif
}
inline int trailing_zeros(std::uint64_t v) { // v != 0
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(v);
#else
    int count = 0;
    while (!(v & 1U)) {
        v >>= 1;
        ++count;
    }
    return count;
#endif
}
} // namespace detail

class Solver {
    const Model &p_;
    Options o_;
    const int n_, m_;
    const bool partition_;
    const std::uint64_t revision_;
    std::uint64_t rng_;
    std::vector<unsigned char> x_, locked_;
    // Enumerate selected IDs in the original ascending order without scanning
    // every candidate byte. This also preserves objective summation/RNG order.
    std::vector<std::uint64_t> selected_bits_;
    int selected_count_ = 0;
    std::vector<int> count_, fixed_count_;
    std::vector<std::vector<int>> active_rows_;
    std::vector<double> c_, down_, up_, delta_, base_weight_, price_, best_price_;
    std::vector<double> greedy_gain_;
    std::vector<int> free_count_;
    std::vector<std::uint64_t> tabu_;
    Result best_;
    Statistics stats_;
    bool impossible_ = false, greedy_ = true;
    int init_phase_ = 0, init_col_ = 0, init_row_ = 0;
    std::size_t fixed_on_count_ = 0;
    std::size_t init_incidence_ = 0;
    double init_price_ = std::numeric_limits<double>::infinity();
    int greedy_col_ = 0, greedy_chosen_ = -1;
    double greedy_score_ = std::numeric_limits<double>::infinity();
    bool random_greedy_ = false;
    std::vector<double> greedy_noise_, greedy_weight_;
    std::vector<double> greedy_cost_;
    int greedy_setup_ = -1;
    double greedy_alpha_ = 1;
    bool greedy_priced_ = false;
    int rebuild_ = 0;
    std::uint64_t last_repair_ = 0, stalls_ = 0;
    Cost bound_ = -std::numeric_limits<Cost>::infinity();
    // A flip's incidence updates may be paused. Never expose this in-flight state.
    std::vector<int> pending_;
    std::size_t pending_at_ = 0, update_row_ = 0, update_col_ = 0;
    bool flip_started_ = false, row_started_ = false;
    double row_add_change_ = 0, row_remove_change_ = 0;
    double row_gain_change_ = 0;
    int row_free_change_ = 0;
    // Persistent 1-/2-flip scan; the state and penalties are fixed until it finishes.
    int scan_ = 0, single_ = -1;
    double single_value_ = 0;
    std::vector<int> origins_, touched_;
    std::vector<double> correction_;
    std::vector<std::uint32_t> touched_stamp_;
    std::uint32_t pair_stamp_ = 1;
    std::size_t origin_ = 0, pair_row_ = 0, pair_col_ = 0, pair_eval_ = 0;
    bool origins_ready_ = false, pair_ready_ = false;
    int swap_out_ = -1, swap_in_ = -1;
    double swap_value_ = 0;
    std::size_t origin_rotation_ = 0;
    std::vector<std::pair<int, int>> exchanges_;
    std::vector<std::pair<double, int>> partners_;
    std::size_t composite_a_ = 0, composite_b_ = 1;
    std::vector<int> net_change_, changed_rows_;
    std::vector<double> penalty_down_change_, penalty_up_change_;
    int penalty_row_ = -1, penalty_scale_col_ = 0;
    std::size_t penalty_incidence_ = 0;
    double penalty_scale_ = 1;
    // Lagrangian subgradient passes are also resumable.
    bool dual_active_ = false;
    int dual_col_ = 0;
    double dual_value_ = 0, dual_abs_ = 0;
    std::vector<int> dual_count_;
    // Nonnegative transformed costs for weighted exact-cover neighborhood search.
    std::vector<double> shift_, search_cost_;
    struct Frame {
        std::vector<std::uint64_t> remaining;
        double cost = 0;
        double lower = 0;
        int incoming = -1, scan_row = 0;
        std::size_t scan_col = 0, next = 0;
        std::vector<int> candidates, trial;
        bool ready = false;
    };
    struct Repair {
        bool active = false, prepared = false;
        std::vector<int> fixed, removed, answer;
        std::vector<unsigned char> fixed_flag, preferred;
        std::vector<std::uint64_t> free;
        std::vector<std::vector<int>> rows;
        std::vector<double> lower_weight;
        std::vector<Frame> stack;
        int prepare_row = 0;
        std::size_t prepare_col = 0;
        std::vector<std::pair<double, int>> ranked;
        double min_share = std::numeric_limits<double>::infinity();
        std::uint64_t nodes = 0;
        double limit = std::numeric_limits<double>::infinity();
        Cost fixed_cost = 0;
    } repair_;
    // Cover-only state. gain[j] counts uncovered row weights; loss[j] counts
    // singly covered row weights. XOR identifies the sole covering column.
    struct CoverSearch {
        std::vector<double> weight, gain, loss;
        std::vector<int> owner, chosen, fixed, position, missing, missing_position;
        std::vector<std::uint64_t> age, banned;
        Cost cost = 0, correction = 0;
        int phase = 0, column = 0, scan = 0, candidate = -1, row = -1;
        int removals = 0, prune = 0, weight_row = 0;
        int remove_limit = 2, add_candidate = -1;
        bool random_remove = true, weight_after_remove = false;
        double value = 0;
        bool fallback = false, weight_started = false;
        std::size_t incidence = 0;
        int pending = -1, sign = 0, next_phase = 0;
        std::size_t flip_row = 0, flip_incidence = 0;
        bool flip_started = false, row_started = false;
        double gain_change = 0;
    } cover_;
    struct Ejection {
        bool active = false;
        int phase = 0, column = 0, scan = 0, row = -1, candidate = -1, weight_row = 0;
        int moves = 0, apply_at = 0;
        std::size_t incidence = 0, apply_row = 0;
        bool weight_started = false, apply_started = false;
        double value = 0;
        Cost cost = 0, initial_cost = 0;
        std::uint64_t stamp = 0;
        std::vector<int> owner, selected, fixed, position, missing, missing_position, move;
        std::vector<unsigned char> chosen;
        std::vector<std::uint64_t> seen, banned, age;
        std::vector<double> base, weight, score;
    } ejection_;
    std::uint64_t ejection_next_ = 256, ejection_wait_ = 256;
    int ejection_successes_ = 0;

    // All internal IDs originate from validated model incidences or hints.
    // The public Model::set() accessor retains its checked access.
    const Model::Set &set(int j) const { return p_.sets_[j]; }

    void start_ejection() {
        auto &s = ejection_;
        s = Ejection{};
        s.active = true;
        s.initial_cost = best_.cost;
        s.owner.assign(n_, -1);
        s.position.assign(m_, -1);
        s.chosen.assign(m_, 0);
        s.missing_position.assign(n_, -1);
        s.seen.assign(m_, 0);
        s.banned.assign(m_, 0);
        s.age.assign(m_, 0);
        s.base.assign(n_, std::numeric_limits<double>::infinity());
        s.weight.resize(n_);
        s.score.resize(m_);
        for (int j : best_.selected) {
            if (locked_[j] == 1) {
                if (set(j).elements.empty()) continue;
                j = p_.cheapest_equivalent(j);
            }
            s.cost += set(j).cost;
            if (locked_[j] == 2) {
                s.fixed.push_back(j);
                continue;
            }
            s.position[j] = static_cast<int>(s.selected.size());
            s.selected.push_back(j);
            s.chosen[j] = 1;
            for (int i : set(j).elements) s.owner[i] = j;
        }
        ++stats_.ejection_batches;
    }
    void ejection_remember() {
        auto &s = ejection_;
        if (!s.missing.empty() || s.cost >= best_.cost) return;
        std::vector<int> selected = s.fixed;
        selected.insert(selected.end(), s.selected.begin(), s.selected.end());
        std::sort(selected.begin(), selected.end());
        Result candidate = p_.evaluate(selected);
        s.cost = candidate.cost;
        if (candidate.feasible && candidate.cost < best_.cost) {
            best_ = std::move(candidate);
            ++stats_.improvements;
            ++stats_.ejection_improvements;
            stalls_ = 0;
        }
    }
    bool ejection_step(detail::Budget &budget) {
        auto &s = ejection_;
        if (s.phase == 0) {
            while (s.column < m_) {
                const auto end = budget.batch_end(s.column, m_);
                if (end == s.column) return false;
                for (auto at = s.column; at < end; ++at) {
                    const int j = at;
                    if (locked_[j]) continue;
                    const double share =
                        search_cost_[j] / std::max<std::size_t>(1, set(j).elements.size());
                    for (int i : set(j).elements) s.base[i] = std::min(s.base[i], share);
                }
                s.column = end;
            }
            for (int i = 0; i < n_; ++i) {
                s.base[i] = std::max(1e-9, std::isfinite(s.base[i]) ? s.base[i] : 1.0);
                s.weight[i] = 2 * s.base[i];
                if (s.owner[i] < 0) {
                    s.missing_position[i] = static_cast<int>(s.missing.size());
                    s.missing.push_back(i);
                }
            }
            s.column = 0;
            s.phase = 1;
        }
        if (s.phase == 1) {
            while (s.column < m_) {
                const auto end = budget.batch_end(s.column, m_);
                if (end == s.column) return false;
                for (auto at = s.column; at < end; ++at) {
                    const int j = at;
                    s.score[j] = search_cost_[j];
                    for (int i : set(j).elements) s.score[j] -= s.weight[i];
                }
                s.column = end;
            }
            s.phase = 2;
        }
        if (s.phase == 2) {
            if (s.moves >= o_.ejection_moves || s.selected.empty()) {
                const bool improved = best_.cost < s.initial_cost;
                s.active = false;
                if (improved) {
                    load(best_.selected);
                    ejection_wait_ = 256;
                    ++ejection_successes_;
                } else {
                    ejection_wait_ = std::min<std::uint64_t>(8192, ejection_wait_ * 2);
                    ejection_successes_ = 0;
                }
                ejection_next_ =
                    stats_.iterations + (improved && ejection_successes_ < 4 ? 0 : ejection_wait_);
                if (ejection_successes_ >= 4) ejection_successes_ = 0;
                return true;
            }
            if (s.missing.empty()) {
                ejection_remember();
                for (double &value : s.weight) value *= 0.9;
                s.column = 0;
                s.phase = 3;
            } else {
                s.phase = 4;
                s.row = -1;
            }
        }
        if (s.phase == 3) {
            while (s.column < m_) {
                const auto end = budget.batch_end(s.column, m_, 1023U);
                if (end == s.column) return false;
                for (auto at = s.column; at < end; ++at) {
                    const int j = at;
                    s.score[j] = search_cost_[j] - 0.9 * (search_cost_[j] - s.score[j]);
                }
                s.column = end;
            }
            int out = -1;
            double value = std::numeric_limits<double>::infinity();
            for (int j : s.selected) {
                const double v = -s.score[j] * (0.8 + 0.4 * uniform());
                if (v < value) {
                    value = v;
                    out = j;
                }
            }
            if (out >= 0) {
                s.move = {out};
                s.apply_at = 0;
                s.apply_row = 0;
                s.apply_started = false;
                s.phase = 5;
                s.candidate = -1;
            }
            return true;
        }
        if (s.phase == 4) {
            if (s.missing.empty()) {
                s.phase = 2;
                return true;
            }
            if (s.row < 0) {
                s.row = s.missing[random() % s.missing.size()];
                s.scan = 0;
                s.candidate = -1;
                s.value = std::numeric_limits<double>::infinity();
            }
            const auto &candidates = row(s.row);
            while (s.scan < static_cast<int>(candidates.size())) {
                const auto end = budget.batch_end(s.scan, static_cast<int>(candidates.size()));
                if (end == s.scan) return false;
                for (auto at = s.scan; at < end; ++at) {
                    const int j = candidates[at];
                    if (s.banned[j] > stats_.ejection_steps) continue;
                    double value = s.score[j];
                    if (++s.stamp == 0) {
                        std::fill(s.seen.begin(), s.seen.end(), 0);
                        ++s.stamp;
                    }
                    for (int i : set(j).elements) {
                        const int k = s.owner[i];
                        if (k >= 0 && s.seen[k] != s.stamp) {
                            s.seen[k] = s.stamp;
                            value -= s.score[k];
                        }
                    }
                    if (value < s.value ||
                        (value == s.value && s.candidate >= 0 && s.age[j] < s.age[s.candidate])) {
                        s.value = value;
                        s.candidate = j;
                    }
                }
                s.scan = end;
            }
            s.move.clear();
            if (s.candidate >= 0) {
                for (int i : set(s.candidate).elements) {
                    const int k = s.owner[i];
                    if (k >= 0 && std::find(s.move.begin(), s.move.end(), k) == s.move.end())
                        s.move.push_back(k);
                }
                s.move.push_back(s.candidate);
            }
            s.apply_at = 0;
            s.apply_row = 0;
            s.apply_started = false;
            s.phase = 5;
        }
        if (s.phase == 5) {
            while (s.apply_at < static_cast<int>(s.move.size())) {
                const int j = s.move[s.apply_at];
                const bool adding = j == s.candidate;
                if (!s.apply_started) {
                    s.chosen[j] = adding;
                    s.cost += (adding ? 1 : -1) * set(j).cost;
                    if (adding) {
                        s.position[j] = static_cast<int>(s.selected.size());
                        s.selected.push_back(j);
                        s.age[j] = stats_.ejection_steps;
                    } else {
                        const int last = s.selected.back();
                        s.selected[s.position[j]] = last;
                        s.position[last] = s.position[j];
                        s.selected.pop_back();
                        s.position[j] = -1;
                        s.banned[j] = stats_.ejection_steps + 5;
                    }
                    s.apply_started = true;
                }
                const auto &elements = set(j).elements;
                while (s.apply_row < elements.size()) {
                    if (budget.poll()) return false;
                    const int i = elements[s.apply_row++];
                    s.owner[i] = adding ? j : -1;
                    if (adding) {
                        const int last = s.missing.back();
                        s.missing[s.missing_position[i]] = last;
                        s.missing_position[last] = s.missing_position[i];
                        s.missing.pop_back();
                        s.missing_position[i] = -1;
                    } else {
                        s.missing_position[i] = static_cast<int>(s.missing.size());
                        s.missing.push_back(i);
                    }
                }
                ++s.apply_at;
                s.apply_row = 0;
                s.apply_started = false;
            }
            if (s.candidate < 0 && !s.move.empty()) {
                s.phase = s.missing.empty() ? 2 : 4;
                s.row = -1;
                return true;
            }
            ++stats_.iterations;
            ++stats_.ejection_steps;
            ++s.moves;
            ejection_remember();
            s.weight_row = 0;
            s.incidence = 0;
            s.weight_started = false;
            s.phase = 6;
        }
        if (s.phase == 6) {
            while (s.weight_row < static_cast<int>(s.missing.size())) {
                const int i = s.missing[s.weight_row];
                if (!s.weight_started) {
                    s.weight[i] += s.base[i];
                    s.weight_started = true;
                }
                const auto &candidates = row(i);
                const double change = s.base[i];
                while (s.incidence < candidates.size()) {
                    const auto end = budget.batch_end(
                        s.incidence, static_cast<decltype(s.incidence)>(candidates.size()), 1023U);
                    if (end == s.incidence) return false;
                    for (auto at = s.incidence; at < end; ++at) {
                        s.score[candidates[at]] -= change;
                    }
                    s.incidence = end;
                }
                ++s.weight_row;
                s.incidence = 0;
                s.weight_started = false;
            }
            s.phase = 2;
            return true;
        }
        return true;
    }

    void cover_cost(Cost change) {
        auto &s = cover_;
        const Cost y = change - s.correction, next = s.cost + y;
        s.correction = (next - s.cost) - y;
        s.cost = next;
    }
    void cover_remember() {
        auto &s = cover_;
        if (s.missing.size() > static_cast<std::size_t>(best_.violation) ||
            (s.missing.size() == static_cast<std::size_t>(best_.violation) && s.cost >= best_.cost))
            return;
        Result r;
        r.selected = s.fixed;
        r.selected.insert(r.selected.end(), s.chosen.begin(), s.chosen.end());
        std::sort(r.selected.begin(), r.selected.end());
        r.coverage = count_;
        r.uncovered = static_cast<int>(s.missing.size());
        r.violation = r.uncovered;
        r.feasible = !r.uncovered;
        Cost correction = 0;
        for (int j : r.selected) {
            const Cost y = set(j).cost - correction, next = r.cost + y;
            correction = (next - r.cost) - y;
            r.cost = next;
        }
        s.cost = r.cost;
        s.correction = 0;
        if (r.violation < best_.violation || r.cost < best_.cost) {
            best_ = std::move(r);
            ++stats_.improvements;
        }
    }
    void cover_queue(int j, int next_phase) {
        auto &s = cover_;
        s.pending = j;
        s.next_phase = next_phase;
        s.phase = 8;
        s.flip_row = s.flip_incidence = 0;
        s.flip_started = s.row_started = false;
    }
    bool cover_step(detail::Budget &budget) {
        auto &s = cover_;
        if (s.phase == 0) {
            s.weight.assign(n_, 1);
            s.gain.assign(m_, 0);
            s.loss.assign(m_, 0);
            s.owner.assign(n_, 0);
            s.position.assign(m_, -1);
            s.missing_position.assign(n_, -1);
            s.age.assign(m_, 0);
            s.banned.assign(m_, 0);
            for (int i = 0; i < n_; ++i)
                if (!count_[i]) {
                    s.missing_position[i] = static_cast<int>(s.missing.size());
                    s.missing.push_back(i);
                }
            s.phase = 1;
        }
        if (s.phase == 1) {
            while (s.column < m_) {
                const auto end = budget.batch_end(s.column, m_);
                if (end == s.column) return false;
                for (auto at = s.column; at < end; ++at) {
                    const int j = at;
                    if (x_[j]) {
                        cover_cost(set(j).cost);
                        if (locked_[j] == 2)
                            s.fixed.push_back(j);
                        else {
                            s.position[j] = static_cast<int>(s.chosen.size());
                            s.chosen.push_back(j);
                        }
                    }
                    for (int i : set(j).elements) {
                        if (!count_[i]) s.gain[j] += 1;
                        if (x_[j]) {
                            s.owner[i] ^= j + 1;
                            if (count_[i] == 1) s.loss[j] += 1;
                        }
                    }
                }
                s.column = end;
            }
            s.phase = 2;
            s.scan = 0;
            s.candidate = -1;
            s.value = -1;
        }
        if (s.phase == 2) {
            if (s.missing.empty()) {
                s.phase = 3;
                s.prune = 0;
                return true;
            }
            while (s.scan < m_) {
                const auto end = budget.batch_end(s.scan, m_, 1023U);
                if (end == s.scan) return false;
                for (auto at = s.scan; at < end; ++at) {
                    const int j = at;
                    if (locked_[j] || x_[j] || s.gain[j] <= 0) continue;
                    const double value =
                        c_[j] > 0 ? s.gain[j] / c_[j] : std::numeric_limits<double>::infinity();
                    if (value > s.value) {
                        s.value = value;
                        s.candidate = j;
                    }
                }
                s.scan = end;
            }
            if (s.candidate >= 0) cover_queue(s.candidate, 2);
            s.scan = 0;
            s.candidate = -1;
            s.value = -1;
            return true;
        }
        if (s.phase == 3) {
            while (s.prune < static_cast<int>(s.chosen.size())) {
                if (budget.sparse_poll()) return false;
                const int j = s.chosen[s.prune];
                if (s.loss[j] <= 0) {
                    cover_queue(j, 3);
                    return true;
                }
                ++s.prune;
            }
            cover_remember();
            s.phase = 4;
            return true;
        }
        if (s.phase == 4) {
            ++stats_.cover_cycles;
            s.removals = 0;
            s.scan = 0;
            s.remove_limit = o_.cover_cost_bound ? 1 : 2;
            s.random_remove = !o_.cover_cost_bound;
            s.weight_after_remove = false;
            s.candidate = -1;
            s.value = std::numeric_limits<double>::infinity();
            s.phase = 5;
        }
        if (s.phase == 5) {
            if (s.removals >= s.remove_limit || s.chosen.empty()) {
                if (s.chosen.empty() && s.removals == 0) ++stats_.iterations;
                s.phase = s.weight_after_remove ? 7 : 6;
                s.row = -1;
                s.weight_row = 0;
                s.incidence = 0;
                s.weight_started = false;
                return true;
            }
            while (s.scan < static_cast<int>(s.chosen.size())) {
                if (budget.sparse_poll()) return false;
                const int j = s.chosen[s.scan++];
                if (c_[j] <= 0) continue; // Free columns never need destruction.
                const double value =
                    s.loss[j] / c_[j] * (s.random_remove ? 0.8 + 0.4 * uniform() : 1.0);
                if (value < s.value ||
                    (value == s.value && s.candidate >= 0 && s.age[j] < s.age[s.candidate])) {
                    s.value = value;
                    s.candidate = j;
                }
            }
            if (s.candidate < 0) {
                s.phase = 6;
                s.row = -1;
                ++stats_.iterations;
                return true;
            }
            s.banned[s.candidate] = stats_.cover_cycles + (o_.cover_cost_bound ? 3 : 0);
            ++s.removals;
            cover_queue(s.candidate, 5);
            s.scan = 0;
            s.candidate = -1;
            s.value = std::numeric_limits<double>::infinity();
            return true;
        }
        if (s.phase == 6) {
            if (s.missing.empty()) {
                s.phase = 3;
                s.prune = 0;
                return true;
            }
            if (s.row < 0) {
                if (o_.cover_cost_bound) ++stats_.cover_cycles;
                s.row = s.missing[random() % s.missing.size()];
                s.scan = 0;
                s.candidate = -1;
                s.value = -1;
                s.fallback = false;
            }
            const auto &candidates = row(s.row);
            while (s.scan < static_cast<int>(candidates.size())) {
                const auto end =
                    budget.batch_end(s.scan, static_cast<int>(candidates.size()), 1023U);
                if (end == s.scan) return false;
                for (auto at = s.scan; at < end; ++at) {
                    const int j = candidates[at];
                    if (x_[j] ||
                        (!s.fallback && (o_.cover_cost_bound ? s.banned[j] > stats_.cover_cycles
                                                             : s.banned[j] == stats_.cover_cycles)))
                        continue;
                    const double value =
                        c_[j] > 0 ? s.gain[j] / c_[j] : std::numeric_limits<double>::infinity();
                    if (value > s.value ||
                        (value == s.value && s.candidate >= 0 && s.age[j] < s.age[s.candidate])) {
                        s.value = value;
                        s.candidate = j;
                    }
                }
                s.scan = end;
            }
            if (o_.cover_cost_bound &&
                (s.candidate < 0 || s.cost + set(s.candidate).cost >= best_.cost)) {
                s.add_candidate = s.candidate;
                s.phase = 9;
                s.scan = 0;
                s.candidate = -1;
                s.value = std::numeric_limits<double>::infinity();
                return true;
            }
            if (s.candidate < 0 && !s.fallback) {
                s.fallback = true;
                s.scan = 0;
                return true;
            }
            if (s.candidate >= 0) {
                cover_queue(s.candidate, 7);
                s.weight_row = 0;
                s.incidence = 0;
                s.weight_started = false;
            }
            return true;
        }
        if (s.phase == 7) {
            while (s.weight_row < static_cast<int>(s.missing.size())) {
                const int i = s.missing[s.weight_row];
                if (!s.weight_started) {
                    s.weight[i] += 1;
                    s.weight_started = true;
                }
                const auto &candidates = row(i);
                while (s.incidence < candidates.size()) {
                    const auto end = budget.batch_end(
                        s.incidence, static_cast<decltype(s.incidence)>(candidates.size()), 1023U);
                    if (end == s.incidence) return false;
                    for (auto at = s.incidence; at < end; ++at) { s.gain[candidates[at]] += 1; }
                    s.incidence = end;
                }
                ++s.weight_row;
                s.incidence = 0;
                s.weight_started = false;
            }
            s.phase = 6;
            s.row = -1;
            return true;
        }
        if (s.phase == 8) {
            const int j = s.pending;
            if (!s.flip_started) {
                s.sign = x_[j] ? -1 : 1;
                flip_selected(j);
                cover_cost(s.sign * set(j).cost);
                s.age[j] = stats_.cover_cycles;
                if (s.sign > 0) {
                    s.position[j] = static_cast<int>(s.chosen.size());
                    s.chosen.push_back(j);
                } else {
                    const int last = s.chosen.back();
                    s.chosen[s.position[j]] = last;
                    s.position[last] = s.position[j];
                    s.chosen.pop_back();
                    s.position[j] = -1;
                }
                s.flip_started = true;
            }
            const auto &elements = set(j).elements;
            while (s.flip_row < elements.size()) {
                const int i = elements[s.flip_row];
                if (!s.row_started) {
                    const int old = count_[i];
                    if (old == 1) s.loss[s.owner[i] - 1] -= s.weight[i];
                    s.owner[i] ^= j + 1;
                    count_[i] += s.sign;
                    if (count_[i] == 1) s.loss[s.owner[i] - 1] += s.weight[i];
                    s.gain_change = 0;
                    if (old == 0) {
                        const int last = s.missing.back();
                        s.missing[s.missing_position[i]] = last;
                        s.missing_position[last] = s.missing_position[i];
                        s.missing.pop_back();
                        s.missing_position[i] = -1;
                        s.gain_change = -s.weight[i];
                    } else if (!count_[i]) {
                        s.missing_position[i] = static_cast<int>(s.missing.size());
                        s.missing.push_back(i);
                        s.gain_change = s.weight[i];
                    }
                    s.row_started = true;
                }
                if (s.gain_change != 0) {
                    const auto &candidates = row(i);
                    const double change = s.gain_change;
                    while (s.flip_incidence < candidates.size()) {
                        const auto end = budget.batch_end(
                            s.flip_incidence,
                            static_cast<decltype(s.flip_incidence)>(candidates.size()), 1023U);
                        if (end == s.flip_incidence) return false;
                        for (auto at = s.flip_incidence; at < end; ++at) {
                            s.gain[candidates[at]] += change;
                        }
                        s.flip_incidence = end;
                    }
                }
                ++s.flip_row;
                s.flip_incidence = 0;
                s.row_started = false;
            }
            s.pending = -1;
            s.phase = s.next_phase;
            ++stats_.iterations;
            ++stats_.flips;
            cover_remember();
            return true;
        }
        if (s.phase == 9) {
            while (s.scan < static_cast<int>(s.chosen.size())) {
                if (budget.sparse_poll()) return false;
                const int j = s.chosen[s.scan++];
                if (c_[j] <= 0) continue;
                const double value = s.loss[j] / c_[j];
                if (value < s.value ||
                    (value == s.value && s.candidate >= 0 && s.age[j] < s.age[s.candidate])) {
                    s.value = value;
                    s.candidate = j;
                }
            }
            const int in = s.add_candidate, out = s.candidate;
            const double add_score = in < 0       ? -1
                                     : c_[in] > 0 ? s.gain[in] / c_[in]
                                                  : std::numeric_limits<double>::infinity();
            if (in >= 0 && out >= 0 && add_score > s.value &&
                s.cost + set(in).cost - set(out).cost < best_.cost) {
                s.banned[out] = stats_.cover_cycles + 3;
                cover_queue(out, 10);
                return true;
            }
            s.phase = 5;
            s.removals = 0;
            s.remove_limit = 2;
            s.random_remove = true;
            s.weight_after_remove = true;
            s.scan = 0;
            s.candidate = -1;
            s.value = std::numeric_limits<double>::infinity();
            return true;
        }
        if (s.phase == 10) {
            cover_queue(s.add_candidate, 7);
            s.weight_row = 0;
            s.incidence = 0;
            s.weight_started = false;
        }
        return true;
    }

    std::uint64_t random() {
        std::uint64_t z = (rng_ += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(random() >> 11) * 0x1.0p-53; }
    const std::vector<int> &row(int i) const { return active_rows_[i]; }
    bool initialize(detail::Budget &budget) {
        while (init_phase_ < 2) {
            while (init_col_ < m_) {
                const auto end = budget.batch_end(init_col_, m_);
                if (end == init_col_) return false;
                for (auto at = init_col_; at < end; ++at) {
                    const int j = at;
                    if (init_phase_ == 0) {
                        c_[j] = static_cast<double>(set(j).cost / p_.scale());
                        if (set(j).elements.empty() || (!partition_ && set(j).cost < 0)) {
                            x_[j] = set(j).cost < 0;
                            locked_[j] = x_[j] ? 2 : 1;
                            fixed_on_count_ += x_[j] != 0;
                        } else if (p_.cheapest_equivalent(j) != j) {
                            if (x_[j]) x_[p_.cheapest_equivalent(j)] = 1;
                            locked_[j] = 1;
                            x_[j] = 0;
                        }
                        if (partition_ && c_[j] < 0)
                            for (int i : set(j).elements)
                                shift_[i] = std::min(shift_[i], c_[j] / set(j).elements.size());
                    } else {
                        double v = c_[j];
                        for (int i : set(j).elements) v -= shift_[i];
                        search_cost_[j] = std::max(0.0, v);
                        if (x_[j])
                            for (int i : set(j).elements) {
                                ++count_[i];
                                if (locked_[j] == 2) ++fixed_count_[i];
                            }
                        if (!locked_[j])
                            for (int i : set(j).elements) active_rows_[i].push_back(j);
                    }
                }
                init_col_ = end;
            }
            ++init_phase_;
            init_col_ = 0;
        }
        while (init_row_ < n_) {
            const int i = init_row_;
            const auto &row = this->row(i);
            while (init_incidence_ < row.size()) {
                const auto end = budget.batch_end(
                    init_incidence_, static_cast<decltype(init_incidence_)>(row.size()));
                if (end == init_incidence_) return false;
                for (auto at = init_incidence_; at < end; ++at) {
                    const int j = row[at];
                    init_price_ = std::min(init_price_, search_cost_[j] / set(j).elements.size());
                }
                init_incidence_ = end;
            }
            base_weight_[i] =
                std::isfinite(init_price_) ? std::max(init_price_, 1.0 / std::max(1, n_)) : 1;
            down_[i] = 2 * base_weight_[i];
            up_[i] = partition_ ? down_[i] : 0;
            price_[i] = std::isfinite(init_price_) ? init_price_ + shift_[i] : 0;
            if (!partition_) price_[i] = std::max(0.0, price_[i]);
            ++init_row_;
            init_incidence_ = 0;
            init_price_ = std::numeric_limits<double>::infinity();
        }
        init_phase_ = 3;
        index_selected();
        remember();
        return true;
    }
    double add_penalty(int i, int k) const { return k == 0 ? -down_[i] : up_[i]; }
    double remove_penalty(int i, int k) const { return k == 1 ? down_[i] : -up_[i]; }
    bool admissible(int j) const { return !locked_[j] && tabu_[j] <= stats_.iterations; }
    void clear_touched() {
        touched_.clear();
        if (++pair_stamp_ == 0) {
            std::fill(touched_stamp_.begin(), touched_stamp_.end(), 0);
            ++pair_stamp_;
        }
    }
    void reset_scan() {
        scan_ = 0;
        single_ = -1;
        single_value_ = 0;
        origins_ready_ = false;
        origins_.clear();
        origin_ = pair_row_ = pair_col_ = pair_eval_ = 0;
        clear_touched();
        pair_ready_ = false;
        swap_out_ = swap_in_ = -1;
        swap_value_ = 0;
        exchanges_.clear();
        partners_.clear();
        composite_a_ = 0;
        composite_b_ = 1;
    }
    void index_selected() {
        std::fill(selected_bits_.begin(), selected_bits_.end(), 0);
        selected_count_ = 0;
        for (int j = 0; j < m_; ++j)
            if (x_[j]) {
                selected_bits_[j / 64] |= UINT64_C(1) << (j % 64);
                ++selected_count_;
            }
    }
    void flip_selected(int j) {
        x_[j] ^= 1;
        selected_bits_[j / 64] ^= UINT64_C(1) << (j % 64);
        selected_count_ += x_[j] ? 1 : -1;
    }
    template <class Function> void for_selected(Function &&function) const {
        for (std::size_t w = 0; w < selected_bits_.size(); ++w)
            for (std::uint64_t bits = selected_bits_[w]; bits; bits &= bits - 1)
                function(static_cast<int>(w * 64) + detail::trailing_zeros(bits));
    }
    std::vector<int> selected() const {
        std::vector<int> v;
        v.reserve(selected_count_);
        for_selected([&](int j) { v.push_back(j); });
        return v;
    }
    void remember() {
        std::int64_t violation = 0;
        for (int k : count_) violation += !k ? 1 : (partition_ ? k - 1 : 0);
        if (violation > best_.violation) {
            ++stalls_;
            return;
        }
        // Coverage is already exact at a transaction boundary. Rewalking every
        // incidence here would make a dense, redundant seed defeat short deadlines.
        Result r;
        r.selected = selected();
        r.coverage = count_;
        r.violation = violation;
        r.feasible = violation == 0;
        for (int k : count_) r.uncovered += k == 0;
        Cost correction = 0;
        for (int j : r.selected) {
            const Cost y = set(j).cost - correction, sum = r.cost + y;
            correction = (sum - r.cost) - y;
            r.cost = sum;
        }
        if (r.violation < best_.violation || r.cost < best_.cost) {
            best_ = std::move(r);
            ++stats_.improvements;
            stalls_ = 0;
        } else
            ++stalls_;
    }
    void load(const std::vector<int> &v) {
        std::fill(x_.begin(), x_.end(), 0);
        count_ = fixed_count_;
        for (int j : v)
            x_[locked_[j] == 1 && !set(j).elements.empty() ? p_.cheapest_equivalent(j) : j] = 1;
        for (int j = 0; j < m_; ++j)
            if (locked_[j]) x_[j] = locked_[j] == 2;
        for (int j = 0; j < m_; ++j)
            if (x_[j] && locked_[j] != 2)
                for (int i : set(j).elements) ++count_[i];
        index_selected();
        std::fill(tabu_.begin(), tabu_.end(), 0);
        rebuild_ = 0;
        reset_scan();
    }
    void queue(int a, int b = -1) {
        pending_ = {a};
        if (b >= 0) {
            pending_.push_back(b);
            ++stats_.swaps;
        }
        pending_at_ = update_row_ = update_col_ = 0;
        flip_started_ = row_started_ = false;
        reset_scan();
    }
    void queue_composite(std::vector<int> move) {
        pending_ = std::move(move);
        ++stats_.composite_moves;
        pending_at_ = update_row_ = update_col_ = 0;
        flip_started_ = row_started_ = false;
        reset_scan();
    }
    double composite_delta(const int *move, int size) {
        double value = 0;
        for (int k = 0; k < size; ++k) {
            const int j = move[k];
            const int sign = x_[j] ? -1 : 1;
            value += sign * c_[j];
            for (int i : set(j).elements) {
                if (net_change_[i] == 0) changed_rows_.push_back(i);
                net_change_[i] += sign;
            }
        }
        // A row may occur twice in changed_rows_ after its net change returns to
        // zero. Clearing it after evaluation makes subsequent occurrences no-ops.
        for (int i : changed_rows_) {
            const int old = count_[i], now = old + net_change_[i];
            value += down_[i] * ((now == 0) - (old == 0)) +
                     up_[i] * (std::max(0, now - 1) - std::max(0, old - 1));
            net_change_[i] = 0;
        }
        changed_rows_.clear();
        return value;
    }
    void queue_penalties(double factor) {
        reset_scan();
        if (!o_.incremental_penalties) {
            rebuild_ = 0;
            return;
        }
        penalty_row_ = 0;
        penalty_incidence_ = 0;
        penalty_scale_col_ = 0;
        penalty_scale_ = factor;
    }
    bool finish_penalties(detail::Budget &budget) {
        while (penalty_row_ < n_) {
            const int i = penalty_row_;
            if (penalty_down_change_[i] || penalty_up_change_[i]) {
                const auto &candidates = row(i);
                while (penalty_incidence_ < candidates.size()) {
                    const auto end = budget.batch_end(
                        penalty_incidence_,
                        static_cast<decltype(penalty_incidence_)>(candidates.size()));
                    if (end == penalty_incidence_) return false;
                    for (auto at = penalty_incidence_; at < end; ++at) {
                        const int j = candidates[at];
                        delta_[j] += x_[j] ? (count_[i] == 1 ? penalty_down_change_[i]
                                                             : -penalty_up_change_[i])
                                           : (count_[i] == 0 ? -penalty_down_change_[i]
                                                             : penalty_up_change_[i]);
                    }
                    penalty_incidence_ = end;
                }
            }
            ++penalty_row_;
            penalty_incidence_ = 0;
        }
        if (penalty_scale_ != 1)
            while (penalty_scale_col_ < m_) {
                if (budget.poll()) return false;
                const int j = penalty_scale_col_++;
                const double base = x_[j] ? -c_[j] : c_[j];
                delta_[j] = base + penalty_scale_ * (delta_[j] - base);
            }
        penalty_row_ = -1;
        // Limit accumulated roundoff without rebuilding every weighted incidence.
        if (stats_.penalty_updates % 256 == 0) rebuild_ = 0;
        return true;
    }
    bool adapt_penalties() {
        if (!o_.adaptive_penalties || !best_.feasible) return false;
        penalty_down_change_.assign(n_, 0);
        penalty_up_change_.assign(n_, 0);
        double factor = 1;
        double objective = 0, norm = 0;
        for_selected([&](int j) { objective += c_[j]; });
        for (int i = 0; i < n_; ++i) {
            const int missing = count_[i] == 0;
            const int excess = partition_ ? std::max(0, count_[i] - 1) : 0;
            objective += down_[i] * missing + up_[i] * excess;
            norm += missing + static_cast<double>(excess) * excess;
        }
        const double target = static_cast<double>(best_.cost / p_.scale());
        if (objective >= target - 1e-12) {
            // Removing a selected positive-cost column has delta=-c+beta*p.
            // Calibrate beta so roughly 10% of removable columns become attractive.
            std::vector<double> thresholds;
            for_selected([&](int j) {
                if (!locked_[j] && c_[j] > 0) {
                    const double penalty = delta_[j] + c_[j];
                    if (penalty > 0) thresholds.push_back(c_[j] / penalty);
                }
            });
            double beta = 0.9;
            if (!thresholds.empty()) {
                const auto at = thresholds.begin() + (thresholds.size() - 1) / 10;
                std::nth_element(thresholds.begin(), at, thresholds.end(), std::greater<double>());
                beta = std::clamp(0.99 * *at, 0.05, 0.99);
            }
            for (int i = 0; i < n_; ++i) {
                down_[i] *= beta;
                up_[i] *= beta;
            }
            factor = beta;
            ++stats_.penalty_decreases;
        } else if (norm > 0) {
            const double step = (target - objective) / norm;
            for (int i = 0; i < n_; ++i) {
                if (!count_[i]) {
                    down_[i] += step;
                    penalty_down_change_[i] = step;
                }
                if (partition_ && count_[i] > 1) {
                    penalty_up_change_[i] = step * (count_[i] - 1);
                    up_[i] += penalty_up_change_[i];
                }
            }
        } else
            return false;
        ++stats_.penalty_updates;
        ++stats_.iterations;
        ++stalls_;
        queue_penalties(factor);
        return true;
    }
    bool finish_move(detail::Budget &budget) {
        while (pending_at_ < pending_.size()) {
            const int j = pending_[pending_at_];
            if (!flip_started_) {
                flip_selected(j);
                flip_started_ = true;
            }
            const int direction = x_[j] ? 1 : -1;
            const auto &elements = set(j).elements;
            while (update_row_ < elements.size()) {
                const int i = elements[update_row_];
                if (!row_started_) {
                    const int old = count_[i], now = old + direction;
                    row_add_change_ = add_penalty(i, now) - add_penalty(i, old);
                    row_remove_change_ = remove_penalty(i, now) - remove_penalty(i, old);
                    row_free_change_ = (now == 0) - (old == 0);
                    row_gain_change_ = row_free_change_ * (random_greedy_ ? greedy_weight_[i] : 1.0);
                    count_[i] = now;
                    row_started_ = true;
                }
                if (row_add_change_ != 0 || row_remove_change_ != 0) {
                    const auto &row = this->row(i);
                    const double add = row_add_change_, remove = row_remove_change_,
                                 gain = row_gain_change_;
                    const int free = row_free_change_;
                    while (update_col_ < row.size()) {
                        const auto end = budget.batch_end(
                            update_col_, static_cast<decltype(update_col_)>(row.size()));
                        if (end == update_col_) return false;
                        for (auto at = update_col_; at < end; ++at) {
                            const int k = row[at];
                            if (k != j) delta_[k] += x_[k] ? remove : add;
                            greedy_gain_[k] += gain;
                            free_count_[k] += free;
                        }
                        update_col_ = end;
                    }
                }
                ++update_row_;
                update_col_ = 0;
                row_started_ = false;
            }
            double d = x_[j] ? -c_[j] : c_[j];
            for (int i : elements)
                d += x_[j] ? remove_penalty(i, count_[i]) : add_penalty(i, count_[i]);
            delta_[j] = d;
            tabu_[j] = stats_.iterations + 3 + random() % 7;
            ++stats_.flips;
            ++pending_at_;
            update_row_ = 0;
            flip_started_ = false;
        }
        pending_.clear();
        ++stats_.iterations;
        remember();
        return true;
    }
    bool rebuild(detail::Budget &budget) {
        while (rebuild_ < m_) {
            const auto end = budget.batch_end(rebuild_, m_);
            if (end == rebuild_) return false;
            for (auto at = rebuild_; at < end; ++at) {
                int j = at;
                double d = x_[j] ? -c_[j] : c_[j];
                greedy_gain_[j] = 0;
                free_count_[j] = 0;
                for (int i : set(j).elements) {
                    d += x_[j] ? remove_penalty(i, count_[i]) : add_penalty(i, count_[i]);
                    if (!count_[i]) {
                        ++free_count_[j];
                        greedy_gain_[j] += random_greedy_ ? greedy_weight_[i] : 1.0;
                    }
                }
                delta_[j] = d;
            }
            rebuild_ = end;
        }
        return true;
    }
    bool greedy_step(detail::Budget &budget) {
        while (greedy_col_ < m_) {
            const auto end = budget.batch_end(greedy_col_, m_);
            if (end == greedy_col_) return false;
            for (auto at = greedy_col_; at < end; ++at) {
                const int j = at;
                if (x_[j] || locked_[j]) continue;
                const double gain = greedy_gain_[j];
                if (!free_count_[j] || gain <= 0 ||
                    (partition_ && free_count_[j] != static_cast<int>(set(j).elements.size())))
                    continue;
                const double s = (random_greedy_ && !greedy_cost_.empty() ? greedy_cost_[j]
                                  : partition_                            ? search_cost_[j]
                                                                          : c_[j]) /
                                 gain * (random_greedy_ ? greedy_noise_[j] : 1.0);
                if (s < greedy_score_) {
                    greedy_score_ = s;
                    greedy_chosen_ = j;
                }
            }
            greedy_col_ = end;
        }
        const int chosen = greedy_chosen_;
        greedy_col_ = 0;
        greedy_chosen_ = -1;
        greedy_score_ = std::numeric_limits<double>::infinity();
        if (chosen < 0) {
            greedy_ = false;
            return true;
        }
        queue(chosen);
        return true;
    }
    bool scan_move(detail::Budget &budget) {
        while (scan_ < m_) {
            const auto end = budget.batch_end(scan_, m_);
            if (end == scan_) return false;
            double value = single_value_;
            int chosen = single_;
            for (auto at = scan_; at < end; ++at) {
                const int j = at;
                if (admissible(j) && delta_[j] < value - 1e-12) {
                    chosen = j;
                    value = delta_[j];
                }
            }
            single_value_ = value;
            single_ = chosen;
            scan_ = end;
        }
        if (single_ >= 0) {
            queue(single_);
            return true;
        }
        if (!origins_ready_) {
            std::vector<int> all;
            for_selected([&](int j) {
                if (admissible(j)) all.push_back(j);
            });
            if (!all.empty() && o_.swap_origins) {
                std::sort(
                    all.begin(), all.end(), [&](int a, int b) { return delta_[a] < delta_[b]; });
                const int k = std::min<int>(o_.swap_origins, static_cast<int>(all.size()));
                for (int i = 0; i < k; ++i)
                    origins_.push_back(all[(origin_rotation_ + i) % all.size()]);
                origin_rotation_ += k;
            }
            origins_ready_ = true;
        }
        while (origin_ < origins_.size()) {
            int out = origins_[origin_];
            const auto &elements = set(out).elements;
            while (!pair_ready_ && pair_row_ < elements.size()) {
                int i = elements[pair_row_];
                if (count_[i] == 1) {
                    const auto &row = this->row(i);
                    const double change = down_[i] + up_[i];
                    while (pair_col_ < row.size()) {
                        const auto end = budget.batch_end(
                            pair_col_, static_cast<decltype(pair_col_)>(row.size()));
                        if (end == pair_col_) return false;
                        for (auto at = pair_col_; at < end; ++at) {
                            const int in = row[at];
                            if (x_[in] || !admissible(in)) continue;
                            if (touched_stamp_[in] != pair_stamp_) {
                                touched_stamp_[in] = pair_stamp_;
                                touched_.push_back(in);
                                correction_[in] = 0;
                            }
                            correction_[in] += change;
                        }
                        pair_col_ = end;
                    }
                }
                ++pair_row_;
                pair_col_ = 0;
            }
            pair_ready_ = true;
            while (pair_eval_ < touched_.size()) {
                const auto end = budget.batch_end(
                    pair_eval_, static_cast<decltype(pair_eval_)>(touched_.size()));
                if (end == pair_eval_) return false;
                for (auto at = pair_eval_; at < end; ++at) {
                    const int in = touched_[at];
                    const double d = delta_[out] + delta_[in] - correction_[in];
                    if (d < swap_value_ - 1e-12) {
                        swap_value_ = d;
                        swap_out_ = out;
                        swap_in_ = in;
                    }
                    if (o_.composite_moves && best_.feasible) {
                        const std::pair<double, int> partner(d, in);
                        // Keep the same two lexicographically smallest entries.
                        // Most candidates cannot enter this already sorted pair.
                        if (partners_.size() < 2) {
                            partners_.push_back(partner);
                            if (partners_.size() == 2 && partners_[1] < partners_[0])
                                std::swap(partners_[0], partners_[1]);
                        } else if (partner < partners_[1]) {
                            partners_[1] = partner;
                            if (partners_[1] < partners_[0]) std::swap(partners_[0], partners_[1]);
                        }
                    }
                }
                pair_eval_ = end;
            }
            if (swap_out_ >= 0) {
                queue(swap_out_, swap_in_);
                return true;
            }
            for (auto partner : partners_) exchanges_.emplace_back(out, partner.second);
            partners_.clear();
            clear_touched();
            ++origin_;
            pair_row_ = pair_col_ = pair_eval_ = 0;
            pair_ready_ = false;
        }
        while (composite_a_ < exchanges_.size()) {
            while (composite_b_ < exchanges_.size()) {
                if (budget.expired()) return false;
                const auto a = exchanges_[composite_a_], b = exchanges_[composite_b_++];
                int move[4] = {a.first, a.second}, size = 2;
                if (b.first != a.first) move[size++] = b.first;
                if (b.second != a.second) move[size++] = b.second;
                if (size >= 3 && composite_delta(move, size) < -1e-12) {
                    queue_composite(std::vector<int>(move, move + size));
                    return true;
                }
            }
            ++composite_a_;
            composite_b_ = composite_a_ + 1;
        }
        if (adapt_penalties()) return true;
        // Strategic oscillation: raise prices on violated rows, or destroy a feasible solution.
        bool violated = false;
        penalty_down_change_.assign(n_, 0);
        penalty_up_change_.assign(n_, 0);
        for (int i = 0; i < n_; ++i) {
            if (!count_[i]) {
                penalty_down_change_[i] = base_weight_[i] * (0.5 + uniform());
                down_[i] += penalty_down_change_[i];
                violated = true;
            }
            if (partition_ && count_[i] > 1) {
                penalty_up_change_[i] = base_weight_[i] * (0.5 + uniform());
                up_[i] += penalty_up_change_[i];
                violated = true;
            }
        }
        if (violated) {
            // Keep penalties bounded without using them as dual certificates.
            double max_weight = 0;
            for (int i = 0; i < n_; ++i) max_weight = std::max({max_weight, down_[i], up_[i]});
            if (max_weight > 1e6)
                for (int i = 0; i < n_; ++i) {
                    down_[i] *= 0.01;
                    up_[i] *= 0.01;
                }
            ++stats_.penalty_updates;
            ++stats_.iterations;
            ++stalls_;
            queue_penalties(max_weight > 1e6 ? 0.01 : 1);
            // Tabu is temporary; do not freeze a violated state behind it.
        } else {
            int out = -1;
            double value = std::numeric_limits<double>::infinity();
            for_selected([&](int j) {
                if (!locked_[j]) {
                    const double v = (delta_[j] + 1e-8) * (0.5 + uniform());
                    if (v < value) {
                        value = v;
                        out = j;
                    }
                }
            });
            if (out >= 0)
                queue(out);
            else {
                ++stats_.iterations;
                reset_scan();
            }
        }
        return true;
    }
    void start_dual() {
        dual_active_ = true;
        dual_col_ = 0;
        dual_value_ = dual_abs_ = 0;
        std::fill(dual_count_.begin(), dual_count_.end(), 0);
        for (double u : price_) {
            dual_value_ += u;
            dual_abs_ += std::abs(u);
        }
    }
    bool dual_pass(detail::Budget &budget) {
        while (dual_col_ < m_) {
            const auto end = budget.batch_end(dual_col_, m_);
            if (end == dual_col_) return false;
            for (auto at = dual_col_; at < end; ++at) {
                int j = at;
                double reduced = c_[j];
                for (int i : set(j).elements) reduced -= price_[i];
                dual_abs_ += std::abs(c_[j]);
                if (reduced < 0) {
                    dual_value_ += reduced;
                    for (int i : set(j).elements) ++dual_count_[i];
                }
            }
            dual_col_ = end;
        }
        // A conservative floating-point allowance, not an exact rational certificate.
        const Cost value =
            (static_cast<Cost>(dual_value_) - 1e-10L * (1 + dual_abs_) * (n_ + 1)) * p_.scale();
        if (value > bound_) {
            bound_ = value;
            best_price_ = price_;
        }
        double norm = 0;
        for (int k : dual_count_) norm += (1.0 - k) * (1.0 - k);
        if (norm > 0) {
            const double target = best_.feasible
                                      ? static_cast<double>(best_.cost / p_.scale())
                                      : dual_value_ + std::max(1.0, std::abs(dual_value_) * 0.2);
            const double step = (1.4 / (1 + 0.05 * stats_.dual_iterations)) *
                                std::max(0.0, target - dual_value_) / norm;
            for (int i = 0; i < n_; ++i)
                price_[i] = std::clamp(
                    price_[i] + step * (1 - dual_count_[i]), partition_ ? -1e6 : 0.0, 1e6);
        }
        ++stats_.dual_iterations;
        dual_active_ = false;
        if (o_.priced_restarts && !partition_ && !o_.fast_cover &&
            stats_.dual_iterations == static_cast<std::uint64_t>(o_.dual_iterations))
            restart_greedy(true);
        return true;
    }
    bool contains(const std::vector<std::uint64_t> &remaining, int j) const {
        const auto &b = set(j).bits;
        for (std::size_t w = 0; w < b.size(); ++w)
            if (b[w] & ~remaining[w]) return false;
        return true;
    }
    int gain(const std::vector<std::uint64_t> &remaining, int j) const {
        int k = 0;
        const auto &b = set(j).bits;
        for (std::size_t w = 0; w < b.size(); ++w) k += detail::popcount(b[w] & remaining[w]);
        return k;
    }
    bool empty(const std::vector<std::uint64_t> &b) const {
        for (auto w : b)
            if (w) return false;
        return true;
    }
    double rank(int j) const {
        double reduced = c_[j];
        for (int i : set(j).elements) reduced -= best_price_[i];
        return reduced + 0.15 * search_cost_[j] / std::max<std::size_t>(1, set(j).elements.size());
    }
    void restart_greedy(bool priced) {
        load({});
        greedy_ = true;
        random_greedy_ = true;
        greedy_noise_.resize(m_);
        greedy_weight_.assign(n_, 1);
        greedy_cost_.resize(m_);
        greedy_alpha_ = 0.2 + 0.8 * uniform();
        greedy_priced_ = priced;
        greedy_setup_ = 0;
        for (int i = 0; i < n_; ++i) {
            down_[i] = 2 * base_weight_[i];
            up_[i] = partition_ ? down_[i] : 0;
        }
    }
    bool prepare_greedy(detail::Budget &budget) {
        while (greedy_setup_ < m_) {
            const auto end = budget.batch_end(greedy_setup_, m_);
            if (end == greedy_setup_) return false;
            for (auto at = greedy_setup_; at < end; ++at) {
                const int j = at;
                greedy_noise_[j] = 0.7 + uniform() * 0.6;
                double reduced = c_[j];
                if (greedy_priced_)
                    for (int i : set(j).elements) reduced -= best_price_[i];
                greedy_cost_[j] = greedy_priced_
                                      ? std::max(0.0, reduced) + greedy_alpha_ * search_cost_[j]
                                      : search_cost_[j];
            }
            greedy_setup_ = end;
        }
        greedy_setup_ = -1;
        return true;
    }
    bool start_repair(detail::Budget &budget) {
        repair_ = Repair{};
        auto &r = repair_;
        r.active = true;
        r.fixed_flag.assign(m_, 0);
        r.preferred.assign(m_, 0);
        r.free.assign((static_cast<std::size_t>(n_) + 63) / 64, 0);
        std::vector<int> source = best_.selected;
        // For an infeasible partition, keep a randomized disjoint subset as the anchor.
        if (!best_.feasible)
            for (std::size_t k = source.size(); k > 1; --k) {
                if (budget.poll()) return false;
                std::swap(source[k - 1], source[random() % k]);
            }
        std::vector<int> covered(n_, 0), removable;
        for (int j : source) {
            if (budget.poll()) return false;
            bool ok = true;
            if (partition_)
                for (int i : set(j).elements)
                    if (covered[i]) {
                        ok = false;
                        break;
                    }
            if (!ok) continue;
            r.fixed.push_back(j);
            r.fixed_flag[j] = 1;
            for (int i : set(j).elements) ++covered[i];
            if (!locked_[j]) removable.push_back(j);
        }
        // Remove the conflicts of a promising entering column. In an exact partition
        // selected columns never overlap one another, so selected-column overlap alone
        // would provide no useful clustering (e.g. expensive singleton anchors).
        int pivot = -1;
        if (partition_ && n_ > 48 && !removable.empty() && (stats_.repairs % 4 != 3)) {
            double best_rank = -std::numeric_limits<double>::infinity();
            for (int t = 0; t < 64; ++t) {
                const int anchor = removable[random() % removable.size()];
                if (set(anchor).elements.empty()) continue;
                const auto &elements = set(anchor).elements;
                const int row = elements[random() % elements.size()];
                const auto &candidates = this->row(row);
                if (candidates.empty()) continue;
                const int j = candidates[random() % candidates.size()];
                if (r.fixed_flag[j] || locked_[j]) continue;
                double replaced = 0;
                int conflicts = 0;
                for (int k : removable) {
                    if (budget.poll()) return false;
                    bool overlap = false;
                    for (int i : set(k).elements)
                        if (set(j).bits[i / 64] & (UINT64_C(1) << (i % 64))) {
                            overlap = true;
                            break;
                        }
                    if (overlap) {
                        replaced += c_[k];
                        ++conflicts;
                    }
                }
                const double score = (replaced - c_[j]) * (0.5 + uniform());
                if (conflicts && conflicts <= o_.repair_max_remove && score > best_rank) {
                    best_rank = score;
                    pivot = j;
                }
            }
        }
        if (pivot >= 0) {
            r.preferred[pivot] = 1;
            for (std::size_t k = 0; k < removable.size();) {
                if (budget.poll()) return false;
                int j = removable[k];
                bool overlap = false;
                for (int i : set(j).elements)
                    if (set(pivot).bits[i / 64] & (UINT64_C(1) << (i % 64))) {
                        overlap = true;
                        break;
                    }
                if (!overlap) {
                    ++k;
                    continue;
                }
                r.fixed_flag[j] = 0;
                r.removed.push_back(j);
                r.preferred[j] = 1;
                for (int i : set(j).elements) --covered[i];
                removable[k] = removable.back();
                removable.pop_back();
            }
        }
        // Some full and random neighborhoods retain diversity beyond one pivot.
        const bool full_repair = n_ <= (partition_ ? 48 : 24) ||
                                 (partition_ && !best_.feasible && stats_.repairs % 4 == 0);
        const int amount =
            full_repair ? static_cast<int>(removable.size())
            : pivot >= 0
                ? 0
                : std::min<int>(
                      static_cast<int>(removable.size()),
                      2 + static_cast<int>(
                              random() %
                              static_cast<std::uint64_t>(std::max(
                                  1, std::min(o_.repair_max_remove,
                                              2 + static_cast<int>(stats_.repairs % 11))))));
        std::vector<unsigned char> near(n_, 0);
        for (int t = 0; t < amount; ++t) {
            std::size_t pos = random() % removable.size();
            if (t && (random() & 3U)) {
                int best_overlap = -1;
                for (std::size_t k = 0; k < removable.size(); ++k) {
                    if (budget.poll()) return false;
                    int overlap = 0;
                    for (int i : set(removable[k]).elements) overlap += near[i];
                    if (overlap > best_overlap) {
                        best_overlap = overlap;
                        pos = k;
                    }
                }
            }
            int j = removable[pos];
            removable[pos] = removable.back();
            removable.pop_back();
            r.fixed_flag[j] = 0;
            r.removed.push_back(j);
            r.preferred[j] = 1;
            for (int i : set(j).elements) {
                --covered[i];
                near[i] = 1;
                // Incidence itself defines overlap; no global column-pair graph is stored.
            }
        }
        std::vector<int> fixed;
        for (int j : r.fixed)
            if (r.fixed_flag[j]) {
                fixed.push_back(j);
                r.fixed_cost += set(j).cost;
            }
        r.fixed = std::move(fixed);
        for (int i = 0; i < n_; ++i)
            if (!covered[i]) r.free[i / 64] |= UINT64_C(1) << (i % 64);
        if (best_.feasible) {
            r.limit = 0;
            for (int j : r.removed) r.limit += search_cost_[j];
        }
        r.rows.resize(n_);
        r.lower_weight.assign(n_, 0);
        ++stats_.repairs;
        return true;
    }
    bool prepare_repair(detail::Budget &budget) {
        auto &r = repair_;
        while (r.prepare_row < n_) {
            if (budget.expired()) return false;
            const int i = r.prepare_row;
            if (!(r.free[i / 64] & (UINT64_C(1) << (i % 64)))) {
                ++r.prepare_row;
                continue;
            }
            // Preserve a partially priced row across arbitrarily short run() slices.
            auto &ranked = r.ranked;
            while (r.prepare_col < row(i).size()) {
                const auto end = budget.batch_end(
                    r.prepare_col, static_cast<decltype(r.prepare_col)>(row(i).size()));
                if (end == r.prepare_col) return false;
                for (auto at = r.prepare_col; at < end; ++at) {
                    const int j = row(i)[at];
                    if (r.fixed_flag[j] || locked_[j] || (partition_ && !contains(r.free, j)))
                        continue;
                    r.min_share =
                        std::min(r.min_share, search_cost_[j] / std::max(1, gain(r.free, j)));
                    ranked.emplace_back(r.preferred[j] ? -1e100 : rank(j), j);
                }
                r.prepare_col = end;
            }
            const std::size_t cap =
                o_.repair_columns_per_row == 0
                    ? ranked.size()
                    : std::min<std::size_t>(ranked.size(), o_.repair_columns_per_row);
            if (cap < ranked.size())
                std::nth_element(ranked.begin(), ranked.begin() + cap, ranked.end());
            auto &row = r.rows[i];
            row.clear();
            for (std::size_t k = 0; k < cap; ++k) row.push_back(ranked[k].second);
            for (int j : r.removed)
                if (std::binary_search(set(j).elements.begin(), set(j).elements.end(), i) &&
                    std::find(row.begin(), row.end(), j) == row.end())
                    row.push_back(j);
            r.lower_weight[i] = std::isfinite(r.min_share) ? r.min_share : 0;
            ++r.prepare_row;
            r.prepare_col = 0;
            ranked.clear();
            r.min_share = std::numeric_limits<double>::infinity();
        }
        Frame root;
        root.remaining = r.free;
        for (double w : r.lower_weight) root.lower += w;
        r.stack.push_back(std::move(root));
        r.prepared = true;
        return true;
    }
    bool build_frame(Frame &f, detail::Budget &budget) {
        auto &r = repair_;
        while (f.scan_row < n_) {
            int i = f.scan_row;
            if (!(f.remaining[i / 64] & (UINT64_C(1) << (i % 64)))) {
                ++f.scan_row;
                continue;
            }
            const auto &row = r.rows[i];
            while (f.scan_col < row.size()) {
                if (budget.poll()) return false;
                const int j = row[f.scan_col++];
                if ((!partition_ || contains(f.remaining, j)) &&
                    f.cost + search_cost_[j] < r.limit - 1e-12)
                    f.trial.push_back(j);
                if (!f.candidates.empty() && f.trial.size() >= f.candidates.size()) break;
            }
            if (f.trial.empty()) {
                f.candidates.clear();
                f.ready = true;
                return true;
            }
            if (f.candidates.empty() || f.trial.size() < f.candidates.size())
                f.candidates = f.trial;
            f.trial.clear();
            f.scan_col = 0;
            ++f.scan_row;
            if (f.candidates.size() == 1) break;
        }
        std::sort(f.candidates.begin(), f.candidates.end(), [&](int a, int b) {
            const double ra = search_cost_[a] / std::max(1, gain(f.remaining, a));
            const double rb = search_cost_[b] / std::max(1, gain(f.remaining, b));
            return ra != rb ? ra < rb : a < b;
        });
        f.ready = true;
        ++r.nodes;
        ++stats_.repair_nodes;
        return true;
    }
    bool repair_search(detail::Budget &budget) {
        auto &r = repair_;
        if (!r.prepared && !prepare_repair(budget)) return false;
        while (!r.stack.empty() && r.nodes < static_cast<std::uint64_t>(o_.repair_nodes)) {
            if (budget.expired()) return false;
            Frame &f = r.stack.back();
            if (empty(f.remaining)) {
                if (f.cost < r.limit - 1e-12) {
                    r.limit = f.cost;
                    r.answer = r.fixed;
                    for (const Frame &a : r.stack)
                        if (a.incoming >= 0) r.answer.push_back(a.incoming);
                    const Result candidate = p_.evaluate(r.answer);
                    if (candidate.feasible && (!best_.feasible || candidate.cost < best_.cost)) {
                        best_ = candidate;
                        ++stats_.improvements;
                        stalls_ = 0;
                    }
                }
                r.stack.pop_back();
                continue;
            }
            if (f.cost + f.lower >= r.limit - 1e-12) {
                r.stack.pop_back();
                continue;
            }
            if (!f.ready && !build_frame(f, budget)) return false;
            if (f.next == f.candidates.size()) {
                r.stack.pop_back();
                continue;
            }
            const int j = f.candidates[f.next++];
            if (f.cost + search_cost_[j] >= r.limit - 1e-12) continue;
            Frame next;
            next.remaining = f.remaining;
            next.cost = f.cost + search_cost_[j];
            next.incoming = j;
            next.lower = f.lower;
            for (int i : set(j).elements)
                if (f.remaining[i / 64] & (UINT64_C(1) << (i % 64)))
                    next.lower -= r.lower_weight[i];
            const auto &bits = set(j).bits;
            for (std::size_t w = 0; w < bits.size(); ++w) next.remaining[w] &= ~bits[w];
            r.stack.push_back(std::move(next));
        }
        r.active = false;
        if (stats_.repairs % 16 == 0 && stalls_ > static_cast<std::uint64_t>(o_.repair_interval)) {
            restart_greedy(o_.priced_restarts && !partition_);
        } else if (!r.answer.empty())
            load(r.answer);
        else if (best_.feasible && stats_.repairs % 8 == 0 &&
                 stalls_ > 4 * static_cast<std::uint64_t>(std::max(1, o_.repair_interval)))
            load(best_.selected);
        return true;
    }

public:
    explicit Solver(const Model &model, Options options = {})
        : p_(model), o_(std::move(options)), n_(p_.elements()), m_(p_.sets()),
          partition_(p_.mode() == Mode::Partition), revision_(p_.revision()), rng_(o_.seed),
          x_(m_, 0), locked_(m_, 0), selected_bits_((static_cast<std::size_t>(m_) + 63) / 64),
          count_(n_, 0), fixed_count_(n_, 0), active_rows_(n_), c_(m_), down_(n_), up_(n_, 0),
          delta_(m_), base_weight_(n_, 1), price_(n_, 0), best_price_(n_, 0), greedy_gain_(m_, 0),
          free_count_(m_, 0), tabu_(m_, 0), correction_(m_, 0), touched_stamp_(m_, 0),
          net_change_(n_, 0), dual_count_(n_, 0), shift_(n_, 0), search_cost_(m_) {
        if (!std::isfinite(o_.safety_margin) || o_.safety_margin < 0 || o_.swap_origins < 0 ||
            o_.repair_interval < 0 || o_.repair_nodes < 0 || o_.repair_columns_per_row < 0 ||
            o_.repair_max_remove < 1 || o_.dual_iterations < 0 || o_.ejection_moves < 0 ||
            std::isnan(o_.time_limit) || o_.time_limit < 0)
            throw std::invalid_argument("set_cover_partition: invalid options");
        best_ = p_.evaluate(o_.initial_sets);
        for (int j : o_.initial_sets) x_[j] = 1;
        for (int i = 0; i < n_; ++i)
            if (p_.row(i).empty()) impossible_ = true;
    }
    Solver(Model &&, Options = {}) = delete; // Avoid a dangling model reference.

    Result run(double seconds) { return run_from(seconds, Clock::now()); }
    Result run() { return run(o_.time_limit); }
    // Allows solve() to include Solver initialization in the same wall-time budget.
    Result run_from(double seconds, Clock::time_point start) {
        if (p_.revision() != revision_)
            throw std::logic_error("set_cover_partition: model changed during session");
        detail::Budget budget(start, seconds, o_.deadline, o_.safety_margin);
        const std::uint64_t before = stats_.iterations;
        StopReason reason = StopReason::TimeLimit;
        if (impossible_)
            reason = StopReason::NoCoveringSet;
        else
            while (!budget.expired()) {
                if (stats_.iterations - before >= o_.iteration_limit) {
                    reason = StopReason::IterationLimit;
                    break;
                }
                if (init_phase_ < 3) {
                    if (!initialize(budget)) break;
                    continue;
                }
                if (!n_) {
                    reason = StopReason::Exhausted;
                    break;
                }
                if (!partition_ && o_.fast_cover) {
                    if (dual_active_) {
                        if (!dual_pass(budget)) break;
                        continue;
                    }
                    if ((cover_.phase == 4 || (cover_.phase == 6 && cover_.row < 0)) &&
                        stats_.dual_iterations < static_cast<std::uint64_t>(o_.dual_iterations) &&
                        stats_.dual_iterations * 2048 <= stats_.cover_cycles) {
                        start_dual();
                        continue;
                    }
                    if (!cover_step(budget)) break;
                    continue;
                }
                if (!pending_.empty()) {
                    if (!finish_move(budget)) break;
                    continue;
                }
                if (penalty_row_ >= 0) {
                    if (!finish_penalties(budget)) break;
                    continue;
                }
                if (greedy_setup_ >= 0) {
                    if (!prepare_greedy(budget)) break;
                    continue;
                }
                if (rebuild_ < m_) {
                    if (!rebuild(budget)) break;
                    continue;
                }
                if (ejection_.active) {
                    if (!ejection_step(budget)) break;
                    continue;
                }
                if (repair_.active) {
                    if (!repair_search(budget)) break;
                    continue;
                }
                if (dual_active_) {
                    if (!dual_pass(budget)) break;
                    continue;
                }
                if (greedy_) {
                    if (!greedy_step(budget)) break;
                    continue;
                }
                if (partition_ && best_.feasible && o_.ejection_moves &&
                    stats_.iterations >= ejection_next_) {
                    start_ejection();
                    continue;
                }
                // First reach feasibility with inexpensive neighborhoods. Once feasible,
                // amortize longer pricing over more local moves instead of monopolizing
                // the short budget with successive full-column passes.
                const int dual_limit =
                    best_.feasible ? o_.dual_iterations : std::min(32, o_.dual_iterations);
                const std::uint64_t dual_spacing = best_.feasible ? 16 : 4;
                if (stats_.dual_iterations < static_cast<std::uint64_t>(dual_limit) &&
                    stats_.dual_iterations * dual_spacing <= stats_.iterations) {
                    start_dual();
                    continue;
                }
                if (o_.repair_interval && o_.repair_nodes &&
                    best_.selected.size() <= static_cast<std::size_t>(n_) + fixed_on_count_ &&
                    (!partition_ || best_.feasible || n_ <= 48 || stats_.iterations >= 8192) &&
                    stats_.iterations >=
                        last_repair_ + static_cast<std::uint64_t>(o_.repair_interval)) {
                    last_repair_ = stats_.iterations;
                    if (!start_repair(budget)) {
                        repair_.active = false;
                        break;
                    }
                    continue;
                }
                if (!scan_move(budget)) break;
            }
        Result result = best_;
        result.lower_bound = bound_;
        result.reason = reason;
        result.stats = stats_;
        result.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }
};

inline Result solve(const Model &model, Options options = {}) {
    const auto start = Clock::now();
    const double seconds = options.time_limit;
    Solver solver(model, std::move(options));
    return solver.run_from(seconds, start);
}
} // namespace set_cover_partition
#endif
