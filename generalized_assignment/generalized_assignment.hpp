#ifndef GENERALIZED_ASSIGNMENT_HPP
#define GENERALIZED_ASSIGNMENT_HPP

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

namespace cpopt {
namespace gap {

enum class StopReason { TimeLimit, IterationLimit, NoMoves, ImpossibleJob };

struct Options {
    double time_limit = 1.0; // Wall seconds, per solve() call.
    std::uint64_t iteration_limit = std::numeric_limits<std::uint64_t>::max();
    double safety_margin = 0.0002; // Also reserve 1% (at most 2 ms).
    int swap_candidates = 12;      // Per job; all pairs for n <= 160.
    int chain_trials = 64;
    int chain_depth = 12;    // 0 disables chains; maximum 16.
    int repair_nodes = 8192; // Subset search work budget; 0 disables.
    bool path_relinking = true;
    bool restarts = true;
    bool initial_prices = true;     // UB-free dual initialization before construction.
    bool incremental_shifts = true; // Strict tabu row cache for n >= 256.
    int double_shift_trials = 64;   // Sampled two-job paths every eight iterations.
};

struct Tolerance {
    long double absolute = 1e-10L;
    long double relative = 1e-12L;
};

template <class Cost = std::int64_t, class Resource = Cost> class Solver {
    static_assert(std::is_arithmetic<Cost>::value && !std::is_same<Cost, bool>::value,
                  "Cost must be a real arithmetic type");
    static_assert(std::is_arithmetic<Resource>::value && !std::is_same<Resource, bool>::value,
                  "Resource must be a real arithmetic type");
    static_assert(!std::is_integral<Cost>::value ||
                      (std::is_signed<Cost>::value && sizeof(Cost) <= 8),
                  "Use signed integers up to int64_t");
    static_assert(!std::is_integral<Resource>::value ||
                      (std::is_signed<Resource>::value && sizeof(Resource) <= 8),
                  "Use signed integers up to int64_t");

public:
    using cost_type =
        typename std::conditional<std::is_integral<Cost>::value, std::int64_t, long double>::type;
    using resource_type =
        typename std::conditional<std::is_integral<Resource>::value, std::int64_t, long double>::type;
    struct Result {
        std::vector<int> assignment; // job -> agent; -1 means unassigned.
        std::vector<resource_type> loads;
        cost_type objective = 0; // Sum over assigned jobs, even if incomplete.
        int unassigned = 0;
        long double violation = 0; // Sum of capacity excess beyond tolerance.
        bool feasible = false;
        StopReason reason = StopReason::TimeLimit;
        double elapsed = 0;           // This call, including all search work.
        std::uint64_t iterations = 0; // Cumulative committed search moves.
        std::uint64_t moves_evaluated = 0;
        std::uint64_t restarts = 0;
    };

    // Matrices are job-major [job][agent]. Empty allowed means every pair is allowed.
    // Construction uses O(n*m) storage and O(m*n*log(n)) preprocessing, outside solve()'s budget.
    Solver(const std::vector<std::vector<Cost>> &costs,
           const std::vector<std::vector<Resource>> &resources,
           const std::vector<Resource> &capacities,
           const std::vector<std::vector<unsigned char>> &allowed = {}, std::uint64_t seed = 1,
           Tolerance tolerance = {})
        : tolerance_(tolerance), rng_(seed) {
        if (!std::isfinite(tolerance.absolute) || !std::isfinite(tolerance.relative) ||
            tolerance.absolute < 0 || tolerance.relative < 0)
            throw std::invalid_argument("gap: invalid tolerance");
        if (costs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            capacities.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("gap: dimensions exceed int range");
        n_ = static_cast<int>(costs.size());
        m_ = static_cast<int>(capacities.size());
        if (resources.size() != costs.size() || (!allowed.empty() && allowed.size() != costs.size()))
            throw std::invalid_argument("gap: inconsistent matrix rows");
        if (m_ && costs.size() > std::numeric_limits<std::size_t>::max() / capacities.size())
            throw std::length_error("gap: matrix size overflow");
        const std::size_t size = costs.size() * capacities.size();
        c_.resize(size);
        a_.resize(size);
        allowed_.resize(size);
        usable_.resize(size);
        cap_.resize(m_);
        tolerance_by_agent_.resize(m_);
        scale_.resize(m_);
        penalty_.resize(m_);
        price_.assign(m_, 0);
        price_load_.resize(m_);
        std::vector<long double> resource_sum(m_, 0);
        std::vector<std::uint64_t> integer_resource_sum(m_, 0);
        const auto integer_limit =
            static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() / 8);
        const long double real_limit = std::numeric_limits<long double>::max() / 1024;
        for (int i = 0; i < m_; ++i) {
            if (!std::isfinite(static_cast<long double>(capacities[i])) || capacities[i] < 0)
                throw std::invalid_argument("gap: capacity must be finite and nonnegative");
            if constexpr (std::is_integral<Resource>::value) {
                if (static_cast<std::uint64_t>(capacities[i]) > integer_limit)
                    throw std::overflow_error("gap: integer capacity exceeds int64 safety range");
                tolerance_by_agent_[i] = 0;
            } else {
                if (capacities[i] > real_limit)
                    throw std::overflow_error("gap: capacity too large");
                tolerance_by_agent_[i] = tolerance.absolute + tolerance.relative * capacities[i];
                if (!std::isfinite(tolerance_by_agent_[i]))
                    throw std::overflow_error("gap: tolerance overflow");
            }
            cap_[i] = capacities[i];
        }
        long double cost_range_sum = 0, real_cost_sum = 0;
        std::uint64_t integer_cost_sum = 0;
        order_.resize(n_);
        alternatives_.assign(n_, 0);
        min_cost_.resize(n_);
        std::vector<long double> regret(n_, 0);
        for (int j = 0; j < n_; ++j) {
            if (costs[j].size() != capacities.size() || resources[j].size() != capacities.size() ||
                (!allowed.empty() && allowed[j].size() != capacities.size()))
                throw std::invalid_argument("gap: inconsistent matrix columns");
            long double minc = infinity(), second = infinity(), maxc = -infinity(), maxabs = 0;
            std::uint64_t maxabs_integer = 0;
            for (int i = 0; i < m_; ++i) {
                const std::size_t p = index(j, i);
                const auto cv = costs[j][i];
                const auto av = resources[j][i];
                if (!std::isfinite(static_cast<long double>(cv)) ||
                    !std::isfinite(static_cast<long double>(av)) || av < 0)
                    throw std::invalid_argument(
                        "gap: finite costs and nonnegative finite resources required");
                if constexpr (std::is_integral<Cost>::value) {
                    const std::int64_t v = cv;
                    const std::uint64_t ab = v < 0 ? static_cast<std::uint64_t>(-(v + 1)) + 1
                                                   : static_cast<std::uint64_t>(v);
                    maxabs_integer = std::max(maxabs_integer, ab);
                } else
                    maxabs = std::max(maxabs, std::fabs(static_cast<long double>(cv)));
                if constexpr (std::is_integral<Resource>::value) {
                    const auto v = static_cast<std::uint64_t>(av);
                    if (v > integer_limit - integer_resource_sum[i])
                        throw std::overflow_error(
                            "gap: resource column sum exceeds int64 safety range");
                    integer_resource_sum[i] += v;
                }
                resource_sum[i] += static_cast<long double>(av);
                if (!std::isfinite(resource_sum[i]) || resource_sum[i] > real_limit)
                    throw std::overflow_error("gap: resource sum too large");
                c_[p] = cv;
                a_[p] = av;
                allowed_[p] = allowed.empty() || allowed[j][i] != 0;
                usable_[p] = allowed_[p] && excess(i, a_[p]) == 0;
                if (usable_[p]) {
                    ++alternatives_[j];
                    const long double v = cv;
                    if (v < minc) {
                        second = minc;
                        minc = v;
                        min_cost_[j] = cv;
                    } else if (v < second)
                        second = v;
                    maxc = std::max(maxc, v);
                }
            }
            if constexpr (std::is_integral<Cost>::value) {
                if (maxabs_integer > integer_limit - integer_cost_sum)
                    throw std::overflow_error(
                        "gap: cost magnitude sum exceeds int64 safety range");
                integer_cost_sum += maxabs_integer;
            } else {
                real_cost_sum += maxabs;
                if (!std::isfinite(real_cost_sum) || real_cost_sum > real_limit)
                    throw std::overflow_error("gap: cost magnitude sum too large");
            }
            if (!alternatives_[j])
                impossible_ = true;
            else {
                cost_range_sum += maxc - minc;
                regret[j] = std::isfinite(second) ? second - minc : infinity();
            }
            order_[j] = j;
        }
        cost_scale_ = cost_range_sum > 0 ? std::max(std::numeric_limits<long double>::min(),
                                                    cost_range_sum / std::max(1, n_))
                                         : 1;
        for (int i = 0; i < m_; ++i) {
            scale_[i] =
                std::max(static_cast<long double>(cap_[i]), resource_sum[i] / std::max(1, n_));
            if (!(scale_[i] > 0)) scale_[i] = 1;
            const long double avg = resource_sum[i] / std::max(1, n_);
            penalty_[i] = avg > 0 ? std::clamp(2 * scale_[i] / avg, 0.1L, 1e8L) : 1;
        }
        normalized_cost_.resize(size);
        normalized_resource_.resize(size);
        for (int j = 0; j < n_; ++j)
            for (int i = 0; i < m_; ++i) {
                const auto p = index(j, i);
                normalized_cost_[p] = static_cast<long double>(c_[p] - min_cost_[j]) / cost_scale_;
                normalized_resource_[p] = a_[p] / scale_[i];
            }
        std::stable_sort(order_.begin(), order_.end(), [&](int x, int y) {
            if (alternatives_[x] != alternatives_[y]) return alternatives_[x] < alternatives_[y];
            return regret[x] > regret[y];
        });
        x_.assign(n_, -1);
        load_.assign(m_, 0);
        tabu_.assign(size, 0);
        refresh_load_.resize(m_);
        load_excess_.assign(m_, 0);
        refresh_excess_.resize(m_);
        shift_dirty_.assign(m_, 1);
        shift_changed_.resize(m_);
        escape_.resize(n_);
        coordinate_agent_.resize(n_);
        coordinate_value_.resize(n_);
        breakpoints_.resize(n_);
        shift_target_.assign(n_, -1);
        shift_value_.resize(n_);
        shift_expiry_.assign(n_, 0);
        head_.assign(m_, -1);
        next_.assign(n_, -1);
        prev_.assign(n_, -1);
        scratch_.assignment.resize(n_);
        scratch_.loads.resize(m_);
        best_.assignment.assign(n_, -1);
        best_.loads.assign(m_, 0);
        best_.unassigned = n_;
        best_.feasible = n_ == 0;
        unassigned_ = n_;
        for (auto &elite : elite_) elite.resize(n_);
        elite_staging_.resize(n_);
        for (auto &records : pair_records_) records.resize(2048);
        pair_coordinates_.resize(2048);
        pair_fenwick_.resize(2049);
        chain_candidates_.resize(n_);
        chain_count_.resize(n_);
        chain_available_.resize(n_);
        chain_eligible_.resize(n_);
        dual_job_.resize(n_);
        dual_job_staging_.resize(n_);
        dual_best_price_.resize(m_);
        chain_query_head_.resize(n_);
        chain_query_next_.resize(n_);
        chain_data_.resize(size);
        std::vector<int> resource_order(n_);
        for (int i = 0; i < m_; ++i) {
            auto first = resource_order.begin();
            std::iota(first, first + n_, 0);
            std::sort(first, first + n_,
                      [&](int j, int k) { return a_[index(j, i)] < a_[index(k, i)]; });
            for (int k = 0; k < n_; ++k) {
                const auto p = index(first[k], i);
                chain_data_[static_cast<std::size_t>(i) * n_ + k] = {
                    a_[p], 0.0L - normalized_cost_[p], first[k], usable_[p]};
            }
        }
    }

    // Warm starts may be partial (-1) or capacity-infeasible. They are hints, not fixed jobs.
    // Validation and rebuilding are outside solve()'s time budget.
    void set_start(const std::vector<int> &assignment) {
        if (assignment.size() != static_cast<std::size_t>(n_))
            throw std::invalid_argument("gap: wrong start size");
        for (int j = 0; j < n_; ++j)
            if (assignment[j] < -1 || assignment[j] >= m_ ||
                (assignment[j] >= 0 && !allowed_[index(j, assignment[j])]))
                throw std::invalid_argument("gap: start contains invalid or forbidden agent");
        x_ = assignment;
        rebuild_unlimited();
        std::fill(tabu_.begin(), tabu_.end(), 0);
        best_.assignment = x_;
        best_.loads = load_;
        best_.objective = objective_;
        best_.unassigned = unassigned_;
        best_.violation = total_excess();
        best_.feasible = !unassigned_ && best_.violation == 0;
        stalled_ = 0;
        filling_cursor_ = 0;
        elite_size_ = 0;
    }

    // Returned reference remains owned by this Solver. Copy it if a durable snapshot is needed.
    const Result &solve(const Options &options = {}) {
        const auto begin = Clock::now();
        if (!std::isfinite(options.time_limit) || options.time_limit < 0 ||
            !std::isfinite(options.safety_margin) || options.safety_margin < 0 ||
            options.swap_candidates < 0 || options.swap_candidates > 100000 ||
            options.chain_trials < 0 || options.chain_trials > 100000 ||
            options.double_shift_trials < 0 || options.double_shift_trials > 100000 ||
            options.repair_nodes < 0 || options.chain_depth < 0 || options.chain_depth > max_chain)
            throw std::invalid_argument("gap: invalid search options");
        const double margin =
            std::max(options.safety_margin, std::min(0.002, options.time_limit * 0.01));
        const double usable_time = std::max(0.0, options.time_limit - margin);
        // Clamp before converting a floating duration to an integral clock duration.
        const double safe_time = std::min(usable_time, 315360000.0);
        deadline_ = begin + std::chrono::duration_cast<Clock::duration>(
                                std::chrono::duration<double>(safe_time));
        timed_out_ = false;
        ticks_ = 0;
        use_shift_cache_ = options.incremental_shifts;
        const std::uint64_t initial_iterations = iteration_;
        StopReason reason = StopReason::TimeLimit;
        if (impossible_)
            reason = StopReason::ImpossibleJob;
        else if (n_ == 0 || m_ == 0)
            reason = StopReason::NoMoves;
        else if (options.iteration_limit == 0)
            reason = StopReason::IterationLimit;
        else if (!expired(true)) {
            if (options.initial_prices && !initial_prices_ready_ && !best_.feasible)
                initialize_prices();
            if (!timed_out_ && fill_missing()) {
                save_best();
                if (best_.feasible && !dual_ready_) train_prices();
                while (!expired(true)) {
                    if (iteration_ - initial_iterations >= options.iteration_limit) {
                        reason = StopReason::IterationLimit;
                        break;
                    }
                    Move best_move;
                    search_shifts(best_move);
                    if (!timed_out_ && options.swap_candidates)
                        search_swaps(best_move, options.swap_candidates);
                    if (!timed_out_ && options.double_shift_trials && iteration_ % 8 == 0)
                        search_double_shifts(best_move, options.double_shift_trials);
                    if (!timed_out_ && options.chain_depth >= 2 && options.chain_trials &&
                        (best_move.score >= -1e-12L || stalled_ > 20))
                        search_chains(best_move, options.chain_trials, options.chain_depth);
                    if (!timed_out_ && !best_move.length) search_shifts(best_move, true);
                    if (timed_out_) break;
                    if (!best_move.length) {
                        reason = StopReason::NoMoves;
                        break;
                    }
                    for (int t = 0; t < best_move.length; ++t) {
                        const int j = best_move.jobs[t], old = x_[j];
                        if (old >= 0)
                            tabu_[index(j, old)] = iteration_ + 4 + random_index(9) +
                                                   static_cast<std::uint64_t>(std::sqrt(n_) / 2);
                        relocate(j, best_move.agents[t]);
                    }
                    ++iteration_;
                    ++stalled_;
                    if constexpr (!std::is_integral<Resource>::value ||
                                  !std::is_integral<Cost>::value)
                        if (iteration_ % 128 == 0) rebase();
                    save_best();
                    if (best_.feasible && !dual_ready_) train_prices();
                    if (timed_out_) break;
                    if (iteration_ % 32 == 0 && options.repair_nodes > 0) {
                        repair_subset(options.repair_nodes);
                        if (timed_out_) break;
                        save_best();
                    }
                    if (iteration_ % 16 == 0) adapt_penalties();
                    if (!dual_ready_ && iteration_ % 64 == 0) update_prices();
                    if (timed_out_) break;
                    if (options.restarts &&
                        stalled_ >= static_cast<std::uint64_t>(std::max(150, n_ / 2))) {
                        perturb(options.path_relinking, options.repair_nodes);
                        if (timed_out_ || !fill_missing()) break;
                        save_best();
                    }
                }
            }
            if (unassigned_ && !timed_out_) save_best();
        }
        best_.reason = reason;
        best_.iterations = iteration_;
        best_.moves_evaluated = evaluated_;
        best_.restarts = restarts_;
        best_.elapsed = std::chrono::duration<double>(Clock::now() - begin).count();
        return best_;
    }

    const Result &best() const { return best_; }
    const Result &solve_until(std::chrono::steady_clock::time_point deadline, Options options = {}) {
        const auto now = Clock::now();
        options.time_limit =
            deadline <= now ? 0.0 : std::chrono::duration<double>(deadline - now).count();
        return solve(options);
    }
    const std::vector<int> &current_assignment() const { return x_; }
    int jobs() const { return n_; }
    int agents() const { return m_; }

    // Independent recomputation; does not trust the search's incremental totals.
    Result evaluate(const std::vector<int> &assignment) const {
        if (assignment.size() != static_cast<std::size_t>(n_))
            throw std::invalid_argument("gap: wrong assignment size");
        Result out;
        out.assignment = assignment;
        out.loads.assign(m_, 0);
        for (int j = 0; j < n_; ++j) {
            const int i = assignment[j];
            if (i < -1 || i >= m_ || (i >= 0 && !allowed_[index(j, i)]))
                throw std::invalid_argument("gap: invalid or forbidden assignment");
            if (i < 0)
                ++out.unassigned;
            else {
                out.objective += c_[index(j, i)];
                out.loads[i] += a_[index(j, i)];
            }
        }
        for (int i = 0; i < m_; ++i) out.violation += excess(i, out.loads[i]);
        out.feasible = !out.unassigned && out.violation == 0;
        return out;
    }

private:
    using Clock = std::chrono::steady_clock;
    static constexpr int max_chain = 16;
    struct Move {
        std::array<int, max_chain> jobs{}, agents{};
        int length = 0;
        long double score = std::numeric_limits<long double>::infinity();
    };
    struct Breakpoint {
        long double value = 0, weight = 0;
    };
    std::vector<Breakpoint> breakpoints_;
    std::vector<int> coordinate_agent_, escape_;
    std::vector<long double> coordinate_value_;
    bool initial_prices_ready_ = false, use_shift_cache_ = true;
    struct PairRecord {
        resource_type left = 0, right = 0;
        cost_type cost = 0;
        unsigned mask = 0;
    };
    std::array<std::vector<PairRecord>, 2> pair_records_;
    std::vector<resource_type> pair_coordinates_;
    std::vector<int> pair_fenwick_;
    int n_ = 0, m_ = 0, unassigned_ = 0, filling_cursor_ = 0;
    bool impossible_ = false, timed_out_ = false, visited_feasible_ = false;
    Tolerance tolerance_;
    std::vector<cost_type> c_, min_cost_;
    std::vector<long double> normalized_cost_, normalized_resource_;
    std::vector<resource_type> a_, cap_, load_, refresh_load_;
    std::vector<unsigned char> allowed_, usable_;
    std::vector<long double> tolerance_by_agent_, scale_, penalty_, price_, price_load_,
        load_excess_, refresh_excess_;
    std::vector<int> x_, order_, alternatives_, head_, next_, prev_;
    struct ChainData {
        resource_type resource;
        long double score;
        int job;
        unsigned char usable;
    };
    std::vector<ChainData> chain_data_;
    bool chain_scores_ready_ = true;
    std::vector<int> chain_count_, chain_query_head_, chain_query_next_;
    std::vector<std::array<int, 4>> chain_candidates_;
    std::vector<resource_type> chain_available_;
    std::vector<unsigned char> chain_eligible_;
    std::vector<long double> dual_job_, dual_job_staging_, dual_best_price_;
    bool dual_ready_ = false;
    std::vector<std::uint64_t> tabu_, shift_expiry_;
    std::vector<unsigned char> shift_dirty_;
    std::vector<int> shift_changed_, shift_target_;
    std::vector<long double> shift_value_;
    std::array<std::vector<int>, 8> elite_;
    std::array<cost_type, 8> elite_cost_{};
    std::vector<int> elite_staging_;
    int elite_size_ = 0;
    Result best_, scratch_;
    cost_type objective_ = 0;
    long double cost_scale_ = 1;
    std::uint64_t rng_, iteration_ = 0, evaluated_ = 0, stalled_ = 0, restarts_ = 0,
                        price_iterations_ = 0;
    unsigned ticks_ = 0;
    Clock::time_point deadline_{};

    static long double infinity() { return std::numeric_limits<long double>::infinity(); }
    std::size_t index(int j, int i) const {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(m_) +
               static_cast<std::size_t>(i);
    }
    std::uint64_t random_bits() {
        // SplitMix64's published arithmetic construction, independently expressed.
        rng_ += UINT64_C(0x9e3779b97f4a7c15);
        auto z = rng_;
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    int random_index(int size) {
        return static_cast<int>(random_bits() % static_cast<std::uint64_t>(size));
    }
    long double random_unit() {
        return static_cast<long double>(random_bits() >> 11) / 9007199254740992.0L;
    }
    bool expired(bool force = false) {
        if (timed_out_) return true;
        if (force || (++ticks_ & 63U) == 0) timed_out_ = Clock::now() >= deadline_;
        return timed_out_;
    }
    bool expired_batch(unsigned &remaining) {
        if (timed_out_) return true;
        if (--remaining == 0) {
            remaining = 64;
            timed_out_ = Clock::now() >= deadline_;
        }
        return timed_out_;
    }
    long double excess(int i, resource_type value) const {
        if constexpr (std::is_integral<Resource>::value)
            return value > cap_[i] ? static_cast<long double>(value - cap_[i]) : 0;
        else
            return std::max(0.0L, (value - cap_[i]) - tolerance_by_agent_[i]);
    }
    long double penalty_delta(int i, resource_type after) const {
        return penalty_[i] * ((excess(i, after) - load_excess_[i]) / scale_[i]);
    }
    long double total_excess() const {
        long double v = 0;
        for (int i = 0; i < m_; ++i) v += excess(i, load_[i]);
        return v;
    }
    void unlink(int j, int i) {
        if (prev_[j] < 0)
            head_[i] = next_[j];
        else
            next_[prev_[j]] = next_[j];
        if (next_[j] >= 0) prev_[next_[j]] = prev_[j];
    }
    void relocate(int j, int to, bool track = true) {
        const int from = x_[j];
        if (from == to) return;
        if (track) {
            if (from >= 0) shift_dirty_[from] = 1;
            if (to >= 0) shift_dirty_[to] = 1;
        }
        if (from >= 0) {
            const auto p = index(j, from);
            objective_ -= c_[p];
            load_[from] -= a_[p];
            unlink(j, from);
            load_excess_[from] = excess(from, load_[from]);
        } else
            --unassigned_;
        x_[j] = to;
        if (to >= 0) {
            const auto p = index(j, to);
            objective_ += c_[p];
            load_[to] += a_[p];
            load_excess_[to] = excess(to, load_[to]);
            prev_[j] = -1;
            next_[j] = head_[to];
            if (head_[to] >= 0) prev_[head_[to]] = j;
            head_[to] = j;
        } else {
            ++unassigned_;
            next_[j] = prev_[j] = -1;
        }
    }
    void rebuild_unlimited() {
        std::fill(shift_dirty_.begin(), shift_dirty_.end(), 1);
        std::fill(load_.begin(), load_.end(), 0);
        std::fill(head_.begin(), head_.end(), -1);
        objective_ = 0;
        unassigned_ = 0;
        for (int j = 0; j < n_; ++j) {
            const int i = x_[j];
            next_[j] = prev_[j] = -1;
            if (i < 0) {
                ++unassigned_;
                continue;
            }
            load_[i] += a_[index(j, i)];
            objective_ += c_[index(j, i)];
            next_[j] = head_[i];
            if (head_[i] >= 0) prev_[head_[i]] = j;
            head_[i] = j;
        }
        for (int i = 0; i < m_; ++i) load_excess_[i] = excess(i, load_[i]);
    }
    bool fill_missing() {
        for (; filling_cursor_ < n_; ++filling_cursor_) {
            const int j = order_[filling_cursor_];
            if (expired()) return false;
            if (x_[j] >= 0) continue;
            int chosen = -1;
            long double best_value = infinity();
            bool found_feasible = false;
            for (int i = 0; i < m_; ++i) {
                if (expired()) return false;
                const auto p = index(j, i);
                if (!usable_[p]) continue;
                const auto after = load_[i] + a_[p];
                const bool fits = excess(i, after) == 0;
                const long double value =
                    normalized_cost_[p] + price_[i] * (normalized_resource_[p]) +
                    penalty_delta(i, after) + (restarts_ ? 0.2L * random_unit() : 0);
                if ((fits && !found_feasible) || (fits == found_feasible && value < best_value)) {
                    chosen = i;
                    best_value = value;
                    found_feasible = fits;
                }
            }
            if (chosen >= 0) relocate(j, chosen);
        }
        filling_cursor_ = 0;
        return unassigned_ == 0;
    }
    // Stage the entire result, then publish with O(1) swaps. A timeout cannot corrupt the incumbent.
    void save_best() {
        // A cheap filter; floating-point candidates are re-evaluated before publication.
        long double v = 0;
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            v += excess(i, load_[i]);
        }
        visited_feasible_ = visited_feasible_ || (!unassigned_ && v == 0);
        if (best_.feasible && (unassigned_ || v > 0 || objective_ >= best_.objective)) return;
        if (!best_.feasible && unassigned_ > best_.unassigned) return;
        if (!best_.feasible && unassigned_ == best_.unassigned && v > best_.violation) return;
        scratch_.objective = 0;
        scratch_.unassigned = 0;
        scratch_.violation = 0;
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            scratch_.loads[i] = 0;
        }
        for (int j = 0; j < n_; ++j) {
            if (expired()) return;
            const int i = x_[j];
            scratch_.assignment[j] = i;
            if (i < 0)
                ++scratch_.unassigned;
            else {
                scratch_.loads[i] += a_[index(j, i)];
                scratch_.objective += c_[index(j, i)];
            }
        }
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            scratch_.violation += excess(i, scratch_.loads[i]);
        }
        scratch_.feasible = !scratch_.unassigned && scratch_.violation == 0;
        const bool improve =
            scratch_.feasible ? (!best_.feasible || scratch_.objective < best_.objective)
                              : (!best_.feasible && (scratch_.unassigned < best_.unassigned ||
                                                     (scratch_.unassigned == best_.unassigned &&
                                                      (scratch_.violation < best_.violation ||
                                                       (scratch_.violation == best_.violation &&
                                                        scratch_.objective < best_.objective)))));
        if (improve) {
            best_.assignment.swap(scratch_.assignment);
            best_.loads.swap(scratch_.loads);
            best_.objective = scratch_.objective;
            best_.unassigned = scratch_.unassigned;
            best_.violation = scratch_.violation;
            best_.feasible = scratch_.feasible;
            stalled_ = 0;
            if (best_.feasible) remember_elite();
        }
    }
    void remember_elite() {
        int slot = -1, closest = n_ + 1;
        for (int e = 0; e < elite_size_; ++e) {
            int distance = 0;
            for (int j = 0; j < n_; ++j) {
                if (expired()) return;
                distance += elite_[e][j] != best_.assignment[j];
            }
            if (distance < closest) {
                closest = distance;
                slot = e;
            }
        }
        if (closest == 0) return;
        if (closest >= std::max(2, n_ / 25)) {
            if (elite_size_ < 8)
                slot = elite_size_;
            else
                slot = static_cast<int>(std::max_element(elite_cost_.begin(), elite_cost_.end()) -
                                        elite_cost_.begin());
        }
        if (slot < 0) slot = 0;
        for (int j = 0; j < n_; ++j) {
            if (expired()) return;
            elite_staging_[j] = best_.assignment[j];
        }
        elite_[slot].swap(elite_staging_);
        elite_cost_[slot] = best_.objective;
        if (slot == elite_size_) ++elite_size_;
    }
    long double shift_delta(int j, int to, long double removal) const {
        const int from = x_[j];
        const auto p = index(j, from), q = index(j, to);
        return static_cast<long double>(c_[q] - c_[p]) / cost_scale_ + removal +
               penalty_delta(to, load_[to] + a_[q]);
    }
    long double shift_delta(int j, int to) const {
        const int from = x_[j];
        return shift_delta(j, to, penalty_delta(from, load_[from] - a_[index(j, from)]));
    }
    bool aspiration(cost_type candidate_cost, int a, resource_type la, int b, resource_type lb) {
        if (best_.feasible && candidate_cost >= best_.objective) return false;
        if (excess(a, la) > 0 || excess(b, lb) > 0) return false;
        for (int i = 0; i < m_; ++i) {
            if (expired()) return false;
            if (i != a && i != b && excess(i, load_[i]) > 0) return false;
        }
        return true;
    }
    void search_shifts(Move &move, bool ignore_tabu = false) {
        unsigned clock_budget = 1;
        if (n_ >= 256 && use_shift_cache_ && !ignore_tabu) {
            search_shifts_cached(move);
            return;
        }
        const int start = random_index(n_);
        for (int k = 0; k < n_; ++k) {
            const int j = k < n_ - start ? start + k : k - (n_ - start);
            const int from = x_[j];
            const auto removal = penalty_delta(from, load_[from] - a_[index(j, from)]);
            for (int to = 0; to < m_; ++to) {
                if (expired_batch(clock_budget)) return;
                const auto q = index(j, to);
                if (to == from || !usable_[q]) continue;
                ++evaluated_;
                const long double score = shift_delta(j, to, removal);
                if (score >= move.score) continue;
                const auto p = index(j, from);
                if (!ignore_tabu && tabu_[q] > iteration_ &&
                    !aspiration(objective_ + c_[q] - c_[p], from, load_[from] - a_[p], to,
                                load_[to] + a_[q]))
                    continue;
                move.score = score;
                move.length = 1;
                move.jobs[0] = j;
                move.agents[0] = to;
            }
        }
    }
    // A shift depends only on its source and destination loads/penalties.
    // Tabu expiration invalidates the corresponding row even without a load change.
    void search_shifts_cached(Move &move) {
        unsigned clock_budget = 1;
        int changed = 0;
        for (int i = 0; i < m_; ++i) {
            if (expired_batch(clock_budget)) return;
            if (shift_dirty_[i]) shift_changed_[changed++] = i;
        }
        const int start = random_index(n_);
        for (int k = 0; k < n_; ++k) {
            if (expired_batch(clock_budget)) return;
            const int j = k < n_ - start ? start + k : k - (n_ - start), from = x_[j];
            int target = shift_target_[j];
            const bool full = target < 0 || shift_dirty_[from] || shift_dirty_[target] ||
                              shift_expiry_[j] <= iteration_;
            const auto removal = penalty_delta(from, load_[from] - a_[index(j, from)]);
            auto best = full ? infinity() : shift_value_[j];
            auto next_expiry = full ? std::numeric_limits<std::uint64_t>::max() : shift_expiry_[j];
            if (full) target = -1;
            for (int q = 0; q < (full ? m_ : changed); ++q) {
                if (expired_batch(clock_budget)) return;
                const int to = full ? q : shift_changed_[q];
                const auto p = index(j, to);
                if (to == from || !usable_[p]) continue;
                ++evaluated_;
                if (tabu_[p] > iteration_) {
                    next_expiry = std::min(next_expiry, tabu_[p]);
                    continue;
                }
                const auto value = shift_delta(j, to, removal);
                if (value < best) {
                    best = value;
                    target = to;
                }
            }
            shift_target_[j] = target;
            shift_value_[j] = best;
            shift_expiry_[j] = next_expiry;
            if (target >= 0 && best < move.score) {
                move.score = best;
                move.length = 1;
                move.jobs[0] = j;
                move.agents[0] = target;
            }
        }
        // Only a complete refresh can acknowledge dirty agents.
        for (int i = 0; i < m_; ++i) {
            if (expired_batch(clock_budget)) return;
            shift_dirty_[i] = 0;
        }
    }
    struct SwapSource {
        int job, agent;
        std::size_t row;
        cost_type cost;
        resource_type remaining;
    };
    SwapSource swap_source(int j) const {
        const auto row = static_cast<std::size_t>(j) * m_;
        const int from = x_[j];
        const auto p = row + from;
        return {j, from, row, c_[p], static_cast<resource_type>(load_[from] - a_[p])};
    }
    void consider_swap(const SwapSource &source, int k, Move &move) {
        const int j = source.job, from = source.agent, to = x_[k];
        if (from == to || !usable_[source.row + to] || !usable_[index(k, from)]) return;
        ++evaluated_;
        const auto jq = source.row + to, kp = index(k, from), kq = index(k, to);
        const resource_type lf = source.remaining + a_[kp], lt = load_[to] - a_[kq] + a_[jq];
        const cost_type dc = (c_[jq] - source.cost) + (c_[kp] - c_[kq]);
        const long double score = static_cast<long double>(dc) / cost_scale_ +
                                  penalty_delta(from, lf) + penalty_delta(to, lt);
        if (score >= move.score) return;
        if ((tabu_[jq] > iteration_ || tabu_[kp] > iteration_) &&
            !aspiration(objective_ + dc, from, lf, to, lt))
            return;
        move.length = 2;
        move.jobs[0] = j;
        move.jobs[1] = k;
        move.agents[0] = to;
        move.agents[1] = from;
        move.score = score;
    }
    void consider_swap(int j, int k, Move &move) { consider_swap(swap_source(j), k, move); }
    void search_swaps(Move &move, int candidates) {
        unsigned clock_budget = 1;
        const int start = random_index(n_);
        for (int t = 0; t < n_; ++t) {
            const int j = t < n_ - start ? start + t : t - (n_ - start);
            const auto source = swap_source(j);
            if (n_ <= 160) {
                for (int k = j + 1; k < n_; ++k) {
                    if (expired_batch(clock_budget)) return;
                    consider_swap(source, k, move);
                }
            } else {
                for (int s = 0; s < candidates; ++s) {
                    if (expired_batch(clock_budget)) return;
                    consider_swap(source, random_index(n_), move);
                }
            }
        }
    }
    // Exact one-coordinate optimization of the capacity-relaxed dual, without an incumbent.
    // The weighted threshold is selected with an interruptible partition, not a library sort.
    void initialize_prices() {
        const auto now = Clock::now();
        if (deadline_ - now < std::chrono::milliseconds(5)) return;
        initial_prices_ready_ = true;
        const auto limit = now + (deadline_ - now) / 8;
        const auto adjusted = [&](int j, int i) {
            const auto p = index(j, i);
            return normalized_cost_[p] + price_[i] * (normalized_resource_[p]);
        };
        // First optimize a common resource price. This moves along a coupled direction
        // that one-coordinate ascent cannot traverse efficiently on nearly flat dual faces.
        long double capacity = 0;
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            capacity += cap_[i];
        }
        if (capacity > 0 && std::isfinite(capacity)) {
            const auto unit =
                std::max(std::numeric_limits<long double>::min(), capacity / std::max(1, n_));
            long double lo = 0, hi = 1, t = 0, best_t = 0, best_bound = -infinity();
            bool bracketed = false;
            for (int pass = 0; pass < 36; ++pass) {
                if (pass && Clock::now() >= limit) break;
                long double bound = -t * (capacity / unit), resource = 0;
                for (int j = 0; j < n_; ++j) {
                    long double value = infinity(), consumption = 0;
                    for (int i = 0; i < m_; ++i) {
                        if (expired()) return;
                        const auto p = index(j, i);
                        if (!usable_[p]) continue;
                        const auto av = a_[p] / unit;
                        const auto v = normalized_cost_[p] + (t > 0 ? t * av : 0);
                        if (v < value || (v == value && av < consumption)) {
                            value = v;
                            consumption = av;
                        }
                    }
                    bound += value;
                    resource += consumption;
                }
                if (!std::isfinite(bound) || !std::isfinite(resource)) break;
                if (bound > best_bound) {
                    best_bound = bound;
                    best_t = t;
                }
                if (resource <= capacity / unit) {
                    if (t == 0) break;
                    hi = t;
                    bracketed = true;
                } else
                    lo = t;
                if (bracketed)
                    t = (lo + hi) * .5L;
                else {
                    t = hi;
                    hi *= 2;
                }
            }
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                price_[i] = best_t > 0 ? std::min(best_t * (scale_[i] / unit), 1e10L) : 0;
            }
        }
        if (Clock::now() >= limit) return;
        const auto refresh = [&](int j) {
            int agent = -1;
            long double value = infinity();
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                if (!usable_[index(j, i)]) continue;
                const auto v = adjusted(j, i);
                if (v < value) {
                    value = v;
                    agent = i;
                }
            }
            coordinate_agent_[j] = agent;
            coordinate_value_[j] = value;
        };
        for (int j = 0; j < n_; ++j) {
            refresh(j);
            if (timed_out_) return;
        }
        for (int sweep = 0; sweep < 4; ++sweep) {
            for (int i = 0; i < m_; ++i) {
                if (Clock::now() >= limit) return;
                int count = 0;
                long double remaining = cap_[i] / scale_[i], weight = 0;
                for (int j = 0; j < n_; ++j) {
                    if (expired()) return;
                    const auto p = index(j, i);
                    if (!usable_[p] || a_[p] == 0) continue;
                    long double other = coordinate_value_[j];
                    if (coordinate_agent_[j] == i) {
                        other = infinity();
                        for (int h = 0; h < m_; ++h) {
                            if (expired()) return;
                            if (h != i && usable_[index(j, h)])
                                other = std::min(other, adjusted(j, h));
                        }
                    }
                    const auto resource = normalized_resource_[p];
                    if (!std::isfinite(other)) {
                        remaining -= resource;
                        continue;
                    }
                    const auto threshold = (other - normalized_cost_[p]) / resource;
                    if (threshold > 0) {
                        breakpoints_[count++] = {threshold, resource};
                        weight += resource;
                    }
                }
                long double value = 0;
                if (remaining < 0)
                    value = 1e10L;
                else if (weight > remaining) {
                    int first = 0, last = count;
                    while (first < last) {
                        if (expired()) return;
                        const auto pivot = breakpoints_[first + (last - first) / 2].value;
                        int high = first, scan = first, low = last;
                        long double greater_weight = 0, equal_weight = 0;
                        while (scan < low) {
                            if (expired()) return;
                            const auto v = breakpoints_[scan];
                            if (v.value > pivot) {
                                greater_weight += v.weight;
                                std::swap(breakpoints_[high++], breakpoints_[scan++]);
                            } else if (v.value < pivot)
                                std::swap(breakpoints_[scan], breakpoints_[--low]);
                            else {
                                equal_weight += v.weight;
                                ++scan;
                            }
                        }
                        if (remaining < greater_weight)
                            last = high;
                        else if (remaining <= greater_weight + equal_weight) {
                            value = pivot;
                            break;
                        } else {
                            remaining -= greater_weight + equal_weight;
                            first = low;
                        }
                    }
                }
                price_[i] = std::clamp(value, 0.0L, 1e10L);
                for (int j = 0; j < n_; ++j) {
                    if (expired()) return;
                    if (coordinate_agent_[j] == i)
                        refresh(j);
                    else if (usable_[index(j, i)]) {
                        const auto v = adjusted(j, i);
                        if (v < coordinate_value_[j]) {
                            coordinate_value_[j] = v;
                            coordinate_agent_[j] = i;
                        }
                    }
                    if (timed_out_) return;
                }
            }
        }
    }
    // Capacity-relaxed LP dual. Prices are normalized by resource and cost scales.
    // The best feasible dual iterate supplies job potentials for chain selection.
    void train_prices() {
        const auto now = Clock::now();
        // Tiny continuation slices must remain available for ordinary search.
        if (deadline_ - now < std::chrono::milliseconds(5)) return;
        const auto training_end = now + (deadline_ - now) / 5;
        long double upper = 0;
        for (int j = 0; j < n_; ++j) {
            if (expired()) return;
            upper += normalized_cost_[index(j, best_.assignment[j])];
        }
        long double best_bound = -infinity(), factor = 1.5L;
        int stale = 0;
        for (int pass = 0; pass < 120; ++pass) {
            if (pass && Clock::now() >= training_end) break;
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                price_load_[i] = 0;
            }
            long double bound = 0;
            for (int j = 0; j < n_; ++j) {
                int chosen = -1;
                long double value = infinity();
                for (int i = 0; i < m_; ++i) {
                    if (expired()) return;
                    const auto p = index(j, i);
                    if (!usable_[p]) continue;
                    const auto v = normalized_cost_[p] + price_[i] * (normalized_resource_[p]);
                    if (v < value) {
                        value = v;
                        chosen = i;
                    }
                }
                bound += value;
                price_load_[chosen] += normalized_resource_[index(j, chosen)];
            }
            long double norm = 0;
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                bound -= price_[i] * (cap_[i] / scale_[i]);
                price_load_[i] -= cap_[i] / scale_[i];
                const auto g = price_load_[i];
                if (price_[i] > 0 || g > 0) norm += g * g;
            }
            if (bound > best_bound) {
                best_bound = bound;
                stale = 0;
                for (int i = 0; i < m_; ++i) {
                    if (expired()) return;
                    dual_best_price_[i] = price_[i];
                }
            } else if (++stale % 15 == 0)
                factor *= .5L;
            if (norm < 1e-20L || factor < .005L) break;
            const auto step = factor * std::max(0.0L, upper - bound) / norm;
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                price_[i] = std::clamp(price_[i] + step * price_load_[i], 0.0L, 1e10L);
            }
        }
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            price_[i] = dual_best_price_[i];
        }
        for (int j = 0; j < n_; ++j) {
            long double value = infinity();
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                const auto p = index(j, i);
                if (!usable_[p]) continue;
                value = std::min(value, normalized_cost_[p] + price_[i] * (normalized_resource_[p]));
            }
            dual_job_staging_[j] = value;
        }
        dual_job_.swap(dual_job_staging_);
        chain_scores_ready_ = false;
        dual_ready_ = true;
    }
    void prepare_chains() {
        unsigned clock_budget = 1;
        // Job potentials change only when a complete dual training result is published.
        // An interrupted refresh stays invalid and is rebuilt before any query uses it.
        if (!chain_scores_ready_) {
            for (int i = 0; i < m_; ++i) {
                const auto offset = static_cast<std::size_t>(i) * n_;
                for (int k = 0; k < n_; ++k) {
                    if (expired_batch(clock_budget)) return;
                    auto &data = chain_data_[offset + k];
                    data.score = dual_job_[data.job] - normalized_cost_[index(data.job, i)];
                }
            }
            chain_scores_ready_ = true;
        }
        // Each assigned job asks for one prefix. Bucket those queries instead of
        // materializing all n prefixes for each agent. Clear also after an interruption.
        for (int j = 0; j < n_; ++j) {
            if (expired_batch(clock_budget)) return;
            chain_query_head_[j] = -1;
        }
        for (int i = 0; i < m_; ++i) {
            if (expired_batch(clock_budget)) return;
            resource_type least = std::numeric_limits<resource_type>::max();
            const auto offset = static_cast<std::size_t>(i) * n_;
            for (int k = 0; k < n_; ++k) {
                if (expired_batch(clock_budget)) return;
                const auto &data = chain_data_[offset + k];
                const int j = data.job;
                if (x_[j] != i && data.usable) {
                    least = data.resource;
                    break;
                }
            }
            const resource_type over = load_[i] > cap_[i] ? load_[i] - cap_[i] : 0;
            for (int j = head_[i]; j >= 0; j = next_[j]) {
                if (expired_batch(clock_budget)) return;
                const resource_type used = a_[index(j, i)];
                chain_available_[j] = used > over ? used - over : used;
                chain_eligible_[j] = least <= chain_available_[j];
            }
        }
        // Prefix maxima take O(n*m); the job threshold queries take O(n*log(n)).
        // Retain up to four equal maxima; no fallback to a lower score on a cycle.
        for (int i = 0; i < m_; ++i) {
            const auto offset = static_cast<std::size_t>(i) * n_;
            for (int j = head_[i]; j >= 0; j = next_[j]) {
                if (expired_batch(clock_budget)) return;
                int lo = 0, hi = n_;
                while (lo < hi) {
                    const int mid = lo + (hi - lo) / 2;
                    if (chain_data_[offset + mid].resource <= chain_available_[j])
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                chain_count_[j] = 0;
                if (lo) {
                    chain_query_next_[j] = chain_query_head_[lo - 1];
                    chain_query_head_[lo - 1] = j;
                }
            }
            std::array<int, 4> candidates{};
            int count = 0, ties = 0;
            long double highest = -infinity();
            for (int k = 0; k < n_; ++k) {
                if (expired_batch(clock_budget)) return;
                const auto &data = chain_data_[offset + k];
                const int j = data.job;
                if (x_[j] != i && chain_eligible_[j] && data.usable) {
                    const auto score = data.score;
                    if (score > highest) {
                        highest = score;
                        count = ties = 0;
                    }
                    if (score == highest) {
                        ++ties;
                        if (count < 4)
                            candidates[count++] = j;
                        else {
                            const int q = random_index(ties);
                            if (q < 4) candidates[q] = j;
                        }
                    }
                }
                if (chain_query_head_[k] >= 0) {
                    for (int job = chain_query_head_[k]; job >= 0; job = chain_query_next_[job]) {
                        if (expired_batch(clock_budget)) return;
                        chain_candidates_[job] = candidates;
                        chain_count_[job] = count;
                    }
                    chain_query_head_[k] = -1;
                }
            }
        }
    }
    void search_double_shifts(Move &move, int trials) {
        unsigned clock_budget = 1;
        if (m_ < 3) return;
        const bool cached = n_ >= 256 && use_shift_cache_;
        if (!cached)
            for (int j = 0; j < n_; ++j) {
                int target = -1;
                long double value = infinity();
                for (int i = 0; i < m_; ++i) {
                    if (expired_batch(clock_budget)) return;
                    if (i == x_[j] || !usable_[index(j, i)] || tabu_[index(j, i)] > iteration_)
                        continue;
                    const auto v = shift_delta(j, i);
                    if (v < value) {
                        value = v;
                        target = i;
                    }
                }
                escape_[j] = target;
            }
        const int start = random_index(n_);
        for (int q = 0; q < std::min(n_, trials); ++q) {
            const int j = q < n_ - start ? start + q : q - (n_ - start), from = x_[j];
            std::array<int, 3> targets{{-1, -1, -1}};
            std::array<long double, 3> values{{infinity(), infinity(), infinity()}};
            for (int i = 0; i < m_; ++i) {
                if (expired_batch(clock_budget)) return;
                const auto p = index(j, i);
                if (i == from || !usable_[p] || tabu_[p] > iteration_) continue;
                const auto v = normalized_cost_[p] + price_[i] * (normalized_resource_[p]);
                for (int t = 0; t < 3; ++t)
                    if (v < values[t]) {
                        for (int h = 2; h > t; --h) {
                            values[h] = values[h - 1];
                            targets[h] = targets[h - 1];
                        }
                        values[t] = v;
                        targets[t] = i;
                        break;
                    }
            }
            for (int to : targets)
                if (to >= 0)
                    for (int k = head_[to]; k >= 0; k = next_[k]) {
                        if (expired_batch(clock_budget)) return;
                        consider_swap(j, k, move);
                        const int last = cached ? shift_target_[k] : escape_[k];
                        if (last < 0 || last == from) continue;
                        const auto jp = index(j, from), jq = index(j, to), kp = index(k, to),
                                   kq = index(k, last);
                        if (tabu_[kq] > iteration_) continue;
                        const cost_type delta = c_[jq] - c_[jp] + c_[kq] - c_[kp];
                        const auto score = static_cast<long double>(delta) / cost_scale_ +
                                           penalty_delta(from, load_[from] - a_[jp]) +
                                           penalty_delta(to, load_[to] - a_[kp] + a_[jq]) +
                                           penalty_delta(last, load_[last] + a_[kq]);
                        ++evaluated_;
                        if (score < move.score) {
                            move.score = score;
                            move.length = 2;
                            move.jobs[0] = j;
                            move.agents[0] = to;
                            move.jobs[1] = k;
                            move.agents[1] = last;
                        }
                    }
        }
    }
    void search_chains(Move &best_move, int trials, int depth_limit) {
        prepare_chains();
        const int start = random_index(n_);
        for (int trial = 0; trial < std::min(trials, n_); ++trial) {
            if (expired()) return;
            std::array<int, max_chain> jobs{}, from{}, to{};
            const int root = trial < n_ - start ? start + trial : trial - (n_ - start);
            const int origin = x_[root];
            if (!chain_eligible_[root]) continue;
            const auto original_cost = objective_;
            std::array<int, max_chain> touched{};
            std::array<resource_type, max_chain> saved_loads{};
            int touched_count = 1;
            touched[0] = origin;
            saved_loads[0] = load_[origin];
            long double penalty_change =
                penalty_delta(origin, load_[origin] - a_[index(root, origin)]);
            relocate(root, -1, false);
            int hole = origin, length = 0, last = root;
            for (int depth = 1; depth < depth_limit; ++depth) {
                int chosen = -1, ties = 0;
                for (int t = 0; t < chain_count_[last]; ++t) {
                    const int j = chain_candidates_[last][t];
                    bool seen = j == root;
                    for (int k = 0; k < length; ++k) seen = seen || jobs[k] == j;
                    if (!seen && tabu_[index(j, hole)] <= iteration_ && random_index(++ties) == 0)
                        chosen = j;
                }
                if (chosen < 0) break;
                const int source = x_[chosen];
                bool known = false;
                for (int q = 0; q < touched_count; ++q) known = known || touched[q] == source;
                if (!known) {
                    touched[touched_count] = source;
                    saved_loads[touched_count++] = load_[source];
                }
                penalty_change += penalty_delta(source, load_[source] - a_[index(chosen, source)]) +
                                  penalty_delta(hole, load_[hole] + a_[index(chosen, hole)]);
                jobs[length] = chosen;
                from[length] = source;
                to[length] = hole;
                relocate(chosen, hole, false);
                ++length;
                hole = source;
                last = chosen;
                for (int i = 0; i < m_; ++i) {
                    if (expired()) break;
                    const auto rp = index(root, i);
                    if (!usable_[rp] || tabu_[rp] > iteration_) continue;
                    ++evaluated_;
                    const auto delta = (objective_ - original_cost) + c_[rp];
                    const long double score = static_cast<long double>(delta) / cost_scale_ +
                                              penalty_change + penalty_delta(i, load_[i] + a_[rp]);
                    if (score < best_move.score) {
                        best_move.score = score;
                        best_move.length = length + 1;
                        for (int k = 0; k < length; ++k) {
                            best_move.jobs[k] = jobs[k];
                            best_move.agents[k] = to[k];
                        }
                        best_move.jobs[length] = root;
                        best_move.agents[length] = i;
                    }
                }
                if (timed_out_) break;
            }
            for (int k = length - 1; k >= 0; --k) relocate(jobs[k], from[k], false);
            relocate(root, origin, false);
            for (int q = 0; q < touched_count; ++q) {
                const int i = touched[q];
                load_[i] = saved_loads[q];
                load_excess_[i] = excess(i, load_[i]);
            }
            objective_ = original_cost;
            if (timed_out_) return;
        }
    }
    void adapt_penalties() {
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            shift_dirty_[i] = 1;
            if (visited_feasible_)
                penalty_[i] = std::max(0.05L, penalty_[i] * 0.9L);
            else if (excess(i, load_[i]) > 0)
                penalty_[i] = std::min(1e12L, penalty_[i] * 1.2L);
        }
        visited_feasible_ = false;
    }
    void rebase() {
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            refresh_load_[i] = 0;
        }
        cost_type cost = 0;
        for (int j = 0; j < n_; ++j) {
            if (expired()) return;
            if (x_[j] >= 0) {
                const auto p = index(j, x_[j]);
                cost += c_[p];
                refresh_load_[x_[j]] += a_[p];
            }
        }
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            refresh_excess_[i] = excess(i, refresh_load_[i]);
        }
        load_.swap(refresh_load_);
        load_excess_.swap(refresh_excess_);
        objective_ = cost;
        std::fill(shift_dirty_.begin(), shift_dirty_.end(), 1);
    }
    void update_prices() {
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            price_load_[i] = 0;
        }
        for (int j = 0; j < n_; ++j) {
            int chosen = -1;
            long double best_value = infinity();
            for (int i = 0; i < m_; ++i) {
                if (expired()) return;
                const auto p = index(j, i);
                if (!usable_[p]) continue;
                const long double v = normalized_cost_[p] + price_[i] * (normalized_resource_[p]);
                if (v < best_value) {
                    best_value = v;
                    chosen = i;
                }
            }
            if (chosen >= 0) price_load_[chosen] += normalized_resource_[index(j, chosen)];
        }
        ++price_iterations_;
        const long double step = 2 / std::sqrt(static_cast<long double>(price_iterations_));
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            price_[i] =
                std::clamp(price_[i] + step * (price_load_[i] - cap_[i] / scale_[i]), 0.0L, 1e10L);
        }
    }
    void repair_subset(int node_limit) {
        if (m_ < 2 || unassigned_) return;
        for (int i = 0; i < m_; ++i) {
            if (expired()) return;
            if (excess(i, load_[i]) > 0) return;
        }
        const int width = 2 + random_index(std::min(4, m_) - 1);
        std::array<int, 4> agents{};
        agents[0] = random_index(m_);
        for (int k = 1; k < width; ++k) {
            int selected = -1;
            long double best = infinity();
            // Grow a group around inexpensive moves currently blocked by capacity.
            if (random_index(4)) {
                const int root = agents[random_index(k)];
                int sampled = -1, count = 0;
                for (int j = head_[root]; j >= 0; j = next_[j]) {
                    if (expired()) return;
                    if (random_index(++count) == 0) sampled = j;
                }
                if (sampled >= 0)
                    for (int i = 0; i < m_; ++i) {
                        if (expired()) return;
                        bool used = false;
                        for (int t = 0; t < k; ++t) used = used || agents[t] == i;
                        if (used || !usable_[index(sampled, i)]) continue;
                        const auto p = index(sampled, i);
                        const auto value =
                            normalized_cost_[p] + price_[i] * (normalized_resource_[p]);
                        if (value < best) {
                            best = value;
                            selected = i;
                        }
                    }
            }
            if (selected < 0) {
                int remaining = m_ - k, choice = random_index(remaining);
                for (int i = 0; i < m_; ++i) {
                    if (expired()) return;
                    bool used = false;
                    for (int t = 0; t < k; ++t) used = used || agents[t] == i;
                    if (!used && choice-- == 0) {
                        selected = i;
                        break;
                    }
                }
            }
            agents[k] = selected;
        }
        constexpr int maximum = 40;
        const int limit = width == 2 ? 22 : maximum;
        std::array<int, maximum> jobs{}, answer{}, path{};
        std::array<long double, maximum> priority{};
        int count = 0;
        for (int k = 0; k < width; ++k)
            for (int j = head_[agents[k]]; j >= 0; j = next_[j]) {
                if (expired()) return;
                int choices = 0;
                long double alternate = infinity();
                const auto old = index(j, x_[j]);
                const auto current =
                    normalized_cost_[old] + price_[x_[j]] * (normalized_resource_[old]);
                for (int t = 0; t < width; ++t)
                    if (usable_[index(j, agents[t])]) {
                        ++choices;
                        const int i = agents[t];
                        const auto p = index(j, i);
                        if (i != x_[j])
                            alternate =
                                std::min(alternate, normalized_cost_[p] +
                                                        price_[i] * (normalized_resource_[p]));
                    }
                if (choices < 2) continue;
                // Ambiguous assignments and expensive current choices receive higher priority.
                const auto key = alternate - current + random_unit() * .8L;
                int slot = count;
                if (count == limit) {
                    slot = static_cast<int>(
                        std::max_element(priority.begin(), priority.begin() + count) -
                        priority.begin());
                    if (key >= priority[slot]) continue;
                } else
                    ++count;
                jobs[slot] = j;
                priority[slot] = key;
            }
        if (count < 2) return;
        std::array<resource_type, 4> loads{};
        for (int t = 0; t < width; ++t) loads[t] = load_[agents[t]];
        cost_type bound = 0;
        for (int k = 0; k < count; ++k) {
            const int j = jobs[k], i = x_[j];
            answer[k] = i;
            bound += c_[index(j, i)];
            for (int t = 0; t < width; ++t)
                if (i == agents[t]) loads[t] -= a_[index(j, i)];
        }
        const auto original = bound;
        if (width == 2 && count <= 22 &&
            node_limit >= (1 << (count / 2)) + (1 << (count - count / 2))) {
            std::array<int, 22> small_jobs{}, small_answer{};
            std::copy_n(jobs.begin(), count, small_jobs.begin());
            std::copy_n(answer.begin(), count, small_answer.begin());
            repair_pair_mitm(small_jobs, count, {agents[0], agents[1], 0}, {loads[0], loads[1], 0},
                             small_answer, bound);
            if (!timed_out_ && bound < original)
                for (int k = 0; k < count; ++k) relocate(jobs[k], small_answer[k]);
            return;
        }
        std::array<long double, 4> lambda{}, best_lambda{}, remaining{}, assigned{};
        for (int t = 0; t < width; ++t) {
            lambda[t] = price_[agents[t]];
            remaining[t] = (static_cast<long double>(cap_[agents[t]] - loads[t]) +
                            tolerance_by_agent_[agents[t]]) /
                           scale_[agents[t]];
        }
        long double upper = 0;
        for (int k = 0; k < count; ++k) upper += normalized_cost_[index(jobs[k], x_[jobs[k]])];
        long double best_dual = -infinity(), factor = 1.0L;
        for (int pass = 0; pass < 32; ++pass) {
            assigned.fill(0);
            long double dual = 0;
            for (int k = 0; k < count; ++k) {
                if (expired()) return;
                const int j = jobs[k];
                int chosen = -1;
                long double low = infinity();
                for (int t = 0; t < width; ++t) {
                    const int i = agents[t];
                    const auto p = index(j, i);
                    if (!usable_[p]) continue;
                    const auto v = normalized_cost_[p] + lambda[t] * (normalized_resource_[p]);
                    if (v < low) {
                        low = v;
                        chosen = t;
                    }
                }
                dual += low;
                assigned[chosen] += normalized_resource_[index(j, agents[chosen])];
            }
            long double norm = 0;
            for (int t = 0; t < width; ++t) {
                dual -= lambda[t] * remaining[t];
                assigned[t] -= remaining[t];
                if (lambda[t] > 0 || assigned[t] > 0) norm += assigned[t] * assigned[t];
            }
            if (dual > best_dual) {
                best_dual = dual;
                best_lambda = lambda;
            }
            if (norm < 1e-20L) break;
            const auto step = factor * std::max(0.0L, upper - dual) / norm;
            for (int t = 0; t < width; ++t)
                lambda[t] = std::clamp(lambda[t] + step * assigned[t], 0.0L, 1e10L);
            if (pass % 8 == 7) factor *= .6L;
        }
        lambda = best_lambda;
        const auto adjusted = [&](int j, int t) {
            const int i = agents[t];
            const auto p = index(j, i);
            return normalized_cost_[p] + lambda[t] * (normalized_resource_[p]);
        };
        const auto regret = [&](int j) {
            long double low = infinity(), second = infinity();
            for (int t = 0; t < width; ++t)
                if (usable_[index(j, agents[t])]) {
                    const auto v = adjusted(j, t);
                    if (v < low) {
                        second = low;
                        low = v;
                    } else
                        second = std::min(second, v);
                }
            return second - low;
        };
        std::sort(jobs.begin(), jobs.begin() + count,
                  [&](int j, int k) { return regret(j) > regret(k); });
        std::array<std::array<int, 4>, maximum> choices{};
        std::array<int, maximum> choice_count{};
        std::array<cost_type, maximum + 1> suffix{}, min_suffix{};
        std::array<long double, maximum + 1> dual_suffix{};
        for (int k = count - 1; k >= 0; --k) {
            const int j = jobs[k];
            answer[k] = x_[j];
            cost_type low = std::numeric_limits<cost_type>::max();
            for (int t = 0; t < width; ++t)
                if (usable_[index(j, agents[t])]) {
                    choices[k][choice_count[k]++] = t;
                    low = std::min(low, c_[index(j, agents[t])]);
                }
            std::sort(choices[k].begin(), choices[k].begin() + choice_count[k],
                      [&](int a, int b) { return adjusted(j, a) < adjusted(j, b); });
            suffix[k] = suffix[k + 1] + low;
            min_suffix[k] = min_suffix[k + 1] + min_cost_[j];
            dual_suffix[k] = dual_suffix[k + 1] + adjusted(j, choices[k][0]);
        }
        long double initial_price = 0;
        for (int t = 0; t < width; ++t) initial_price -= lambda[t] * remaining[t];
        int nodes = 0;
        const auto dfs = [&](auto &&self, int k, cost_type cost, long double used_price) -> void {
            if (++nodes > node_limit || expired()) return;
            if (cost + suffix[k] >= bound) return;
            // Compare after subtracting integral constants, to preserve small integer cost changes.
            const auto dual_gap =
                static_cast<long double>(cost + min_suffix[k] - bound) / cost_scale_ +
                dual_suffix[k] + used_price;
            const auto roundoff = 1e-10L * (1 + std::fabs(dual_suffix[k]) + std::fabs(used_price));
            if (dual_gap > roundoff) return;
            if (k == count) {
                bound = cost;
                answer = path;
                return;
            }
            const int j = jobs[k];
            for (int q = 0; q < choice_count[k]; ++q) {
                const int t = choices[k][q], i = agents[t];
                const auto p = index(j, i);
                const auto previous = loads[t];
                loads[t] += a_[p];
                if (excess(i, loads[t]) == 0) {
                    path[k] = i;
                    self(self, k + 1, cost + c_[p],
                         used_price + lambda[t] * (normalized_resource_[p]));
                }
                loads[t] = previous;
                if (nodes >= node_limit || timed_out_) return;
            }
        };
        dfs(dfs, 0, 0, initial_price);
        if (timed_out_ || bound >= original) return;
        for (int k = 0; k < count; ++k) relocate(jobs[k], answer[k]);
    }
    // Independent meet-in-the-middle solution of a two-agent subset.
    // Two resource constraints become an offline 2D dominance minimum query.
    void repair_pair_mitm(const std::array<int, 22> &jobs, int count,
                          const std::array<int, 3> &agents, const std::array<resource_type, 3> &base,
                          std::array<int, 22> &answer, cost_type &bound) {
        // The bounded library sorts below are not interruptible. Leave ample room before starting.
        if (deadline_ - Clock::now() < std::chrono::milliseconds(2)) return;
        const int split = count / 2;
        std::array<int, 2> sizes{};
        for (int half = 0; half < 2; ++half) {
            const int begin = half == 0 ? 0 : split, end = half == 0 ? split : count;
            auto &records = pair_records_[half];
            records[0] = {};
            // Reuse assignment prefixes. Each sum still adds jobs in the original
            // order (also for floating point), and final masks stay in ascending order.
            for (int k = begin; k < end; ++k) {
                const unsigned bit = 1U << (k - begin);
                const auto left = index(jobs[k], agents[0]), right = index(jobs[k], agents[1]);
                for (unsigned mask = 0; mask < bit; ++mask) {
                    if (expired()) return;
                    const auto prefix = records[mask];
                    auto l = prefix, r = prefix;
                    l.cost += c_[left];
                    l.left += a_[left];
                    r.cost += c_[right];
                    r.right += a_[right];
                    r.mask = mask | bit;
                    records[mask] = l;
                    records[mask | bit] = r;
                }
            }
            for (unsigned mask = 0; mask < (1U << (end - begin)); ++mask) {
                if (expired()) return;
                const auto r = records[mask];
                if (excess(agents[0], base[0] + r.left) == 0 &&
                    excess(agents[1], base[1] + r.right) == 0)
                    records[sizes[half]++] = r;
            }
        }
        auto &first = pair_records_[0];
        auto &second = pair_records_[1];
        std::sort(first.begin(), first.begin() + sizes[0],
                  [](const PairRecord &a, const PairRecord &b) { return a.left < b.left; });
        if (expired(true)) return;
        std::sort(second.begin(), second.begin() + sizes[1],
                  [](const PairRecord &a, const PairRecord &b) { return a.left > b.left; });
        if (expired(true)) return;
        for (int k = 0; k < sizes[0]; ++k) pair_coordinates_[k] = first[k].right;
        std::sort(pair_coordinates_.begin(), pair_coordinates_.begin() + sizes[0]);
        if (expired(true)) return;
        std::fill(pair_fenwick_.begin(), pair_fenwick_.end(), -1);
        int inserted = 0;
        for (int t = 0; t < sizes[1]; ++t) {
            if (expired()) return;
            const auto &r = second[t];
            while (inserted < sizes[0] &&
                   excess(agents[0], base[0] + first[inserted].left + r.left) == 0) {
                if (expired()) return;
                int k = static_cast<int>(std::lower_bound(pair_coordinates_.begin(),
                                                          pair_coordinates_.begin() + sizes[0],
                                                          first[inserted].right) -
                                         pair_coordinates_.begin()) +
                        1;
                for (; k <= sizes[0]; k += k & -k) {
                    const int old = pair_fenwick_[k];
                    if (old < 0 || first[inserted].cost < first[old].cost)
                        pair_fenwick_[k] = inserted;
                }
                ++inserted;
            }
            int lo = 0, hi = sizes[0];
            while (lo < hi) {
                const int mid = lo + (hi - lo) / 2;
                if (excess(agents[1], base[1] + pair_coordinates_[mid] + r.right) == 0)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            int best = -1;
            for (int k = lo; k > 0; k -= k & -k) {
                const int candidate = pair_fenwick_[k];
                if (candidate >= 0 && (best < 0 || first[candidate].cost < first[best].cost))
                    best = candidate;
            }
            if (best >= 0 && first[best].cost + r.cost < bound) {
                bound = first[best].cost + r.cost;
                for (int k = 0; k < split; ++k) answer[k] = agents[(first[best].mask >> k) & 1U];
                for (int k = split; k < count; ++k)
                    answer[k] = agents[(r.mask >> (k - split)) & 1U];
            }
        }
    }
    void relink() {
        if (elite_size_ < 2) return;
        const int target = random_index(elite_size_);
        int distance = 0;
        for (int j = 0; j < n_; ++j) {
            if (expired()) return;
            distance += x_[j] != elite_[target][j];
        }
        for (int step = 0; step < (distance / 2 + distance % 2); ++step) {
            int job = -1;
            long double delta = infinity();
            for (int j = 0; j < n_; ++j) {
                if (expired()) return;
                if (x_[j] == elite_[target][j]) continue;
                const auto value = shift_delta(j, elite_[target][j]);
                if (value < delta) {
                    delta = value;
                    job = j;
                }
            }
            if (job < 0) return;
            relocate(job, elite_[target][job]);
            // Do not update the elite set while its target is in use.
        }
    }
    void perturb(bool use_relinking, int repair_nodes) {
        ++restarts_;
        stalled_ = 0;
        // Return toward the incumbent through individual moves; interruption leaves a valid state.
        if (best_.feasible && restarts_ % 3 != 0) {
            for (int j = 0; j < n_; ++j) {
                if (expired()) return;
                relocate(j, best_.assignment[j]);
            }
            for (int k = 0; k < 16 && repair_nodes > 0; ++k) {
                repair_subset(repair_nodes);
                if (timed_out_) return;
                save_best();
            }
        }
        if (use_relinking && restarts_ % 2 == 0) {
            relink();
            if (timed_out_) return;
        }
        for (int j = n_ - 1; j > 0; --j) {
            if (expired()) return;
            std::swap(order_[j], order_[random_index(j + 1)]);
        }
        const int count =
            restarts_ % 7 == 0 ? n_ : std::min(n_, std::max(3, n_ / (restarts_ % 3 == 0 ? 3 : 8)));
        for (int k = 0; k < count; ++k) {
            if (expired()) return;
            relocate(restarts_ % 7 == 0 ? k : random_index(n_), -1);
        }
        filling_cursor_ = 0;
    }
};

} // namespace gap
} // namespace cpopt
#endif
