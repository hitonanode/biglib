#ifndef QUADRATIC_ASSIGNMENT_HPP
#define QUADRATIC_ASSIGNMENT_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace quadratic_assignment {
using Clock = std::chrono::steady_clock;
using Cost = double;
enum class Search { hybrid, breakout, tabu };
enum class StopReason { time_limit, iteration_limit, exhaustive };

struct Options {
    double time_limit = 1.0;
    Clock::time_point deadline = Clock::time_point::max();
    double safety_margin = 0.001;
    std::uint64_t seed = 1;
    std::uint64_t max_iterations = std::numeric_limits<std::uint64_t>::max(); // Per run.
    Search search = Search::hybrid;
    int exact_threshold = 8; // 0..10; exhaustive search below this size.
    double perturbation_fraction = 0.15;
    int stagnation_optima = 100;  // Strong breakout perturbation threshold.
    int hybrid_epoch = 12;        // n * this many swaps per search mode.
    int restart_stagnation = 100; // n * this many swaps without improvement; 0 disables.
    double relative_tolerance = 1e-12;
    std::vector<int> initial_assignment; // Empty, or n entries; -1 means unassigned.
};

struct Statistics {
    double elapsed = 0; // This call. Other counters are cumulative.
    std::uint64_t iterations = 0, improvements = 0, local_optima = 0;
    std::uint64_t perturbations = 0, perturbation_moves = 0, mode_switches = 0;
    std::uint64_t delta_evaluations = 0, delta_updates = 0, permutations = 0;
    std::uint64_t rebuilds = 0, restarts = 0;
    bool symmetric = false, delta_ready = false, batched_initialization = false;
};

struct Result {
    std::vector<int> assignment;
    Cost cost = 0;
    StopReason reason = StopReason::time_limit;
    Statistics stats;
};

namespace detail {
inline Clock::time_point deadline(Clock::time_point start, double seconds, Clock::time_point outer) {
    if (std::isnan(seconds) || seconds < 0)
        throw std::invalid_argument("quadratic_assignment: invalid time_limit");
    const double range = std::chrono::duration<double>(Clock::time_point::max() - start).count();
    if (seconds < range)
        outer = std::min(outer, start + std::chrono::duration_cast<Clock::duration>(
                                            std::chrono::duration<double>(seconds)));
    return outer;
}
struct Budget {
    Clock::time_point stop;
    Budget(Clock::time_point start, double seconds, Clock::time_point outer, double margin) {
        stop = deadline(start, seconds, outer);
        if (stop <= start) {
            stop = start;
            return;
        }
        if (stop != Clock::time_point::max()) {
            const double remaining = std::chrono::duration<double>(stop - start).count();
            stop -= std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(std::min(margin, remaining * 0.1)));
        }
    }
    bool expired() const { return stop != Clock::time_point::max() && Clock::now() >= stop; }
};
inline int dimension(std::size_t n) {
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        (n && n > std::numeric_limits<std::size_t>::max() / n / sizeof(Cost)))
        throw std::length_error("quadratic_assignment: matrix too large");
    return static_cast<int>(n);
}
} // namespace detail

// Owns its input matrices. Copying Solver is a complete in-memory checkpoint.
class Solver {
    enum class Phase { build, scan, update, verify, reset_history, exact, done };
    static constexpr Cost inf = std::numeric_limits<Cost>::infinity();
    int n_;
    Options options_;
    std::vector<Cost> a_, at_, b_, bt_, delta_, gradient_;
    std::vector<Cost> permuted_b_, permuted_bt_, row_minimum_;
    std::vector<std::int32_t> narrow_a_, narrow_at_, narrow_b_, narrow_bt_;
    std::vector<int> p_, best_, last_optimum_;
    std::vector<std::uint64_t> tabu_, history_;
    std::vector<Cost> outgoing_, incoming_, location_out_, location_in_, cross_sum_;
    Cost current_ = 0, best_cost_ = 0, absolute_sum_ = 0;
    std::uint64_t random_, last_improvement_ = 0, epoch_start_ = 0, last_restart_ = 0;
    Statistics stats_;
    Phase phase_ = Phase::build, after_verify_ = Phase::scan;
    bool symmetric_ = false, tabu_mode_ = false, have_optimum_ = false, cross_update_ = false;
    bool exact_objective_ = false;
    int row_ = 0, column_ = 1, swap_r_ = -1, swap_s_ = -1;
    int gradient_row_ = 0, gradient_k_ = 0;
    int permutation_row_ = 0;
    int candidate_r_ = -1, candidate_s_ = -1, fallback_r_ = -1, fallback_s_ = -1;
    Cost candidate_delta_ = inf, fallback_delta_ = inf;
    std::uint64_t candidate_age_ = 0, candidate_ties_ = 0;
    int perturb_left_ = 0, perturb_kind_ = 0, jump_ = 0, stalled_optima_ = 0;
    std::uint64_t optimum_improvements_ = 0;
    int verify_row_ = 0;
    long double verify_sum_ = 0, verify_correction_ = 0, verify_absolute_ = 0;

    // IEEE binary64 values have an order-preserving integer key after flipping
    // the sign bit (nonnegative) or all bits (negative). Integer min reductions
    // can be vectorized without changing any floating-point arithmetic.
    using Minimum =
        std::conditional_t<sizeof(Cost) == sizeof(std::uint64_t) && std::numeric_limits<Cost>::is_iec559 &&
                               std::numeric_limits<Cost>::digits == 53,
                           std::uint64_t, Cost>;
    static Minimum minimum_key(Cost value) {
        if constexpr (std::is_same_v<Minimum, Cost>)
            return value;
        else {
            std::uint64_t bits;
            std::memcpy(&bits, &value, sizeof bits);
            return bits ^ ((0 - (bits >> 63)) | (UINT64_C(1) << 63));
        }
    }
    static Cost minimum_value(Minimum key) {
        if constexpr (std::is_same_v<Minimum, Cost>)
            return key;
        else {
            std::uint64_t bits = key;
            bits ^= ((bits >> 63) - 1) | (UINT64_C(1) << 63);
            Cost value;
            std::memcpy(&value, &bits, sizeof value);
            return value;
        }
    }
    std::size_t ix(int i, int j) const { return static_cast<std::size_t>(i) * n_ + j; }
    std::uint64_t random() {
        std::uint64_t z = (random_ += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    std::uint64_t bounded(std::uint64_t bound) { return random() % bound; }
    double unit() { return (random() >> 11) * (1.0 / 9007199254740992.0); }
    double tolerance() const { return options_.relative_tolerance * std::max(1.0, absolute_sum_); }
    void validate_options() const {
        if (std::isnan(options_.time_limit) || options_.time_limit < 0 ||
            !std::isfinite(options_.safety_margin) || options_.safety_margin < 0 ||
            !std::isfinite(options_.relative_tolerance) || options_.relative_tolerance < 0 ||
            !std::isfinite(options_.perturbation_fraction) || options_.perturbation_fraction <= 0 ||
            options_.perturbation_fraction > 1 || options_.stagnation_optima < 1 ||
            options_.hybrid_epoch < 1 || options_.restart_stagnation < 0 ||
            options_.exact_threshold < 0 || options_.exact_threshold > 10 ||
            (options_.search != Search::hybrid && options_.search != Search::breakout &&
             options_.search != Search::tabu))
            throw std::invalid_argument("quadratic_assignment: invalid options");
    }
    void validate_permutation(const std::vector<int> &p, bool partial) const {
        if (p.size() != static_cast<std::size_t>(n_))
            throw std::invalid_argument("quadratic_assignment: assignment size mismatch");
        std::vector<unsigned char> seen(n_, 0);
        for (int x : p) {
            if (partial && x == -1) continue;
            if (x < 0 || x >= n_ || seen[x])
                throw std::invalid_argument("quadratic_assignment: invalid or duplicate location");
            seen[x] = 1;
        }
    }
    // Used only when sum |A[i,j] B[p[i],p[j]]| <= 2^48 and all stored
    // entries are integers. Every product and partial sum is then exact.
    template <bool Permuted, class Entry>
    std::pair<Cost, Cost> exact_row(const Cost *a, const Entry *b, const std::vector<int> &p) const {
        Cost s0 = 0, s1 = 0, s2 = 0, s3 = 0, a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        const auto location = [&](int j) { return Permuted ? j : p[j]; };
        int j = 0;
        for (; j <= n_ - 4; j += 4) {
            const Cost t0 = a[j] * b[location(j)], t1 = a[j + 1] * b[location(j + 1)];
            const Cost t2 = a[j + 2] * b[location(j + 2)], t3 = a[j + 3] * b[location(j + 3)];
            s0 += t0;
            s1 += t1;
            s2 += t2;
            s3 += t3;
            a0 += std::abs(t0);
            a1 += std::abs(t1);
            a2 += std::abs(t2);
            a3 += std::abs(t3);
        }
        for (; j < n_; ++j) {
            const Cost t = a[j] * b[location(j)];
            s0 += t;
            a0 += std::abs(t);
        }
        return {((s0 + s1) + s2) + s3, ((a0 + a1) + a2) + a3};
    }
    Cost value(const std::vector<int> &p, Cost *absolute = nullptr) const {
        if (exact_objective_) {
            Cost sum = 0, abs_sum = 0;
            for (int i = 0; i < n_; ++i) {
                const auto row = exact_row<false>(a_.data() + ix(i, 0), b_.data() + ix(p[i], 0), p);
                sum += row.first;
                abs_sum += row.second;
            }
            if (absolute) *absolute = abs_sum;
            return sum;
        }
        long double sum = 0, correction = 0, abs_sum = 0;
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < n_; ++j) {
                const long double term = static_cast<long double>(a_[ix(i, j)]) * b_[ix(p[i], p[j])];
                const long double y = term - correction, next = sum + y;
                correction = (next - sum) - y;
                sum = next;
                abs_sum += std::abs(term);
            }
        if (absolute) *absolute = static_cast<Cost>(abs_sum);
        return static_cast<Cost>(sum);
    }

    // O(n), valid for asymmetric matrices, arbitrary diagonals and signed entries.
    Cost input_delta(int i, int j) const {
        const int x = p_[i], y = p_[j];
        const Cost *ai = a_.data() + ix(i, 0);
        const Cost *aj = a_.data() + ix(j, 0);
        const Cost *bx = b_.data() + ix(x, 0);
        const Cost *by = b_.data() + ix(y, 0);
        Cost d = (ai[i] - aj[j]) * (by[y] - bx[x]);
        if (symmetric_) {
            for (int k = 0; k < n_; ++k)
                if (k != i && k != j) d += 2 * (ai[k] - aj[k]) * (by[p_[k]] - bx[p_[k]]);
        } else {
            d += (ai[j] - aj[i]) * (by[x] - bx[y]);
            const Cost *ati = at_.data() + ix(i, 0);
            const Cost *atj = at_.data() + ix(j, 0);
            const Cost *btx = bt_.data() + ix(x, 0);
            const Cost *bty = bt_.data() + ix(y, 0);
            for (int k = 0; k < n_; ++k)
                if (k != i && k != j)
                    d += (ai[k] - aj[k]) * (by[p_[k]] - bx[p_[k]]) +
                         (ati[k] - atj[k]) * (bty[p_[k]] - btx[p_[k]]);
        }
        return d;
    }
    // Keep the input formula's term order and endpoint branches for floating
    // arithmetic. Only the row layout changes; no reassociation is used here.
    Cost floating_cached_delta(int i, int j) const {
        const int x = i, y = j;
        const Cost *ai = a_.data() + ix(i, 0);
        const Cost *aj = a_.data() + ix(j, 0);
        const Cost *bx = permuted_b_.data() + ix(x, 0);
        const Cost *by = permuted_b_.data() + ix(y, 0);
        Cost d = (ai[i] - aj[j]) * (by[y] - bx[x]);
        if (symmetric_) {
            for (int k = 0; k < n_; ++k)
                if (k != i && k != j) d += 2 * (ai[k] - aj[k]) * (by[k] - bx[k]);
        } else {
            d += (ai[j] - aj[i]) * (by[x] - bx[y]);
            const Cost *ati = at_.data() + ix(i, 0);
            const Cost *atj = at_.data() + ix(j, 0);
            const Cost *btx = permuted_bt_.data() + ix(x, 0);
            const Cost *bty = permuted_bt_.data() + ix(y, 0);
            for (int k = 0; k < n_; ++k)
                if (k != i && k != j)
                    d += (ai[k] - aj[k]) * (by[k] - bx[k]) + (ati[k] - atj[k]) * (bty[k] - btx[k]);
        }
        return d;
    }
    Cost cached_direct_delta(int i, int j) const {
        const int x = i, y = j;
        const Cost *ai = a_.data() + ix(i, 0);
        const Cost *aj = a_.data() + ix(j, 0);
        const Cost *data = permuted_b_.data();
        const Cost *bx = data + ix(x, 0);
        const Cost *by = data + ix(y, 0);
        // This path is restricted to n * max|A| * max|B| <= 2^48 integers.
        // The int64_t reduction permits vectorization without fast-math. Each
        // term/sum is exact; including the endpoint subtraction, intermediate
        // magnitudes are at most 20n max|A| max|B| < 2^53.
        const auto accumulate = [&](auto term) {
            std::int64_t sum = 0;
            for (int k = 0; k < n_; ++k) sum += static_cast<std::int64_t>(term(k));
            return static_cast<Cost>(sum) - term(i) - term(j);
        };
        Cost d = (ai[i] - aj[j]) * (by[y] - bx[x]);
        if (symmetric_) {
            d += accumulate([&](int k) { return 2 * (ai[k] - aj[k]) * (by[k] - bx[k]); });
        } else {
            d += (ai[j] - aj[i]) * (by[x] - bx[y]);
            const Cost *ati = at_.data() + ix(i, 0);
            const Cost *atj = at_.data() + ix(j, 0);
            const Cost *transposed = permuted_bt_.data();
            const Cost *btx = transposed + ix(x, 0);
            const Cost *bty = transposed + ix(y, 0);
            d += accumulate([&](int k) {
                return (ai[k] - aj[k]) * (by[k] - bx[k]) + (ati[k] - atj[k]) * (bty[k] - btx[k]);
            });
        }
        return d;
    }
    Cost narrow_direct_delta(int i, int j) const {
        const auto *ai = narrow_a_.data() + ix(i, 0);
        const auto *aj = narrow_a_.data() + ix(j, 0);
        const auto *bi = narrow_b_.data() + ix(i, 0);
        const auto *bj = narrow_b_.data() + ix(j, 0);
        // Both differences fit int32_t. Widen before multiplying; all sums
        // are bounded by the existing n max|A| max|B| <= 2^48 guard.
        const auto accumulate = [&](auto term) {
            std::int64_t sum = 0;
            for (int k = 0; k < n_; ++k) sum += term(k);
            return sum - term(i) - term(j);
        };
        std::int64_t d = static_cast<std::int64_t>(ai[i] - aj[j]) * (bj[j] - bi[i]);
        if (symmetric_) {
            d += 2 * accumulate([&](int k) {
                     return static_cast<std::int64_t>(ai[k] - aj[k]) * (bj[k] - bi[k]);
                 });
        } else {
            d += static_cast<std::int64_t>(ai[j] - aj[i]) * (bj[i] - bi[j]);
            const auto *ati = narrow_at_.data() + ix(i, 0);
            const auto *atj = narrow_at_.data() + ix(j, 0);
            const auto *bti = narrow_bt_.data() + ix(i, 0);
            const auto *btj = narrow_bt_.data() + ix(j, 0);
            d += accumulate([&](int k) {
                return static_cast<std::int64_t>(ai[k] - aj[k]) * (bj[k] - bi[k]) +
                       static_cast<std::int64_t>(ati[k] - atj[k]) * (btj[k] - bti[k]);
            });
        }
        return static_cast<Cost>(d);
    }
    Cost direct_delta(int i, int j) const {
        if (!narrow_a_.empty()) return narrow_direct_delta(i, j);
        return gradient_.empty() ? floating_cached_delta(i, j) : cached_direct_delta(i, j);
    }
    Cost minimum_row(int i) const {
        const Cost *values = delta_.data() + ix(i, 0);
        Cost result = inf;
        for (int j = i + 1; j < n_; ++j) result = std::fmin(result, values[j]);
        return result;
    }
    Cost cached_delta(int i, int j) const { return i < j ? delta_[ix(i, j)] : delta_[ix(j, i)]; }
    void reset_scan() {
        row_ = 0;
        column_ = 1;
        candidate_r_ = fallback_r_ = -1;
        candidate_delta_ = fallback_delta_ = inf;
        candidate_age_ = 0;
        candidate_ties_ = 0;
    }
    void begin_verify(Phase next) {
        verify_row_ = 0;
        verify_sum_ = verify_correction_ = verify_absolute_ = 0;
        after_verify_ = next;
        phase_ = Phase::verify;
    }
    bool verify(detail::Budget &budget) {
        for (; verify_row_ < n_; ++verify_row_) {
            if (budget.expired()) return false;
            const Cost *a = a_.data() + ix(verify_row_, 0);
            const Cost *b = b_.data() + ix(p_[verify_row_], 0);
            if (exact_objective_) {
                const auto row =
                    permutation_row_ != n_ ? exact_row<false>(a, b, p_)
                    : !narrow_a_.empty()
                        ? exact_row<true>(a, narrow_b_.data() + ix(verify_row_, 0), p_)
                        : exact_row<true>(a, permuted_b_.data() + ix(verify_row_, 0), p_);
                verify_sum_ += row.first;
                verify_absolute_ += row.second;
                continue;
            }
            long double sum = verify_sum_, correction = verify_correction_,
                        absolute = verify_absolute_;
            const auto accumulate = [&](const auto *values, auto location) {
                for (int j = 0; j < n_; ++j) {
                    const long double t = static_cast<long double>(a[j]) * values[location(j)];
                    const long double y = t - correction, next = sum + y;
                    correction = (next - sum) - y;
                    sum = next;
                    absolute += std::abs(t);
                }
            };
            if (permutation_row_ == n_ && !narrow_a_.empty())
                accumulate(narrow_b_.data() + ix(verify_row_, 0), [](int j) { return j; });
            else if (permutation_row_ == n_)
                accumulate(permuted_b_.data() + ix(verify_row_, 0), [](int j) { return j; });
            else
                accumulate(b, [&](int j) { return p_[j]; });
            verify_sum_ = sum;
            verify_correction_ = correction;
            verify_absolute_ = absolute;
        }
        current_ = static_cast<Cost>(verify_sum_);
        absolute_sum_ = static_cast<Cost>(verify_absolute_);
        if (current_ < best_cost_) {
            best_cost_ = current_;
            best_ = p_;
            ++stats_.improvements;
            last_improvement_ = stats_.iterations;
            stalled_optima_ = 0;
        }
        phase_ = after_verify_;
        if (phase_ == Phase::build) gradient_row_ = gradient_k_ = 0;
        return true;
    }
    // Adjacent output rows reuse the same B row. The template also handles the
    // final odd row without accessing beyond either matrix.
    template <int Rows> bool gradient_rows(detail::Budget &budget) {
        Cost *g = gradient_.data() + ix(gradient_row_, 0);
        Cost *h = g + (Rows == 2 ? n_ : 0);
        if (gradient_k_ == 0) std::fill_n(g, Rows * static_cast<std::size_t>(n_), 0);
        for (; gradient_k_ < n_; ++gradient_k_) {
            if ((gradient_k_ & 3) == 0 && budget.expired()) return false;
            const int k = gradient_k_;
            const Cost *out = b_.data() + ix(p_[k], 0);
            const Cost a0 = a_[ix(gradient_row_, k)];
            const Cost a1 = Rows == 2 ? a_[ix(gradient_row_ + 1, k)] : 0;
            if (symmetric_) {
                const Cost c0 = 2 * a0, c1 = 2 * a1;
                if (c1 == 0) {
                    if (c0 != 0)
                        for (int l = 0; l < n_; ++l) g[l] += c0 * out[l];
                } else if (c0 == 0) {
                    for (int l = 0; l < n_; ++l) h[l] += c1 * out[l];
                } else {
                    for (int l = 0; l < n_; ++l) {
                        g[l] += c0 * out[l];
                        h[l] += c1 * out[l];
                    }
                }
            } else {
                const Cost *in = bt_.data() + ix(p_[k], 0);
                const Cost t0 = at_[ix(gradient_row_, k)];
                const Cost t1 = Rows == 2 ? at_[ix(gradient_row_ + 1, k)] : 0;
                if (a1 == 0 && t1 == 0) {
                    if (a0 != 0 || t0 != 0)
                        for (int l = 0; l < n_; ++l) g[l] += a0 * in[l] + t0 * out[l];
                } else if (a0 == 0 && t0 == 0) {
                    for (int l = 0; l < n_; ++l) h[l] += a1 * in[l] + t1 * out[l];
                } else {
                    for (int l = 0; l < n_; ++l) {
                        g[l] += a0 * in[l] + t0 * out[l];
                        h[l] += a1 * in[l] + t1 * out[l];
                    }
                }
            }
        }
        gradient_row_ += Rows;
        gradient_k_ = 0;
        return true;
    }
    bool build(detail::Budget &budget) {
        // B[p[i],p[k]] in facility order makes every direct delta a sequential
        // read. The cache is updated by row/column swaps in O(n) after each move.
        for (; permutation_row_ < n_; ++permutation_row_) {
            if (budget.expired()) return false;
            const int i = permutation_row_;
            for (int k = 0; k < n_; ++k) {
                if (!narrow_a_.empty()) {
                    narrow_b_[ix(i, k)] = static_cast<std::int32_t>(b_[ix(p_[i], p_[k])]);
                    if (!symmetric_)
                        narrow_bt_[ix(i, k)] = static_cast<std::int32_t>(bt_[ix(p_[i], p_[k])]);
                } else {
                    permuted_b_[ix(i, k)] = b_[ix(p_[i], p_[k])];
                    if (!symmetric_) permuted_bt_[ix(i, k)] = bt_[ix(p_[i], p_[k])];
                }
            }
        }
        // For bounded integer entries these sums are exact in binary64. The
        // contiguous inner products avoid a gather for every swap and facility.
        if (!gradient_.empty()) {
            while (gradient_row_ < n_ - 1)
                if (!gradient_rows<2>(budget)) return false;
            if (gradient_row_ < n_ && !gradient_rows<1>(budget)) return false;
        }
        unsigned checks = 0;
        for (; row_ < n_; ++row_, column_ = row_ + 1) {
            for (; column_ < n_; ++column_) {
                if ((checks++ & 7U) == 0 && budget.expired()) return false;
                if (gradient_.empty())
                    delta_[ix(row_, column_)] = direct_delta(row_, column_);
                else {
                    const int i = row_, j = column_, x = p_[i], y = p_[j];
                    delta_[ix(i, j)] =
                        gradient_[ix(i, y)] + gradient_[ix(j, x)] - gradient_[ix(i, x)] -
                        gradient_[ix(j, y)] +
                        (a_[ix(i, i)] + a_[ix(j, j)] - a_[ix(i, j)] - a_[ix(j, i)]) *
                            (b_[ix(x, x)] + b_[ix(y, y)] - b_[ix(x, y)] - b_[ix(y, x)]);
                }
                ++stats_.delta_evaluations;
            }
            row_minimum_[row_] = minimum_row(row_);
        }
        stats_.delta_ready = true;
        phase_ = Phase::scan;
        reset_scan();
        return true;
    }
    void accept(int r, int s, bool perturbation) {
        const int u = p_[r], v = p_[s];
        const Cost change = delta_[ix(r, s)];
        // The correction for every disjoint swap is the sum of two outer products.
        for (int i = 0; i < n_; ++i) {
            outgoing_[i] = a_[ix(i, r)] - a_[ix(i, s)];
            location_out_[i] = b_[ix(p_[i], v)] - b_[ix(p_[i], u)];
            if (!symmetric_) {
                incoming_[i] = a_[ix(r, i)] - a_[ix(s, i)];
                location_in_[i] = b_[ix(v, p_[i])] - b_[ix(u, p_[i])];
            }
            if (i != r && i != s) {
                cross_sum_[i] = cached_delta(r, i) + cached_delta(s, i) - change;
                if (!symmetric_) {
                    const int w = p_[i];
                    const Cost ha = a_[ix(r, s)] - a_[ix(r, i)] - a_[ix(s, r)] + a_[ix(s, i)] +
                                    a_[ix(i, r)] - a_[ix(i, s)];
                    const Cost hb = b_[ix(u, v)] - b_[ix(u, w)] - b_[ix(v, u)] + b_[ix(v, w)] +
                                    b_[ix(w, u)] - b_[ix(w, v)];
                    cross_sum_[i] += ha * hb;
                }
            }
        }
        ++stats_.iterations;
        const auto tenure = static_cast<std::uint64_t>(
            std::max(1.0, n_ * (tabu_mode_ ? 0.2 + 1.6 * unit() : 0.9 + 0.2 * unit())));
        tabu_[ix(r, u)] = stats_.iterations + tenure;
        tabu_[ix(s, v)] = stats_.iterations + tenure;
        history_[ix(r, s)] = stats_.iterations;
        std::swap(p_[r], p_[s]);
        current_ += change;
        for (int k = 0; !permuted_b_.empty() && k < n_; ++k) {
            std::swap(permuted_b_[ix(r, k)], permuted_b_[ix(s, k)]);
            if (!symmetric_) std::swap(permuted_bt_[ix(r, k)], permuted_bt_[ix(s, k)]);
        }
        for (int k = 0; !permuted_b_.empty() && k < n_; ++k) {
            std::swap(permuted_b_[ix(k, r)], permuted_b_[ix(k, s)]);
            if (!symmetric_) std::swap(permuted_bt_[ix(k, r)], permuted_bt_[ix(k, s)]);
        }
        if (!narrow_a_.empty()) {
            for (int k = 0; k < n_; ++k) {
                std::swap(narrow_b_[ix(r, k)], narrow_b_[ix(s, k)]);
                if (!symmetric_) std::swap(narrow_bt_[ix(r, k)], narrow_bt_[ix(s, k)]);
            }
            for (int k = 0; k < n_; ++k) {
                std::swap(narrow_b_[ix(k, r)], narrow_b_[ix(k, s)]);
                if (!symmetric_) std::swap(narrow_bt_[ix(k, r)], narrow_bt_[ix(k, s)]);
            }
        }
        delta_[ix(r, s)] = -change;
        // Rebuild row minima alongside the two stages of the table update.
        // No scan observes this incomplete state, including after a resume.
        std::fill(row_minimum_.begin(), row_minimum_.end(), inf);
        row_minimum_[r] = -change;
        swap_r_ = r;
        swap_s_ = s;
        cross_update_ = false;
        row_ = 0;
        column_ = 1;
        phase_ = Phase::update;
        stats_.delta_ready = false;
        if (perturbation) {
            --perturb_left_;
            ++stats_.perturbation_moves;
        }
        // Validate every prospective record, and periodically remove objective drift.
        if (current_ < best_cost_ + tolerance() || (stats_.iterations & 255U) == 0)
            begin_verify(Phase::update);
    }
    bool update(detail::Budget &budget) {
        // Long loops are interruptible. p_ cannot change until the whole table is ready.
        if (!cross_update_) {
            for (; row_ < n_; ++row_) {
                if (budget.expired()) return false;
                const int i = row_;
                if (i == swap_r_ || i == swap_s_) continue;
                Cost *values = delta_.data() + ix(i, 0);
                const Cost *out = outgoing_.data();
                const Cost *loc = location_out_.data();
                const Cost oi = out[i], li = loc[i];
                Minimum minimum = minimum_key(inf);
                // Cross entries are recomputed below. Masking them first lets
                // the disjoint update and minimum use one contiguous loop.
                if (swap_r_ > i) values[swap_r_] = inf;
                if (swap_s_ > i) values[swap_s_] = inf;
                const auto span = [&](int first, int last) {
                    if (symmetric_) {
                        for (int j = first; j < last; ++j) {
                            values[j] += 2 * ((oi - out[j]) * (loc[j] - li));
                            minimum = std::min(minimum, minimum_key(values[j]));
                        }
                    } else {
                        const Cost *in = incoming_.data();
                        const Cost *lin = location_in_.data();
                        const Cost ii = in[i], lni = lin[i];
                        for (int j = first; j < last; ++j) {
                            values[j] +=
                                (oi - out[j]) * (loc[j] - li) + (ii - in[j]) * (lin[j] - lni);
                            minimum = std::min(minimum, minimum_key(values[j]));
                        }
                    }
                };
                span(i + 1, n_);
                stats_.delta_updates +=
                    static_cast<std::uint64_t>(n_ - i - 1 - (swap_r_ > i) - (swap_s_ > i));
                row_minimum_[i] = minimum_value(minimum);
            }
            cross_update_ = true;
            row_ = 0;
        }
        for (; row_ < n_; ++row_) {
            if ((row_ & 7) == 0 && budget.expired()) return false;
            const int k = row_;
            if (k == swap_r_ || k == swap_s_) continue;
            const Cost d = direct_delta(std::min(swap_r_, k), std::max(swap_r_, k));
            delta_[ix(std::min(swap_r_, k), std::max(swap_r_, k))] = d;
            delta_[ix(std::min(swap_s_, k), std::max(swap_s_, k))] = cross_sum_[k] - d;
            const int first = std::min(swap_r_, k), second = std::min(swap_s_, k);
            row_minimum_[first] = std::fmin(row_minimum_[first], d);
            row_minimum_[second] = std::fmin(row_minimum_[second], cross_sum_[k] - d);
            ++stats_.delta_evaluations;
            ++stats_.delta_updates;
        }
        reset_scan();
        stats_.delta_ready = true;
        phase_ = Phase::scan;
#ifdef QUADRATIC_ASSIGNMENT_CHECK_INVARIANTS
        check_invariants();
#endif
        // Rebuild occasionally to bound accumulated floating-point cache drift.
        if ((stats_.iterations & 4095U) == 0) {
            phase_ = Phase::build;
            stats_.delta_ready = false;
            ++stats_.rebuilds;
            gradient_row_ = gradient_k_ = 0;
        }
        return true;
    }
    void start_perturbation() {
        ++stats_.local_optima;
        ++stats_.perturbations;
        if (stats_.improvements == optimum_improvements_)
            ++stalled_optima_;
        else
            stalled_optima_ = 0;
        optimum_improvements_ = stats_.improvements;
        const int base = std::max(2, static_cast<int>(n_ * options_.perturbation_fraction));
        if (have_optimum_ && p_ == last_optimum_)
            jump_ = std::min(n_, jump_ + 1);
        else
            jump_ = base;
        last_optimum_ = p_;
        have_optimum_ = true;
        if (stalled_optima_ >= options_.stagnation_optima) {
            jump_ = std::max(2, static_cast<int>(n_ * (0.4 + 0.2 * unit())));
            stalled_optima_ = 0;
            perturb_kind_ = 2;
        } else {
            const double directed = std::max(
                0.75, std::exp(-static_cast<double>(stalled_optima_) / options_.stagnation_optima));
            perturb_kind_ = unit() < directed ? 0 : (unit() < 0.7 ? 1 : 2);
        }
        perturb_left_ = jump_;
        reset_scan();
    }
    void restart() {
        // Preserve the incumbent, but leave a basin that has resisted both modes.
        p_ = best_;
        for (int i = n_; i > 1; --i) std::swap(p_[i - 1], p_[bounded(i)]);
        permutation_row_ = 0;
        perturb_left_ = 0;
        have_optimum_ = false;
        tabu_mode_ = options_.search == Search::tabu;
        epoch_start_ = last_restart_ = stats_.iterations;
        reset_scan();
        stats_.delta_ready = false;
        ++stats_.restarts;
        phase_ = Phase::reset_history;
    }
    bool reset_history(detail::Budget &budget) {
        for (; row_ < n_; ++row_) {
            if (budget.expired()) return false;
            std::fill_n(tabu_.data() + ix(row_, 0), n_, 0);
            std::fill_n(history_.data() + ix(row_, 0), n_, 0);
        }
        reset_scan();
        begin_verify(Phase::build);
        return true;
    }
    bool scan(detail::Budget &budget) {
        if (row_ == 0 && column_ == 1 && !perturb_left_ && options_.restart_stagnation &&
            stats_.iterations - std::max(last_improvement_, last_restart_) >=
                static_cast<std::uint64_t>(options_.restart_stagnation) * n_) {
            restart();
            return true;
        }
        if (row_ == 0 && column_ == 1 && options_.search == Search::hybrid && !perturb_left_ &&
            stats_.iterations - epoch_start_ >=
                static_cast<std::uint64_t>(options_.hybrid_epoch) * n_) {
            tabu_mode_ = !tabu_mode_;
            epoch_start_ = stats_.iterations;
            ++stats_.mode_switches;
        }
        if (!tabu_mode_ && perturb_left_ && perturb_kind_ == 2) {
            int r = static_cast<int>(bounded(n_)), s = static_cast<int>(bounded(n_ - 1));
            if (s >= r) ++s;
            if (r > s) std::swap(r, s);
            accept(r, s, true);
            return true;
        }
        const bool recency = !tabu_mode_ && perturb_left_ && perturb_kind_ == 1;
        const bool directed = tabu_mode_ || perturb_left_;
        const Cost eps = tolerance();
        unsigned checks = 0;
        for (; row_ < n_; ++row_, column_ = row_ + 1) {
            const int i = row_;
            // A whole row strictly worse than an allowed candidate cannot
            // affect either candidate or fallback. Equal values still run
            // the original tie-breaking code, in the original order.
            if (row_minimum_[i] > candidate_delta_) {
                if ((checks++ & 255U) == 0 && budget.expired()) return false;
                continue;
            }
            for (; column_ < n_; ++column_) {
                if ((checks++ & 255U) == 0 && budget.expired()) return false;
                const int j = column_;
                const Cost d = delta_[ix(i, j)];
                // fallback_delta_ <= candidate_delta_. A worse delta cannot
                // change either choice or consume tie-breaking randomness.
                // Recency selection keeps candidate_delta_ at infinity.
                if (d > candidate_delta_) continue;
                if (d < fallback_delta_) {
                    fallback_delta_ = d;
                    fallback_r_ = i;
                    fallback_s_ = j;
                }
                if (recency) {
                    const auto age = history_[ix(i, j)];
                    if (candidate_r_ < 0 || age < candidate_age_) {
                        candidate_age_ = age;
                        candidate_r_ = i;
                        candidate_s_ = j;
                        candidate_ties_ = 1;
                    } else if (age == candidate_age_ && bounded(++candidate_ties_) == 0) {
                        candidate_r_ = i;
                        candidate_s_ = j;
                    }
                    continue;
                }
                bool allowed = true;
                if (directed) {
                    if (tabu_mode_)
                        allowed = tabu_[ix(i, p_[j])] <= stats_.iterations ||
                                  tabu_[ix(j, p_[i])] <= stats_.iterations;
                    else
                        allowed = history_[ix(i, j)] == 0 || stats_.iterations - history_[ix(i, j)] >
                                                                 static_cast<std::uint64_t>(n_);
                    allowed = allowed || current_ + d < best_cost_ - eps;
                }
                if (allowed && d < candidate_delta_) {
                    candidate_delta_ = d;
                    candidate_r_ = i;
                    candidate_s_ = j;
                    candidate_ties_ = 1;
                } else if (allowed && d == candidate_delta_ && bounded(++candidate_ties_) == 0) {
                    candidate_r_ = i;
                    candidate_s_ = j;
                }
            }
        }
        if (!directed && candidate_delta_ >= -eps) {
            start_perturbation();
            return true;
        }
        // On tiny/degenerate neighborhoods all moves can be tabu. Relax only then.
        if (candidate_r_ < 0) {
            candidate_r_ = fallback_r_;
            candidate_s_ = fallback_s_;
        }
        accept(candidate_r_, candidate_s_, !tabu_mode_ && perturb_left_ > 0);
        return true;
    }
    bool exact(detail::Budget &budget) {
        if (budget.expired()) return false;
        current_ = value(p_);
        ++stats_.permutations;
        ++stats_.iterations;
        if (current_ < best_cost_) {
            best_cost_ = current_;
            best_ = p_;
            ++stats_.improvements;
        }
        if (!std::next_permutation(p_.begin(), p_.end())) {
            p_ = best_;
            phase_ = Phase::done;
        }
        return true;
    }
    Result result(StopReason reason, Clock::time_point start) const {
        Result r;
        r.assignment = best_;
        r.cost = best_cost_;
        r.reason = reason;
        r.stats = stats_;
        r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        return r;
    }

public:
    template <class Flow, class Distance>
    Solver(int n, Flow flow, Distance distance, const Options &options = {})
        : n_(n), options_(options), random_(options.seed) {
        if (n < 0) throw std::invalid_argument("quadratic_assignment: negative dimension");
        detail::dimension(static_cast<std::size_t>(n));
        validate_options();
        const std::size_t count = static_cast<std::size_t>(n) * n;
        a_.resize(count);
        b_.resize(count);
        at_.resize(count);
        bt_.resize(count);
        delta_.resize(count);
        tabu_.resize(count);
        history_.resize(count);
        p_.resize(n);
        best_.resize(n);
        last_optimum_.resize(n);
        outgoing_.resize(n);
        incoming_.resize(n);
        location_out_.resize(n);
        location_in_.resize(n);
        cross_sum_.resize(n);
        row_minimum_.resize(n);
        const Cost limit = std::sqrt(std::numeric_limits<Cost>::max()) / (32.0 * std::max(n, 1));
        bool integral = true;
        Cost max_a = 0, max_b = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const Cost a = static_cast<Cost>(flow(i, j)), b = static_cast<Cost>(distance(i, j));
                if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a) > limit ||
                    std::abs(b) > limit)
                    throw std::invalid_argument(
                        "quadratic_assignment: non-finite or excessive matrix entry");
                a_[ix(i, j)] = at_[ix(j, i)] = a;
                b_[ix(i, j)] = bt_[ix(j, i)] = b;
                integral = integral && std::trunc(a) == a && std::trunc(b) == b;
                max_a = std::max(max_a, std::abs(a));
                max_b = std::max(max_b, std::abs(b));
            }
        symmetric_ = a_ == at_ && b_ == bt_;
        if (integral && (max_b == 0 || max_a <= (281474976710656.0 / std::max(n, 1)) / max_b)) {
            gradient_.resize(count);
            // A signed difference of any two entries must still fit int32_t.
            constexpr Cost bound = std::numeric_limits<std::int32_t>::max() / 2;
            if (max_a <= bound && max_b <= bound) {
                narrow_a_.assign(a_.begin(), a_.end());
                narrow_b_.resize(count);
                if (!symmetric_) {
                    narrow_at_.assign(at_.begin(), at_.end());
                    narrow_bt_.resize(count);
                }
            }
        }
        if (narrow_a_.empty()) {
            permuted_b_.resize(count);
            if (!symmetric_) permuted_bt_.resize(count);
        }
        exact_objective_ =
            integral && (max_b == 0 ||
                         max_a <= ((281474976710656.0 / std::max(n, 1)) / std::max(n, 1)) / max_b);
        set_start(options_.initial_assignment);
    }
    template <class A, class B>
    Solver(const std::vector<std::vector<A>> &a, const std::vector<std::vector<B>> &b,
           const Options &options = {})
        : Solver(
              checked_size(a, b), [&a](int i, int j) { return a[i][j]; },
              [&b](int i, int j) { return b[i][j]; }, options) {}

    template <class A, class B>
    static int
    checked_size(const std::vector<std::vector<A>> &a, const std::vector<std::vector<B>> &b) {
        if (a.size() != b.size())
            throw std::invalid_argument("quadratic_assignment: matrix size mismatch");
        const int n = detail::dimension(a.size());
        for (int i = 0; i < n; ++i)
            if (a[i].size() != a.size() || b[i].size() != a.size())
                throw std::invalid_argument("quadratic_assignment: matrix must be square");
        return n;
    }
    void set_start(const std::vector<int> &assignment) {
        // Validate before changing any persistent state.
        if (!assignment.empty()) validate_permutation(assignment, true);
        std::vector<int> next(n_, -1), free_locations;
        std::vector<unsigned char> used(n_, 0);
        if (!assignment.empty()) next = assignment;
        for (int x : next)
            if (x >= 0) used[x] = 1;
        for (int x = 0; x < n_; ++x)
            if (!used[x]) free_locations.push_back(x);
        for (std::size_t i = free_locations.size(); i > 1; --i)
            std::swap(free_locations[i - 1], free_locations[bounded(i)]);
        std::size_t k = 0;
        for (int &x : next)
            if (x < 0) x = free_locations[k++];
        p_ = next;
        best_ = next;
        current_ = best_cost_ = value(p_, &absolute_sum_);
        std::fill(tabu_.begin(), tabu_.end(), 0);
        std::fill(history_.begin(), history_.end(), 0);
        stats_ = Statistics{};
        stats_.symmetric = symmetric_;
        stats_.batched_initialization = !gradient_.empty();
        last_improvement_ = epoch_start_ = optimum_improvements_ = last_restart_ = 0;
        perturb_left_ = stalled_optima_ = 0;
        jump_ = 2;
        have_optimum_ = false;
        tabu_mode_ = options_.search == Search::tabu;
        reset_scan();
        phase_ = Phase::build;
        gradient_row_ = gradient_k_ = 0;
        permutation_row_ = 0;
        if (n_ <= 1)
            phase_ = Phase::done;
        else if (n_ <= options_.exact_threshold) {
            std::iota(p_.begin(), p_.end(), 0);
            phase_ = Phase::exact;
        }
    }
    int size() const { return n_; }
    Cost evaluate(const std::vector<int> &assignment) const {
        validate_permutation(assignment, false);
        return value(assignment);
    }
    const std::vector<int> &current_assignment() const { return p_; }
    Result best() const {
        return result(phase_ == Phase::done ? StopReason::exhaustive : StopReason::iteration_limit,
                      Clock::now());
    }
    Result run() { return run(options_.time_limit); }
    Result run(double seconds, Clock::time_point outer_deadline = Clock::time_point::max()) {
        const auto start = Clock::now();
        detail::Budget budget(
            start, seconds, std::min(options_.deadline, outer_deadline), options_.safety_margin);
        const auto initial_iterations = stats_.iterations;
        while (phase_ != Phase::done) {
            if (budget.expired()) return result(StopReason::time_limit, start);
            if ((phase_ == Phase::scan || phase_ == Phase::build || phase_ == Phase::exact) &&
                stats_.iterations - initial_iterations >= options_.max_iterations)
                return result(StopReason::iteration_limit, start);
            bool complete = true;
            switch (phase_) {
            case Phase::build: complete = build(budget); break;
            case Phase::scan: complete = scan(budget); break;
            case Phase::update: complete = update(budget); break;
            case Phase::verify: complete = verify(budget); break;
            case Phase::reset_history: complete = reset_history(budget); break;
            case Phase::exact: complete = exact(budget); break;
            case Phase::done: break;
            }
            if (!complete) return result(StopReason::time_limit, start);
        }
        return result(StopReason::exhaustive, start);
    }
    Result run_until(Clock::time_point deadline) {
        return run(std::numeric_limits<double>::infinity(), deadline);
    }
#ifdef QUADRATIC_ASSIGNMENT_CHECK_INVARIANTS
    Cost debug_delta(int i, int j) const {
        return stats_.delta_ready ? cached_delta(i, j) : input_delta(i, j);
    }
    void check_invariants() const {
        validate_permutation(p_, false);
        validate_permutation(best_, false);
        if (value(best_) != best_cost_)
            throw std::logic_error("quadratic_assignment: best objective mismatch");
        if (stats_.delta_ready)
            for (int i = 0; i < n_; ++i)
                for (int j = i + 1; j < n_; ++j) {
                    const Cost expected = input_delta(i, j), actual = delta_[ix(i, j)];
                    if (row_minimum_[i] > actual)
                        throw std::logic_error("quadratic_assignment: invalid row minimum");
                    if (std::abs(expected - actual) >
                        1e-8 * std::max({1.0, std::abs(expected), absolute_sum_}))
                        throw std::logic_error("quadratic_assignment: delta cache mismatch");
                }
    }
#endif
};

template <class Flow, class Distance>
inline Solver make_solver(int n, Flow flow, Distance distance, const Options &options = {}) {
    return Solver(n, std::move(flow), std::move(distance), options);
}
template <class A, class B>
inline Solver make_solver(const std::vector<std::vector<A>> &a,
                          const std::vector<std::vector<B>> &b, const Options &options = {}) {
    return Solver(a, b, options);
}

template <class Flow, class Distance>
inline Result solve(int n, Flow flow, Distance distance, const Options &options = {}) {
    const auto start = Clock::now();
    const auto stop = detail::deadline(start, options.time_limit, options.deadline);
    auto solver = make_solver(n, std::move(flow), std::move(distance), options);
    auto r = solver.run(options.time_limit, stop);
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    return r;
}
template <class A, class B>
inline Result solve(const std::vector<std::vector<A>> &a, const std::vector<std::vector<B>> &b,
                    const Options &options = {}) {
    const auto start = Clock::now();
    const auto stop = detail::deadline(start, options.time_limit, options.deadline);
    auto solver = make_solver(a, b, options);
    auto r = solver.run(options.time_limit, stop);
    r.stats.elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    return r;
}
} // namespace quadratic_assignment
#endif
