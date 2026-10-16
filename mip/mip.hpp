#ifndef mip_HPP
#define mip_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cpopt {

inline constexpr double inf = std::numeric_limits<double>::infinity();
enum class Status { Optimal, Infeasible, Unbounded, Limit, NumericalError, RelaxationUnbounded };
enum class Sense { Minimize, Maximize };
enum class VarType { Continuous, Integer, Binary };
enum class Relation { LessEqual, Equal, GreaterEqual };
using Terms = std::vector<std::pair<int, double>>;

inline const char *status_name(Status s) {
    switch (s) {
    case Status::Optimal: return "optimal";
    case Status::Infeasible: return "infeasible";
    case Status::Unbounded: return "unbounded";
    case Status::Limit: return "limit";
    case Status::NumericalError: return "numerical_error";
    case Status::RelaxationUnbounded: return "relaxation_unbounded";
    }
    return "unknown";
}

struct Options {
    double time_limit = inf; // 秒。前処理を含む solve() ごとの経過時間。
    std::uint64_t node_limit = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t iteration_limit = std::numeric_limits<std::uint64_t>::max();
    double feasibility_tolerance = 1e-7;
    double integrality_tolerance = 1e-7;
    double pivot_tolerance = 1e-11;
    double absolute_gap = 1e-7;
    double relative_gap = 0;
    bool presolve = true;
    bool scaling = true;
    bool dualize = true; // 連続 LP で有利な形なら双対問題を解く。
    bool warm_start = true;
    bool reuse_root = true; // 同じ Model の再求解で、最適な根基底を保存・再利用。
    bool pseudocost = true;
    bool rounding_heuristic = true;
    bool integer_objective = true;
    bool reduced_cost_fixing = true; // 暫定解を改善できない非基底整数変数を固定・列削除。
    int root_cut_rounds = -1;        // -1: 小数の整数変数が1個のときだけ3巡、0: 無効。
    int cuts_per_round = 8;
    int strong_branching_candidates = 4;
    int strong_branching_depth = 3;
    int strong_branching_iterations = 40;
    std::vector<double> initial_solution; // 元の変数順。実行可能性を検査してから使う。
};

struct Statistics {
    std::uint64_t nodes = 0, pivots = 0, strong_branching_pivots = 0;
    std::uint64_t heuristic_solutions = 0;
    int cuts = 0;
    int reduced_cost_fixings = 0;
    int numerical_restarts = 0;
    bool solved_dual = false;
    bool reused_root = false;
    std::uint64_t root_pivots = 0; // 根の初回求解・再最適化・復旧の総 pivot 数。
    int original_variables = 0, original_constraints = 0;
    int presolved_variables = 0, presolved_constraints = 0;
    double seconds = 0;
};

struct Result {
    Status status = Status::NumericalError;
    bool has_solution = false;
    double objective = std::numeric_limits<double>::quiet_NaN();
    // 最大化なら最適値の上界、最小化なら下界。解がなくても有効。
    double best_bound = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> x;
    Statistics stats;
    double gap() const { return has_solution ? std::abs(best_bound - objective) : inf; }
};

namespace detail {
// 部分的に構築した前処理結果を誤って利用せず、保持済みの暫定解へ戻る。
struct TimeLimitReached {};
struct Variable {
    double lower, upper;
    VarType type;
};
struct Row {
    Terms terms;
    double rhs;
}; // 常に <=

struct Budget {
    const Options &opt;
    Statistics &stats;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    static constexpr std::size_t check_interval = 4096;
    const bool timed = std::isfinite(opt.time_limit);
    std::size_t work_left = 0;
    // 係数ごとには時計を読まない。処理量を積算し、約4096要素ごとに読む。
    void checkpoint(std::size_t work = check_interval) {
        if (!timed) return;
        if (work < work_left) {
            work_left -= work;
            return;
        }
        work_left = check_interval;
        if (elapsed() >= opt.time_limit) throw TimeLimitReached{};
    }
    // 密なループはブロック内に条件分岐を入れず、自動ベクトル化を維持する。
    template <class F> void chunks(std::size_t size, F f) {
        if (!timed) {
            f(0, size);
            return;
        }
        for (std::size_t begin = 0; begin < size; begin += check_interval) {
            const auto end = std::min(size, begin + check_interval);
            checkpoint(end - begin);
            f(begin, end);
        }
    }
    template <class F> static bool scan(Budget *budget, std::size_t size, F f) {
        if (!budget || !budget->timed) return f(0, size);
        for (std::size_t begin = 0; begin < size; begin += check_interval) {
            const auto end = std::min(size, begin + check_interval);
            budget->checkpoint(end - begin);
            if (!f(begin, end)) return false;
        }
        return true;
    }
    template <class T> std::vector<T> copy(const std::vector<T> &source) {
        if (!timed) return source;
        checkpoint(1);
        std::vector<T> result;
        result.reserve(source.size());
        chunks(source.size(), [&](std::size_t begin, std::size_t end) {
            result.insert(result.end(), source.begin() + begin, source.begin() + end);
        });
        return result;
    }
    template <class T> std::vector<T> filled(std::size_t size, const T &value = T{}) {
        if (!timed) return std::vector<T>(size, value);
        checkpoint(1);
        std::vector<T> result;
        result.reserve(size);
        chunks(size, [&](std::size_t begin, std::size_t end) {
            result.insert(result.end(), end - begin, value);
        });
        return result;
    }
    bool stop() const {
        return stats.pivots >= opt.iteration_limit ||
               (std::isfinite(opt.time_limit) && elapsed() >= opt.time_limit);
    }
};

// 辞書: x_B = rhs - D x_N, z = rhs_z - D_z x_N。
// 初期列数は元変数数 + 人工変数 + RHS。固定列は後から削除できる。
// slack の列を別に持たず交換する。
struct Tableau {
    int n = 0, m = 0, stride = 0, next_id = 0, original_n = 0;
    static constexpr int fixed_zero = std::numeric_limits<int>::max();
    bool primal_feasible = false;
    double objective_scale = 1;
    std::vector<double> dense;
    std::vector<std::shared_ptr<std::vector<double>>> data;
    std::vector<int> basic, nonbasic;
    std::vector<int> location; // fixed_zero: 0固定、それ以外は非負: 列、負: -行-1。人工変数を除く。
    std::vector<int> pivot_columns;   // 疎な pivot 行の非ゼロ列。作業領域なので複製不要。
    std::vector<double *> pivot_rows; // detach 済みの行。pivot 中だけ有効。
    const double *row(int i) const {
        return data.empty() ? dense.data() + std::size_t(i) * stride : data[i]->data();
    }
    double *write_row(int i, Budget *budget = nullptr) {
        if (data.empty()) return dense.data() + std::size_t(i) * stride;
        if (data[i].use_count() != 1)
            data[i] =
                std::make_shared<std::vector<double>>(budget ? budget->copy(*data[i]) : *data[i]);
        return data[i]->data();
    }
    double value() const { return row(m)[n + 1] / objective_scale; }
    void scale_objective(const std::vector<double> &c) {
        double largest = 0;
        for (double v : c) largest = std::max(largest, std::abs(v));
        int e = 0;
        if (largest) std::frexp(largest, &e);
        // 小さすぎる目的係数だけ拡大する。大きい目的を縮めると、例えば
        // 1e12*x + y の y の寄与が pivot 許容誤差より小さくなってしまう。
        objective_scale = std::ldexp(1.0, -std::max(-1000, std::min(0, e)));
    }

    Tableau() = default;
    Tableau(const Tableau &other) { *this = other; }
    Tableau &operator=(const Tableau &other) {
        if (this == &other) return *this;
        n = other.n;
        m = other.m;
        stride = other.stride;
        next_id = other.next_id;
        original_n = other.original_n;
        objective_scale = other.objective_scale;
        primal_feasible = other.primal_feasible;
        // 複製直後の分枝行1本で再確保しない。値のコピーと行追加で二重に
        // 配列をコピーするコストを避ける（探索・分枝の順序は同じ）。
        if (!other.data.empty()) data.reserve(other.data.size() + 1);
        data = other.data;
        dense.reserve(other.dense.size() + (other.dense.empty() ? 0 : std::size_t(stride) * 2));
        dense.assign(other.dense.begin(), other.dense.end());
        basic.reserve(other.basic.size() + 1);
        location.reserve(other.location.size() + 1);
        basic = other.basic;
        nonbasic = other.nonbasic;
        location = other.location;
        return *this;
    }
    // 再求解の保存・復元時は、既存の不変スナップショットを壊さず複製する。
    Tableau clone(Budget &budget) const {
        Tableau t;
        t.n = n;
        t.m = m;
        t.stride = stride;
        t.next_id = next_id;
        t.original_n = original_n;
        t.primal_feasible = primal_feasible;
        t.objective_scale = objective_scale;
        budget.checkpoint();
        t.data = budget.copy(data);
        t.dense.reserve(dense.size() + (dense.empty() ? 0 : std::size_t(stride) * 2));
        budget.chunks(dense.size(), [&](std::size_t begin, std::size_t end) {
            t.dense.insert(t.dense.end(), dense.begin() + begin, dense.begin() + end);
        });
        t.basic = budget.copy(basic);
        t.nonbasic = budget.copy(nonbasic);
        t.location = budget.copy(location);
        return t;
    }
    Tableau(Tableau &&) = default;
    Tableau &operator=(Tableau &&) = default;
    static void copy_coefficients(const std::vector<double> &a, double *dest, Budget *budget) {
        auto copy = [&](std::size_t begin, std::size_t end) {
            std::copy(a.begin() + begin, a.begin() + end, dest + begin);
        };
        if (budget)
            budget->chunks(a.size(), copy);
        else
            copy(0, a.size());
    }
    static void copy_coefficients(const Terms &a, double *dest, Budget *budget) {
        auto copy = [&](std::size_t begin, std::size_t end) {
            for (auto k = begin; k < end; ++k) dest[a[k].first] = a[k].second;
        };
        if (budget)
            budget->chunks(a.size(), copy);
        else
            copy(0, a.size());
    }
    template <class Matrix>
    Tableau(const Matrix &a, const std::vector<double> &b, const std::vector<double> &c,
            Budget *budget = nullptr)
        : n(int(c.size())), m(int(b.size())), stride(n + 2), next_id(n + m), original_n(n),
          basic(m), nonbasic(n + 1), location(n + m) {
        if (stride < 128)
            dense = budget ? budget->filled<double>(std::size_t(m + 2) * stride)
                           : std::vector<double>(std::size_t(m + 2) * stride);
        else {
            data.reserve(m + 2);
            for (int i = 0; i < m + 2; ++i)
                data.push_back(std::make_shared<std::vector<double>>(
                    budget ? budget->filled<double>(stride) : std::vector<double>(stride)));
        }
        for (int i = 0; i < m; ++i) {
            double *q = write_row(i);
            copy_coefficients(a[i], q, budget);
            q[n] = -1;
            q[n + 1] = b[i];
            basic[i] = n + i;
            location[n + i] = -i - 1;
        }
        scale_objective(c);
        double *z = write_row(m);
        for (int j = 0; j < n; ++j) {
            nonbasic[j] = j;
            location[j] = j;
            z[j] = -c[j] * objective_scale;
        }
        nonbasic[n] = -1;
        write_row(m + 1)[n] = 1;
    }

    void pivot(int r, int s, Budget &budget) {
        double *p = write_row(r, &budget);
        pivot_rows.clear();
        if (!data.empty())
            for (int i = 0; i < m + 2; ++i)
                if (i != r && row(i)[s] != 0) { pivot_rows.push_back(write_row(i, &budget)); }
        const double inv = 1 / p[s];
        // フロー等では辞書を密に保存していても pivot 行は疎なことが多い。
        // 密な行は早めに走査を打ち切り、従来の連続ループを使う。
        pivot_columns.clear();
        const int sparse_limit = stride >= 64 ? stride / 4 : 0;
        if (sparse_limit) {
            for (int j = 0; j < stride; ++j)
                if (p[j] != 0) {
                    pivot_columns.push_back(j);
                    if (int(pivot_columns.size()) > sparse_limit) break;
                }
        }
        const bool sparse = sparse_limit && int(pivot_columns.size()) <= sparse_limit;
        // 連続メモリでの rank-one 更新。ゼロ行は丸ごと省く。
        auto update = [&](double *q) {
            const double f = q[s] * inv;
            if (f == 0) return;
            if (sparse) {
                for (int j : pivot_columns) q[j] -= f * p[j];
            } else {
                for (int j = 0; j < stride; ++j) q[j] -= f * p[j];
            }
            q[s] = -f;
        };
        // 格納方式の分岐と共有ポインタの参照は、数値更新の内側に持ち込まない。
        if (data.empty()) {
            for (double *q = dense.data(); q != p; q += stride) update(q);
            double *end = dense.data() + dense.size();
            for (double *q = p + stride; q != end; q += stride) update(q);
        } else
            for (double *q : pivot_rows) update(q);
        for (int j = 0; j < stride; ++j) p[j] *= inv;
        p[s] = inv;
        std::swap(basic[r], nonbasic[s]);
        if (basic[r] >= 0) location[basic[r]] = -r - 1;
        if (nonbasic[s] >= 0) location[nonbasic[s]] = s;
        ++budget.stats.pivots;
    }

    Status primal(int objective_row, Budget &budget) {
        const double eps = budget.opt.pivot_tolerance;
        int degenerate = 0;
        bool bland = false;
        while (true) {
            if (budget.stop()) return Status::Limit;
            const double *z = row(objective_row);
            int s = -1;
            for (int j = 0; j <= n; ++j) {
                if (nonbasic[j] == -1 || z[j] >= -eps) continue;
                if (s < 0 || (bland ? nonbasic[j] < nonbasic[s] : z[j] < z[s])) s = j;
            }
            if (s < 0) return Status::Optimal;
            int r = -1;
            double ratio = inf;
            for (int i = 0; i < m; ++i) {
                const double *q = row(i);
                if (q[s] <= eps) continue;
                const double v = std::max(0.0, q[n + 1]) / q[s];
                if (v < ratio - eps ||
                    (std::abs(v - ratio) <= eps && (r < 0 || basic[i] < basic[r]))) {
                    r = i;
                    ratio = v;
                }
            }
            if (r < 0) return Status::Unbounded;
            if (ratio <= eps)
                ++degenerate;
            else {
                degenerate = 0;
                bland = false;
            }
            if (degenerate > 32) bland = true; // 退化が続いたら循環を避ける。
            pivot(r, s, budget);
            if (!std::isfinite(value())) return Status::NumericalError;
        }
    }

    Status solve(Budget &budget) {
        if (budget.stop()) return Status::Limit;
        int r = -1;
        for (int i = 0; i < m; ++i)
            if (r < 0 || row(i)[n + 1] < row(r)[n + 1]) r = i;
        if (r >= 0 && row(r)[n + 1] < -budget.opt.pivot_tolerance) {
            // 最小費用フローや被覆では、slack 基底が既に双対実行可能。
            // 人工変数で全行を混ぜる Phase I を省き、直接実行可能性を回復する。
            bool dual_feasible = true;
            for (int j = 0; j <= n; ++j)
                if (nonbasic[j] != -1 && row(m)[j] < 0) dual_feasible = false;
            if (dual_feasible) return dual(budget);
            pivot(r, n, budget);
            Status s = primal(m + 1, budget);
            if (s == Status::Limit || s == Status::NumericalError) return s;
            if (s != Status::Optimal) return Status::NumericalError;
            if (row(m + 1)[n + 1] < -budget.opt.pivot_tolerance) return Status::Infeasible;
            for (int i = 0; i < m; ++i)
                if (basic[i] == -1) {
                    int jbest = -1;
                    for (int j = 0; j <= n; ++j)
                        if (nonbasic[j] != -1 && std::abs(row(i)[j]) > budget.opt.pivot_tolerance &&
                            (jbest < 0 || std::abs(row(i)[j]) > std::abs(row(i)[jbest])))
                            jbest = j;
                    if (jbest >= 0) pivot(i, jbest, budget);
                }
        }
        // 人工変数が非基底なら選ばない。冗長行の基底に残った場合も値は 0。
        primal_feasible = true;
        return primal(m, budget);
    }

    std::vector<double> solution() const {
        std::vector<double> x(original_n);
        for (int j = 0; j < original_n; ++j)
            if (location[j] < 0) x[j] = row(-location[j] - 1)[n + 1];
        return x;
    }

    // 正の reduced cost と整数の最小移動幅で、改善解では0と証明できる列を除く。
    // 元変数/slackのIDは保ち、辞書の物理的な列数だけを縮める。
    int fix_nonbasic(const std::vector<double> &step, double gap, Budget &budget,
                     int minimum_removed = 1) {
        // 多くのノードでは削除しない。pivot の作業領域に候補だけを記録し、
        // 閾値に達してから location と残す列の配列を作る。
        pivot_columns.clear();
        const double margin = budget.opt.feasibility_tolerance * (1 + std::abs(value()));
        const double *z = row(m);
        for (int j = 0; j <= n; ++j) {
            budget.checkpoint(1);
            const int id = nonbasic[j];
            if (id >= 0 && step[id] > 0 && z[j] / objective_scale * step[id] > gap + margin)
                pivot_columns.push_back(j);
        }
        const int removed = int(pivot_columns.size());
        if (removed < minimum_removed) return 0;
        auto positions = budget.copy(location);
        std::vector<int> keep;
        keep.reserve(n + 1 - removed);
        int k = 0;
        for (int j = 0; j <= n; ++j) {
            budget.checkpoint(1);
            if (k < removed && pivot_columns[k] == j) {
                positions[nonbasic[j]] = fixed_zero;
                ++k;
            } else
                keep.push_back(j);
        }
        const int width = int(keep.size()) + 1;
        std::vector<std::shared_ptr<std::vector<double>>> compact;
        std::vector<double> compact_dense;
        if (width < 128)
            compact_dense = budget.filled<double>(std::size_t(m + 2) * width);
        else
            compact.reserve(m + 2);
        for (int i = 0; i < m + 2; ++i) {
            auto values = width < 128 ? std::vector<double>{} : budget.filled<double>(width);
            double *dest =
                width < 128 ? compact_dense.data() + std::size_t(i) * width : values.data();
            const double *source = row(i);
            budget.chunks(keep.size(), [&](std::size_t begin, std::size_t end) {
                for (auto j = begin; j < end; ++j) dest[j] = source[keep[j]];
            });
            dest[width - 1] = source[n + 1];
            if (width >= 128)
                compact.push_back(std::make_shared<std::vector<double>>(std::move(values)));
        }
        std::vector<int> columns;
        columns.reserve(keep.size());
        for (int j : keep) {
            const int id = nonbasic[j];
            if (id >= 0) positions[id] = int(columns.size());
            columns.push_back(id);
        }
        data = std::move(compact);
        dense = std::move(compact_dense);
        nonbasic = std::move(columns);
        location = std::move(positions);
        n -= removed;
        stride = width;
        return removed;
    }

    // 最適辞書に a*x <= rhs を追加。目的行の双対実行可能性は保存される。
    template <class Coefficients = Terms>
    void add_constraint(const Coefficients &a, double rhs, Budget *budget = nullptr,
                        double multiplier = 1) {
        primal_feasible = false;
        // 時間割込みのない狭い辞書では、追加先に直接組み立てて一時行を省く。
        // budget ありの場合は従来通り、完成した行だけを辞書へ追加する。
        const bool inplace = data.empty() && !budget;
        std::vector<double> buffer;
        double *q;
        if (inplace) {
            dense.resize(std::size_t(m + 3) * stride);
            std::copy_backward(write_row(m), write_row(m + 2), write_row(m + 3));
            q = write_row(m);
            std::fill(q, q + stride, 0);
        } else {
            buffer.resize(stride);
            q = buffer.data();
        }
        q[n + 1] = rhs;
        for (const auto &term : a) {
            if (budget) budget->checkpoint(1);
            const int var = term.first;
            const double coef = term.second * multiplier;
            const int pos = location[var];
            if (pos == fixed_zero) continue;
            if (pos >= 0)
                q[pos] += coef;
            else {
                const double *p = row(-pos - 1);
                auto subtract = [&](std::size_t begin, std::size_t end) {
                    for (auto j = begin; j < end; ++j) q[j] -= coef * p[j];
                };
                if (budget)
                    budget->chunks(stride, subtract);
                else
                    subtract(0, stride);
            }
        }
        if (data.empty() && !inplace) {
            dense.resize(std::size_t(m + 3) * stride);
            std::copy_backward(write_row(m), write_row(m + 2), write_row(m + 3));
            std::copy(buffer.begin(), buffer.end(), write_row(m));
        } else if (!data.empty())
            data.insert(data.begin() + m, std::make_shared<std::vector<double>>(std::move(buffer)));
        basic.push_back(next_id++);
        location.push_back(-m - 1);
        ++m;
    }

    Status dual(Budget &budget, int max_pivots = std::numeric_limits<int>::max()) {
        const double eps = budget.opt.pivot_tolerance;
        int iter = 0;
        int degenerate = 0;
        bool bland = false;
        while (true) {
            if (budget.stop()) return Status::Limit;
            int r = -1;
            for (int i = 0; i < m; ++i)
                if (row(i)[n + 1] < -eps &&
                    (r < 0 || (bland ? basic[i] < basic[r] : row(i)[n + 1] < row(r)[n + 1])))
                    r = i;
            if (r < 0) {
                primal_feasible = true;
                return Status::Optimal;
            }
            if (iter++ >= max_pivots) return Status::Limit;
            int s = -1;
            double ratio = inf;
            const double *q = row(r);
            const double *z = row(m);
            for (int j = 0; j <= n; ++j) {
                if (nonbasic[j] == -1 || q[j] >= -eps) continue;
                const double v = std::max(0.0, z[j]) / -q[j];
                if (v < ratio - eps ||
                    (std::abs(v - ratio) <= eps && (s < 0 || nonbasic[j] < nonbasic[s]))) {
                    s = j;
                    ratio = v;
                }
            }
            if (s < 0) return Status::Infeasible;
            if (ratio <= eps)
                ++degenerate;
            else {
                degenerate = 0;
                bland = false;
            }
            if (degenerate > 32) bland = true;
            pivot(r, s, budget);
            if (!std::isfinite(value())) return Status::NumericalError;
        }
    }
};

// 元変数は固定・平行移動・反転・自由変数の分割のいずれかで、対応は高々2項。
// 可変長配列を変数ごとに確保せず、根の保存・復元でもまとめてコピーする。
struct Mapping {
    struct Columns {
        std::pair<int, double> entries[2];
        int count = 0;
        auto begin() const { return entries; }
        auto end() const { return entries + count; }
        int size() const { return count; }
        bool empty() const { return count == 0; }
        const auto &operator[](int j) const { return entries[j]; }
        void emplace_back(int j, double a) { entries[count++] = {j, a}; }
    } terms;
    double offset = 0;
};
// 変換済みの係数行は追加後に変更しない。根の保存・復元では行を共有し、
// 追加行・カットのリストだけを独立に持つ。Model のコピーにも変更が漏れない。
struct ImmutableRows {
    std::vector<std::shared_ptr<const Terms>> rows;
    std::size_t size() const { return rows.size(); }
    void reserve(std::size_t n) { rows.reserve(n); }
    const Terms &operator[](std::size_t i) const { return *rows[i]; }
    void push_terms(Terms row) { rows.push_back(std::make_shared<const Terms>(std::move(row))); }
    void push_back(const std::vector<double> &row) {
        Terms terms;
        for (int j = 0; j < int(row.size()); ++j)
            if (row[j]) terms.emplace_back(j, row[j]);
        push_terms(std::move(terms));
    }
};
struct Prepared {
    std::vector<Variable> vars;
    std::vector<Row> rows;
    std::vector<Mapping> map;
    ImmutableRows a;
    std::vector<double> b, c;
    double offset = 0;
    bool infeasible = false;
    Prepared clone(Budget &budget) const {
        Prepared p;
        p.vars = budget.copy(vars);
        p.map = budget.copy(map);
        p.b = budget.copy(b);
        p.c = budget.copy(c);
        p.offset = offset;
        p.infeasible = infeasible;
        p.a.rows = budget.copy(a.rows);
        // rows は前処理だけで使う作業領域。保存・再求解には不要。
        return p;
    }
};

// Gomory mixed-integer cut。教科書の式から生成し、非基底 slack を元の行で
// 消去して保持する。これにより cold start でも同じカットを再現できる。
inline int gomory_cuts(Tableau &tab, Prepared &p, const std::vector<bool> &integer, Budget &budget) {
    const auto &opt = budget.opt;
    struct Cut {
        Terms terms;
        double score;
    };
    std::vector<Cut> cuts;
    for (int i = 0; i < tab.m; ++i) {
        budget.checkpoint(1);
        int var = tab.basic[i];
        if (var < 0 || var >= tab.original_n || !integer[var]) continue;
        const double *row = tab.row(i);
        const double f = row[tab.n + 1] - std::floor(row[tab.n + 1]);
        if (f < .01 || f > .99) continue;
        Terms terms;
        double norm = 0, largest = 0;
        for (int j = 0; j <= tab.n; ++j) {
            budget.checkpoint(1);
            const int id = tab.nonbasic[j];
            if (id < 0) continue;
            const double a = row[j];
            double alpha;
            if (id < tab.original_n && integer[id]) {
                const double fj = a - std::floor(a);
                alpha = fj <= f ? fj / f : (1 - fj) / (1 - f);
            } else
                alpha = a >= 0 ? a / f : -a / (1 - f);
            if (alpha == 0) continue;
            terms.emplace_back(id, -alpha);
            norm += alpha * alpha;
            largest = std::max(largest, alpha);
        }
        // 極端な係数の行は切らず、通常の分枝で処理する。
        if (!std::isfinite(norm) || largest > 1e6) continue;
        cuts.push_back({std::move(terms), 1 / std::max(1e-12, std::sqrt(norm))});
    }
    std::stable_sort(
        cuts.begin(), cuts.end(), [](const Cut &a, const Cut &b) { return a.score > b.score; });
    int added = 0;
    for (const auto &cut : cuts) {
        budget.checkpoint();
        if (added >= opt.cuts_per_round) break;
        std::vector<long double> a(tab.original_n);
        // 丸め誤差に対してわずかに外側へ緩める。
        long double rhs = -1 + 1e-8;
        for (const auto &term : cut.terms) {
            const int id = term.first;
            const double coef = term.second;
            if (id < tab.original_n)
                a[id] += coef;
            else {
                const int r = id - tab.original_n;
                rhs -= (long double)coef * p.b[r];
                budget.chunks(p.a[r].size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) {
                        auto [j, v] = p.a[r][k];
                        a[j] -= (long double)coef * v;
                    }
                });
            }
        }
        Terms coefficients;
        double max_abs = 0;
        budget.chunks(a.size(), [&](std::size_t begin, std::size_t end) {
            for (auto j = begin; j < end; ++j) {
                const double v = double(a[j]);
                max_abs = std::max(max_abs, std::abs(v));
                if (v) coefficients.emplace_back(int(j), v);
            }
        });
        if (!std::isfinite(double(rhs)) || !std::isfinite(max_abs) || max_abs > 1e8) continue;
        p.a.push_terms(std::move(coefficients));
        p.b.push_back(double(rhs));
        tab.add_constraint(cut.terms, -1 + 1e-8);
        ++added;
    }
    return added;
}

inline double bound_round(double v, bool lower, double tol) {
    return lower ? std::ceil(v - tol) : std::floor(v + tol);
}

// 各行の区間活動量から境界を伝播する。無限境界も扱う。
inline bool propagate(const std::vector<Row> &rows, std::vector<Variable> &vars, Budget &budget,
                      int rounds = 8) {
    const auto &opt = budget.opt;
    for (int pass = 0; pass < rounds; ++pass) {
        bool changed = false;
        for (const Row &r : rows) {
            budget.checkpoint(1);
            long double sum = 0;
            int infinite = 0;
            budget.chunks(r.terms.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) {
                    const auto [j, a] = r.terms[k];
                    const double b = a > 0 ? vars[j].lower : vars[j].upper;
                    if (std::isfinite(b))
                        sum += (long double)a * b;
                    else
                        ++infinite;
                }
            });
            const long double slack =
                opt.feasibility_tolerance + 1e-12L * (std::abs(sum) + std::abs(r.rhs));
            if (!infinite && sum > (long double)r.rhs + slack) return false;
            if (!Budget::scan(&budget, r.terms.size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) {
                        const auto [j, a] = r.terms[k];
                        Variable &v = vars[j];
                        const double old = a > 0 ? v.lower : v.upper;
                        const bool finite = std::isfinite(old);
                        if (infinite - !finite != 0) continue;
                        const long double other = sum - (finite ? (long double)a * old : 0);
                        double b = double(((long double)r.rhs - other) / a);
                        if (!std::isfinite(b)) continue;
                        const bool lower = a < 0;
                        if (v.type != VarType::Continuous)
                            b = bound_round(b, lower, opt.integrality_tolerance);
                        double &target = lower ? v.lower : v.upper;
                        if (lower ? b > target : b < target) {
                            // 各行内での活動量は更新前の境界に対応する。古い境界からも有効な推論。
                            // 微小な連続境界の更新による収束ループを抑える。
                            if (std::isfinite(target) &&
                                std::abs(b - target) <= 1e-10 * (1 + std::abs(b)))
                                continue;
                            target = b;
                            changed = true;
                        }
                        if (v.lower > v.upper + opt.feasibility_tolerance) return false;
                    }
                    return true;
                }))
                return false;
        }
        if (!changed) break;
    }
    return true;
}
} // namespace detail

class Model {
    std::vector<detail::Variable> vars_;
    std::vector<detail::Row> rows_;
    std::vector<double> objective_;
    Sense sense_ = Sense::Minimize;
    double constant_ = 0;
    struct RootCache {
        detail::Prepared prepared;
        detail::Tableau tableau;
        std::size_t rows = 0;
        bool presolve = false, scaling = false, dualize = false, relax = false;
        double feasibility = 0, integrality = 0, pivot = 0;
        bool matches(const Options &o, bool r) const {
            return presolve == o.presolve && scaling == o.scaling && dualize == o.dualize &&
                   relax == r && feasibility == o.feasibility_tolerance &&
                   integrality == o.integrality_tolerance && pivot == o.pivot_tolerance;
        }
    };
    // Model のコピー間でも読み取り専用。追加行やカットは必ず別の作業用辞書に入れる。
    mutable std::shared_ptr<const RootCache> root_cache_;

    void check_var(int j) const {
        if (j < 0 || j >= int(vars_.size())) throw std::invalid_argument("invalid variable index");
    }
    Terms normalize(Terms terms) const {
        for (auto [j, a] : terms) {
            check_var(j);
            if (!std::isfinite(a)) throw std::invalid_argument("non-finite coefficient");
        }
        std::sort(terms.begin(), terms.end());
        std::size_t size = 0;
        for (auto [j, a] : terms) {
            if (size && terms[size - 1].first == j)
                terms[size - 1].second += a;
            else
                terms[size++] = {j, a};
        }
        terms.resize(size);
        terms.erase(std::remove_if(terms.begin(), terms.end(), [](auto p) { return p.second == 0; }),
                    terms.end());
        for (auto p : terms)
            if (!std::isfinite(p.second)) throw std::invalid_argument("coefficient overflow");
        return terms;
    }

    detail::Prepared prepare(detail::Budget &budget, bool relax) const {
        const auto &opt = budget.opt;
        detail::Prepared p;
        p.vars = budget.copy(vars_);
        if (!budget.timed)
            p.rows = rows_;
        else {
            p.rows.reserve(rows_.size());
            for (const auto &r : rows_) {
                budget.checkpoint(1);
                p.rows.push_back({budget.copy(r.terms), r.rhs});
            }
        }
        for (auto &v : p.vars) {
            budget.checkpoint(1);
            if (relax) v.type = VarType::Continuous;
            if (v.type != VarType::Continuous) {
                v.lower = detail::bound_round(v.lower, true, opt.integrality_tolerance);
                v.upper = detail::bound_round(v.upper, false, opt.integrality_tolerance);
            }
            if (v.lower > v.upper) {
                p.infeasible = true;
                return p;
            }
        }
        if (opt.presolve) {
            // 整数係数だけの行は gcd で RHS を切り下げられる。
            for (auto &r : p.rows) {
                budget.checkpoint(1);
                long long g = 0;
                bool integral = true;
                for (auto [j, a] : r.terms) {
                    budget.checkpoint(1);
                    if (p.vars[j].type == VarType::Continuous || std::abs(a) > 1e9 ||
                        a != std::round(a)) {
                        integral = false;
                        break;
                    }
                    g = std::gcd(g, (long long)std::llround(std::abs(a)));
                }
                if (integral && g > 0) {
                    r.rhs = std::floor(r.rhs / g + opt.integrality_tolerance);
                    budget.chunks(r.terms.size(), [&](std::size_t begin, std::size_t end) {
                        for (auto k = begin; k < end; ++k) r.terms[k].second /= g;
                    });
                }
            }
            if (!detail::propagate(p.rows, p.vars, budget)) {
                p.infeasible = true;
                return p;
            }
        }
        p.map = budget.filled<detail::Mapping>(vars_.size());
        p.c.reserve(vars_.size());
        p.a.reserve(p.rows.size() + p.vars.size());
        p.b.reserve(p.rows.size() + p.vars.size());
        const double sign = sense_ == Sense::Maximize ? 1 : -1;
        p.offset = sign * constant_;
        auto new_column = [&](int j, double coef) {
            p.map[j].terms.emplace_back(int(p.c.size()), coef);
            p.c.push_back(sign * objective_[j] * coef);
        };
        for (int j = 0; j < int(vars_.size()); ++j) {
            budget.checkpoint(1);
            const auto &v = p.vars[j];
            if (v.lower == v.upper)
                p.map[j].offset = v.lower;
            else if (std::isfinite(v.lower)) {
                p.map[j].offset = v.lower;
                new_column(j, 1);
            } else if (std::isfinite(v.upper)) {
                p.map[j].offset = v.upper;
                new_column(j, -1);
            } else {
                new_column(j, 1);
                new_column(j, -1);
            }
            p.offset += sign * objective_[j] * p.map[j].offset;
        }
        auto append = [&](Terms a, double b) {
            double largest = 0;
            budget.chunks(a.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k)
                    largest = std::max(largest, std::abs(a[k].second));
            });
            if (largest == 0) {
                if (b < -opt.feasibility_tolerance) p.infeasible = true;
                return;
            }
            if (opt.scaling) {
                // 2 の冪で行をスケール。係数比を悪化させず丸めも避ける。
                int e = 0;
                std::frexp(largest, &e);
                double scale = std::ldexp(1.0, -std::max(-500, std::min(500, e)));
                budget.chunks(a.size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) a[k].second *= scale;
                });
                b *= scale;
            }
            p.a.push_terms(std::move(a));
            p.b.push_back(b);
        };
        for (const auto &r : p.rows) {
            budget.checkpoint(1);
            long double rhs = r.rhs;
            Terms a;
            a.reserve(r.terms.size());
            long double maximum = 0;
            bool bounded = true;
            budget.chunks(r.terms.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) {
                    const auto [j, coef] = r.terms[k];
                    const auto &v = p.vars[j];
                    double extreme = coef > 0 ? v.upper : v.lower;
                    if (std::isfinite(extreme))
                        maximum += (long double)coef * extreme;
                    else
                        bounded = false;
                    rhs -= (long double)coef * p.map[j].offset;
                    for (auto [column, scale] : p.map[j].terms)
                        a.emplace_back(column, coef * scale);
                }
            });
            if (opt.presolve && bounded && maximum <= (long double)r.rhs) continue;
            append(std::move(a), double(rhs));
        }
        // 非負係数の残存行から直接導ける上限は重複して追加しない。
        // 伝播で得た上限をすべて行にすると、割当 LP 等で n^2 行増えてしまう。
        auto implied_upper = budget.filled<double>(p.c.size(), inf);
        if (opt.presolve)
            for (int i = 0; i < int(p.a.size()); ++i) {
                bool nonnegative = true;
                budget.chunks(p.a[i].size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k)
                        if (p.a[i][k].second < 0) {
                            nonnegative = false;
                            break;
                        }
                });
                if (!nonnegative) continue;
                budget.chunks(p.a[i].size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) {
                        auto [j, v] = p.a[i][k];
                        if (v > 0) implied_upper[j] = std::min(implied_upper[j], p.b[i] / v);
                    }
                });
            }
        for (int j = 0; j < int(p.vars.size()); ++j) {
            budget.checkpoint(1);
            const auto &v = p.vars[j];
            if (v.lower == v.upper || !std::isfinite(v.lower) || !std::isfinite(v.upper)) continue;
            if (implied_upper[p.map[j].terms[0].first] <= v.upper - v.lower) continue;
            Terms a;
            a.emplace_back(p.map[j].terms[0].first, 1);
            append(std::move(a), v.upper - v.lower);
        }
        return p;
    }

    // 古い前処理で得た境界と変数変換は、制約を追加しても有効なまま。
    // 再度の列消去を避け、追加行だけを保存済みの座標に変換する。
    void append_to_root(detail::Prepared &p, detail::Tableau &tab, std::size_t first,
                        detail::Budget &budget) const {
        const auto &opt = budget.opt;
        for (std::size_t i = first; i < rows_.size(); ++i) {
            budget.checkpoint();
            const auto &r = rows_[i];
            long long g = 0;
            bool integral = opt.presolve;
            if (integral)
                for (auto [j, coef] : r.terms) {
                    budget.checkpoint(1);
                    if (p.vars[j].type == VarType::Continuous || std::abs(coef) > 1e9 ||
                        coef != std::round(coef)) {
                        integral = false;
                        break;
                    }
                    g = std::gcd(g, (long long)std::llround(std::abs(coef)));
                }
            const double divisor = integral && g ? double(g) : 1;
            long double rhs =
                integral && g ? std::floor(r.rhs / divisor + opt.integrality_tolerance) : r.rhs;
            Terms a;
            a.reserve(r.terms.size());
            budget.chunks(r.terms.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) {
                    const auto [j, raw] = r.terms[k];
                    const double coef = raw / divisor;
                    rhs -= (long double)coef * p.map[j].offset;
                    for (auto [column, multiplier] : p.map[j].terms)
                        a.emplace_back(column, coef * multiplier);
                }
            });
            double largest = 0;
            budget.chunks(a.size(), [&](std::size_t begin, std::size_t end) {
                for (auto j = begin; j < end; ++j)
                    largest = std::max(largest, std::abs(a[j].second));
            });
            if (!largest) {
                if (rhs < -opt.feasibility_tolerance) {
                    p.infeasible = true;
                    return;
                }
                continue;
            }
            double scale = 1;
            if (opt.scaling) {
                int e = 0;
                std::frexp(largest, &e);
                scale = std::ldexp(1.0, -std::max(-500, std::min(500, e)));
            }
            budget.chunks(a.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) a[k].second *= scale;
            });
            const double b = double(rhs) * scale;
            tab.add_constraint(a, b, &budget);
            p.a.push_terms(std::move(a));
            p.b.push_back(b);
        }
    }

    void save_root(detail::Prepared &p, detail::Tableau &tab, detail::Budget &budget, bool relax,
                   bool consume) const {
        budget.checkpoint();
        auto next = std::make_shared<RootCache>();
        if (consume) {
            // この根で最適性を証明済みなら、直後に返すので辞書の複製は不要。
            next->prepared = std::move(p);
            next->prepared.rows.clear();
            next->tableau = std::move(tab);
        } else {
            next->prepared = p.clone(budget);
            next->tableau = tab.clone(budget);
        }
        next->rows = rows_.size();
        const auto &o = budget.opt;
        next->presolve = o.presolve;
        next->scaling = o.scaling;
        next->dualize = o.dualize;
        next->relax = relax;
        next->feasibility = o.feasibility_tolerance;
        next->integrality = o.integrality_tolerance;
        next->pivot = o.pivot_tolerance;
        root_cache_ = std::move(next); // 完全なスナップショットができてから差し替える。
    }

    std::vector<double> recover(const detail::Prepared &p, const std::vector<double> &y) const {
        std::vector<double> x(vars_.size());
        for (int j = 0; j < int(x.size()); ++j) {
            x[j] = p.map[j].offset;
            for (auto [k, a] : p.map[j].terms) x[j] += a * y[k];
        }
        return x;
    }

    // 元変数を直接復元し、変換後の全列ぶんの一時ベクトルを作らない。
    std::vector<double> recover(const detail::Prepared &p, const detail::Tableau &tab) const {
        std::vector<double> x(vars_.size());
        for (int j = 0; j < int(x.size()); ++j) {
            x[j] = p.map[j].offset;
            for (auto [k, a] : p.map[j].terms) {
                const int pos = tab.location[k];
                const double y = pos < 0 ? tab.row(-pos - 1)[tab.n + 1] : 0;
                x[j] += a * y;
            }
        }
        return x;
    }

    double objective_value(const std::vector<double> &x) const {
        long double value = constant_;
        for (int j = 0; j < int(x.size()); ++j) value += (long double)objective_[j] * x[j];
        return double(value);
    }

    bool
    feasible(const std::vector<double> &x, double tol, bool relax, detail::Budget *budget) const {
        if (x.size() != vars_.size()) return false;
        if (!detail::Budget::scan(budget, x.size(), [&](std::size_t begin, std::size_t end) {
                for (auto j = begin; j < end; ++j) {
                    if (!std::isfinite(x[j]) || x[j] < vars_[j].lower - tol ||
                        x[j] > vars_[j].upper + tol)
                        return false;
                    if (!relax && vars_[j].type != VarType::Continuous &&
                        std::abs(x[j] - std::round(x[j])) > tol)
                        return false;
                }
                return true;
            }))
            return false;
        for (const auto &r : rows_) {
            if (budget) budget->checkpoint(1);
            long double lhs = 0, magnitude = std::abs(r.rhs);
            auto evaluate = [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) {
                    const auto [j, a] = r.terms[k];
                    lhs += (long double)a * x[j];
                    magnitude += std::abs((long double)a * x[j]);
                }
            };
            if (budget)
                budget->chunks(r.terms.size(), evaluate);
            else
                evaluate(0, r.terms.size());
            if (lhs > (long double)r.rhs + tol + 1e-12L * magnitude) return false;
        }
        return true;
    }

public:
    int add_variable(double lower = 0, double upper = inf, VarType type = VarType::Continuous) {
        if (std::isnan(lower) || std::isnan(upper) || lower == inf || upper == -inf)
            throw std::invalid_argument("invalid variable bound");
        if (type == VarType::Binary) {
            lower = std::max(lower, 0.0);
            upper = std::min(upper, 1.0);
        }
        root_cache_.reset();
        vars_.push_back({lower, upper, type});
        objective_.push_back(0);
        return int(vars_.size()) - 1;
    }
    int add_integer(double lower = 0, double upper = inf) {
        return add_variable(lower, upper, VarType::Integer);
    }
    int add_binary() { return add_variable(0, 1, VarType::Binary); }
    int num_variables() const { return int(vars_.size()); }
    void clear_root_cache() const { root_cache_.reset(); }
    void add_constraint(Terms terms, Relation rel, double rhs) {
        if (!std::isfinite(rhs)) throw std::invalid_argument("non-finite rhs");
        terms = normalize(std::move(terms));
        if (rel == Relation::LessEqual) {
            rows_.push_back({std::move(terms), rhs});
            return;
        }
        if (rel != Relation::GreaterEqual) rows_.push_back({terms, rhs});
        if (rel != Relation::LessEqual) {
            for (auto &t : terms) t.second = -t.second;
            rows_.push_back({std::move(terms), -rhs});
        }
    }
    void add_le(Terms a, double b) { add_constraint(std::move(a), Relation::LessEqual, b); }
    void add_ge(Terms a, double b) { add_constraint(std::move(a), Relation::GreaterEqual, b); }
    void add_eq(Terms a, double b) { add_constraint(std::move(a), Relation::Equal, b); }
    void set_objective(Terms terms, Sense sense = Sense::Minimize, double constant = 0) {
        if (!std::isfinite(constant)) throw std::invalid_argument("non-finite objective constant");
        terms = normalize(std::move(terms));
        root_cache_.reset();
        std::fill(objective_.begin(), objective_.end(), 0);
        for (auto [j, a] : terms) objective_[j] = a;
        sense_ = sense;
        constant_ = constant;
    }
    void minimize(Terms terms, double constant = 0) {
        set_objective(std::move(terms), Sense::Minimize, constant);
    }
    void maximize(Terms terms, double constant = 0) {
        set_objective(std::move(terms), Sense::Maximize, constant);
    }

    bool is_feasible(const std::vector<double> &x, double tol = 1e-7, bool relax = false) const {
        return feasible(x, tol, relax, nullptr);
    }

    // relax=true なら整数変数も連続変数として解く。
    Result solve(const Options &opt = {}, bool relax = false) const {
        if (!(opt.time_limit >= 0) || !(opt.feasibility_tolerance > 0) ||
            !(opt.integrality_tolerance > 0) || !(opt.pivot_tolerance > 0) ||
            !(opt.absolute_gap >= 0) || !(opt.relative_gap >= 0) ||
            !std::isfinite(opt.feasibility_tolerance) ||
            !std::isfinite(opt.integrality_tolerance) || !std::isfinite(opt.pivot_tolerance) ||
            !std::isfinite(opt.absolute_gap) || !std::isfinite(opt.relative_gap) ||
            opt.strong_branching_candidates < 0 || opt.strong_branching_iterations < 0 ||
            opt.strong_branching_depth < 0 || opt.root_cut_rounds < -1 || opt.cuts_per_round < 0)
            throw std::invalid_argument("invalid solver options");
        Result ret;
        detail::Budget budget{opt, ret.stats};
        ret.stats.original_variables = int(vars_.size());
        ret.stats.original_constraints = int(rows_.size());
        const double sign = sense_ == Sense::Maximize ? 1 : -1;
        double incumbent = -inf;
        double global_bound = inf;
        auto finish = [&](Status status) {
            if (status == Status::Infeasible && ret.has_solution) status = Status::NumericalError;
            // カットで条件が悪くなった場合は元のモデルへ戻る。暫定解と総予算を継承。
            if (status == Status::NumericalError && ret.stats.cuts) {
                Options retry = opt;
                retry.root_cut_rounds = 0;
                retry.time_limit = std::max(0.0, opt.time_limit - budget.elapsed());
                retry.iteration_limit =
                    opt.iteration_limit - std::min(opt.iteration_limit, ret.stats.pivots);
                retry.node_limit = opt.node_limit - std::min(opt.node_limit, ret.stats.nodes);
                if (ret.has_solution) retry.initial_solution = ret.x;
                Result fresh = solve(retry, relax);
                fresh.stats.nodes += ret.stats.nodes;
                fresh.stats.pivots += ret.stats.pivots;
                fresh.stats.strong_branching_pivots += ret.stats.strong_branching_pivots;
                fresh.stats.root_pivots += ret.stats.root_pivots;
                fresh.stats.reused_root = fresh.stats.reused_root || ret.stats.reused_root;
                fresh.stats.heuristic_solutions += ret.stats.heuristic_solutions;
                fresh.stats.cuts += ret.stats.cuts;
                fresh.stats.reduced_cost_fixings += ret.stats.reduced_cost_fixings;
                fresh.stats.numerical_restarts += ret.stats.numerical_restarts + 1;
                fresh.stats.seconds = budget.elapsed();
                return fresh;
            }
            ret.status = status;
            if (ret.has_solution) global_bound = std::max(global_bound, incumbent);
            ret.best_bound = sign * global_bound;
            ret.stats.seconds = budget.elapsed();
            return ret;
        };
        auto remember = [&](const std::vector<double> &x, bool heuristic) {
            double val = objective_value(x) * sign;
            if (!std::isfinite(val) || val <= incumbent) return false;
            incumbent = val;
            ret.has_solution = true;
            ret.x = x;
            ret.objective = sign * val;
            if (heuristic) ++ret.stats.heuristic_solutions;
            return true;
        };
        auto accept = [&](const std::vector<double> &x, bool heuristic, bool timed_check = true) {
            if (!feasible(x, opt.feasibility_tolerance, relax, timed_check ? &budget : nullptr))
                return false;
            return remember(x, heuristic);
        };
        try {
            if (!opt.initial_solution.empty()) {
                if (!is_feasible(opt.initial_solution, opt.feasibility_tolerance, relax))
                    throw std::invalid_argument("initial_solution is not feasible");
                remember(opt.initial_solution, true);
            }
            if (budget.stop() || opt.node_limit == 0) return finish(Status::Limit);
            auto cached = opt.reuse_root && root_cache_ && root_cache_->matches(opt, relax)
                              ? root_cache_
                              : nullptr;
            bool reuse = bool(cached);
            detail::Prepared p;
            detail::Tableau root;
            std::vector<int> integers;
            std::vector<double> root_x;
            for (;;) {
                if (budget.stop() || ret.stats.nodes >= opt.node_limit)
                    return finish(Status::Limit);
                if (reuse) {
                    ret.stats.reused_root = true;
                    global_bound = cached->tableau.value() + cached->prepared.offset;
                    p = cached->prepared.clone(budget);
                    root = cached->tableau.clone(budget);
                    append_to_root(p, root, cached->rows, budget);
                } else
                    p = prepare(budget, relax);
                budget.checkpoint();
                if (p.infeasible) {
                    global_bound = -inf;
                    return finish(Status::Infeasible);
                }
                if (opt.rounding_heuristic) {
                    // 根 LP が時間切れでも使える簡単な暫定解を先に検査する。
                    for (bool upper : {false, true}) {
                        auto x = budget.filled<double>(p.vars.size());
                        budget.chunks(p.vars.size(), [&](std::size_t begin, std::size_t end) {
                            for (auto j = begin; j < end; ++j) {
                                const auto &v = p.vars[j];
                                double target = upper ? v.upper : v.lower;
                                x[j] = std::isfinite(target)
                                           ? target
                                           : std::max(v.lower, std::min(v.upper, 0.0));
                            }
                        });
                        accept(x, true);
                    }
                }
                ret.stats.presolved_variables = int(p.c.size());
                ret.stats.presolved_constraints = int(p.b.size());
                integers.clear();
                for (int j = 0; j < int(vars_.size()); ++j)
                    if (!relax && vars_[j].type != VarType::Continuous && !p.map[j].terms.empty())
                        integers.push_back(j);
                // 多数の非負変数を持つ covering/assignment LP は、双対のほうが小さい
                // 基底で Phase I も不要になることが多い。MIP は分枝用の主辞書を保持する。
                if (!reuse && opt.dualize && integers.empty() && !p.b.empty() &&
                    p.c.size() > 2 * p.b.size() &&
                    std::all_of(p.c.begin(), p.c.end(), [](double c) { return c <= 0; })) {
                    std::vector<std::vector<double>> a;
                    a.reserve(p.c.size());
                    for (std::size_t j = 0; j < p.c.size(); ++j)
                        a.push_back(budget.filled<double>(p.b.size()));
                    std::vector<double> b = p.c, c = p.b;
                    for (double &v : b) v = -v;
                    for (double &v : c) v = -v;
                    for (int i = 0; i < int(p.b.size()); ++i)
                        budget.chunks(p.a[i].size(), [&](std::size_t begin, std::size_t end) {
                            for (auto k = begin; k < end; ++k) {
                                auto [j, v] = p.a[i][k];
                                a[j][i] = -v;
                            }
                        });
                    detail::Tableau dual(a, b, c, &budget);
                    ++ret.stats.nodes;
                    const auto before = ret.stats.pivots;
                    const Status s = dual.solve(budget);
                    ret.stats.root_pivots += ret.stats.pivots - before;
                    if (s == Status::Limit) return finish(s);
                    if (s == Status::Optimal) {
                        // 双対の各制約の slack に対する reduced cost が元 LP の変数値。
                        std::vector<double> y(p.c.size());
                        for (int j = 0; j < int(y.size()); ++j) {
                            const int pos = dual.location[dual.n + j];
                            if (pos >= 0) y[j] = dual.row(dual.m)[pos] / dual.objective_scale;
                        }
                        auto x = recover(p, y);
                        const double obj = sign * objective_value(x);
                        const double bound = p.offset - dual.value();
                        if (std::isfinite(obj) && std::isfinite(bound) &&
                            is_feasible(x, opt.feasibility_tolerance, relax) &&
                            std::abs(obj - bound) <=
                                opt.feasibility_tolerance * (1 + std::abs(obj))) {
                            remember(x, false);
                            global_bound = incumbent;
                            ret.stats.solved_dual = true;
                            return finish(Status::Optimal);
                        }
                    }
                    // 双対の状態だけでは曖昧な場合や復元誤差がある場合は、主問題で判定。
                }
                if (budget.stop() || ret.stats.nodes >= opt.node_limit)
                    return finish(Status::Limit);
                if (!reuse) root = detail::Tableau(p.a, p.b, p.c, &budget);
                ++ret.stats.nodes;
                const auto before = ret.stats.pivots;
                // 長い再最適化は一度だけ作り直す。総時間・pivot・node の予算は継続。
                const int soft_limit =
                    int(std::min<long long>(std::numeric_limits<int>::max(),
                                            std::max(256LL, 4LL * (root.n + (long long)root.m))));
                Status root_status = reuse ? root.dual(budget, soft_limit) : root.solve(budget);
                ret.stats.root_pivots += ret.stats.pivots - before;
                if (root_status == Status::Optimal) {
                    root_x = recover(p, root);
                    const double obj = sign * objective_value(root_x);
                    const double bound = root.value() + p.offset;
                    bool valid =
                        std::isfinite(obj) && std::isfinite(bound) &&
                        feasible(root_x, opt.feasibility_tolerance, true, &budget) &&
                        std::abs(obj - bound) <= opt.feasibility_tolerance * (1 + std::abs(obj));
                    for (int j = 0; j <= root.n && valid; ++j)
                        if (root.nonbasic[j] != -1)
                            valid = std::isfinite(root.row(root.m)[j]) &&
                                    root.row(root.m)[j] >= -10 * opt.pivot_tolerance;
                    if (!valid) root_status = Status::NumericalError;
                }
                if (reuse && root_status != Status::Optimal &&
                    (root_status != Status::Limit || !budget.stop())) {
                    // 不変の根以外は保存しない。再利用中の不整合・実行不能も新規求解で確認。
                    root_cache_.reset();
                    cached.reset();
                    reuse = false;
                    global_bound = inf;
                    ++ret.stats.numerical_restarts;
                    continue;
                }
                if (root_status != Status::Optimal) {
                    // Phase II の打切りなら現在の基底も実行可能候補。整数性も検査する。
                    if (root_status == Status::Limit && root.primal_feasible)
                        accept(recover(p, root), false, false);
                    if (root_status == Status::Infeasible) global_bound = -inf;
                    if (root_status == Status::Unbounded &&
                        !is_feasible(recover(p, root), opt.feasibility_tolerance, true))
                        root_status = Status::NumericalError;
                    if (root_status == Status::Unbounded && !integers.empty())
                        root_status = Status::RelaxationUnbounded;
                    return finish(root_status);
                }
                break;
            }
            global_bound = root.value() + p.offset;
            // 保存の途中で時間切れになっても、既に確認できた整数解は失わない。
            bool integral_root = true;
            auto rounded = root_x;
            for (int j : integers) {
                if (std::abs(rounded[j] - std::round(rounded[j])) > opt.integrality_tolerance)
                    integral_root = false;
                rounded[j] = std::round(rounded[j]);
            }
            bool root_proven = false;
            if (integral_root && accept(rounded, false))
                root_proven = std::abs(incumbent - global_bound) <=
                              opt.feasibility_tolerance * (1 + std::abs(incumbent));
            // 同じ値の initial_solution が既にあっても根で最適性を証明できる。
            if (integral_root && ret.has_solution &&
                std::abs(incumbent - global_bound) <=
                    opt.feasibility_tolerance * (1 + std::abs(incumbent)))
                root_proven = true;
            if (opt.reuse_root && (!reuse || cached->rows != rows_.size())) {
                try {
                    save_root(p, root, budget, relax, root_proven);
                } catch (const detail::TimeLimitReached &) {
                    if (!root_proven) throw;
                }
            }
            if (root_proven) {
                global_bound = incumbent;
                return finish(Status::Optimal);
            }
            if (integers.empty()) return finish(Status::NumericalError);

            std::vector<double> rounding_buffer;
            auto round_candidate = [&](const std::vector<double> &x) {
                for (int mode = 0; mode < 3; ++mode) {
                    rounding_buffer.assign(x.begin(), x.end());
                    auto &y = rounding_buffer;
                    for (int j : integers)
                        y[j] = mode == 0   ? std::round(y[j])
                               : mode == 1 ? std::floor(y[j])
                                           : std::ceil(y[j]);
                    accept(y, true);
                }
            };
            if (opt.rounding_heuristic) round_candidate(root_x);
            if (opt.rounding_heuristic && ret.has_solution) {
                // 実行可能な丸め解から、目的を改善する整数の1座標移動を試す。
                // 問題の形は仮定せず、各候補を元の全制約で検査する。
                auto order = budget.copy(integers);
                std::sort(order.begin(), order.end(), [&](int a, int b) {
                    budget.checkpoint(1);
                    return std::abs(objective_[a]) > std::abs(objective_[b]);
                });
                auto x = budget.copy(ret.x);
                for (int j : order) {
                    budget.checkpoint(1);
                    const double gain = sign * objective_[j];
                    if (gain == 0) continue;
                    const double old = x[j], next = old + (gain > 0 ? 1 : -1);
                    if (next < p.vars[j].lower || next > p.vars[j].upper) continue;
                    x[j] = next;
                    if (feasible(x, opt.feasibility_tolerance, false, &budget))
                        remember(x, true);
                    else
                        x[j] = old;
                }
            }

            std::vector<bool> integer_column(p.c.size(), false);
            for (int j : integers)
                if (p.map[j].terms.size() == 1 && p.map[j].offset == std::round(p.map[j].offset))
                    integer_column[p.map[j].terms[0].first] = true;
            int cut_rounds = opt.root_cut_rounds;
            if (cut_rounds < 0) {
                const auto x = recover(p, root);
                int fractional = 0;
                for (int j : integers)
                    fractional += std::abs(x[j] - std::round(x[j])) > opt.integrality_tolerance;
                cut_rounds = fractional == 1 ? 3 : 0;
            }
            for (int round = 0; round < cut_rounds; ++round) {
                if (budget.stop()) return finish(Status::Limit);
                const int added = detail::gomory_cuts(root, p, integer_column, budget);
                if (!added) break;
                ret.stats.cuts += added;
                Status s = root.dual(budget);
                if (s == Status::Infeasible) {
                    global_bound = -inf;
                    return finish(s);
                }
                if (s != Status::Optimal) return finish(s);
                global_bound = root.value() + p.offset;
            }

            // 整数目的関数なら LP 上界を格子に切り下げ、無駄な探索を省く。
            double lattice_base = sign * constant_, lattice_step = 0;
            if (opt.integer_objective) {
                long long g = 0;
                bool valid = true;
                for (int j = 0; j < int(vars_.size()); ++j) {
                    if (p.map[j].terms.empty()) {
                        lattice_base += sign * objective_[j] * p.map[j].offset;
                        continue;
                    }
                    double a = objective_[j];
                    if (a == 0) continue;
                    if (vars_[j].type == VarType::Continuous || a != std::round(a) ||
                        std::abs(a) > 1e9) {
                        valid = false;
                        break;
                    }
                    g = std::gcd(g, (long long)std::llround(std::abs(a)));
                }
                if (valid) lattice_step = double(g);
            }
            auto tighten = [&](double bound) {
                if (lattice_step && std::isfinite(bound)) {
                    const double slack =
                        opt.integrality_tolerance +
                        1e-12 * (std::abs(bound) + std::abs(lattice_base)) / lattice_step;
                    return lattice_base +
                           lattice_step * std::floor((bound - lattice_base) / lattice_step + slack);
                }
                return bound;
            };
            auto prunable = [&](double b) {
                return ret.has_solution && tighten(b) <= incumbent + opt.absolute_gap +
                                                             opt.relative_gap * std::abs(incumbent);
            };
            if (opt.rounding_heuristic) round_candidate(recover(p, root));
            if (prunable(global_bound)) {
                global_bound = tighten(global_bound);
                return finish(Status::Optimal);
            }

            auto improvement_target = [&] {
                if (!lattice_step) return incumbent;
                return lattice_base +
                       lattice_step * (std::floor((incumbent - lattice_base) / lattice_step) + 1);
            };
            std::vector<double> step;
            if (opt.reduced_cost_fixing) {
                step = budget.filled<double>(root.next_id);
                for (int j = 0; j < root.original_n; ++j)
                    if (integer_column[j]) step[j] = 1;
                // 整数変数の上下限に対応する slack も、スケーリング後の刻みを持つ。
                for (std::size_t i = 0; i < p.a.size(); ++i) {
                    budget.checkpoint(1);
                    if (p.a[i].size() != 1 || !integer_column[p.a[i][0].first]) continue;
                    const double scale = std::abs(p.a[i][0].second);
                    const double rhs = p.b[i] / scale;
                    if (std::isfinite(rhs) && rhs == std::round(rhs))
                        step[root.original_n + i] = scale;
                }
                if (ret.has_solution)
                    ret.stats.reduced_cost_fixings +=
                        root.fix_nonbasic(step, global_bound - improvement_target(), budget);
            }

            struct Node {
                detail::Tableau tab;
                std::vector<detail::Row> branches;
                double bound;
                int depth;
                bool solved;
                int branch_var = -1;
                bool up = false;
                double distance = 0, parent_bound = 0;
            };
            std::vector<Node> pending;
            pending.push_back({std::move(root), {}, global_bound, 0, true});
            struct Pseudo {
                double down = 0, up = 0;
                int nd = 0, nu = 0;
            };
            std::vector<Pseudo> pc(vars_.size());
            auto add_branch = [&](detail::Tableau &tab, int j, bool up, double bound) {
                budget.checkpoint();
                const double s = up ? -1 : 1;
                tab.add_constraint(p.map[j].terms, s * (bound - p.map[j].offset), nullptr, s);
            };
            auto update_pc = [&](int j, bool up, double distance, double gain) {
                if (!std::isfinite(gain)) return;
                auto &v = pc[j];
                (up ? v.up : v.down) += std::max(0.0, gain) / std::max(distance, 1e-6);
                ++(up ? v.nu : v.nd);
            };
            std::vector<std::pair<double, int>> candidates;
            while (!pending.empty()) {
                global_bound = incumbent;
                for (const auto &node : pending)
                    global_bound = std::max(global_bound, tighten(node.bound));
                if (prunable(global_bound)) return finish(Status::Optimal);
                if (budget.stop()) return finish(Status::Limit);
                Node node = std::move(pending.back());
                pending.pop_back();
                if (prunable(node.bound)) continue;
                if (node.depth) {
                    if (ret.stats.nodes >= opt.node_limit) return finish(Status::Limit);
                    ++ret.stats.nodes;
                }
                if (!node.solved) {
                    Status s;
                    if (opt.warm_start)
                        s = node.tab.dual(budget);
                    else {
                        node.tab = detail::Tableau(p.a, p.b, p.c, &budget);
                        // コールド比較用: 分枝行を元の非基底辞書に追加し Phase I から解く。
                        for (const auto &b : node.branches) {
                            int j = b.terms[0].first;
                            bool up = b.terms[0].second < 0;
                            add_branch(node.tab, j, up, up ? -b.rhs : b.rhs);
                            node.tab.write_row(node.tab.m - 1)[node.tab.n] = -1;
                        }
                        s = node.tab.solve(budget);
                    }
                    if (s == Status::Infeasible) continue;
                    if (s != Status::Optimal) return finish(s);
                    node.bound = node.tab.value() + p.offset;
                    if (node.branch_var >= 0)
                        update_pc(node.branch_var, node.up, node.distance,
                                  node.parent_bound - node.bound);
                    if (prunable(node.bound)) continue;
                }
                auto x = recover(p, node.tab);
                if (!is_feasible(x, opt.feasibility_tolerance, true))
                    return finish(Status::NumericalError);
                if (opt.rounding_heuristic) round_candidate(x);
                if (prunable(node.bound)) continue;
                if (opt.reduced_cost_fixing && ret.has_solution && node.depth) {
                    step.resize(node.tab.next_id, 1);
                    ret.stats.reduced_cost_fixings +=
                        node.tab.fix_nonbasic(step, node.bound - improvement_target(), budget,
                                              std::max(4, node.tab.n / 8));
                }
                candidates.clear();
                for (int j : integers) {
                    const double f = x[j] - std::floor(x[j]);
                    if (std::min(f, 1 - f) <= opt.integrality_tolerance) continue;
                    const auto &v = pc[j];
                    double d = f, u = 1 - f;
                    if (opt.pseudocost) {
                        d *= v.nd ? v.down / v.nd : 1;
                        u *= v.nu ? v.up / v.nu : 1;
                    }
                    candidates.emplace_back(std::min(d, u) + 0.1 * std::max(d, u), j);
                }
                if (candidates.empty()) {
                    for (int j : integers) x[j] = std::round(x[j]);
                    if (!is_feasible(x, opt.feasibility_tolerance, false))
                        return finish(Status::NumericalError);
                    remember(x, false);
                    continue;
                }
                std::sort(
                    candidates.begin(), candidates.end(), std::greater<std::pair<double, int>>());
                int branch = candidates[0].second;
                double best_score = -1;
                detail::Tableau best_down, best_up;
                Status down_status = Status::Limit, up_status = Status::Limit;
                bool probed = false;
                const int probes =
                    node.depth < opt.strong_branching_depth
                        ? std::min(int(candidates.size()), opt.strong_branching_candidates)
                        : 0;
                for (int k = 0; k < probes; ++k) {
                    int j = candidates[k].second;
                    const double f = x[j] - std::floor(x[j]);
                    budget.checkpoint();
                    detail::Tableau down = node.tab, up = node.tab;
                    add_branch(down, j, false, std::floor(x[j]));
                    add_branch(up, j, true, std::ceil(x[j]));
                    auto before = ret.stats.pivots;
                    auto ds = down.dual(budget, opt.strong_branching_iterations);
                    auto us = up.dual(budget, opt.strong_branching_iterations);
                    ret.stats.strong_branching_pivots += ret.stats.pivots - before;
                    if (budget.stop()) return finish(Status::Limit);
                    if (ds == Status::NumericalError || us == Status::NumericalError)
                        return finish(Status::NumericalError);
                    double dg = ds == Status::Infeasible
                                    ? inf
                                    : std::max(0.0, node.bound - (down.value() + p.offset));
                    double ug = us == Status::Infeasible
                                    ? inf
                                    : std::max(0.0, node.bound - (up.value() + p.offset));
                    if (ds == Status::Optimal) update_pc(j, false, f, dg);
                    if (us == Status::Optimal) update_pc(j, true, 1 - f, ug);
                    if (opt.rounding_heuristic) {
                        if (ds == Status::Optimal) round_candidate(recover(p, down));
                        if (us == Status::Optimal) round_candidate(recover(p, up));
                    }
                    const double score = std::min(dg, ug) + 0.1 * std::max(dg, ug);
                    if (!probed || score > best_score) {
                        probed = true;
                        branch = j;
                        best_score = score;
                        best_down = std::move(down);
                        best_up = std::move(up);
                        down_status = ds;
                        up_status = us;
                    }
                    if (ds == Status::Infeasible || us == Status::Infeasible) break;
                }
                if (!probed) {
                    budget.checkpoint();
                    best_down = node.tab;
                    best_up = std::move(node.tab);
                    add_branch(best_down, branch, false, std::floor(x[branch]));
                    add_branch(best_up, branch, true, std::ceil(x[branch]));
                }
                auto push = [&](detail::Tableau tab, bool up, Status s) {
                    if (s == Status::Infeasible) return;
                    double b = std::min(node.bound, tab.value() + p.offset);
                    if (prunable(b)) return;
                    std::vector<detail::Row> branches;
                    if (!opt.warm_start) {
                        branches = node.branches;
                        branches.push_back({{{branch, up ? -1.0 : 1.0}},
                                            up ? -std::ceil(x[branch]) : std::floor(x[branch])});
                    }
                    const double distance =
                        up ? std::ceil(x[branch]) - x[branch] : x[branch] - std::floor(x[branch]);
                    pending.push_back({std::move(tab), std::move(branches), b, node.depth + 1,
                                       s == Status::Optimal, branch, up, distance, node.bound});
                };
                // 有望な側を先に深さ優先探索し、暫定解を早めに得る。
                bool up_first = best_up.value() > best_down.value();
                if (up_first) {
                    push(std::move(best_down), false, down_status);
                    push(std::move(best_up), true, up_status);
                } else {
                    push(std::move(best_up), true, up_status);
                    push(std::move(best_down), false, down_status);
                }
            }
            global_bound = incumbent;
            return finish(ret.has_solution ? Status::Optimal : Status::Infeasible);
        } catch (const detail::TimeLimitReached &) { return finish(Status::Limit); }
    }
    Result solve_relaxation(const Options &opt = {}) const { return solve(opt, true); }
};

// 同じ変数で制約追加・目的関数変更を繰り返す LP 向け。
// maximize c*x, A*x<=b, x>=0。整数制約・変数変換・presolve は Model を使う。
// 制約追加だけなら双対法、目的関数変更だけなら主単体法で基底を再利用する。
class LPSolver {
    enum class State { Cold, Optimal, Primal, Dual };
    State state_ = State::Cold;
    std::vector<std::vector<double>> a_;
    std::vector<double> b_, c_;
    detail::Tableau tab_;
    void check_coefficients(const std::vector<double> &row) const {
        if (row.size() != c_.size()) throw std::invalid_argument("LP dimension mismatch");
        for (double v : row)
            if (!std::isfinite(v)) throw std::invalid_argument("non-finite LP coefficient");
    }

public:
    LPSolver(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<double> c)
        : a_(std::move(a)), b_(std::move(b)), c_(std::move(c)) {
        if (a_.size() != b_.size()) throw std::invalid_argument("LP A/b dimension mismatch");
        check_coefficients(c_);
        for (const auto &row : a_) check_coefficients(row);
        for (double v : b_)
            if (!std::isfinite(v)) throw std::invalid_argument("non-finite LP rhs");
    }
    void add_le(std::vector<double> a, double b) {
        check_coefficients(a);
        if (!std::isfinite(b)) throw std::invalid_argument("non-finite LP rhs");
        if (state_ == State::Optimal || state_ == State::Dual) {
            Terms terms;
            for (int j = 0; j < int(a.size()); ++j)
                if (a[j]) terms.emplace_back(j, a[j]);
            tab_.add_constraint(terms, b);
            state_ = State::Dual;
        } else
            state_ = State::Cold;
        a_.push_back(std::move(a));
        b_.push_back(b);
    }
    void set_objective(std::vector<double> c) {
        check_coefficients(c);
        c_ = std::move(c);
        if (state_ == State::Optimal || state_ == State::Primal) {
            tab_.scale_objective(c_);
            double *z = tab_.write_row(tab_.m);
            std::fill(z, z + tab_.stride, 0);
            for (int j = 0; j <= tab_.n; ++j) {
                int id = tab_.nonbasic[j];
                if (id >= 0 && id < tab_.n) z[j] = -c_[id] * tab_.objective_scale;
            }
            for (int i = 0; i < tab_.m; ++i) {
                int id = tab_.basic[i];
                if (id < 0 || id >= tab_.n || c_[id] == 0) continue;
                const double *row = tab_.row(i);
                for (int j = 0; j < tab_.stride; ++j)
                    z[j] += c_[id] * tab_.objective_scale * row[j];
            }
            state_ = State::Primal;
        } else
            state_ = State::Cold;
    }
    Result solve(const Options &opt = {}) {
        if (!(opt.time_limit >= 0) || !(opt.feasibility_tolerance > 0) ||
            !(opt.pivot_tolerance > 0) || !std::isfinite(opt.feasibility_tolerance) ||
            !std::isfinite(opt.pivot_tolerance))
            throw std::invalid_argument("invalid LP options");
        Result r;
        r.best_bound = inf;
        r.stats.original_variables = r.stats.presolved_variables = int(c_.size());
        r.stats.original_constraints = r.stats.presolved_constraints = int(b_.size());
        detail::Budget budget{opt, r.stats};
        try {
            if (budget.stop() || opt.node_limit == 0)
                r.status = Status::Limit;
            else {
                ++r.stats.nodes;
                if (state_ == State::Cold) {
                    tab_ = detail::Tableau(a_, b_, c_, &budget);
                    state_ = State::Primal;
                    r.status = tab_.solve(budget);
                } else if (state_ == State::Dual)
                    r.status = tab_.dual(budget);
                else if (state_ == State::Primal)
                    r.status = tab_.primal(tab_.m, budget);
                else
                    r.status = Status::Optimal;
            }
        } catch (const detail::TimeLimitReached &) { r.status = Status::Limit; }
        if (r.status == Status::Optimal || (r.status == Status::Limit && tab_.primal_feasible &&
                                            state_ != State::Cold && state_ != State::Dual)) {
            r.x = tab_.solution();
            bool valid = true;
            long double objective = 0;
            for (int j = 0; j < int(c_.size()); ++j) {
                if (!std::isfinite(r.x[j]) || r.x[j] < -opt.feasibility_tolerance) valid = false;
                objective += (long double)c_[j] * r.x[j];
            }
            for (int i = 0; i < int(b_.size()); ++i) {
                long double lhs = 0, magnitude = std::abs(b_[i]);
                for (int j = 0; j < int(c_.size()); ++j) {
                    lhs += (long double)a_[i][j] * r.x[j];
                    magnitude += std::abs((long double)a_[i][j] * r.x[j]);
                }
                if (lhs > (long double)b_[i] + opt.feasibility_tolerance + 1e-12L * magnitude)
                    valid = false;
            }
            if (std::abs(double(objective) - tab_.value()) >
                opt.feasibility_tolerance * (1 + std::abs(double(objective))))
                valid = false;
            if (valid && std::isfinite(double(objective))) {
                r.has_solution = true;
                r.objective = double(objective);
                if (r.status == Status::Optimal) r.best_bound = r.objective;
            } else {
                if (r.status == Status::Optimal) r.status = Status::NumericalError;
                r.x.clear();
            }
        }
        if (r.status == Status::Infeasible) r.best_bound = -inf;
        state_ = r.status == Status::Optimal ? State::Optimal : State::Cold;
        r.stats.seconds = budget.elapsed();
        return r;
    }
};

// 行列形式の簡単な入口。maximize c*x subject to A*x <= b, x >= 0。
inline Result linprog(const std::vector<std::vector<double>> &a, const std::vector<double> &b,
                      const std::vector<double> &c, const Options &options = {}) {
    if (a.size() != b.size()) throw std::invalid_argument("A/b dimension mismatch");
    Model model;
    for (std::size_t j = 0; j < c.size(); ++j) model.add_variable();
    Terms cost;
    for (int j = 0; j < int(c.size()); ++j) cost.emplace_back(j, c[j]);
    model.maximize(std::move(cost));
    for (int i = 0; i < int(a.size()); ++i) {
        if (a[i].size() != c.size()) throw std::invalid_argument("A/c dimension mismatch");
        Terms row;
        for (int j = 0; j < int(c.size()); ++j)
            if (a[i][j] != 0) row.emplace_back(j, a[i][j]);
        model.add_le(std::move(row), b[i]);
    }
    return model.solve(options);
}

} // namespace cpopt
#endif // mip_HPP
