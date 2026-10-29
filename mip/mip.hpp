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
#include <unordered_map>
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
    bool advanced_presolve = true; // 変更行の伝播、整数係数強化、重複行除去。
    bool probing = false; // 0/1仮固定による境界・係数強化。探索経路への影響があるため任意。
    bool integer_substitution =
        false; // 整数等式の係数±1の変数を消去。探索経路への影響があるため任意。
    bool equality_substitution = true; // 等式による変数消去。advanced_presolveで対象を広げる。
    bool scaling = true;
    bool dualize = true; // 連続 LP で有利な形なら双対問題を解く。
    bool warm_start = true;
    bool reuse_root = true; // 同じ Model の再求解で、最適な根基底を保存・再利用。
    bool pseudocost = true;
    bool rounding_heuristic = true;
    bool feasibility_pump = true; // 暫定解がないとき、丸め点へのLP射影と小規模な仮固定探索。
    bool diving_heuristic = true; // 根LPから少数の整数仮固定・双対再最適化を試す。
    bool integer_objective = true;
    bool reduced_cost_fixing = true; // 暫定解を改善できない非基底整数変数を固定・列削除。
    bool bounded_variables = true;   // 有限上限を追加行にせず、単体法で直接扱う。
    int sparse_tableau = -1;         // -1: 密度で自動選択・切替、0: 密、1: 疎辞書を維持。
    bool basis_repair = true;        // 特異基底・原問題の残差を検出した場合の再構築と再最適化。
    bool root_diversification =
        true;               // 暫定解がないとき、別の最適基底からも汎用ヒューリスティックを始める。
    bool dual_devex = true; // 双対単体法の行選択を参照基底の辺長で重み付け。
    bool primal_sweeps = true;                       // 大きいモデルの根LP前に疎な座標改善を試す。
    std::uint64_t tableau_cache_limit = 64ULL << 20; // 待機ノードの辞書の見積りバイト上限。
    int root_cut_rounds = -1; // -1: 小数の整数変数が1個のときだけ3巡、0: 無効。
    int cuts_per_round = 8;
    int strong_branching_candidates = 4;
    int strong_branching_depth = 3;
    int strong_branching_iterations = 40;
    std::vector<double> initial_solution; // 元の変数順。実行可能性を検査してから使う。
};

struct Statistics {
    std::uint64_t nodes = 0, pivots = 0, strong_branching_pivots = 0;
    std::uint64_t bound_flips = 0; // pivots に含まれる、基底交換を伴わない上下限間の移動。
    std::uint64_t tableau_evictions = 0, tableau_rebuilds = 0;
    std::uint64_t peak_pending_tableau_bytes = 0;
    std::uint64_t sparse_pivots = 0, sparse_tableaux = 0, dense_conversions = 0;
    std::uint64_t heuristic_solutions = 0;
    std::uint64_t basis_refactorizations = 0; // 元係数からの補償演算による基底再構築。
    std::uint64_t objective_perturbations = 0;
    std::uint64_t pump_pivots = 0, pump_lp_solves = 0;
    std::uint64_t diving_pivots = 0, diving_lp_solves = 0;
    int cuts = 0;
    int reduced_cost_fixings = 0;
    int numerical_restarts = 0;
    int substituted_variables = 0;
    std::uint64_t presolve_bound_changes = 0, strengthened_coefficients = 0;
    std::uint64_t probing_trials = 0, probing_bound_changes = 0;
    int duplicate_rows = 0;
    int implied_integers = 0;
    std::uint64_t presolved_nonzeros = 0;
    double presolve_seconds = 0;
    bool solved_dual = false;
    bool reused_root = false;
    std::uint64_t root_storage_bytes =
        0; // 根LP開始時の配列サイズ（ハッシュの空き枠を含む）。RSSではない。
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

inline bool
propagate(const std::vector<Row> &rows, std::vector<Variable> &vars, Budget &budget, int rounds);

// 実行不能の証拠の再構成だけに使う、上限付きの符号付き整数。
// 基数2^30なら積・桁上がりを標準uint64_tだけで処理できる。
// 大きな式は中止する。LP本体の演算や入力値の丸めには使わない。
struct ExactOverflow {};
struct ExactInteger {
    static constexpr std::uint32_t mask = (1U << 30) - 1;
    std::vector<std::uint32_t> digits;
    bool negative = false;
    ExactInteger(long long value = 0) {
        negative = value < 0;
        std::uint64_t a = negative ? std::uint64_t(-(value + 1)) + 1 : std::uint64_t(value);
        while (a) {
            digits.push_back(std::uint32_t(a & mask));
            a >>= 30;
        }
    }
    void normalize() {
        while (!digits.empty() && !digits.back()) digits.pop_back();
        if (digits.empty()) negative = false;
        if (digits.size() > 256) throw ExactOverflow{};
    }
    bool zero() const { return digits.empty(); }
    ExactInteger absolute() const {
        auto x = *this;
        x.negative = false;
        return x;
    }
    ExactInteger operator-() const {
        auto x = *this;
        if (!x.zero()) x.negative = !x.negative;
        return x;
    }
    static int compare(const ExactInteger &a, const ExactInteger &b) {
        if (a.digits.size() != b.digits.size()) return a.digits.size() < b.digits.size() ? -1 : 1;
        for (std::size_t i = a.digits.size(); i--;)
            if (a.digits[i] != b.digits[i]) return a.digits[i] < b.digits[i] ? -1 : 1;
        return 0;
    }
    friend ExactInteger operator+(const ExactInteger &a, const ExactInteger &b) {
        ExactInteger c;
        if (a.negative == b.negative) {
            c.negative = a.negative;
            c.digits.resize(std::max(a.digits.size(), b.digits.size()) + 1);
            std::uint64_t carry = 0;
            for (std::size_t i = 0; i < c.digits.size(); ++i) {
                carry += (i < a.digits.size() ? a.digits[i] : 0) +
                         std::uint64_t(i < b.digits.size() ? b.digits[i] : 0);
                c.digits[i] = std::uint32_t(carry & mask);
                carry >>= 30;
            }
        } else {
            const bool swapped = compare(a, b) < 0;
            const auto &x = swapped ? b : a;
            const auto &y = swapped ? a : b;
            c.negative = x.negative;
            c.digits.resize(x.digits.size());
            long long borrow = 0;
            for (std::size_t i = 0; i < x.digits.size(); ++i) {
                const long long value =
                    (long long)x.digits[i] - (i < y.digits.size() ? y.digits[i] : 0) - borrow;
                borrow = value < 0;
                c.digits[i] = std::uint32_t(value + (borrow << 30));
            }
        }
        c.normalize();
        return c;
    }
    friend ExactInteger operator-(const ExactInteger &a, const ExactInteger &b) { return a + -b; }
    friend ExactInteger operator*(const ExactInteger &a, const ExactInteger &b) {
        if (a.zero() || b.zero()) return {};
        if (a.digits.size() + b.digits.size() > 257) throw ExactOverflow{};
        ExactInteger c;
        c.negative = a.negative != b.negative;
        c.digits.resize(a.digits.size() + b.digits.size());
        for (std::size_t i = 0; i < a.digits.size(); ++i) {
            std::uint64_t carry = 0;
            for (std::size_t j = 0; j < b.digits.size(); ++j) {
                carry += std::uint64_t(a.digits[i]) * b.digits[j] + c.digits[i + j];
                c.digits[i + j] = std::uint32_t(carry & mask);
                carry >>= 30;
            }
            c.digits[i + b.digits.size()] = std::uint32_t(carry);
        }
        c.normalize();
        return c;
    }
    ExactInteger shifted(int bits) const {
        if (zero()) return {};
        if (bits < 0 || digits.size() + std::size_t(bits / 30) > 256) throw ExactOverflow{};
        ExactInteger c;
        c.negative = negative;
        c.digits.resize(digits.size() + bits / 30 + 1);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < digits.size(); ++i) {
            carry += std::uint64_t(digits[i]) << (bits % 30);
            c.digits[i + bits / 30] = std::uint32_t(carry & mask);
            carry >>= 30;
        }
        c.digits.back() = std::uint32_t(carry);
        c.normalize();
        return c;
    }
    ExactInteger divided(const ExactInteger &divisor) const {
        if (divisor.zero()) throw ExactOverflow{};
        ExactInteger q, remainder;
        q.digits.resize(digits.size());
        q.negative = negative != divisor.negative;
        if (divisor.digits.size() == 1) {
            std::uint64_t carry = 0;
            for (std::size_t i = digits.size(); i--;) {
                carry = (carry << 30) + digits[i];
                q.digits[i] = std::uint32_t(carry / divisor.digits[0]);
                carry %= divisor.digits[0];
            }
            if (carry) throw ExactOverflow{};
        } else {
            const auto d = divisor.absolute();
            for (std::size_t bit = digits.size() * 30; bit--;) {
                remainder = remainder.shifted(1);
                if ((digits[bit / 30] >> (bit % 30)) & 1U) remainder = remainder + ExactInteger(1);
                if (compare(remainder, d) >= 0) {
                    remainder = remainder - d;
                    q.digits[bit / 30] |= 1U << (bit % 30);
                }
            }
            if (!remainder.zero()) throw ExactOverflow{};
        }
        q.normalize();
        return q;
    }
    static std::pair<long long, int> dyadic(double value) {
        if (!std::isfinite(value)) throw ExactOverflow{};
        if (!value) return {0, 0};
        int exponent = 0;
        long long mantissa = (long long)std::ldexp(std::frexp(value, &exponent), 53);
        exponent -= 53;
        while (mantissa % 2 == 0) {
            mantissa /= 2;
            ++exponent;
        }
        return {mantissa, exponent};
    }
};

// 辞書: x_B = rhs - D x_N, z = rhs_z - D_z x_N。
// 初期列数は元変数数 + 人工変数 + RHS。固定列は後から削除できる。
// slack の列を別に持たず交換する。
// 誤差なし変換を用いた2成分の演算。基底再構築時だけ使用する。
// long double と double が同じ精度の処理系でも追加の有効桁を保持する。
struct Accurate {
    double hi = 0, lo = 0;
    Accurate() = default;
    Accurate(double x) : hi(x) {}
    static Accurate parts(double h, double l) {
        Accurate r;
        r.hi = h + l;
        r.lo = l - (r.hi - h);
        return r;
    }
    friend Accurate operator+(Accurate a, Accurate b) {
        const double h = a.hi + b.hi, v = h - a.hi;
        return parts(h, ((a.hi - (h - v)) + (b.hi - v)) + (a.lo + b.lo));
    }
    friend Accurate operator-(Accurate a, Accurate b) {
        b.hi = -b.hi;
        b.lo = -b.lo;
        return a + b;
    }
    friend Accurate operator*(Accurate a, Accurate b) {
        const double h = a.hi * b.hi;
        return parts(h, std::fma(a.hi, b.hi, -h) + a.hi * b.lo + a.lo * b.hi + a.lo * b.lo);
    }
    friend Accurate operator/(Accurate a, Accurate b) {
        const double q = a.hi / b.hi;
        const auto rem = a - b * Accurate(q);
        return parts(q, (rem.hi + rem.lo) / b.hi);
    }
    double value() const { return hi + lo; }
};

struct Tableau {
    struct SparseRow {
        struct Cell {
            int column = -1;
            double value = 0;
        };
        std::vector<Cell> table;
        std::vector<double> full;
        int count = 0, used = 0;
        int rhs_column = -1, rhs_slot = -1;
        static unsigned hash(unsigned x) {
            const std::uint64_t product = std::uint64_t(x) * 2654435761U;
            return unsigned(product ^ (product >> 32));
        }
        std::size_t position(int j) const {
            const std::size_t mask = table.size() - 1;
            std::size_t k = hash(unsigned(j)) & mask, first = table.size();
            while (table[k].column != -1) {
                if (table[k].column == j) return k;
                if (table[k].column == -2 && first == table.size()) first = k;
                k = (k + 1) & mask;
            }
            return first == table.size() ? k : first;
        }
        double get(int j) const {
            if (!full.empty()) return full[j];
            if (j == rhs_column) return rhs_slot < 0 ? 0 : table[rhs_slot].value;
            if (table.empty()) return 0;
            const auto &cell = table[position(j)];
            return cell.column == j ? cell.value : 0;
        }
        void put(std::size_t k, int j, double value) {
            auto &cell = table[k];
            if (cell.column == j) {
                if (value != 0)
                    cell.value = value;
                else {
                    cell.column = -2;
                    --count;
                }
            } else if (value != 0) {
                used += cell.column == -1;
                ++count;
                cell = {j, value};
            }
            if (j == rhs_column) rhs_slot = value != 0 ? int(k) : -1;
        }
        void reserve(int extra, Budget *budget = nullptr) {
            if (!full.empty()) return;
            if (!table.empty() && 2 * (count + extra) <= int(table.size()) &&
                4 * (used + extra) <= 3 * int(table.size()))
                return;
            std::size_t size = 8;
            while (size < std::size_t(2 * (count + extra))) size *= 2;
            SparseRow next;
            next.rhs_column = rhs_column;
            next.table = budget ? budget->filled<Cell>(size) : std::vector<Cell>(size);
            for (const auto &cell : table) {
                if (budget) budget->checkpoint(1);
                if (cell.column >= 0)
                    next.put(next.position(cell.column), cell.column, cell.value);
            }
            table = std::move(next.table);
            count = next.count;
            used = next.used;
            rhs_slot = next.rhs_slot;
        }
        void set(int j, double value) {
            if (!full.empty()) {
                full[j] = value;
                return;
            }
            reserve(1);
            put(position(j), j, value);
        }
        void subtract_product(int j, double f, double a) {
            if (!full.empty()) {
                full[j] -= f * a;
                return;
            }
            const auto k = position(j);
            put(k, j, (table[k].column == j ? table[k].value : 0) - f * a);
        }
        std::uint64_t bytes() const {
            return table.size() * sizeof(Cell) + full.size() * sizeof(double);
        }
    };
    int n = 0, m = 0, stride = 0, next_id = 0, original_n = 0;
    static constexpr int fixed_zero = std::numeric_limits<int>::max();
    bool primal_feasible = false;
    bool slack_fixing = true;  // 再構築で追加した行の整数刻みが不明ならslackを固定しない。
    bool needs_repair = false; // 復旧を中断した辞書の境界を枝刈りに使わない。
    double objective_scale = 1;
    // 基底交換で更新しない原係数。分枝・カットの追加行も元の座標で保持する。
    std::vector<std::shared_ptr<const Row>> source_rows;
    std::shared_ptr<const std::vector<double>> source_cost;
    bool perturbed = false;
    bool stable_ratio = false;
    std::vector<double> dense;
    std::vector<std::shared_ptr<std::vector<double>>> data;
    std::vector<std::shared_ptr<SparseRow>> sparse_rows;
    std::uint64_t sparse_bytes = 0;
    mutable std::vector<double> read_buffer; // row() の疎行読み出し用。次の疎行読み出しまで有効。
    std::vector<int> basic, nonbasic;
    std::vector<int> location; // fixed_zero: 0固定、それ以外は非負: 列、負: -行-1。人工変数を除く。
    std::vector<double> upper;          // 元の構造変数だけ。空なら全て上限なし。
    std::vector<unsigned char> flipped; // 1: 辞書の座標は upper[id] - x[id]。
    std::vector<int> pivot_columns;     // 疎な pivot 行の非ゼロ列。作業領域なので複製不要。
    std::vector<double *> pivot_rows;   // detach 済みの行。pivot 中だけ有効。
    Terms sparse_entries;               // 疎pivotの作業領域。反復ごとの再確保を避ける。
    std::vector<int> sparse_affected;
    const double *dense_row(int i) const {
        return data.empty() ? dense.data() + std::size_t(i) * stride : data[i]->data();
    }
    const double *row(int i) const {
        if (!sparse_rows.empty()) {
            const auto &r = *sparse_rows[i];
            if (!r.full.empty()) return r.full.data();
            read_buffer.assign(stride, 0);
            for (const auto &cell : r.table)
                if (cell.column >= 0) read_buffer[cell.column] = cell.value;
            return read_buffer.data();
        }
        return dense_row(i);
    }
    double coefficient(int i, int j) const {
        return sparse_rows.empty() ? dense_row(i)[j] : sparse_rows[i]->get(j);
    }
    std::uint64_t stored_bytes() const {
        return sparse_rows.empty() ? std::uint64_t(m + 2) * stride * sizeof(double) : sparse_bytes;
    }
    template <class F> void nonzeros(int i, F fn) const {
        if (!sparse_rows.empty() && sparse_rows[i]->full.empty()) {
            for (const auto &cell : sparse_rows[i]->table)
                if (cell.column >= 0) fn(cell.column, cell.value);
        } else {
            const double *q = row(i);
            for (int j = 0; j < stride; ++j)
                if (q[j] != 0) fn(j, q[j]);
        }
    }
    void make_full(SparseRow &r, Budget *budget = nullptr) const {
        if (!r.full.empty()) return;
        auto full = budget ? budget->filled<double>(stride) : std::vector<double>(stride);
        for (const auto &cell : r.table) {
            if (budget) budget->checkpoint(1);
            if (cell.column >= 0) full[cell.column] = cell.value;
        }
        r.full = std::move(full);
        std::vector<SparseRow::Cell>{}.swap(r.table);
        r.count = r.used = 0;
    }
    SparseRow &prepare_row(int i, int extra, Budget &budget) {
        auto &ptr = sparse_rows[i];
        const auto before = ptr->bytes();
        if (ptr.use_count() != 1) {
            auto next = std::make_shared<SparseRow>();
            next->table = budget.copy(ptr->table);
            next->full = budget.copy(ptr->full);
            next->count = ptr->count;
            next->used = ptr->used;
            next->rhs_column = ptr->rhs_column;
            next->rhs_slot = ptr->rhs_slot;
            ptr = std::move(next);
        }
        if (ptr->full.empty() && (ptr->count + extra) * 8 > stride)
            make_full(*ptr, &budget);
        else
            ptr->reserve(extra, &budget);
        sparse_bytes = sparse_bytes - before + ptr->bytes();
        return *ptr;
    }
    double *write_row(int i, Budget *budget = nullptr) {
        if (!sparse_rows.empty()) {
            Options options;
            Statistics stats;
            Budget unlimited{options, stats};
            Budget &b = budget ? *budget : unlimited;
            auto &r = prepare_row(i, 0, b);
            const auto before = r.bytes();
            make_full(r, &b);
            sparse_bytes = sparse_bytes - before + r.bytes();
            return r.full.data();
        }
        return write_dense_row(i, budget);
    }
    double *write_dense_row(int i, Budget *budget = nullptr) {
        if (data.empty()) return dense.data() + std::size_t(i) * stride;
        if (data[i].use_count() != 1)
            data[i] =
                std::make_shared<std::vector<double>>(budget ? budget->copy(*data[i]) : *data[i]);
        return data[i]->data();
    }
    double value() const { return coefficient(m, n + 1) / objective_scale; }
    double bound(int id) const {
        return !upper.empty() && id >= 0 && id < original_n ? upper[id] : inf;
    }
    bool reversed(int id) const {
        return !flipped.empty() && id >= 0 && id < original_n && flipped[id];
    }
    double coordinate(int id) const {
        const int pos = location[id];
        const double v = pos < 0 ? coefficient(-pos - 1, n + 1) : 0;
        return reversed(id) ? upper[id] - v : v;
    }
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
        needs_repair = other.needs_repair;
        slack_fixing = other.slack_fixing;
        // 複製直後の分枝行1本で再確保しない。値のコピーと行追加で二重に
        // 配列をコピーするコストを避ける（探索・分枝の順序は同じ）。
        if (!other.data.empty()) data.reserve(other.data.size() + 1);
        data = other.data;
        sparse_rows = other.sparse_rows;
        sparse_bytes = other.sparse_bytes;
        dense.reserve(other.dense.size() + (other.dense.empty() ? 0 : std::size_t(stride) * 2));
        dense.assign(other.dense.begin(), other.dense.end());
        basic.reserve(other.basic.size() + 1);
        location.reserve(other.location.size() + 1);
        basic = other.basic;
        nonbasic = other.nonbasic;
        location = other.location;
        upper = other.upper;
        flipped = other.flipped;
        source_rows = other.source_rows;
        source_cost = other.source_cost;
        perturbed = other.perturbed;
        stable_ratio = other.stable_ratio;
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
        t.needs_repair = needs_repair;
        t.slack_fixing = slack_fixing;
        t.objective_scale = objective_scale;
        budget.checkpoint();
        t.data = budget.copy(data);
        t.sparse_rows = budget.copy(sparse_rows);
        t.sparse_bytes = sparse_bytes;
        t.dense.reserve(dense.size() + (dense.empty() ? 0 : std::size_t(stride) * 2));
        budget.chunks(dense.size(), [&](std::size_t begin, std::size_t end) {
            t.dense.insert(t.dense.end(), dense.begin() + begin, dense.begin() + end);
        });
        t.basic = budget.copy(basic);
        t.nonbasic = budget.copy(nonbasic);
        t.location = budget.copy(location);
        t.upper = budget.copy(upper);
        t.flipped = budget.copy(flipped);
        t.source_rows = budget.copy(source_rows);
        t.source_cost = source_cost;
        t.perturbed = perturbed;
        t.stable_ratio = stable_ratio;
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
    static Terms sparse_coefficients(const Terms &a, Budget *budget) {
        Terms result = budget ? budget->copy(a) : a;
        if (!std::is_sorted(
                result.begin(), result.end(), [](auto a, auto b) { return a.first < b.first; }))
            std::stable_sort(result.begin(), result.end(), [budget](auto a, auto b) {
                if (budget) budget->checkpoint(1);
                return a.first < b.first;
            });
        std::size_t out = 0;
        for (std::size_t k = 0; k < result.size(); ++k) {
            if (budget) budget->checkpoint(1);
            if (k + 1 < result.size() && result[k].first == result[k + 1].first) continue;
            if (result[k].second) result[out++] = result[k];
        }
        result.resize(out);
        return result;
    }
    static Terms sparse_coefficients(const std::vector<double> &a, Budget *budget) {
        Terms result;
        for (int j = 0; j < int(a.size()); ++j) {
            if (budget) budget->checkpoint(1);
            if (a[j]) result.emplace_back(j, a[j]);
        }
        return result;
    }
    static std::size_t count_nonzero(const Terms &a, Budget *) { return a.size(); }
    static std::size_t count_nonzero(const std::vector<double> &a, Budget *budget) {
        std::size_t result = 0;
        for (double v : a) {
            if (budget) budget->checkpoint(1);
            result += v != 0;
        }
        return result;
    }
    void initialize_row(SparseRow &r, const Terms &entries, Budget *budget = nullptr) const {
        r.rhs_column = n + 1;
        r.reserve(int(entries.size()), budget);
        for (auto [j, a] : entries) {
            if (budget) budget->checkpoint(1);
            r.set(j, a);
        }
        if (r.count * 8 > stride) make_full(r, budget);
    }
    template <class Matrix>
    Tableau(const Matrix &a, const std::vector<double> &b, const std::vector<double> &c,
            Budget *budget = nullptr, const std::vector<double> &bounds = {})
        : n(int(c.size())), m(int(b.size())), stride(n + 2), next_id(n + m), original_n(n),
          basic(m), nonbasic(n + 1), location(n + m) {
        upper = budget ? budget->copy(bounds) : bounds;
        if (!upper.empty())
            flipped = budget ? budget->filled<unsigned char>(n) : std::vector<unsigned char>(n);
        source_cost = std::make_shared<const std::vector<double>>(budget ? budget->copy(c) : c);
        source_rows.reserve(m);
        for (int i = 0; i < m; ++i)
            source_rows.push_back(
                std::make_shared<const Row>(Row{sparse_coefficients(a[i], budget), b[i]}));
        double smallest = inf, largest = 0;
        for (const auto &r : source_rows)
            for (auto [j, value] : r->terms) {
                if (budget) budget->checkpoint(1);
                (void)j;
                if (value) {
                    smallest = std::min(smallest, std::abs(value));
                    largest = std::max(largest, std::abs(value));
                }
            }
        // 係数の桁差が大きい辞書では、初めから二段階の比率検定を使う。
        // それ以外でも小pivot・残差による再構築を契機に切り替える。
        stable_ratio = largest > 1000 * smallest;
        const int sparse_mode = budget ? budget->opt.sparse_tableau : -1;
        bool use_sparse = sparse_mode == 1;
        if (sparse_mode == -1 && std::uint64_t(n) * m >= (1ULL << 20)) {
            std::uint64_t nnz = 0;
            for (int i = 0; i < m; ++i) {
                if (budget) budget->checkpoint(1);
                nnz += count_nonzero(a[i], budget);
            }
            use_sparse = nnz * 128 <= std::uint64_t(n) * m;
        }
        if (use_sparse) {
            sparse_rows.reserve(m + 2);
            for (int i = 0; i < m; ++i) {
                auto r = std::make_shared<SparseRow>();
                auto entries = sparse_coefficients(a[i], budget);
                entries.emplace_back(n, -1);
                if (b[i]) entries.emplace_back(n + 1, b[i]);
                initialize_row(*r, entries, budget);
                sparse_bytes += r->bytes();
                sparse_rows.push_back(std::move(r));
                basic[i] = n + i;
                location[n + i] = -i - 1;
            }
            scale_objective(c);
            for (int i = 0; i < 2; ++i) {
                auto r = std::make_shared<SparseRow>();
                r->full = budget ? budget->filled<double>(stride) : std::vector<double>(stride);
                sparse_bytes += r->bytes();
                sparse_rows.push_back(std::move(r));
            }
            for (int j = 0; j < n; ++j) {
                nonbasic[j] = j;
                location[j] = j;
                sparse_rows[m]->full[j] = -c[j] * objective_scale;
            }
            nonbasic[n] = -1;
            sparse_rows[m + 1]->full[n] = 1;
            if (budget) ++budget->stats.sparse_tableaux;
            return;
        }
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

    void to_dense(Budget &budget) {
        if (sparse_rows.empty()) return;
        std::vector<std::shared_ptr<std::vector<double>>> next;
        next.reserve(m + 2);
        for (int i = 0; i < m + 2; ++i) {
            auto values = budget.filled<double>(stride);
            nonzeros(i, [&](int j, double v) {
                budget.checkpoint(1);
                values[j] = v;
            });
            next.push_back(std::make_shared<std::vector<double>>(std::move(values)));
        }
        data = std::move(next);
        sparse_rows.clear();
        sparse_bytes = 0;
        ++budget.stats.dense_conversions;
    }

    // 共有行の分離・ハッシュ容量の確保を全て済ませてから数値を更新する。
    // 数値更新は通常の密pivot同様に一括。確保中の期限切れでも基底は一貫している。
    void sparse_pivot(int r, int s, Budget &budget) {
        const double inv = 1 / coefficient(r, s);
        auto &entries = sparse_entries;
        entries.clear();
        nonzeros(r, [&](int j, double a) {
            budget.checkpoint(1);
            entries.emplace_back(j, a);
        });
        auto &affected = sparse_affected;
        affected.clear();
        for (int i = 0; i < m + 2; ++i) {
            budget.checkpoint(1);
            if (i == r || coefficient(i, s) == 0) continue;
            prepare_row(i, int(entries.size()), budget);
            affected.push_back(i);
        }
        prepare_row(r, 1, budget);
        for (int i : affected) {
            auto &q = *sparse_rows[i];
            const double f = q.get(s) * inv;
            for (auto [j, a] : entries)
                if (j != s) q.subtract_product(j, f, a);
            // 予約済み容量のためsetの追加確保を避ける。
            if (!q.full.empty())
                q.full[s] = -f;
            else
                q.put(q.position(s), s, -f);
        }
        auto &p = *sparse_rows[r];
        if (!p.full.empty()) {
            for (double &a : p.full) a *= inv;
            p.full[s] = inv;
        } else {
            for (auto &cell : p.table)
                if (cell.column >= 0) cell.value *= inv;
            p.put(p.position(s), s, inv);
        }
        std::swap(basic[r], nonbasic[s]);
        if (basic[r] >= 0) location[basic[r]] = -r - 1;
        if (nonbasic[s] >= 0) location[nonbasic[s]] = s;
        ++budget.stats.pivots;
        ++budget.stats.sparse_pivots;
    }

    void pivot(int r, int s, Budget &budget) {
        if (!sparse_rows.empty()) {
            if (budget.opt.sparse_tableau != 1 &&
                sparse_bytes * 16 > std::uint64_t(m + 2) * stride * sizeof(double))
                to_dense(budget);
            else {
                sparse_pivot(r, s, budget);
                return;
            }
        }
        double *p = write_dense_row(r, &budget);
        pivot_rows.clear();
        if (!data.empty())
            for (int i = 0; i < m + 2; ++i)
                if (i != r && dense_row(i)[s] != 0) {
                    pivot_rows.push_back(write_dense_row(i, &budget));
                }
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

    // 基底変数の上限を下限0に読み替える。実際の点・目的値は変わらない。
    void reverse_basic(int r, Budget &budget) {
        const int id = basic[r];
        if (!sparse_rows.empty()) {
            auto &q = prepare_row(r, 1, budget);
            const double rhs = upper[id] - q.get(n + 1);
            if (!q.full.empty()) {
                for (double &a : q.full) a = -a;
                q.full[n + 1] = rhs;
            } else {
                for (auto &cell : q.table)
                    if (cell.column >= 0) cell.value = -cell.value;
                q.put(q.position(n + 1), n + 1, rhs);
            }
            flipped[id] ^= 1;
            return;
        }
        double *q = write_dense_row(r, &budget);
        for (int j = 0; j <= n; ++j) q[j] = -q[j];
        q[n + 1] = upper[id] - q[n + 1];
        flipped[id] ^= 1;
    }

    // 非基底変数が反対側の上限に達した場合、基底交換なしに座標を反転。
    // 共有行の確保中は数値を変更しない。中断後も辞書は一貫したまま。
    void flip_nonbasic(int s, Budget &budget) {
        const int id = nonbasic[s];
        const double u = upper[id];
        if (!sparse_rows.empty()) {
            auto &affected = sparse_affected;
            affected.clear();
            for (int i = 0; i < m + 2; ++i) {
                budget.checkpoint(1);
                if (coefficient(i, s) == 0) continue;
                prepare_row(i, 1, budget);
                affected.push_back(i);
            }
            for (int i : affected) {
                auto &q = *sparse_rows[i];
                const double a = q.get(s);
                const double rhs = q.get(n + 1) - a * u;
                if (!q.full.empty()) {
                    q.full[n + 1] = rhs;
                    q.full[s] = -a;
                } else {
                    q.put(q.position(n + 1), n + 1, rhs);
                    q.put(q.position(s), s, -a);
                }
            }
            flipped[id] ^= 1;
            ++budget.stats.pivots;
            ++budget.stats.bound_flips;
            return;
        }
        pivot_rows.clear();
        for (int i = 0; i < m + 2; ++i)
            if (dense_row(i)[s] != 0) pivot_rows.push_back(write_dense_row(i, &budget));
        for (double *q : pivot_rows) {
            q[n + 1] -= q[s] * u;
            q[s] = -q[s];
        }
        flipped[id] ^= 1;
        ++budget.stats.pivots;
        ++budget.stats.bound_flips;
    }

    bool consistent(Budget &budget) const {
        if (!source_cost || perturbed) return true;
        // 同じ変数を各非ゼロ係数から何度もたどらず、基底の値を一度だけ読む。
        std::vector<double> values;
        values.reserve(original_n);
        budget.chunks(original_n, [&](std::size_t first, std::size_t last) {
            for (auto j = first; j < last; ++j) values.push_back(coordinate(int(j)));
        });
        for (int i = 0; i < m; ++i) {
            budget.checkpoint(1);
            const auto &src = *source_rows[i];
            double residual = src.rhs - coordinate(original_n + i);
            double magnitude = 1 + std::abs(src.rhs) + std::abs(coordinate(original_n + i));
            for (auto [j, a] : src.terms) {
                budget.checkpoint(1);
                const double x = values[j];
                residual = std::fma(-a, x, residual);
                magnitude += std::abs(a * x);
            }
            if (!std::isfinite(residual) || std::abs(residual) > 2e-12 * magnitude) return false;
        }
        double residual = value(), magnitude = 1 + std::abs(residual);
        for (int j = 0; j < original_n; ++j) {
            budget.checkpoint(1);
            const double x = values[j], c = (*source_cost)[j];
            residual = std::fma(-c, x, residual);
            magnitude += std::abs(c * x);
        }
        return std::isfinite(residual) && std::abs(residual) <= 2e-12 * magnitude;
    }
    bool dual_feasible(double eps) const {
        if (needs_repair) return false;
        const double *z = row(m);
        for (int j = 0; j <= n; ++j)
            if (nonbasic[j] >= 0 && (!std::isfinite(z[j]) || z[j] < -eps)) return false;
        return true;
    }

    // 制約の残差が正常なら、目的行だけを補償加算で作り直す。
    // 目的行の誤差のためにO(m*n*min(m,n))の基底再構築を行う費用を避ける。
    bool restore_objective(Budget &budget) {
        if (!source_cost || perturbed) return false;
        auto z = budget.filled<Accurate>(stride);
        for (int id = 0; id < original_n; ++id) {
            budget.checkpoint(1);
            double cost = (*source_cost)[id] * objective_scale;
            if (cost == 0) continue;
            if (reversed(id)) {
                z[n + 1] = z[n + 1] + Accurate(cost) * Accurate(upper[id]);
                cost = -cost;
            }
            const int pos = location[id];
            if (pos == fixed_zero) continue;
            if (pos >= 0)
                z[pos] = z[pos] - Accurate(cost);
            else
                nonzeros(-pos - 1, [&](int j, double a) {
                    budget.checkpoint(1);
                    z[j] = z[j] + Accurate(cost) * Accurate(a);
                });
        }
        auto values = budget.filled<double>(stride);
        for (int j = 0; j < stride; ++j) {
            budget.checkpoint(1);
            values[j] = z[j].value();
            if (!std::isfinite(values[j])) return false;
        }
        // 全て準備してから短い配列を一括交換。途中で目的行を半更新にしない。
        double *target = write_row(m, &budget);
        std::copy(values.begin(), values.end(), target);
        return true;
    }

    // 元の係数から現在の基底を再構築。別領域で完成させてから交換するため、
    // 途中の時間切れは元の辞書・暫定解を変更しない。
    bool refactor(Budget &budget, bool repair_basis = false) {
        if (!source_cost || perturbed) return false;
        bool repaired = false;
        for (int id : basic)
            if (id < 0) {
                if (!repair_basis) return false;
                repaired = true;
            }
        const int width = original_n + 2;
        auto work = budget.filled<Accurate>(std::size_t(m + 1) * width);
        auto get = [&](int i) { return work.data() + std::size_t(i) * width; };
        auto rows = budget.filled<int>(m);
        auto cols = budget.filled<int>(original_n + 1);
        auto wanted = budget.filled<unsigned char>(next_id);
        for (int id : basic)
            if (id >= 0) wanted[id] = 1;
        for (int i = 0; i < m; ++i) {
            rows[i] = original_n + i;
            auto *q = get(i);
            q[original_n] = -1;
            q[original_n + 1] = source_rows[i]->rhs;
            for (auto [j, a] : source_rows[i]->terms) {
                budget.checkpoint(1);
                q[j] = q[j] + Accurate(a);
            }
        }
        for (int j = 0; j < original_n; ++j) {
            cols[j] = j;
            get(m)[j] = -(*source_cost)[j] * objective_scale;
            if (reversed(j))
                for (int i = 0; i <= m; ++i) {
                    budget.checkpoint(1);
                    auto *q = get(i);
                    q[original_n + 1] = q[original_n + 1] - q[j] * Accurate(upper[j]);
                    q[j] = Accurate(0) - q[j];
                }
        }
        cols[original_n] = -1;
        while (true) {
            int r = -1, s = -1;
            double largest = 0;
            bool remaining = false;
            for (int j = 0; j < original_n; ++j)
                if (cols[j] >= 0 && wanted[cols[j]]) {
                    remaining = true;
                    for (int i = 0; i < m; ++i)
                        if (!wanted[rows[i]]) {
                            budget.checkpoint(1);
                            const double a = std::abs(get(i)[j].value());
                            if (a > largest) {
                                largest = a;
                                r = i;
                                s = j;
                            }
                        }
                }
            if (!remaining) break;
            if (r < 0 || largest < 1e-25 || !std::isfinite(largest)) {
                if (!repair_basis) return false;
                // 丸め誤差で選ばれた特異な基底は、独立な列と残ったslackで修復する。
                // 固定済みの非基底変数は復活させない。
                for (int id : rows)
                    if (location[id] == fixed_zero) return false;
                repaired = true;
                break;
            }
            auto *p = get(r);
            const Accurate inverse = Accurate(1) / p[s];
            for (int i = 0; i <= m; ++i)
                if (i != r) {
                    auto *q = get(i);
                    if (q[s].value() == 0) continue;
                    const auto f = q[s] * inverse;
                    budget.chunks(width, [&](std::size_t first, std::size_t last) {
                        for (auto j = first; j < last; ++j)
                            if (int(j) != s) q[j] = q[j] - f * p[j];
                    });
                    q[s] = Accurate(0) - f;
                }
            for (int j = 0; j < width; ++j) p[j] = p[j] * inverse;
            p[s] = inverse;
            std::swap(rows[r], cols[s]);
        }
        auto row_position = budget.filled<int>(next_id, -1);
        auto col_position = budget.filled<int>(next_id, -1);
        for (int i = 0; i < m; ++i) row_position[rows[i]] = i;
        for (int j = 0; j <= original_n; ++j)
            if (cols[j] >= 0) col_position[cols[j]] = j;
        Tableau next = clone(budget);
        if (repaired) {
            next.basic = rows;
            next.nonbasic.clear();
            for (int id : cols)
                if (id < 0 || location[id] != fixed_zero) next.nonbasic.push_back(id);
            if (next.nonbasic.size() != std::size_t(n + 1)) return false;
            for (int i = 0; i < m; ++i) next.location[next.basic[i]] = -i - 1;
            for (int j = 0; j <= n; ++j)
                if (next.nonbasic[j] >= 0) next.location[next.nonbasic[j]] = j;
        }
        for (int i = 0; i <= m; ++i) {
            double *dest = next.write_row(i, &budget);
            const auto *src = get(i == m ? m : row_position[next.basic[i]]);
            for (int j = 0; j <= n; ++j) {
                budget.checkpoint(1);
                dest[j] =
                    src[next.nonbasic[j] < 0 ? original_n : col_position[next.nonbasic[j]]].value();
                if (!std::isfinite(dest[j])) return false;
            }
            dest[n + 1] = src[original_n + 1].value();
            if (!std::isfinite(dest[n + 1])) return false;
        }
        *this = std::move(next);
        ++budget.stats.basis_refactorizations;
        return true;
    }

    Status primal(int objective_row, Budget &budget, bool allow_recovery = true) {
        const double eps = budget.opt.pivot_tolerance;
        int degenerate = 0;
        bool bland = false;
        while (true) {
            if (budget.stop()) return Status::Limit;
            if (allow_recovery && objective_row == m && budget.opt.basis_repair && !perturbed &&
                degenerate > std::max(64, m + n))
                return recover_basis(budget);
            const double *z = row(objective_row);
            int s = -1;
            for (int j = 0; j <= n; ++j) {
                if (nonbasic[j] == -1 || z[j] >= -eps) continue;
                if (s < 0 || (bland ? nonbasic[j] < nonbasic[s] : z[j] < z[s])) s = j;
            }
            if (s < 0) return Status::Optimal;
            int r = -1;
            double ratio = inf, largest = 0;
            bool upper_leaves = false;
            const double cap = bound(nonbasic[s]);
            double relaxed = cap;
            auto scan = [&](auto consider) {
                if (sparse_rows.empty()) {
                    for (int i = 0; i < m; ++i) {
                        const double *q = dense_row(i);
                        consider(i, q[s], q[n + 1]);
                    }
                } else
                    for (int i = 0; i < m; ++i)
                        consider(i, coefficient(i, s), coefficient(i, n + 1));
            };
            if (stable_ratio && !bland)
                scan([&](int i, double a, double rhs) {
                    const bool at_upper = a < -eps && std::isfinite(bound(basic[i]));
                    if (a <= eps && !at_upper) return;
                    const double slack = at_upper ? bound(basic[i]) - rhs : rhs;
                    relaxed = std::min(relaxed, (std::max(0.0, slack) + eps) / std::abs(a));
                });
            scan([&](int i, double a, double rhs) {
                const bool at_upper = a < -eps && std::isfinite(bound(basic[i]));
                if (a <= eps && !at_upper) return;
                const double v = std::max(0.0, at_upper ? bound(basic[i]) - rhs : rhs) / std::abs(a);
                if ((!stable_ratio || bland) ? ((bland ? v < ratio : v < ratio - eps) ||
                                                ((bland ? v == ratio : std::abs(v - ratio) <= eps) &&
                                                 (r < 0 || basic[i] < basic[r])))
                                             : (v <= relaxed && (r < 0 || std::abs(a) > largest))) {
                    r = i;
                    ratio = v;
                    largest = std::abs(a);
                    upper_leaves = at_upper;
                }
            });
            const bool flip = std::isfinite(cap) && (r < 0 || cap < (bland ? ratio : ratio - eps));
            if (r < 0 && !flip) return Status::Unbounded;
            if (flip) ratio = cap;
            if (ratio <= eps)
                ++degenerate;
            else {
                degenerate = 0;
                bland = false;
            }
            if (degenerate > 32) bland = true; // 退化が続いたら循環を避ける。
            if (flip)
                flip_nonbasic(s, budget);
            else {
                if (upper_leaves) reverse_basic(r, budget);
                pivot(r, s, budget);
            }
            if (!std::isfinite(value())) return Status::NumericalError;
        }
    }

    // 浮動小数点の乗数で相殺しきれない無上限変数がある場合だけ、
    // 小さな支持集合の零空間を分数なし消去で求め直す。得られた証拠は
    // 入力doubleをそのまま二進有理数として、全列・全境界で再検査する。
    // 支持・桁数の上限を超えたら通常の基底復旧へ戻り、証拠を推測しない。
    bool
    exact_infeasible(int objective_row, Budget &budget, const std::vector<Variable> &vars) const {
        try {
            std::vector<int> support;
            std::vector<double> approximate;
            const double *z = row(objective_row);
            for (int i = 0; i < m; ++i) {
                budget.checkpoint(1);
                const int pos = location[original_n + i];
                const double weight = pos >= 0 && pos != fixed_zero          ? z[pos]
                                      : pos < 0 && -pos - 1 == objective_row ? 1
                                                                             : 0;
                if (std::isfinite(weight) && weight > 0) {
                    support.push_back(i);
                    approximate.push_back(weight);
                }
            }
            const int k = int(support.size());
            if (!k || k > 96 || std::uint64_t(k) * original_n > 65536) return false;
            // 微小な乗数を従属変数に選ぶと、僅かな修正で符号が反転する。
            // 大きい乗数から消去し、小さい乗数は自由変数側に残す。
            std::vector<int> order(k);
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(),
                             [&](int a, int b) { return approximate[a] > approximate[b]; });
            auto old_support = support;
            auto old_approximate = approximate;
            for (int i = 0; i < k; ++i) {
                support[i] = old_support[order[i]];
                approximate[i] = old_approximate[order[i]];
            }
            auto exponents = budget.filled<int>(k);
            auto coefficients = budget.filled<ExactInteger>(std::size_t(original_n) * k);
            auto at = [&](int column, int i) -> ExactInteger & {
                return coefficients[std::size_t(column) * k + i];
            };
            for (int i = 0; i < k; ++i) {
                const auto &src = *source_rows[support[i]];
                for (auto [j, a] : src.terms) {
                    budget.checkpoint(1);
                    if (a) exponents[i] = std::min(exponents[i], ExactInteger::dyadic(a).second);
                }
                for (auto [j, a] : src.terms) {
                    budget.checkpoint(1);
                    const auto part = ExactInteger::dyadic(a);
                    at(j, i) =
                        at(j, i) + ExactInteger(part.first).shifted(part.second - exponents[i]);
                }
            }
            auto sums = budget.filled<double>(original_n);
            auto magnitudes = budget.filled<double>(original_n);
            for (int i = 0; i < k; ++i)
                for (auto [j, a] : source_rows[support[i]]->terms) {
                    budget.checkpoint(1);
                    sums[j] += approximate[i] * a;
                    magnitudes[j] += std::abs(approximate[i] * a);
                }
            std::vector<std::vector<ExactInteger>> equations;
            for (int j = 0; j < original_n; ++j)
                if (!std::isfinite(vars[j].upper) &&
                    std::abs(sums[j]) <= 1e-7 * std::max(1.0, magnitudes[j])) {
                    budget.checkpoint();
                    equations.emplace_back(coefficients.begin() + std::size_t(j) * k,
                                           coefficients.begin() + std::size_t(j + 1) * k);
                }
            std::vector<int> pivots;
            ExactInteger previous(1);
            for (int column = 0; column < k && pivots.size() < equations.size(); ++column) {
                const int rank = int(pivots.size());
                int chosen = rank;
                while (chosen < int(equations.size()) && equations[chosen][column].zero())
                    ++chosen;
                if (chosen == int(equations.size())) continue;
                std::swap(equations[rank], equations[chosen]);
                const auto pivot_value = equations[rank][column];
                for (int i = rank + 1; i < int(equations.size()); ++i) {
                    const auto factor = equations[i][column];
                    for (int j = column + 1; j < k; ++j) {
                        budget.checkpoint(pivot_value.digits.size() > 16 ? Budget::check_interval
                                                                         : 64);
                        equations[i][j] =
                            (pivot_value * equations[i][j] - factor * equations[rank][j])
                                .divided(previous);
                    }
                    equations[i][column] = ExactInteger{};
                }
                previous = pivot_value;
                pivots.push_back(column);
            }
            if (int(pivots.size()) == k) return false;
            auto weights = budget.filled<ExactInteger>(k);
            int weight_exponent = 0;
            for (int i = 0; i < k; ++i)
                weight_exponent = std::min(
                    weight_exponent, ExactInteger::dyadic(approximate[i]).second + exponents[i]);
            for (int i = 0; i < k; ++i)
                if (!std::binary_search(pivots.begin(), pivots.end(), i)) {
                    const auto part = ExactInteger::dyadic(approximate[i]);
                    weights[i] =
                        previous *
                        ExactInteger(part.first).shifted(part.second + exponents[i] - weight_exponent);
                }
            for (int r = int(pivots.size()); r--;) {
                const int column = pivots[r];
                ExactInteger sum;
                for (int j = column + 1; j < k; ++j) {
                    budget.checkpoint(64);
                    sum = sum - equations[r][j] * weights[j];
                }
                weights[column] = sum.divided(equations[r][column]);
            }
            if (previous.negative)
                for (auto &w : weights) w = -w;
            for (const auto &w : weights)
                if (w.negative) return false;
            int exponent = 0;
            for (int i = 0; i < k; ++i)
                exponent =
                    std::min(exponent, ExactInteger::dyadic(source_rows[support[i]]->rhs).second -
                                           exponents[i]);
            for (const auto &v : vars) {
                budget.checkpoint(1);
                if (std::isfinite(v.lower))
                    exponent = std::min(exponent, ExactInteger::dyadic(v.lower).second);
                if (std::isfinite(v.upper))
                    exponent = std::min(exponent, ExactInteger::dyadic(v.upper).second);
            }
            ExactInteger rhs, magnitude = ExactInteger(1).shifted(-exponent);
            for (int i = 0; i < k; ++i) {
                budget.checkpoint(64);
                const auto part = ExactInteger::dyadic(source_rows[support[i]]->rhs);
                const auto term = (weights[i] * ExactInteger(part.first))
                                      .shifted(part.second - exponents[i] - exponent);
                rhs = rhs + term;
                magnitude = magnitude + term.absolute();
            }
            for (int j = 0; j < original_n; ++j) {
                ExactInteger coefficient;
                for (int i = 0; i < k; ++i) {
                    budget.checkpoint(64);
                    coefficient = coefficient + weights[i] * at(j, i);
                }
                if (coefficient.zero()) continue;
                const double endpoint = coefficient.negative ? vars[j].upper : vars[j].lower;
                if (!std::isfinite(endpoint)) return false;
                const auto part = ExactInteger::dyadic(endpoint);
                const auto term =
                    (coefficient * ExactInteger(part.first)).shifted(part.second - exponent);
                rhs = rhs - term;
                magnitude = magnitude + term.absolute();
            }
            return rhs.negative &&
                   ExactInteger::compare((-rhs) * ExactInteger(1000000000000LL), magnitude) > 0;
        } catch (const ExactOverflow &) { return false; }
    }

    // Phase I の負の目的値だけで実行不能としない。非負の行乗数を
    // 元係数に掛け、区間境界上の最小活動量より右辺が小さいことを確認する。
    // 辞書の係数誤差で作られた見かけの矛盾を、探索の証拠に使わない。
    bool verify_infeasible(int objective_row, Budget &budget,
                           const std::vector<Variable> *enclosure = nullptr) const {
        auto endpoint = [&](int j, double a) {
            if (enclosure) return a < 0 ? (*enclosure)[j].upper : (*enclosure)[j].lower;
            return location[j] == fixed_zero ? (reversed(j) ? upper[j] : 0) : a < 0 ? bound(j) : 0;
        };
        auto propagate_certificate = [&] {
            if (enclosure) return exact_infeasible(objective_row, budget, *enclosure);
            auto vars = budget.filled<Variable>(original_n, Variable{0, inf, VarType::Continuous});
            for (int j = 0; j < original_n; ++j) {
                vars[j].upper = bound(j);
                if (location[j] == fixed_zero)
                    vars[j].lower = vars[j].upper = reversed(j) ? upper[j] : 0;
            }
            std::vector<Row> rows;
            rows.reserve(m);
            for (int i = 0; i < m; ++i) {
                rows.push_back(Row{budget.copy(source_rows[i]->terms), source_rows[i]->rhs});
                if (location[original_n + i] == fixed_zero) {
                    auto terms = budget.copy(source_rows[i]->terms);
                    for (auto &term : terms) term.second = -term.second;
                    rows.push_back(Row{std::move(terms), -source_rows[i]->rhs});
                }
            }
            if (!propagate(rows, vars, budget, 8)) return true;
            return verify_infeasible(objective_row, budget, &vars);
        };
        // 単一行と区間境界だけの矛盾も、正当な非負乗数の証拠になる。
        for (const auto &row : source_rows) {
            Accurate minimum;
            double magnitude = 1 + std::abs(row->rhs);
            bool finite = true;
            for (auto [j, a] : row->terms) {
                budget.checkpoint(1);
                const double bound_value = endpoint(j, a);
                if (!std::isfinite(bound_value)) {
                    finite = false;
                    break;
                }
                minimum = minimum + Accurate(a) * Accurate(bound_value);
                magnitude += std::abs(a * bound_value);
            }
            if (finite && std::isfinite(minimum.value()) &&
                (minimum - Accurate(row->rhs)).value() >
                    budget.opt.feasibility_tolerance + 1e-12 * magnitude)
                return true;
        }
        const double *z = row(objective_row);
        double largest = 0;
        auto weights = budget.filled<double>(m);
        for (int i = 0; i < m; ++i) {
            budget.checkpoint(1);
            const int pos = location[original_n + i];
            if (pos >= 0 && pos != fixed_zero && std::isfinite(z[pos]))
                weights[i] = std::max(0.0, z[pos]);
            else if (pos < 0 && -pos - 1 == objective_row)
                weights[i] = 1;
            largest = std::max(largest, weights[i]);
        }
        if (largest == 0) return propagate_certificate();
        // 乗数の消去誤差が無上限変数に微小な負係数を作る場合がある。
        // 微小な乗数を落とした別の証拠を試すが、毎回元の全係数で検証する。
        // 合成後の負係数を0と見なしたり、変数の上限を仮定したりはしない。
        for (int attempt = 0; attempt < 5; ++attempt) {
            const double cutoffs[] = {0., 1e-14, 1e-12, 1e-10, 0.};
            const double cutoff = cutoffs[attempt];
            if (attempt == 4) {
                // 有理な乗数をdoubleに丸めると、相殺する列に負の残差が残る。
                // 連分数で共通の整数倍率を提案し、丸めた証拠を元係数で再検証する。
                auto numerators = budget.filled<long long>(m);
                auto denominators = budget.filled<long long>(m, 1);
                long long common = 1;
                bool representable = true;
                for (int i = 0; i < m; ++i) {
                    budget.checkpoint(1);
                    const double target = weights[i] / largest;
                    if (target < 1e-12) continue;
                    double value = target;
                    long long h0 = 0, h1 = 1, k0 = 1, k1 = 0;
                    bool found = false;
                    for (int k = 0; k < 32 && value <= 1000000; ++k) {
                        const auto integer = (long long)std::floor(value);
                        const long long h = integer * h1 + h0, d = integer * k1 + k0;
                        if (d > 1000000 || !d) break;
                        if (std::abs(target - double(h) / d) < 1e-10) {
                            numerators[i] = h;
                            denominators[i] = d;
                            found = true;
                            break;
                        }
                        const double fraction = value - integer;
                        if (fraction <= 0) break;
                        value = 1 / fraction;
                        h0 = h1;
                        h1 = h;
                        k0 = k1;
                        k1 = d;
                    }
                    if (!found) {
                        representable = false;
                        break;
                    }
                    const auto factor = denominators[i] / std::gcd(common, denominators[i]);
                    if (common > 1000000000000LL / factor) {
                        representable = false;
                        break;
                    }
                    common *= factor;
                }
                if (!representable) break;
                for (int i = 0; i < m; ++i)
                    weights[i] = double(numerators[i] * (common / denominators[i]));
            }
            auto coefficients = budget.filled<Accurate>(original_n);
            auto absolute = budget.filled<double>(original_n);
            Accurate rhs;
            double magnitude = 1;
            for (int i = 0; i < m; ++i) {
                budget.checkpoint(1);
                const double y = weights[i];
                if (!std::isfinite(y)) return false;
                if (y == 0 || y < cutoff * largest) continue;
                const auto &source = *source_rows[i];
                rhs = rhs + Accurate(y) * Accurate(source.rhs);
                magnitude += std::abs(y * source.rhs);
                for (auto [j, coefficient] : source.terms) {
                    budget.checkpoint(1);
                    coefficients[j] = coefficients[j] + Accurate(y) * Accurate(coefficient);
                    absolute[j] += std::abs(y * coefficient);
                }
            }
            bool bounded = true;
            for (int j = 0; j < original_n; ++j) {
                budget.checkpoint(1);
                const double coefficient = coefficients[j].value();
                if (!std::isfinite(coefficient)) return false;
                const double bound_value = endpoint(j, coefficient);
                if (!std::isfinite(bound_value)) {
                    bounded = false;
                    break;
                }
                rhs = rhs - coefficients[j] * Accurate(bound_value);
                magnitude += absolute[j] * std::abs(bound_value);
            }
            if (bounded && std::isfinite(rhs.value()) && std::isfinite(magnitude) &&
                rhs.value() < -1e-12 * magnitude)
                return true;
        }
        return propagate_certificate();
    }

    // 任意の正しい辞書から実行可能性を回復する。違反行だけに人工列を置き、
    // 人工変数を最大違反量にすると全ての基本変数が境界上になるようにする。
    // 有限上限がある基底でも、Phase I の開始点が実行可能となる。
    Status reoptimize(Budget &budget, int repair_attempts = 0) {
        if (budget.stop()) return Status::Limit;
        primal_feasible = false;
        int r = -1, artificial = -1;
        double violation = 0;
        for (int j = 0; j <= n; ++j)
            if (nonbasic[j] == -1) artificial = j;
        if (artificial < 0) return Status::NumericalError;
        for (int i = 0; i < m; ++i) {
            budget.checkpoint(1);
            if (coefficient(i, n + 1) > bound(basic[i])) reverse_basic(i, budget);
            if (coefficient(i, n + 1) < -violation) {
                violation = -coefficient(i, n + 1);
                r = i;
            }
        }
        if (r >= 0 && violation > budget.opt.pivot_tolerance) {
            for (int i = 0; i < m; ++i) {
                double *q = write_row(i, &budget);
                q[artificial] = std::min(0.0, q[n + 1]) / violation;
            }
            double *z = write_row(m + 1, &budget);
            std::fill(z, z + stride, 0);
            z[artificial] = 1;
            write_row(m, &budget)[artificial] = 0;
            pivot(r, artificial, budget);
            Status status = primal(m + 1, budget, false);
            if (status != Status::Optimal) return status;
            if (coefficient(m + 1, n + 1) < -budget.opt.feasibility_tolerance) {
                if (verify_infeasible(m + 1, budget)) return Status::Infeasible;
                if (budget.opt.basis_repair && repair_attempts < 2 && refactor(budget, true)) {
                    Options reliable = budget.opt;
                    reliable.dual_devex = false;
                    Budget retry{reliable, budget.stats, budget.start};
                    return reoptimize(retry, repair_attempts + 1);
                }
                return Status::NumericalError;
            }
            for (int i = 0; i < m; ++i)
                if (basic[i] == -1) {
                    int s = -1;
                    double largest = 0;
                    nonzeros(i, [&](int j, double a) {
                        if (j <= n && nonbasic[j] >= 0 && std::abs(a) > largest) {
                            s = j;
                            largest = std::abs(a);
                        }
                    });
                    if (s < 0 || largest <= budget.opt.pivot_tolerance) {
                        // 冗長行に残った人工変数は微小係数で無理にpivotせず、
                        // 元の列とslackから独立な基底を作る。
                        if (repair_attempts >= 2 || !refactor(budget, true))
                            return Status::NumericalError;
                        return reoptimize(budget, repair_attempts + 1);
                    }
                    pivot(i, s, budget);
                }
        }
        if (!restore_objective(budget)) return Status::NumericalError;
        primal_feasible = true;
        return primal(m, budget, false);
    }

    // 修復した基底でも残差が消えないときは、同じ制約・固定を初期基底から解く。
    // この再試行では再帰的な修復を無効にし、総時間と反復数をそのまま引き継ぐ。
    Status restart_basis(Budget &budget) {
        if (!source_cost) return Status::NumericalError;
        needs_repair = true;
        if (budget.stop()) return Status::Limit;
        std::vector<Terms> a;
        std::vector<double> b;
        a.reserve(m);
        b.reserve(m);
        for (int i = 0; i < m; ++i) {
            a.push_back(budget.copy(source_rows[i]->terms));
            b.push_back(source_rows[i]->rhs);
        }
        for (int i = 0; i < m; ++i)
            if (location[original_n + i] == fixed_zero) {
                auto terms = budget.copy(source_rows[i]->terms);
                for (auto &term : terms) term.second = -term.second;
                a.push_back(std::move(terms));
                b.push_back(-source_rows[i]->rhs);
            }
        auto bounds = upper.empty() ? budget.filled<double>(original_n, inf) : budget.copy(upper);
        for (int j = 0; j < original_n; ++j)
            if (location[j] == fixed_zero) {
                if (reversed(j)) {
                    a.push_back(Terms{{j, -1}});
                    b.push_back(-upper[j]);
                } else
                    bounds[j] = 0;
            }
        Options reliable = budget.opt;
        reliable.basis_repair = false;
        reliable.dual_devex = false;
        Budget retry{reliable, budget.stats, budget.start};
        Tableau fresh(a, b, *source_cost, &retry, bounds);
        // 追加した等式のslackについて、外部のstep配列には刻みがない。
        fresh.slack_fixing = slack_fixing && a.size() == std::size_t(m);
        ++budget.stats.numerical_restarts;
        Status status = fresh.solve(retry);
        if (status == Status::Optimal && !fresh.consistent(retry)) status = Status::NumericalError;
        fresh.needs_repair = status != Status::Optimal;
        *this = std::move(fresh);
        return status;
    }

    Status recover_basis(Budget &budget, int max_pivots = std::numeric_limits<int>::max()) {
        if (!budget.opt.basis_repair) return Status::NumericalError;
        needs_repair = true;
        if (max_pivots <= 0 || budget.stop()) return Status::Limit;
        stable_ratio = true;
        Options local = budget.opt;
        // 行の縮尺で増幅された違反を、復旧時だけ細かいpivot・Phase I判定で扱う。
        local.pivot_tolerance = std::min(local.pivot_tolerance, 1e-13);
        local.feasibility_tolerance = std::min(local.feasibility_tolerance, local.pivot_tolerance);
        local.iteration_limit =
            budget.stats.pivots +
            std::min(std::uint64_t(std::max(0, max_pivots)),
                     local.iteration_limit - std::min(local.iteration_limit, budget.stats.pivots));
        Budget recovery{local, budget.stats, budget.start};
        if (!refactor(recovery, true)) return restart_basis(recovery);
        Status status = reoptimize(recovery);
        if (status == Status::Optimal && !consistent(recovery)) status = Status::NumericalError;
        needs_repair = status != Status::Optimal;
        return status == Status::NumericalError ? restart_basis(recovery) : status;
    }

    Status solve(Budget &budget) {
        if (budget.stop()) return Status::Limit;
        int r = -1;
        for (int i = 0; i < m; ++i)
            if (r < 0 || coefficient(i, n + 1) < coefficient(r, n + 1)) r = i;
        if (r >= 0 && coefficient(r, n + 1) < -budget.opt.pivot_tolerance) {
            // 最小費用フローや被覆では、slack 基底が既に双対実行可能。
            // 人工変数で全行を混ぜる Phase I を省き、直接実行可能性を回復する。
            bool dual_feasible = true;
            for (int j = 0; j <= n; ++j)
                if (nonbasic[j] != -1 && row(m)[j] < 0) dual_feasible = false;
            if (dual_feasible) {
                auto violation = [&] {
                    double sum = 0;
                    for (int i = 0; i < m; ++i) {
                        budget.checkpoint(1);
                        const double rhs = coefficient(i, n + 1);
                        sum += std::max(0.0, std::max(-rhs, rhs - bound(basic[i])));
                    }
                    return sum;
                };
                double initial = violation();
                Status status = dual(budget, 256, true);
                if (status != Status::Limit || budget.stop()) return status;
                if (!stable_ratio || violation() < .5 * initial)
                    return dual(budget, std::numeric_limits<int>::max(), true);
                // 退化した初期双対法だけに小さい正の摂動を加える。
                // 終了後は原係数から真の目的を復元し、主単体法で最適化し直す。
                // 摂動した目的の境界を証明用の探索へ渡してはならない。
                double *z = write_row(m, &budget);
                for (int j = 0; j <= n; ++j)
                    if (nonbasic[j] >= 0)
                        z[j] += 1e-6 * (1 + ((unsigned(nonbasic[j]) * 2654435761U) >> 16) / 65536.0);
                perturbed = true;
                ++budget.stats.objective_perturbations;
                status = dual(budget, std::numeric_limits<int>::max(), true);
                perturbed = false;
                if (status != Status::Optimal) return status;
                if (!refactor(budget)) return recover_basis(budget);
                primal_feasible = true;
                return primal(m, budget);
            }
            pivot(r, n, budget);
            Status s = primal(m + 1, budget);
            if (s == Status::Limit || s == Status::NumericalError) return s;
            if (s != Status::Optimal) return Status::NumericalError;
            // 人工目的の零判定は実行可能性の許容誤差で行う。pivot の選択用の
            // 小さい閾値では、消去の丸め残差だけで可解な問題を実行不能にしてしまう。
            if (row(m + 1)[n + 1] < -budget.opt.feasibility_tolerance) {
                if (verify_infeasible(m + 1, budget)) return Status::Infeasible;
                return recover_basis(budget);
            }
            for (int i = 0; i < m; ++i)
                if (basic[i] == -1) {
                    int jbest = -1;
                    nonzeros(i, [&](int j, double a) {
                        if (j <= n && nonbasic[j] != -1 && std::abs(a) > budget.opt.pivot_tolerance &&
                            (jbest < 0 || std::abs(a) > std::abs(coefficient(i, jbest)) ||
                             (std::abs(a) == std::abs(coefficient(i, jbest)) && j < jbest)))
                            jbest = j;
                    });
                    if (jbest >= 0) pivot(i, jbest, budget);
                }
        }
        // 人工変数が非基底なら選ばない。冗長行の基底に残った場合も値は 0。
        primal_feasible = true;
        return primal(m, budget);
    }

    std::vector<double> solution() const {
        std::vector<double> x(original_n);
        for (int j = 0; j < original_n; ++j) x[j] = coordinate(j);
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
            if (id >= 0 && (slack_fixing || id < original_n) && step[id] > 0 &&
                z[j] / objective_scale * step[id] > gap + margin)
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
        std::vector<std::shared_ptr<SparseRow>> compact_sparse;
        std::uint64_t compact_bytes = 0;
        std::vector<std::shared_ptr<std::vector<double>>> compact;
        std::vector<double> compact_dense;
        if (!sparse_rows.empty()) {
            auto remap = budget.filled<int>(stride, -1);
            for (int j = 0; j < int(keep.size()); ++j) remap[keep[j]] = j;
            remap[n + 1] = width - 1;
            compact_sparse.reserve(m + 2);
            for (int i = 0; i < m + 2; ++i) {
                auto next = std::make_shared<SparseRow>();
                Terms entries;
                nonzeros(i, [&](int j, double a) {
                    budget.checkpoint(1);
                    if (remap[j] >= 0) entries.emplace_back(remap[j], a);
                });
                if (i >= m || entries.size() * 8 > std::size_t(width)) {
                    next->full = budget.filled<double>(width);
                    for (auto [j, a] : entries) next->full[j] = a;
                } else {
                    next->rhs_column = width - 1;
                    next->reserve(int(entries.size()), &budget);
                    for (auto [j, a] : entries) next->set(j, a);
                }
                compact_bytes += next->bytes();
                compact_sparse.push_back(std::move(next));
            }
        } else {
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
        }
        std::vector<int> columns;
        columns.reserve(keep.size());
        for (int j : keep) {
            const int id = nonbasic[j];
            if (id >= 0) positions[id] = int(columns.size());
            columns.push_back(id);
        }
        sparse_rows = std::move(compact_sparse);
        sparse_bytes = compact_bytes;
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
                        double multiplier = 1, bool oriented = false) {
        Row origin;
        origin.rhs = rhs;
        bool slacks = false;
        for (const auto &term : a) slacks |= term.first >= original_n;
        if (!slacks) {
            for (const auto &term : a) {
                if (budget) budget->checkpoint(1);
                double v = term.second * multiplier;
                if (oriented && reversed(term.first)) {
                    origin.rhs -= v * upper[term.first];
                    v = -v;
                }
                origin.terms.emplace_back(term.first, v);
            }
        } else {
            std::vector<Accurate> coefficients(original_n);
            Accurate b(rhs);
            for (const auto &term : a) {
                if (budget) budget->checkpoint(1);
                const int id = term.first;
                const double v = term.second * multiplier;
                if (id < original_n) {
                    if (oriented && reversed(id)) {
                        b = b - Accurate(v) * Accurate(upper[id]);
                        coefficients[id] = coefficients[id] - Accurate(v);
                    } else
                        coefficients[id] = coefficients[id] + Accurate(v);
                } else {
                    const auto &row = *source_rows[id - original_n];
                    b = b - Accurate(v) * Accurate(row.rhs);
                    for (auto [j, c] : row.terms) {
                        if (budget) budget->checkpoint(1);
                        coefficients[j] = coefficients[j] - Accurate(v) * Accurate(c);
                    }
                }
            }
            origin.rhs = b.value();
            for (int j = 0; j < original_n; ++j)
                if (coefficients[j].value()) origin.terms.emplace_back(j, coefficients[j].value());
        }
        auto saved = std::make_shared<const Row>(std::move(origin));
        primal_feasible = false;
        // 時間割込みのない狭い辞書では、追加先に直接組み立てて一時行を省く。
        // budget ありの場合は従来通り、完成した行だけを辞書へ追加する。
        const bool inplace = sparse_rows.empty() && data.empty() && !budget;
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
            double coef = term.second * multiplier;
            if (!oriented && reversed(var)) {
                q[n + 1] -= coef * upper[var];
                coef = -coef;
            }
            const int pos = location[var];
            if (pos == fixed_zero) continue;
            if (pos >= 0)
                q[pos] += coef;
            else if (!sparse_rows.empty()) {
                nonzeros(-pos - 1, [&](int j, double v) {
                    if (budget) budget->checkpoint(1);
                    q[j] -= coef * v;
                });
            } else {
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
        if (!sparse_rows.empty()) {
            auto next = std::make_shared<SparseRow>();
            initialize_row(*next, sparse_coefficients(buffer, budget), budget);
            sparse_bytes += next->bytes();
            sparse_rows.insert(sparse_rows.begin() + m, std::move(next));
        } else if (data.empty() && !inplace) {
            dense.resize(std::size_t(m + 3) * stride);
            std::copy_backward(write_row(m), write_row(m + 2), write_row(m + 3));
            std::copy(buffer.begin(), buffer.end(), write_row(m));
        } else if (!data.empty())
            data.insert(data.begin() + m, std::make_shared<std::vector<double>>(std::move(buffer)));
        source_rows.push_back(std::move(saved));
        basic.push_back(next_id++);
        location.push_back(-m - 1);
        ++m;
    }

    Status
    dual(Budget &budget, int max_pivots = std::numeric_limits<int>::max(), bool initial = false) {
        // 中断したPhase I/主単体法の辞書は、双対実行可能とは限らない。
        if (!perturbed && !dual_feasible(10 * budget.opt.pivot_tolerance)) {
            if (needs_repair || !restore_objective(budget) ||
                !dual_feasible(10 * budget.opt.pivot_tolerance))
                return recover_basis(budget, max_pivots);
        }
        const double eps = budget.opt.pivot_tolerance;
        const double feasibility = budget.opt.feasibility_tolerance;
        int iter = 0, repairs = 0, stagnant = 0;
        double best_violation = inf;
        const bool weighted = initial && budget.opt.dual_devex && m >= 32 && n >= 16;
        std::vector<double> weights;
        std::vector<unsigned char> reference;
        int framework_age = 0;
        auto reset_weights = [&] {
            weights = budget.filled<double>(m, 1);
            reference = budget.filled<unsigned char>(next_id);
            for (int id : basic)
                if (id >= 0) reference[id] = 1;
            framework_age = 0;
        };
        if (weighted) reset_weights();
        bool fresh_basis = false;
        int degenerate = 0;
        bool bland = false;
        while (true) {
            if (budget.stop()) return Status::Limit;
            int r = -1;
            double best_merit = 0, total_violation = 0;
            bool upper_leaves = false;
            auto consider = [&](int i, double rhs) {
                const double upper_slack = bound(basic[i]) - rhs;
                const double violation = std::min(rhs, upper_slack);
                total_violation += std::max(0.0, -violation);
                const double merit = weighted ? -violation / std::sqrt(weights[i]) : -violation;
                if (violation < -eps &&
                    (r < 0 || (bland ? basic[i] < basic[r] : merit > best_merit))) {
                    // 行の縮尺を元に戻すと微小な違反も増幅するため、修復できるなら続ける。
                    // 入る列がない場合だけ、許容幅内の丸め残差を実行不能の証拠にしない。
                    if (violation >= -feasibility) {
                        bool repairable = false;
                        budget.checkpoint(std::size_t(n) + 2);
                        nonzeros(i, [&](int j, double a) {
                            if (j <= n && nonbasic[j] != -1 &&
                                (upper_slack < rhs ? a > eps : a < -eps))
                                repairable = true;
                        });
                        if (!repairable) return;
                    }
                    r = i;
                    best_merit = merit;
                    upper_leaves = upper_slack < rhs;
                }
            };
            if (sparse_rows.empty()) {
                for (int i = 0; i < m; ++i) consider(i, dense_row(i)[n + 1]);
            } else
                for (int i = 0; i < m; ++i) consider(i, coefficient(i, n + 1));
            if (r < 0) {
                if (!consistent(budget)) {
                    stable_ratio = true;
                    if (++repairs > 3 || !refactor(budget))
                        return recover_basis(budget, max_pivots - iter);
                    continue;
                }
                primal_feasible = true;
                if (!dual_feasible(10 * eps)) {
                    // RHSだけが一致していてもreduced costには誤差が残り得る。
                    // 誤差を含む目的を最適化せず、先に真の目的を再構築する。
                    if (!perturbed && !restore_objective(budget)) return Status::NumericalError;
                    Options local;
                    local.dual_devex = budget.opt.dual_devex;
                    local.basis_repair = budget.opt.basis_repair;
                    local.time_limit = budget.opt.time_limit;
                    local.pivot_tolerance = eps;
                    local.sparse_tableau = budget.opt.sparse_tableau;
                    local.iteration_limit =
                        budget.stats.pivots +
                        std::min(std::uint64_t(std::max(0, max_pivots - iter)),
                                 budget.opt.iteration_limit -
                                     std::min(budget.opt.iteration_limit, budget.stats.pivots));
                    Budget cleanup{local, budget.stats, budget.start};
                    const Status status = primal(m, cleanup);
                    if (status == Status::Optimal && !consistent(budget))
                        return recover_basis(cleanup);
                    return status;
                }
                return Status::Optimal;
            }
            if (iter++ >= max_pivots) return Status::Limit;
            if (total_violation < best_violation * (1 - 1e-8)) {
                best_violation = total_violation;
                stagnant = 0;
            } else
                ++stagnant;
            if (budget.opt.basis_repair && !perturbed && stagnant > std::max(64, m + n))
                return recover_basis(budget, max_pivots - iter);
            if (upper_leaves) reverse_basic(r, budget);
            int s = -1;
            double ratio = inf, largest = 0, relaxed = inf;
            const double *z = row(m);
            const double *q = row(r);
            if (stable_ratio && !bland)
                for (int j = 0; j <= n; ++j) {
                    if (nonbasic[j] != -1 && q[j] < -eps)
                        relaxed = std::min(relaxed, (std::max(0.0, z[j]) + eps) / -q[j]);
                }
            for (int j = 0; j <= n; ++j) {
                const double a = q[j];
                if (nonbasic[j] == -1 || a >= -eps) continue;
                const double v = std::max(0.0, z[j]) / -a;
                if ((!stable_ratio || bland) ? ((bland ? v < ratio : v < ratio - eps) ||
                                                ((bland ? v == ratio : std::abs(v - ratio) <= eps) &&
                                                 (s < 0 || nonbasic[j] < nonbasic[s])))
                                             : (v <= relaxed && (s < 0 || -a > largest))) {
                    s = j;
                    ratio = v;
                    largest = -a;
                }
            }
            if (s < 0) {
                if (verify_infeasible(r, budget)) return Status::Infeasible;
                if (++repairs > 3 || !refactor(budget))
                    return recover_basis(budget, max_pivots - iter);
                fresh_basis = true;
                stable_ratio = true;
                continue;
            }
            if (!perturbed && !fresh_basis && std::abs(coefficient(r, s)) < 1e-7) {
                fresh_basis = true;
                stable_ratio = true;
                if (!refactor(budget)) return recover_basis(budget, max_pivots - iter);
                continue;
            }
            if (ratio <= eps)
                ++degenerate;
            else {
                degenerate = 0;
                bland = false;
            }
            if (degenerate > 32) bland = true;
            if (weighted) {
                // 参照基底の列だけで選択行の正確な辺長を求める。
                // 他の行は既にあるpivot列からmax更新し、追加の線形方程式を解かない。
                double exact = basic[r] >= 0 ? reference[basic[r]] : 0;
                nonzeros(r, [&](int j, double a) {
                    budget.checkpoint(1);
                    if (j <= n && nonbasic[j] >= 0 && reference[nonbasic[j]]) exact += a * a;
                });
                exact = std::max(1.0, exact);
                if (!std::isfinite(exact) || exact > 4 * weights[r] || weights[r] > 4 * exact ||
                    framework_age >= m) {
                    reset_weights();
                    exact = 1;
                }
                const double pivot_value = coefficient(r, s);
                const double next_weight = exact / (pivot_value * pivot_value);
                for (int i = 0; i < m; ++i) {
                    budget.checkpoint(1);
                    const double a = coefficient(i, s);
                    weights[i] = std::max(weights[i], next_weight * a * a);
                }
                weights[r] = std::max(1.0, next_weight);
                ++framework_age;
            }
            pivot(r, s, budget);
            fresh_basis = false;
            if (!std::isfinite(value())) return Status::NumericalError;
        }
    }
};

// 通常の変数変換は高々2項を直接保持。等式代入で長くなる式だけ共有配列に置く。
struct Mapping {
    struct Columns {
        std::pair<int, double> entries[2];
        int count = 0;
        std::shared_ptr<const Terms> expanded;
        const std::pair<int, double> *begin() const {
            return expanded ? expanded->data() : entries;
        }
        const std::pair<int, double> *end() const { return begin() + count; }
        int size() const { return count; }
        bool empty() const { return count == 0; }
        const auto &operator[](int j) const { return begin()[j]; }
        void emplace_back(int j, double a) { entries[count++] = {j, a}; }
        void assign(Terms a) {
            count = int(a.size());
            if (count <= 2) {
                expanded.reset();
                std::copy(a.begin(), a.end(), entries);
            } else
                expanded = std::make_shared<const Terms>(std::move(a));
        }
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
    std::vector<double> b, c, upper;
    double offset = 0;
    bool infeasible = false;
    int substitutions = 0;
    Prepared clone(Budget &budget) const {
        Prepared p;
        p.vars = budget.copy(vars);
        p.map = budget.copy(map);
        p.b = budget.copy(b);
        p.c = budget.copy(c);
        p.upper = budget.copy(upper);
        p.offset = offset;
        p.infeasible = infeasible;
        p.substitutions = substitutions;
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
            if (id < tab.original_n) {
                if (tab.reversed(id)) {
                    rhs -= (long double)coef * tab.upper[id];
                    a[id] -= coef;
                } else
                    a[id] += coef;
            } else {
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
        tab.add_constraint(cut.terms, -1 + 1e-8, nullptr, 1, true);
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

// 変わった境界を含む行だけを再検査する。作業上限で止めても、得た境界は有効。
// 仮固定では同じ索引を再利用し、元の境界は変更しない。
struct Propagation {
    const std::vector<Row> &rows;
    std::vector<std::vector<int>> incidence;
    std::vector<int> queue;
    std::vector<unsigned char> queued;
    std::uint64_t size = 0, changes = 0;
    Propagation(const std::vector<Row> &r, std::size_t n, Budget &budget) : rows(r) {
        incidence = budget.filled<std::vector<int>>(n);
        queued = budget.filled<unsigned char>(r.size());
        for (int i = 0; i < int(r.size()); ++i) {
            budget.checkpoint(1);
            size += 1 + r[i].terms.size();
            for (auto [j, a] : r[i].terms) {
                budget.checkpoint(1);
                (void)a;
                incidence[j].push_back(i);
            }
        }
    }
    bool run(std::vector<Variable> &vars, Budget &budget, std::uint64_t &work, int seed = -1) {
        queue.clear();
        budget.chunks(queued.size(), [&](std::size_t a, std::size_t b) {
            std::fill(queued.begin() + a, queued.begin() + b, 0);
        });
        auto enqueue = [&](int i) {
            if (!queued[i]) {
                queued[i] = 1;
                queue.push_back(i);
            }
        };
        if (seed < 0)
            for (int i = 0; i < int(rows.size()); ++i) {
                budget.checkpoint(1);
                enqueue(i);
            }
        else
            for (int i : incidence[seed]) {
                budget.checkpoint(1);
                enqueue(i);
            }
        const auto &opt = budget.opt;
        for (std::size_t head = 0; head < queue.size() && work; ++head) {
            const int i = queue[head];
            queued[i] = 0;
            const Row &row = rows[i];
            const auto cost = 1 + row.terms.size();
            if (cost > work) break;
            work -= cost;
            budget.checkpoint(1);
            long double sum = 0, magnitude = std::abs(row.rhs);
            int infinite = 0;
            for (auto [j, a] : row.terms) {
                budget.checkpoint(1);
                if (a == 0) continue;
                const double b = a > 0 ? vars[j].lower : vars[j].upper;
                if (std::isfinite(b)) {
                    sum += (long double)a * b;
                    magnitude += std::abs((long double)a * b);
                } else
                    ++infinite;
            }
            if (!infinite &&
                sum > (long double)row.rhs + opt.feasibility_tolerance + 1e-12L * magnitude)
                return false;
            for (auto [j, a] : row.terms) {
                budget.checkpoint(1);
                if (a == 0) continue;
                auto &v = vars[j];
                const double old = a > 0 ? v.lower : v.upper;
                const bool finite = std::isfinite(old);
                if (infinite - !finite) continue;
                const long double other = sum - (finite ? (long double)a * old : 0);
                double b = double(((long double)row.rhs - other) / a);
                if (!std::isfinite(b)) continue;
                const bool lower = a < 0;
                if (v.type != VarType::Continuous)
                    b = bound_round(b, lower, opt.integrality_tolerance);
                double &target = lower ? v.lower : v.upper;
                if (!(lower ? b > target : b < target)) continue;
                if (std::isfinite(target) && std::abs(b - target) <= 1e-10 * (1 + std::abs(b)))
                    continue;
                target = b;
                ++changes;
                if (v.lower > v.upper + opt.feasibility_tolerance) return false;
                for (int next : incidence[j]) {
                    budget.checkpoint(1);
                    if (next != i) enqueue(next);
                }
            }
        }
        return true;
    }
};

// 上側活動量Uを使う整数係数強化。delta=U-rhsより大きい整数変数の
// 係数をdeltaへ縮める。最大側の端点では同じ式、それ以外の整数値では
// その行が常に成立するので、整数実行可能集合を保つ。連続変数には適用しない。
inline bool opposite_rows(const Row &a, const Row &b, Budget &budget) {
    if (a.rhs != -b.rhs || a.terms.size() != b.terms.size()) return false;
    for (std::size_t k = 0; k < a.terms.size(); ++k) {
        budget.checkpoint(1);
        if (a.terms[k].first != b.terms[k].first || a.terms[k].second != -b.terms[k].second)
            return false;
    }
    return true;
}

inline bool strengthen(std::vector<Row> &rows, const std::vector<Variable> &vars, Budget &budget) {
    bool changed = false;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        budget.checkpoint(1);
        // 等式を片側ずつ強化すると、後段での変数消去の機会を失う。
        if (i + 1 < rows.size() && opposite_rows(rows[i], rows[i + 1], budget)) {
            ++i;
            continue;
        }
        auto &row = rows[i];
        Accurate maximum;
        bool finite = true;
        for (auto [j, a] : row.terms) {
            budget.checkpoint(1);
            const double b = a > 0 ? vars[j].upper : vars[j].lower;
            if (!std::isfinite(b)) {
                finite = false;
                break;
            }
            maximum = maximum + Accurate(a) * Accurate(b);
        }
        if (!finite) continue;
        const Accurate excess = maximum - Accurate(row.rhs);
        double delta = excess.value();
        if ((Accurate(delta) - excess).value() < 0) delta = std::nextafter(delta, inf);
        if (!(delta > 0) || !std::isfinite(delta)) continue;
        Accurate rhs(row.rhs);
        bool row_changed = false;
        for (auto &term : row.terms) {
            budget.checkpoint(1);
            const auto &v = vars[term.first];
            const double a = term.second;
            if (v.type == VarType::Continuous || v.lower == v.upper ||
                std::abs(a) <= delta + 1e-10 * std::max(1.0, std::abs(a)))
                continue;
            const double next = std::copysign(delta, a), bound = a > 0 ? v.upper : v.lower;
            rhs = rhs - (Accurate(a) - Accurate(next)) * Accurate(bound);
            term.second = next;
            row_changed = true;
            ++budget.stats.strengthened_coefficients;
        }
        if (row_changed) {
            row.rhs = rhs.value();
            if ((Accurate(row.rhs) - rhs).value() < 0) row.rhs = std::nextafter(row.rhs, inf);
            changed = true;
        }
    }
    return changed;
}

inline bool
remove_fixed_terms(std::vector<Row> &rows, const std::vector<Variable> &vars, Budget &budget) {
    bool changed = false;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        budget.checkpoint(1);
        // 右辺を外側へ丸めて等式を微小な区間に変えない。固定項は代入時に扱う。
        if (i + 1 < rows.size() && opposite_rows(rows[i], rows[i + 1], budget)) {
            ++i;
            continue;
        }
        auto &row = rows[i];
        Accurate rhs(row.rhs);
        std::size_t count = 0;
        for (auto [j, a] : row.terms) {
            budget.checkpoint(1);
            if (vars[j].lower == vars[j].upper) {
                rhs = rhs - Accurate(a) * Accurate(vars[j].lower);
                changed = true;
            } else if (a)
                row.terms[count++] = {j, a};
        }
        row.terms.resize(count);
        row.rhs = rhs.value();
        if ((Accurate(row.rhs) - rhs).value() < 0) row.rhs = std::nextafter(row.rhs, inf);
        long long divisor = 0;
        bool integral = true;
        for (auto [j, a] : row.terms) {
            budget.checkpoint(1);
            if (vars[j].type == VarType::Continuous || std::abs(a) > 1e9 || a != std::round(a)) {
                integral = false;
                break;
            }
            divisor = std::gcd(divisor, (long long)std::llround(std::abs(a)));
        }
        if (integral && divisor) {
            const double b = std::floor(row.rhs / divisor + budget.opt.integrality_tolerance);
            changed |= b * divisor != row.rhs || divisor != 1;
            row.rhs = b;
            for (auto &term : row.terms) {
                budget.checkpoint(1);
                term.second /= divisor;
            }
        }
    }
    return changed;
}

// 0/1変数の両側の伝播結果の包絡だけを採用。片側の境界を大域境界にしない。
inline bool probe(Propagation &propagation, std::vector<Row> &rows, std::vector<Variable> &vars,
                  Budget &budget) {
    Options local = budget.opt;
    local.time_limit =
        std::min(local.time_limit, budget.elapsed() + std::min(.02, local.time_limit * .05));
    Budget trial_budget{local, budget.stats, budget.start};
    try {
        std::vector<int> candidates;
        for (int j = 0; j < int(vars.size()); ++j) {
            trial_budget.checkpoint(1);
            if (vars[j].type != VarType::Continuous && vars[j].lower == 0 && vars[j].upper == 1 &&
                !propagation.incidence[j].empty())
                candidates.push_back(j);
        }
        std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
            trial_budget.checkpoint(1);
            const auto x = propagation.incidence[a].size(), y = propagation.incidence[b].size();
            return x != y ? x > y : a < b;
        });
        std::uint64_t work = std::min<std::uint64_t>(2000000, 32 * propagation.size);
        int failed = 0, tried = 0;
        for (int j : candidates) {
            if (!work || tried >= 128 || failed >= 16 || trial_budget.stop()) break;
            if (vars[j].lower == vars[j].upper) continue;
            ++tried;
            auto low = trial_budget.copy(vars), high = trial_budget.copy(vars);
            low[j].upper = 0;
            high[j].lower = 1;
            ++budget.stats.probing_trials;
            bool ok0 = propagation.run(low, trial_budget, work, j);
            ++budget.stats.probing_trials;
            bool ok1 = propagation.run(high, trial_budget, work, j);
            if (!ok0 && !ok1) return false;
            bool changed = false;
            // 片側が不可能なら、もう片側の推論をすべて使える。
            for (int k = 0; k < int(vars.size()); ++k) {
                trial_budget.checkpoint(1);
                double l = !ok0   ? high[k].lower
                           : !ok1 ? low[k].lower
                                  : std::min(low[k].lower, high[k].lower);
                double u = !ok0   ? high[k].upper
                           : !ok1 ? low[k].upper
                                  : std::max(low[k].upper, high[k].upper);
                if (l > vars[k].lower || u < vars[k].upper) {
                    vars[k].lower = std::max(vars[k].lower, l);
                    vars[k].upper = std::min(vars[k].upper, u);
                    ++budget.stats.probing_bound_changes;
                    changed = true;
                }
            }
            if (ok0 && ok1)
                for (int i : propagation.incidence[j]) {
                    trial_budget.checkpoint(1);
                    if ((i > 0 && opposite_rows(rows[i - 1], rows[i], trial_budget)) ||
                        (i + 1 < int(rows.size()) &&
                         opposite_rows(rows[i], rows[i + 1], trial_budget)))
                        continue;
                    auto &row = rows[i];
                    auto term = std::lower_bound(row.terms.begin(), row.terms.end(), j,
                                                 [](const auto &t, int v) { return t.first < v; });
                    if (term == row.terms.end() || term->first != j || term->second == 0) continue;
                    const double a = term->second;
                    const auto &implied = a > 0 ? low : high;
                    Accurate activity;
                    bool finite = true;
                    for (auto [k, c] : row.terms) {
                        trial_budget.checkpoint(1);
                        if (k == j || c == 0) continue;
                        double b = c > 0 ? implied[k].upper : implied[k].lower;
                        if (!std::isfinite(b)) {
                            finite = false;
                            break;
                        }
                        activity = activity + Accurate(c) * Accurate(b);
                    }
                    if (!finite) continue;
                    // x=0側（負係数ならx=1側）の含意境界で、行の余裕を証明する。
                    // 反対側の式は不変に保ち、その余裕だけ係数を縮める。
                    const Accurate slack = Accurate(row.rhs) - activity - Accurate(a < 0 ? a : 0);
                    if (!(slack.value() > 1e-8 * std::max(1.0, std::abs(a)))) continue;
                    Accurate reduced = a > 0 ? Accurate(a) - slack : Accurate(a) + slack;
                    double next = reduced.value();
                    if (a > 0) {
                        if ((Accurate(next) - reduced).value() < 0)
                            next = std::nextafter(next, inf);
                        next = std::max(0.0, next);
                    } else {
                        if ((Accurate(next) - reduced).value() > 0)
                            next = std::nextafter(next, -inf);
                        next = std::min(0.0, next);
                    }
                    if (next == a) continue;
                    if (a > 0) {
                        const Accurate rhs = Accurate(row.rhs) - (Accurate(a) - Accurate(next));
                        row.rhs = rhs.value();
                        if ((Accurate(row.rhs) - rhs).value() < 0)
                            row.rhs = std::nextafter(row.rhs, inf);
                    }
                    term->second = next;
                    ++budget.stats.strengthened_coefficients;
                    changed = true;
                }
            if (changed) {
                failed = 0;
                if (!propagation.run(vars, trial_budget, work)) return false;
            } else
                ++failed;
        }
    } catch (const TimeLimitReached &) { /* 局所予算が尽きても、確定した境界は使える。 */
    }
    for (auto &row : rows) {
        budget.checkpoint(1);
        row.terms.erase(std::remove_if(row.terms.begin(), row.terms.end(),
                                       [&](auto t) {
                                           budget.checkpoint(1);
                                           return t.second == 0;
                                       }),
                        row.terms.end());
    }
    budget.checkpoint();
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
        bool presolve = false, advanced = false, probing = false, integer_substitution = false,
             substitution = false;
        bool scaling = false, dualize = false, relax = false, bounded = false;
        int sparse_mode = -1;
        double feasibility = 0, integrality = 0, pivot = 0;
        bool matches(const Options &o, bool r) const {
            return presolve == o.presolve && advanced == o.advanced_presolve &&
                   probing == o.probing && integer_substitution == o.integer_substitution &&
                   substitution == o.equality_substitution && scaling == o.scaling &&
                   dualize == o.dualize && bounded == o.bounded_variables &&
                   sparse_mode == o.sparse_tableau && relax == r &&
                   feasibility == o.feasibility_tolerance &&
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

    static void combine(Terms &terms, detail::Budget &budget) {
        std::sort(terms.begin(), terms.end(), [&](const auto &a, const auto &b) {
            budget.checkpoint(1);
            return a.first < b.first;
        });
        std::size_t out = 0;
        for (std::size_t begin = 0; begin < terms.size();) {
            const int j = terms[begin].first;
            long double value = 0;
            do {
                budget.checkpoint(1);
                value += terms[begin++].second;
            } while (begin < terms.size() && terms[begin].first == j);
            if (value != 0) terms[out++] = {j, double(value)};
        }
        terms.resize(out);
    }

    struct Substitution {
        int variable;
        double offset;
        Terms terms;
    };
    // q*x + a*y = b から x = (b-a*y)/q を代入する。x の上下限も行として残す。
    // 目的係数が0の連続変数、または整数等式の係数±1の整数変数を選ぶ。
    // 後者は整数格子を保つ。
    std::vector<Substitution> substitute(detail::Prepared &p, detail::Budget &budget) const {
        std::vector<Substitution> eliminated;
        // 代入の連鎖で復元式だけが二次規模に膨らむことも避ける。
        // 未展開の参照数を保守的に伝播し、式の総項数の増加を4n以内に抑える。
        auto references = budget.filled<std::uint64_t>(p.vars.size(), 1);
        std::uint64_t growth = 0;
        const std::uint64_t growth_limit = 4 * std::uint64_t(p.vars.size());
        for (std::size_t i = 0; i + 1 < p.rows.size(); ++i) {
            budget.checkpoint(1);
            const auto &first = p.rows[i];
            const auto &second = p.rows[i + 1];
            if (first.terms.size() < 2 ||
                first.terms.size() > (budget.opt.advanced_presolve ? 16U : 4U) ||
                first.terms.size() != second.terms.size() || first.rhs != -second.rhs)
                continue;
            bool equal = true;
            bool lattice = budget.opt.advanced_presolve && budget.opt.integer_substitution &&
                           first.rhs == std::round(first.rhs) && std::abs(first.rhs) <= 1e12;
            for (auto [j, a] : first.terms)
                lattice = lattice && p.vars[j].type != VarType::Continuous && a == std::round(a) &&
                          std::abs(a) <= 1e9;
            double largest = 0, pivot = 0;
            int variable = -1;
            for (std::size_t k = 0; k < first.terms.size(); ++k) {
                const auto [j, a] = first.terms[k];
                equal = equal && j == second.terms[k].first && a == -second.terms[k].second;
                largest = std::max(largest, std::abs(a));
                if ((vars_[j].type == VarType::Continuous ||
                     p.vars[j].type == VarType::Continuous || (lattice && std::abs(a) == 1)) &&
                    objective_[j] == 0 && p.vars[j].lower != p.vars[j].upper &&
                    (variable < 0 || std::abs(a) > std::abs(pivot))) {
                    variable = j;
                    pivot = a;
                }
            }
            if (!equal || variable < 0 || std::abs(pivot) < largest * 1e-6) continue;
            Substitution rule{variable, first.rhs / pivot, {}};
            for (auto [j, a] : first.terms)
                if (j != variable) rule.terms.emplace_back(j, -a / pivot);
            if (!std::isfinite(rule.offset)) continue;
            const auto added = references[variable] * (rule.terms.size() - 1);
            if (added > growth_limit - growth) continue;
            std::vector<std::pair<std::size_t, detail::Row>> replacements;
            bool valid = true;
            // 生成する係数をすべて確認してから適用。途中の時間切れでは結果を使わない。
            for (std::size_t k = 0; k < p.rows.size(); ++k) {
                budget.checkpoint(1);
                if (k == i || k == i + 1) continue;
                const auto &row = p.rows[k];
                const auto pos = std::lower_bound(row.terms.begin(), row.terms.end(), variable,
                                                  [](const auto &t, int j) { return t.first < j; });
                if (pos == row.terms.end() || pos->first != variable) continue;
                if (replacements.size() >= 128) {
                    valid = false;
                    break;
                } // fill-in の上限。
                const double factor = pos->second;
                detail::Row next{budget.copy(row.terms),
                                 double((long double)row.rhs - (long double)factor * rule.offset)};
                next.terms.erase(next.terms.begin() + (pos - row.terms.begin()));
                for (auto [j, a] : rule.terms) next.terms.emplace_back(j, factor * a);
                combine(next.terms, budget);
                if (!std::isfinite(next.rhs)) valid = false;
                for (auto t : next.terms)
                    if (!std::isfinite(t.second)) valid = false;
                if (!valid) break;
                replacements.emplace_back(k, std::move(next));
            }
            const auto &v = p.vars[variable];
            detail::Row lower, upper;
            lower.rhs = upper.rhs = 0;
            if (std::isfinite(v.lower)) {
                lower.rhs = rule.offset - v.lower;
                for (auto [j, a] : rule.terms) lower.terms.emplace_back(j, -a);
            }
            if (std::isfinite(v.upper)) {
                upper.rhs = v.upper - rule.offset;
                upper.terms = rule.terms;
            }
            if (!valid || !std::isfinite(lower.rhs) || !std::isfinite(upper.rhs)) continue;
            // 新しく対象を広げる等式では、変数数だけでなく非ゼロ数も減る場合に限る。
            // 特に疎なLPを密にするfill-inを防ぐ。従来の短いゼロ目的の消去は維持。
            if (budget.opt.advanced_presolve &&
                (first.terms.size() > 4 || vars_[variable].type != VarType::Continuous)) {
                std::size_t before = first.terms.size() + second.terms.size(),
                            after = lower.terms.size() + upper.terms.size();
                for (const auto &replacement : replacements) {
                    before += p.rows[replacement.first].terms.size();
                    after += replacement.second.terms.size();
                }
                if (after > before) continue;
            }
            for (auto &replacement : replacements)
                p.rows[replacement.first] = std::move(replacement.second);
            p.rows[i] = std::move(lower);
            p.rows[i + 1] = std::move(upper);
            growth += added;
            for (auto term : rule.terms) references[term.first] += references[variable];
            eliminated.push_back(std::move(rule));
            ++budget.stats.substituted_variables;
            ++i;
        }
        return eliminated;
    }

    detail::Prepared prepare(detail::Budget &budget, bool relax) const {
        const auto &opt = budget.opt;
        const double started = budget.elapsed();
        struct Timing {
            detail::Budget &budget;
            double start;
            ~Timing() { budget.stats.presolve_seconds += budget.elapsed() - start; }
        } timing{budget, started};
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
        bool has_integer = false, has_continuous = false;
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
            has_integer |= v.type != VarType::Continuous;
            has_continuous |= v.type == VarType::Continuous;
        }
        if (opt.presolve && opt.advanced_presolve && has_integer && has_continuous) {
            // 整数の右辺・整数変数と係数±1で結ばれた連続変数は、既に整数。
            // 等式の両向きが厳密に一致する場合だけ使い、LP緩和には適用しない。
            for (int pass = 0; pass < 8; ++pass) {
                bool changed = false;
                for (std::size_t i = 0; i + 1 < p.rows.size(); ++i) {
                    budget.checkpoint(1);
                    const auto &a = p.rows[i];
                    const auto &b = p.rows[i + 1];
                    if (a.terms.empty() || a.terms.size() != b.terms.size() || a.rhs != -b.rhs ||
                        a.rhs != std::round(a.rhs) || std::abs(a.rhs) > 1e12)
                        continue;
                    bool valid = true;
                    int variable = -1;
                    for (std::size_t k = 0; k < a.terms.size(); ++k) {
                        budget.checkpoint(1);
                        auto [j, c] = a.terms[k];
                        if (j != b.terms[k].first || c != -b.terms[k].second ||
                            c != std::round(c) || std::abs(c) > 1e9) {
                            valid = false;
                            break;
                        }
                        if (p.vars[j].type == VarType::Continuous) {
                            if (variable >= 0 || std::abs(c) != 1) {
                                valid = false;
                                break;
                            }
                            variable = j;
                        }
                    }
                    if (!valid || variable < 0) continue;
                    auto &v = p.vars[variable];
                    v.type = VarType::Integer;
                    v.lower = detail::bound_round(v.lower, true, opt.integrality_tolerance);
                    v.upper = detail::bound_round(v.upper, false, opt.integrality_tolerance);
                    ++budget.stats.implied_integers;
                    changed = true;
                    if (v.lower > v.upper) {
                        p.infeasible = true;
                        return p;
                    }
                }
                if (!changed) break;
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
            if (opt.advanced_presolve) {
                detail::Propagation propagation(p.rows, p.vars.size(), budget);
                auto work = 16 * propagation.size;
                if (!propagation.run(p.vars, budget, work)) {
                    p.infeasible = true;
                    return p;
                }
                for (int pass = 0; has_integer && pass < 3; ++pass) {
                    if (!detail::strengthen(p.rows, p.vars, budget)) break;
                    work = 8 * propagation.size;
                    if (!propagation.run(p.vars, budget, work)) {
                        p.infeasible = true;
                        return p;
                    }
                }
                budget.stats.presolve_bound_changes += propagation.changes;
                if (opt.probing && has_integer &&
                    !detail::probe(propagation, p.rows, p.vars, budget)) {
                    p.infeasible = true;
                    return p;
                }
                const auto before_final = propagation.changes;
                if (has_integer && detail::strengthen(p.rows, p.vars, budget)) {
                    work = 8 * propagation.size;
                    if (!propagation.run(p.vars, budget, work)) {
                        p.infeasible = true;
                        return p;
                    }
                }
                budget.stats.presolve_bound_changes += propagation.changes - before_final;
            } else if (!detail::propagate(p.rows, p.vars, budget)) {
                p.infeasible = true;
                return p;
            }
            bool fixed = false;
            for (const auto &v : p.vars) {
                budget.checkpoint(1);
                fixed |= v.lower == v.upper;
            }
            if (opt.advanced_presolve && (fixed || budget.stats.strengthened_coefficients))
                for (int pass = 0; pass < 3; ++pass) {
                    if (!detail::remove_fixed_terms(p.rows, p.vars, budget)) break;
                    detail::Propagation propagation(p.rows, p.vars.size(), budget);
                    auto work = 8 * propagation.size;
                    if (!propagation.run(p.vars, budget, work)) {
                        p.infeasible = true;
                        return p;
                    }
                    budget.stats.presolve_bound_changes += propagation.changes;
                }
        }
        std::vector<Substitution> substitutions;
        if (opt.presolve && opt.equality_substitution) substitutions = substitute(p, budget);
        p.substitutions = int(substitutions.size());
        if (!substitutions.empty() && !detail::propagate(p.rows, p.vars, budget)) {
            p.infeasible = true;
            return p;
        }
        auto removed = budget.filled<unsigned char>(vars_.size());
        for (const auto &rule : substitutions) removed[rule.variable] = 1;
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
            if (removed[j]) continue;
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
        // 逆順に合成すると、復元式と追加制約の変換は残った列だけで表せる。
        for (auto it = substitutions.rbegin(); it != substitutions.rend(); ++it) {
            auto &map = p.map[it->variable];
            map.offset = it->offset;
            Terms expression;
            for (auto [j, a] : it->terms) {
                budget.checkpoint(1);
                map.offset += a * p.map[j].offset;
                for (auto [k, b] : p.map[j].terms) {
                    budget.checkpoint(1);
                    expression.emplace_back(k, a * b);
                }
            }
            combine(expression, budget);
            map.terms.assign(std::move(expression));
        }
        std::unordered_map<std::size_t, std::vector<int>> identical;
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
            if (opt.presolve && opt.advanced_presolve) {
                std::size_t hash = a.size();
                for (auto [j, v] : a) {
                    budget.checkpoint(1);
                    hash ^= std::hash<double>{}(v) + std::size_t(j) + 0x9e3779b9U + (hash << 6) +
                            (hash >> 2);
                }
                auto &bucket = identical[hash];
                for (int i : bucket) {
                    bool equal = p.a[i].size() == a.size();
                    for (std::size_t k = 0; equal && k < a.size(); ++k) {
                        budget.checkpoint(1);
                        equal = p.a[i][k] == a[k];
                    }
                    if (equal) {
                        p.b[i] = std::min(p.b[i], b);
                        ++budget.stats.duplicate_rows;
                        return;
                    }
                }
                bucket.push_back(int(p.b.size()));
            }
            budget.stats.presolved_nonzeros += a.size();
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
        if (opt.bounded_variables) p.upper = budget.filled<double>(p.c.size(), inf);
        bool finite_upper = false;
        for (int j = 0; j < int(p.vars.size()); ++j) {
            budget.checkpoint(1);
            if (removed[j]) continue;
            const auto &v = p.vars[j];
            if (v.lower == v.upper || !std::isfinite(v.lower) || !std::isfinite(v.upper)) continue;
            if (implied_upper[p.map[j].terms[0].first] <= v.upper - v.lower) continue;
            if (opt.bounded_variables) {
                p.upper[p.map[j].terms[0].first] = v.upper - v.lower;
                finite_upper = true;
                continue;
            }
            Terms a;
            a.emplace_back(p.map[j].terms[0].first, 1);
            append(std::move(a), v.upper - v.lower);
        }
        if (!finite_upper) p.upper.clear();
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
            // 消去した変数と残った変数が同じ列を参照する場合は係数を加算する。
            if (p.substitutions) combine(a, budget);
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
        next->substitution = o.equality_substitution;
        next->advanced = o.advanced_presolve;
        next->probing = o.probing;
        next->integer_substitution = o.integer_substitution;
        next->scaling = o.scaling;
        next->dualize = o.dualize;
        next->bounded = o.bounded_variables;
        next->sparse_mode = o.sparse_tableau;
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
                const double y = tab.coordinate(k);
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

    // 実行可能な点から目的を改善する方向へ動かす。制約種別は仮定しない。
    // 列の隣接行だけから許される移動量を求め、その行の活動量だけを更新。
    // 中断時も x は完成した座標移動の結果なので、元モデルで検査して回収できる。
    void improve_feasible(std::vector<double> &x, const std::vector<detail::Variable> &vars,
                          detail::Budget &budget) const {
        std::vector<Terms> columns(vars.size());
        auto activity = budget.filled<long double>(rows_.size());
        for (int i = 0; i < int(rows_.size()); ++i) {
            budget.chunks(rows_[i].terms.size(), [&](std::size_t begin, std::size_t end) {
                for (auto k = begin; k < end; ++k) {
                    const auto [j, a] = rows_[i].terms[k];
                    activity[i] += (long double)a * x[j];
                    columns[j].emplace_back(i, a);
                }
            });
        }
        std::vector<int> order(vars.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            budget.checkpoint(1);
            const double va = std::abs(objective_[a]), vb = std::abs(objective_[b]);
            return va == vb ? a < b : va > vb;
        });
        const double sign = sense_ == Sense::Maximize ? 1 : -1;
        for (int pass = 0; pass < 2; ++pass) {
            bool changed = false;
            for (int j : order) {
                budget.checkpoint(1);
                const double gain = sign * objective_[j];
                if (gain == 0) continue;
                const double direction = gain > 0 ? 1 : -1;
                long double distance = gain > 0 ? vars[j].upper - x[j] : x[j] - vars[j].lower;
                budget.chunks(columns[j].size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) {
                        const auto [i, a] = columns[j][k];
                        const double effect = direction * a;
                        if (effect > 0)
                            distance = std::min(
                                distance, ((long double)rows_[i].rhs - activity[i]) / effect);
                    }
                });
                if (vars[j].type != VarType::Continuous) distance = std::floor(distance);
                if (!(distance > 0) || !std::isfinite(double(distance))) continue;
                const double next = x[j] + direction * double(distance);
                const double move = next - x[j];
                if (move == 0) continue;
                x[j] = next;
                budget.chunks(columns[j].size(), [&](std::size_t begin, std::size_t end) {
                    for (auto k = begin; k < end; ++k) {
                        const auto [i, a] = columns[j][k];
                        activity[i] += (long double)a * move;
                    }
                });
                changed = true;
            }
            if (!changed) break;
        }
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
            opt.strong_branching_depth < 0 || opt.sparse_tableau < -1 || opt.sparse_tableau > 1 ||
            opt.root_cut_rounds < -1 || opt.cuts_per_round < 0)
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
            // 変換・カット後の数値不整合は元のモデルで再確認。暫定解と総予算を継承。
            if (status == Status::NumericalError &&
                ((opt.presolve && opt.advanced_presolve) || ret.stats.cuts ||
                 (opt.equality_substitution && ret.stats.substituted_variables) ||
                 (opt.basis_repair && opt.scaling))) {
                Options retry = opt;
                retry.root_cut_rounds = 0;
                retry.advanced_presolve = false;
                // スケーリング後の絶対許容幅が元の行で増幅する場合は、元の単位で再確認。
                if (opt.basis_repair && !opt.advanced_presolve) retry.scaling = false;
                if (ret.stats.substituted_variables) retry.equality_substitution = false;
                retry.time_limit = std::max(0.0, opt.time_limit - budget.elapsed());
                retry.iteration_limit =
                    opt.iteration_limit - std::min(opt.iteration_limit, ret.stats.pivots);
                retry.node_limit = opt.node_limit - std::min(opt.node_limit, ret.stats.nodes);
                if (ret.has_solution) retry.initial_solution = ret.x;
                Result fresh = solve(retry, relax);
                fresh.stats.nodes += ret.stats.nodes;
                fresh.stats.pivots += ret.stats.pivots;
                fresh.stats.strong_branching_pivots += ret.stats.strong_branching_pivots;
                fresh.stats.bound_flips += ret.stats.bound_flips;
                fresh.stats.basis_refactorizations += ret.stats.basis_refactorizations;
                fresh.stats.objective_perturbations += ret.stats.objective_perturbations;
                fresh.stats.pump_pivots += ret.stats.pump_pivots;
                fresh.stats.pump_lp_solves += ret.stats.pump_lp_solves;
                fresh.stats.sparse_pivots += ret.stats.sparse_pivots;
                fresh.stats.sparse_tableaux += ret.stats.sparse_tableaux;
                fresh.stats.dense_conversions += ret.stats.dense_conversions;
                fresh.stats.tableau_evictions += ret.stats.tableau_evictions;
                fresh.stats.tableau_rebuilds += ret.stats.tableau_rebuilds;
                fresh.stats.peak_pending_tableau_bytes = std::max(
                    fresh.stats.peak_pending_tableau_bytes, ret.stats.peak_pending_tableau_bytes);
                fresh.stats.root_pivots += ret.stats.root_pivots;
                fresh.stats.reused_root = fresh.stats.reused_root || ret.stats.reused_root;
                fresh.stats.heuristic_solutions += ret.stats.heuristic_solutions;
                fresh.stats.diving_pivots += ret.stats.diving_pivots;
                fresh.stats.diving_lp_solves += ret.stats.diving_lp_solves;
                fresh.stats.substituted_variables += ret.stats.substituted_variables;
                fresh.stats.presolve_bound_changes += ret.stats.presolve_bound_changes;
                fresh.stats.strengthened_coefficients += ret.stats.strengthened_coefficients;
                fresh.stats.probing_trials += ret.stats.probing_trials;
                fresh.stats.probing_bound_changes += ret.stats.probing_bound_changes;
                fresh.stats.duplicate_rows += ret.stats.duplicate_rows;
                fresh.stats.implied_integers += ret.stats.implied_integers;
                fresh.stats.presolve_seconds += ret.stats.presolve_seconds;
                fresh.stats.cuts += ret.stats.cuts;
                fresh.stats.reduced_cost_fixings += ret.stats.reduced_cost_fixings;
                fresh.stats.numerical_restarts += ret.stats.numerical_restarts + 1;
                fresh.stats.seconds = budget.elapsed();
                return fresh;
            }
            ret.status = status;
            // 不整合のある辞書や実行不能判定から得た境界は証明に使えない。
            // 元の制約で検査済みの暫定解は保持し、境界だけを未確定に戻す。
            if (status == Status::NumericalError) global_bound = inf;
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
                    ret.stats.substituted_variables = p.substitutions;
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
                    // 保存済み根からの短い再最適化では、索引作成の固定費を避ける。
                    if (!reuse && opt.primal_sweeps && vars_.size() >= 128 && ret.has_solution) {
                        auto x = budget.copy(ret.x);
                        try {
                            improve_feasible(x, p.vars, budget);
                        } catch (const detail::TimeLimitReached &) {
                            accept(x, true, false);
                            throw;
                        }
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
                if (!reuse && opt.dualize && p.upper.empty() && integers.empty() && !p.b.empty() &&
                    p.c.size() > 2 * p.b.size() &&
                    std::all_of(p.c.begin(), p.c.end(), [](double c) { return c <= 0; })) {
                    auto a = budget.filled<Terms>(p.c.size());
                    auto b = budget.copy(p.c), c = budget.copy(p.b);
                    for (double &v : b) v = -v;
                    for (double &v : c) v = -v;
                    for (int i = 0; i < int(p.b.size()); ++i)
                        budget.chunks(p.a[i].size(), [&](std::size_t begin, std::size_t end) {
                            for (auto k = begin; k < end; ++k) {
                                auto [j, v] = p.a[i][k];
                                a[j].emplace_back(i, -v);
                            }
                        });
                    detail::Tableau dual(a, b, c, &budget);
                    ret.stats.root_storage_bytes = dual.stored_bytes();
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
                if (!reuse) root = detail::Tableau(p.a, p.b, p.c, &budget, p.upper);
                ret.stats.root_storage_bytes = root.stored_bytes();
                ++ret.stats.nodes;
                const auto before = ret.stats.pivots;
                // 長い再最適化は一度だけ作り直す。総時間・pivot・node の予算は継続。
                const int soft_limit =
                    int(std::min<long long>(std::numeric_limits<int>::max(),
                                            std::max(256LL, 4LL * (root.n + (long long)root.m))));
                Status root_status = reuse ? root.dual(budget, soft_limit) : root.solve(budget);
                ret.stats.root_pivots += ret.stats.pivots - before;
                for (int repair = 0; repair < 2 && root_status == Status::Optimal; ++repair) {
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
                    if (valid) break;
                    const auto repair_start = ret.stats.pivots;
                    root_status = repair == 0 ? root.recover_basis(budget) : Status::NumericalError;
                    ret.stats.root_pivots += ret.stats.pivots - repair_start;
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

            detail::Tableau alternative_root;
            std::vector<double> alternative_point;
            auto diversify_root = [&] {
                if (!opt.root_diversification || ret.has_solution || !opt.dual_devex ||
                    std::uint64_t(p.a.size() + 2) * (p.c.size() + 2) * sizeof(double) >
                        opt.tableau_cache_limit)
                    return;
                Options local = opt;
                local.dual_devex = false;
                local.time_limit = std::min(
                    opt.time_limit, budget.elapsed() + std::min(.02, opt.time_limit * .03));
                local.iteration_limit =
                    ret.stats.pivots +
                    std::min<std::uint64_t>(
                        4096, opt.iteration_limit - std::min(opt.iteration_limit, ret.stats.pivots));
                detail::Budget alternate_budget{local, ret.stats, budget.start};
                const auto before = ret.stats.pivots;
                try {
                    detail::Tableau candidate(p.a, p.b, p.c, &alternate_budget, p.upper);
                    ++ret.stats.diving_lp_solves;
                    if (candidate.solve(alternate_budget) == Status::Optimal) {
                        auto point = recover(p, candidate);
                        if (feasible(point, opt.feasibility_tolerance, true, &alternate_budget)) {
                            alternative_point = std::move(point);
                            alternative_root = std::move(candidate);
                        }
                    }
                } catch (const detail::TimeLimitReached &) {}
                ret.stats.diving_pivots += ret.stats.pivots - before;
            };

            // 丸め点へのL1距離を最小化するLPと、距離目的での小さな分枝探索。
            // ヒューリスティックの境界は証明に使わず、元モデルで検査した解だけ採用。
            auto run_pump = [&] {
                if (!opt.feasibility_pump || ret.has_solution) return;
                Options local = opt;
                local.time_limit = std::min(
                    opt.time_limit, budget.elapsed() + std::min(.25, opt.time_limit * .25));
                local.iteration_limit =
                    ret.stats.pivots +
                    std::min<std::uint64_t>(
                        20000,
                        opt.iteration_limit - std::min(opt.iteration_limit, ret.stats.pivots));
                detail::Budget pump_budget{local, ret.stats, budget.start};
                const auto before = ret.stats.pivots;
                auto account = [&] { ret.stats.pump_pivots += ret.stats.pivots - before; };
                try {
                    int solves = 0;
                    auto x = alternative_point.empty() ? root_x : alternative_point;
                    std::vector<double> previous;
                    std::uint32_t random = 123456789;
                    auto draw = [&] {
                        random ^= random << 13;
                        random ^= random >> 17;
                        random ^= random << 5;
                        return random / 4294967296.0;
                    };
                    for (int round = 0; round < 40 && solves < 256 && !pump_budget.stop(); ++round) {
                        std::vector<double> target;
                        target.reserve(integers.size());
                        for (int j : integers) target.push_back(std::round(x[j]));
                        bool same = target == previous;
                        if (same || (round && round % 7 == 0)) {
                            for (int k = 0; k < int(integers.size()); ++k) {
                                int j = integers[k];
                                double f = x[j] - std::floor(x[j]);
                                if (draw() < .1 + .5 * std::min(f, 1 - f)) {
                                    double direction = draw() < .5 ? -1 : 1;
                                    target[k] = std::max(std::ceil(p.vars[j].lower),
                                                         std::min(std::floor(p.vars[j].upper),
                                                                  target[k] + direction));
                                }
                            }
                        }
                        previous = target;
                        detail::ImmutableRows a;
                        a.rows = pump_budget.copy(p.a.rows);
                        auto b = pump_budget.copy(p.b);
                        auto c = pump_budget.filled<double>(p.c.size());
                        auto upper = pump_budget.copy(p.upper);
                        if (upper.empty()) upper.resize(c.size(), inf);
                        for (int k = 0; k < int(integers.size()); ++k) {
                            pump_budget.checkpoint(1);
                            const int j = integers[k];
                            if (p.map[j].terms.empty()) continue;
                            double direction = target[k] <= p.vars[j].lower   ? -1
                                               : target[k] >= p.vars[j].upper ? 1
                                                                              : 0;
                            if (direction) {
                                for (auto [v, coef] : p.map[j].terms) c[v] += direction * coef;
                            } else {
                                int t = int(c.size());
                                c.push_back(-1);
                                upper.push_back(inf);
                                for (double sign : {1.0, -1.0}) {
                                    Terms row;
                                    for (auto [v, coef] : p.map[j].terms)
                                        row.emplace_back(v, sign * coef);
                                    row.emplace_back(t, -1);
                                    a.push_terms(std::move(row));
                                    b.push_back(sign * (target[k] - p.map[j].offset));
                                }
                            }
                        }
                        detail::Tableau tab(a, b, c, &pump_budget, upper);
                        ++solves;
                        ++ret.stats.pump_lp_solves;
                        Status status = tab.solve(pump_budget);
                        if (status != Status::Optimal) break;
                        x = recover(p, tab.solution());
                        if (!feasible(x, opt.feasibility_tolerance, true, &pump_budget)) break;
                        auto rounded = x;
                        for (int j : integers) rounded[j] = std::round(rounded[j]);
                        if (accept(rounded, true)) break;
                        std::vector<detail::Tableau> stack;
                        stack.push_back(std::move(tab));
                        for (int node = 0;
                             node < 128 && solves < 256 && !stack.empty() && !pump_budget.stop();
                             ++node) {
                            auto near = std::move(stack.back());
                            stack.pop_back();
                            ++solves;
                            ++ret.stats.pump_lp_solves;
                            Status solved = near.dual(pump_budget, 256);
                            if (solved != Status::Optimal) continue;
                            auto y = recover(p, near.solution());
                            int chosen = -1;
                            double score = 0;
                            for (int j : integers) {
                                double f = y[j] - std::floor(y[j]);
                                if (std::min(f, 1 - f) > std::max(score, opt.integrality_tolerance)) {
                                    score = std::min(f, 1 - f);
                                    chosen = j;
                                }
                            }
                            if (chosen < 0) {
                                for (int j : integers) y[j] = std::round(y[j]);
                                if (accept(y, true)) {
                                    account();
                                    return;
                                }
                                continue;
                            }
                            if (stack.size() >= 32) continue;
                            bool up = y[chosen] - std::floor(y[chosen]) >= .5;
                            auto add = [&](detail::Tableau &t, bool direction) {
                                double sign = direction ? -1 : 1;
                                double b = direction ? std::ceil(y[chosen]) : std::floor(y[chosen]);
                                t.add_constraint(p.map[chosen].terms,
                                                 sign * (b - p.map[chosen].offset), &pump_budget,
                                                 sign);
                            };
                            auto other = near.clone(pump_budget);
                            add(other, !up);
                            add(near, up);
                            stack.push_back(std::move(other));
                            stack.push_back(std::move(near));
                        }
                    }
                } catch (const detail::TimeLimitReached &) {}
                account();
                budget.checkpoint();
            };

            // 証明用の根と独立した少数の dive。連続変数は各LPで再調整する。
            // 総時間・総pivot予算を共有し、局所打切りでは検査済みの解だけ残す。
            bool dive_pending = (opt.diving_heuristic || opt.feasibility_pump) && !reuse;
            detail::Tableau dive_root;
            if (dive_pending && opt.diving_heuristic) dive_root = root.clone(budget);
            auto run_dive = [&] {
                dive_pending = false;
                if (!opt.diving_heuristic) return;
                // 暫定解が根の境界に近ければ、追加探索より証明側へ予算を残す。
                const double urgency =
                    ret.has_solution ? std::min(1.0, 10 * std::max(0.0, global_bound - incumbent) /
                                                         std::max(1.0, std::abs(incumbent)))
                                     : 1.0;
                const double allowance =
                    std::min(0.08, opt.time_limit * 0.15) * std::max(.05, urgency);
                Options local;
                local.dual_devex = opt.dual_devex;
                local.basis_repair = opt.basis_repair;
                local.pivot_tolerance = opt.pivot_tolerance;
                local.time_limit = std::min(opt.time_limit, budget.elapsed() + allowance);
                local.iteration_limit =
                    ret.stats.pivots +
                    std::min<std::uint64_t>(
                        4096, opt.iteration_limit - std::min(opt.iteration_limit, ret.stats.pivots));
                detail::Budget dive_budget{local, ret.stats, budget.start};
                const auto before = ret.stats.pivots;
                try {
                    auto down_locks = dive_budget.filled<int>(vars_.size());
                    auto up_locks = dive_budget.filled<int>(vars_.size());
                    for (const auto &row : rows_)
                        for (auto [j, a] : row.terms) {
                            dive_budget.checkpoint(1);
                            ++(a > 0 ? up_locks[j] : down_locks[j]);
                        }
                    for (int attempt = 0; attempt < 3 && !dive_budget.stop(); ++attempt) {
                        auto tab = (attempt > 0 && !alternative_point.empty() ? alternative_root
                                                                              : dive_root)
                                       .clone(dive_budget);
                        for (int depth = 0; depth < 512 && !dive_budget.stop(); ++depth) {
                            auto x = recover(p, tab);
                            int chosen = -1;
                            bool up = false;
                            double best = inf;
                            for (int j : integers) {
                                dive_budget.checkpoint(1);
                                const double f = x[j] - std::floor(x[j]);
                                if (std::min(f, 1 - f) <= opt.integrality_tolerance) continue;
                                // 近い整数、制約を壊しにくい方向、目的改善方向の順に試す。
                                const bool direction = attempt == 0   ? f >= .5
                                                       : attempt == 1 ? up_locks[j] < down_locks[j]
                                                                      : sign * objective_[j] > 0;
                                const double score =
                                    (direction ? 1 - f : f) *
                                    (attempt == 1 ? 1 + (direction ? up_locks[j] : down_locks[j])
                                                  : 1);
                                if (score < best) {
                                    best = score;
                                    chosen = j;
                                    up = direction;
                                }
                            }
                            if (chosen < 0) {
                                for (int j : integers) x[j] = std::round(x[j]);
                                accept(x, true);
                                break;
                            }
                            auto previous = tab.clone(dive_budget);
                            auto branch = [&](bool direction) {
                                const double multiplier = direction ? -1 : 1;
                                const double bound =
                                    direction ? std::ceil(x[chosen]) : std::floor(x[chosen]);
                                tab.add_constraint(p.map[chosen].terms,
                                                   multiplier * (bound - p.map[chosen].offset),
                                                   &dive_budget, multiplier);
                                ++ret.stats.diving_lp_solves;
                                return tab.dual(dive_budget, 256);
                            };
                            Status status = branch(up);
                            if (status == Status::Infeasible) {
                                tab = std::move(previous);
                                status = branch(!up);
                            }
                            if (status != Status::Optimal) break;
                            if (ret.has_solution &&
                                tab.value() + p.offset <= incumbent + opt.absolute_gap)
                                break;
                        }
                    }
                } catch (const detail::TimeLimitReached &) {}
                ret.stats.diving_pivots += ret.stats.pivots - before;
                dive_root = detail::Tableau{};
                budget.checkpoint();
            };

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

            struct BranchPath {
                std::shared_ptr<const BranchPath> parent;
                int variable;
                bool up;
                double bound;
            };
            struct Node {
                detail::Tableau tab;
                std::shared_ptr<const BranchPath> path;
                double bound;
                int depth;
                bool solved;
                int branch_var = -1;
                bool up = false;
                double distance = 0, parent_bound = 0;
                std::uint64_t bytes = 0;
            };
            auto tableau_bytes = [](const detail::Tableau &tab) {
                // 共有行も全て数える保守的な値。辞書を復元可能なキャッシュとして扱う。
                return tab.stored_bytes();
            };
            auto replay_root = root.clone(budget);
            std::vector<Node> pending;
            pending.push_back({std::move(root), {}, global_bound, 0, true});
            pending.back().bytes = tableau_bytes(pending.back().tab);
            std::uint64_t pending_bytes = pending.back().bytes;
            std::size_t oldest_resident = 0;
            auto trim_cache = [&] {
                while (pending_bytes > opt.tableau_cache_limit && oldest_resident < pending.size()) {
                    budget.checkpoint();
                    auto &old = pending[oldest_resident++];
                    if (!old.bytes) continue;
                    pending_bytes -= old.bytes;
                    old.bytes = 0;
                    old.tab = detail::Tableau{};
                    if (old.solved)
                        old.branch_var = -1; // 既に更新した擬似コストを二重計上しない。
                    old.solved = false;
                    ++ret.stats.tableau_evictions;
                }
                ret.stats.peak_pending_tableau_bytes =
                    std::max(ret.stats.peak_pending_tableau_bytes, pending_bytes);
            };
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
            trim_cache();
            const auto search_start = ret.stats.pivots;
            while (!pending.empty()) {
                global_bound = incumbent;
                for (const auto &node : pending)
                    global_bound = std::max(global_bound, tighten(node.bound));
                if (prunable(global_bound)) return finish(Status::Optimal);
                if (budget.stop()) return finish(Status::Limit);
                // 容易なMIPは先に通常の探索で解く。ノード予算が尽きる前には試行する。
                if (dive_pending && (ret.stats.pivots - search_start >= 4096 ||
                                     ret.stats.nodes >= opt.node_limit)) {
                    diversify_root();
                    run_dive();
                    run_pump();
                    alternative_root = detail::Tableau{};
                    alternative_point.clear();
                    if (prunable(global_bound)) return finish(Status::Optimal);
                    if (budget.stop()) return finish(Status::Limit);
                }
                Node node = std::move(pending.back());
                pending.pop_back();
                pending_bytes -= node.bytes;
                oldest_resident = std::min(oldest_resident, pending.size());
                if (prunable(node.bound)) continue;
                if (node.depth) {
                    if (ret.stats.nodes >= opt.node_limit) return finish(Status::Limit);
                    ++ret.stats.nodes;
                }
                if (!node.solved) {
                    Status s;
                    if (!opt.warm_start || !node.bytes) {
                        node.tab = opt.warm_start
                                       ? replay_root.clone(budget)
                                       : detail::Tableau(p.a, p.b, p.c, &budget, p.upper);
                        std::vector<const BranchPath *> path;
                        for (auto b = node.path.get(); b; b = b->parent.get()) {
                            budget.checkpoint(1);
                            path.push_back(b);
                        }
                        for (auto it = path.rbegin(); it != path.rend(); ++it) {
                            const auto &b = **it;
                            add_branch(node.tab, b.variable, b.up, b.bound);
                            if (!opt.warm_start)
                                node.tab.write_row(node.tab.m - 1)[node.tab.n] = -1;
                        }
                        ++ret.stats.tableau_rebuilds;
                    }
                    s = opt.warm_start ? node.tab.dual(budget) : node.tab.solve(budget);
                    if (s == Status::Infeasible) continue;
                    if (s == Status::Unbounded) return finish(Status::NumericalError);
                    if (s != Status::Optimal) return finish(s);
                    node.bound = std::min(node.bound, node.tab.value() + p.offset);
                    if (node.branch_var >= 0)
                        update_pc(node.branch_var, node.up, node.distance,
                                  node.parent_bound - node.bound);
                }
                auto x = recover(p, node.tab);
                if (!feasible(x, opt.feasibility_tolerance, true, &budget)) {
                    const Status recovered = node.tab.recover_basis(budget);
                    if (recovered == Status::Infeasible) continue;
                    if (recovered == Status::Unbounded) return finish(Status::NumericalError);
                    if (recovered != Status::Optimal) return finish(recovered);
                    x = recover(p, node.tab);
                    if (!feasible(x, opt.feasibility_tolerance, true, &budget))
                        return finish(Status::NumericalError);
                    node.bound =
                        std::min(node.depth ? node.parent_bound : inf, node.tab.value() + p.offset);
                }
                if (opt.rounding_heuristic) round_candidate(x);
                if (prunable(node.bound)) continue;
                if (opt.reduced_cost_fixing && ret.has_solution && node.depth) {
                    step.resize(node.tab.next_id, 1);
                    // 再構築した辞書には祖先の列固定を戻さない。保存済みの node.bound が
                    // より強くても、reduced cost は現在の辞書の目的値と組にして使う。
                    ret.stats.reduced_cost_fixings += node.tab.fix_nonbasic(
                        step, node.tab.value() + p.offset - improvement_target(), budget,
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
                    auto checked_probe = [&](detail::Tableau &tab) {
                        const auto start = ret.stats.pivots;
                        Status result = tab.dual(budget, opt.strong_branching_iterations);
                        if (result == Status::Optimal &&
                            !feasible(recover(p, tab), opt.feasibility_tolerance, true, &budget)) {
                            const int left = std::max(
                                0, opt.strong_branching_iterations - int(ret.stats.pivots - start));
                            result = tab.recover_basis(budget, left);
                        }
                        return result == Status::Unbounded ? Status::NumericalError : result;
                    };
                    auto ds = checked_probe(down);
                    auto us = checked_probe(up);
                    ret.stats.strong_branching_pivots += ret.stats.pivots - before;
                    if (budget.stop()) return finish(Status::Limit);
                    if ((ds == Status::Optimal &&
                         !feasible(recover(p, down), opt.feasibility_tolerance, true, &budget)) ||
                        (us == Status::Optimal &&
                         !feasible(recover(p, up), opt.feasibility_tolerance, true, &budget)))
                        return finish(Status::NumericalError);
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
                    double b = node.bound;
                    if (tab.dual_feasible(10 * opt.pivot_tolerance))
                        b = std::min(b, tab.value() + p.offset);
                    if (prunable(b)) return;
                    auto path = std::make_shared<BranchPath>(BranchPath{
                        node.path, branch, up, up ? std::ceil(x[branch]) : std::floor(x[branch])});
                    const double distance =
                        up ? std::ceil(x[branch]) - x[branch] : x[branch] - std::floor(x[branch]);
                    const auto bytes = tableau_bytes(tab);
                    pending_bytes += bytes;
                    pending.push_back({std::move(tab), std::move(path), b, node.depth + 1,
                                       s == Status::Optimal, branch, up, distance, node.bound,
                                       bytes});
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
                trim_cache();
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
            tab_.source_cost = std::make_shared<const std::vector<double>>(c_);
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
            !std::isfinite(opt.pivot_tolerance) || opt.sparse_tableau < -1 || opt.sparse_tableau > 1)
            throw std::invalid_argument("invalid LP options");
        Result r;
        r.best_bound = inf;
        r.stats.original_variables = r.stats.presolved_variables = int(c_.size());
        r.stats.original_constraints = r.stats.presolved_constraints = int(b_.size());
        detail::Budget budget{opt, r.stats};
        if (state_ != State::Cold) r.stats.root_storage_bytes = tab_.stored_bytes();
        try {
            if (budget.stop() || opt.node_limit == 0)
                r.status = Status::Limit;
            else {
                ++r.stats.nodes;
                if (state_ == State::Cold) {
                    tab_ = detail::Tableau(a_, b_, c_, &budget);
                    r.stats.root_storage_bytes = tab_.stored_bytes();
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
