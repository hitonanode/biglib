# LP・MIP ソルバー

**[mip.hpp](mip.hpp) だけで動く C++17 のヘッダオンリー実装。** 標準ライブラリだけで利用できる。ファイル全体を提出コードの先頭に貼り付けられる。

連続・整数・0/1 変数、変数の上下限、自由変数、`<=` / `>=` / `==`、最大化・最小化に対応する。小〜中規模の LP・MIP を対象にした密行列ソルバーで、AHC での時間制限・暫定解の利用、LP・MIP の繰り返し再最適化もできる。

## 最初の例

```cpp
#include "mip.hpp"
#include <iostream>

int main() {
    cpopt::Model m;
    int x = m.add_integer();  // 0 <= x、整数
    int y = m.add_integer();
    m.add_le({{x, 2}, {y, 3}}, 17);
    m.add_le({{x, 5}, {y, 1}}, 19);
    m.maximize({{x, 1}, {y, 2}});
    auto r = m.solve();
    if (r.status == cpopt::Status::Optimal) {
        std::cout << std::llround(r.objective) << '\n';  // 11
        std::cout << r.x[x] << ' ' << r.x[y] << '\n';
    }
}
```

`cpopt::Terms` は `vector<pair<int, double>>`。同じ変数を複数回入れても係数を加算する。係数0の項は省略できる。変数は制約や目的関数に使う前に追加する。

```cpp
int a = m.add_variable();                         // 連続、[0, +inf)
int b = m.add_variable(-5, 8);                    // 連続、[-5, 8]
int c = m.add_variable(-cpopt::inf, cpopt::inf);    // 自由変数
int d = m.add_integer(-3, 10);                    // 整数、[-3, 10]
int e = m.add_binary();                          // 0 または 1
m.add_eq({{a, 1}, {b, -2}}, 3);
m.add_ge({{c, 1}, {d, 1}}, -4);
m.minimize({{a, 5}, {e, -7}}, 12);                // 定数項 +12
auto relaxation = m.solve_relaxation();           // 整数条件だけを外す
```

`add_variable` の第3引数には `VarType::Continuous` / `Integer` / `Binary` も指定できる。`add_constraint(terms, Relation::LessEqual/Equal/GreaterEqual, rhs)` は `add_le/ge/eq` の共通版。目的関数の再設定は、それ以前の目的関数を置き換える。

## 行列形式と繰り返し解く LP

`linprog(A, b, c)` は **最大化** `c*x`、制約 `A*x <= b, x >= 0` を解く。

```cpp
auto r = cpopt::linprog({{1, 1}, {2, 1}}, {4, 5}, {3, 2});
// 最適値9、x=(1, 3)
```

同じ変数のまま制約や目的関数を変える場合は `LPSolver` を使う。

```cpp
cpopt::LPSolver lp({{1, 1}, {2, 1}}, {4, 5}, {3, 2});
auto first = lp.solve();        // 最適値9
lp.add_le({1, 1}, 3);
auto second = lp.solve();       // 最適値8。親の基底から双対法で再最適化
lp.set_objective({0, 10});
auto third = lp.solve();        // 最適値30。主単体法で再最適化
```

制約追加だけ、または目的関数変更だけなら基底を再利用する。両方を変更してから `solve()` すると初期状態から解き直す。制約・変数の削除や既存行の RHS 変更には対応しない。`Model` の MIP 再求解については次節を参照。

`LPSolver` は非負の連続変数の最大化専用。前処理と行スケーリングは行わない。`Options` の時間・ノード・反復制限と LP 許容誤差が有効で、MIP 用オプションや初期解は使わない。

## 制約を追加しながら MIP を繰り返し解く

同じ `Model` の `solve()` は、既定で最適な根 LP の基底を保存する。次回は追加制約を保存済みの変数変換に通し、双対シンプレックス法で根を再最適化する。MIP の途中で時間切れになっても、根 LP が最適に解けて保存が完了していれば次回に利用できる。

```cpp
cpopt::Options o;
o.time_limit = 0.02;             // 各 solve() に20 ms
cpopt::Result r = m.solve(o);
m.add_le({{x, 1}, {y, 1}}, 8);
// 以前の暫定解が追加制約を満たすなら初期解として渡せる。
if (r.has_solution && m.is_feasible(r.x)) o.initial_solution = r.x;
else o.initial_solution.clear();
r = m.solve(o);                  // 保存した根基底から再最適化
// r.stats.reused_root: 再利用を試みたか
// r.stats.root_pivots: 根の求解・再最適化・復旧で使った pivot 数
```

保存対象はカット・分枝を加える前の最適な根。探索木、擬似コスト、カット、整数暫定解は次の呼び出しへ自動では持ち越さない。暫定解を再利用する場合は上のように検査して `initial_solution` を設定する。無効な初期解は例外となる。

変数追加・目的関数再設定は保存基底を破棄する。前処理・スケーリング・自動双対化・数値許容誤差・整数緩和の設定が変わった場合も作り直す。`o.reuse_root=false` で保存・再利用を無効化でき、`m.clear_root_cache()` で保存用メモリを解放できる。連続 LP でも主問題の根は保存するが、自動双対化で解いた結果はこのキャッシュの対象外。

再利用時は古い前処理の変数対応を維持するため、追加制約による新たな固定変数消去や境界伝播は行わない。根の pivot が減っても探索全体が速くなるとは限らず、制約が大きく変わる場合は `clear_root_cache()` して再度前処理する選択肢もある。

保存・復元・追加行の変換にも同じ時間予算を適用する。途中で停止しても古い完全なスナップショットを保持し、不完全な基底は次回に使わない。数値不整合や再最適化の停滞時には、残りの総予算で一度だけ前処理から作り直す。`Model` のコピー間では不変の保存データを共有するが、作業用基底は独立する。同じ `Model` を複数スレッドから同時に `solve()` する使い方には対応しない。

## 実行時間を制限する

```cpp
cpopt::Options o;
o.time_limit = 1.8;                  // solve() 内の経過秒
o.node_limit = 100000;
o.iteration_limit = 1000000;          // strong branching やカットも含む総 pivot 数
o.initial_solution = known_solution; // 任意。元の変数順の実行可能解
auto r = m.solve(o);
if (r.has_solution) {
    // Optimal でなく Limit の場合も、この暫定解は使える。
    use(r.x);
}
```

`initial_solution` が実行不可能なら `invalid_argument` を投げる。初期解を与えなくても、簡単な境界上の解や丸めた解を検査して暫定解を探す。前処理や根 LP の途中で止まると、暫定解がない場合もある。検査済みの解は前処理中の時間切れでも保持する。主単体法の Phase II の途中で止まった場合は、現在の基底からも実行可能解を回収する。MIP の整数性も検査するため、小数の LP 解を整数解として返すことはない。

時間は `solve()` 呼び出しごとに計測し、前処理、根 LP の求解、探索、保存基底の復元と追加制約の変換を含む。

厳密な実時間上限ではなく、メモリ確保・解放、1回の pivot、探索中の辞書コピーには割り込まない。指定した初期解の検査と、打切り時に現在の LP 解を返すための最終検査も完了させる。初期解は `time_limit=0` でも検査して返すため、大きなモデルではその検査時間も必要になる。提出全体の制限から余裕を引いて指定する。

| `Result` の値 | 意味 |
|---|---|
| `status` | 終了理由。下表参照 |
| `has_solution` | `x` と `objective` に実行可能な暫定解があるか |
| `x` | 元の変数番号順の解 |
| `objective` | 定数項を含む元の目的関数値。解なしなら NaN |
| `best_bound` | 最大化なら最適値の上界、最小化なら下界。未計算ならそれぞれ +inf / -inf |
| `gap()` | 暫定解と境界の絶対差。暫定解なしなら +inf |
| `stats` | ノード数、pivot 数、カット数、経過秒など |

| `Status` | 意味 |
|---|---|
| `Optimal` | 指定した数値許容誤差・ギャップ内で最適 |
| `Infeasible` | 実行可能解なし |
| `Unbounded` | LP の目的値が非有界 |
| `Limit` | 時間・ノード・反復の上限。暫定解の有無は別に確認 |
| `NumericalError` | 数値的な不整合を検出。この場合の `best_bound` は証明に使わない |
| `RelaxationUnbounded` | MIP の根 LP が非有界。MIP 自体の実行可能性・非有界性は未判定 |

`stats.nodes` は根の復旧・双対化の試行も含む処理した LP ノード数。strong branching の試行 LP は含めず、その pivot は `stats.strong_branching_pivots` と `stats.pivots` に含む。`node_limit=1` でも根での試行分枝は行う。`stats.numerical_restarts` はカットを取り除いた再試行と、保存基底からの再最適化を断念した復旧の回数。

## 調整できる項目

既定値は [Options](mip.hpp) を参照。

| 項目 | 既定値 | 用途 |
|---|---:|---|
| `presolve` | true | 境界伝播、固定変数・冗長行の除去、整数行の gcd 簡約 |
| `scaling` | true | 行を2の冪でスケーリング。小さすぎる目的係数の拡大は別途実行 |
| `dualize` | true | 変数が多い連続の被覆 LP 等を双対問題へ変換して解く |
| `warm_start` | true | 分枝先の LP を親の基底から再最適化 |
| `reuse_root` | true | 同じ `Model` の呼び出し間で最適な根基底を保存・再利用 |
| `pseudocost` | true | 過去の目的値の改善量から分枝変数を選ぶ |
| `rounding_heuristic` | true | 境界上・丸め後の解を検査して暫定解を得る |
| `integer_objective` | true | 整数係数の整数目的関数に対する上界の切り下げと固定の強化 |
| `reduced_cost_fixing` | true | 根・探索ノードの reduced cost と暫定解から、改善解で0と証明できる非基底整数列・slackを固定して辞書を縮小 |
| `strong_branching_candidates` | 4 | 浅いノードで試す分枝変数数。0で無効 |
| `strong_branching_depth` | 3 | 試行分枝を行う深さ |
| `strong_branching_iterations` | 40 | 各試行 LP の pivot 上限 |
| `root_cut_rounds` | -1 | -1: 小数値の整数変数が1個の根で3巡、0: 無効、正数: 指定巡数 |
| `cuts_per_round` | 8 | 各巡で加える GMI カットの上限 |
| `absolute_gap` / `relative_gap` | 1e-7 / 0 | 暫定解から許容するギャップ |
| `feasibility_tolerance` | 1e-7 | 元の変数・制約での実行可能性判定 |
| `integrality_tolerance` | 1e-7 | 整数とみなす絶対誤差 |
| `pivot_tolerance` | 1e-11 | simplex 内部の係数・可行性の判定 |

カットを常に増やせば速くなるわけではない。多次元ナップサックや安定集合では、探索ノードが減っても各 LP が重くなる例がある。`root_cut_rounds=0` と比較するとよい。小さい LP は前処理の時間が支配的なこともある。

`r.stats.reduced_cost_fixings` は、探索中に reduced cost と暫定解を使って固定した列数の全ノード合計。固定は現在の暫定解に依存するため、次回用の根基底には持ち越さない。

## 数値計算と規模の限界

`double` による浮動小数点ソルバーで、厳密な有理数・整数証明器ではない。係数の桁が極端に違う、ほぼ平行な制約が多い、巨大な整数を1の差まで区別する場合には誤判定の可能性がある。`is_feasible` は絶対許容誤差に加えて、行の活動量に `1e-12` を掛けた誤差を許す。整数の出力は、`Optimal` または実行可能な暫定解を確認してから丸める。

根の tableau の係数は概ね `8 * (制約数+2) * (変換後変数数+2)` バイトで、行を管理する領域も必要。元の行列・変換情報も別に保持する。根キャッシュは根 tableau を保持し、変換後の不変な係数行は作業領域と共有する。幅の広い辞書では変更されない辞書行も共有する。不要なら `reuse_root=false` とし、既に保存済みなら `clear_root_cache()` で解放する。等式は2本の不等式、自由変数は2変数に変換し、有限上限やカット・分枝が行を増やす。分枝限定では待機ノードごとに変更した行や狭い辞書を持つため、LP より多くのメモリを使う。固定列を除いた探索中の辞書はその分だけ小さくなる。

## 使用例

- [abc224_h.cpp](abc224_h.cpp): 二部グラフの被覆
- [abc326_g.cpp](abc326_g.cpp): 連続変数と0/1変数を組み合わせた定式化
- [abc338_c.cpp](abc338_c.cpp): 整数の資源配分
- [abc358_c.cpp](abc358_c.cpp): 集合被覆
- [abc364_e.cpp](abc364_e.cpp): 2制約のナップサック
- [abc354_g.cpp](abc354_g.cpp): 部分文字列の包含関係を制約にした文字列選択
- [abc374_e.cpp](abc374_e.cpp): 予算制約下での整数の機械購入と生産量最大化
- [abc397_g.cpp](abc397_g.cpp): 辺の重み選択と最短距離の最大化
- [arc137_e.cpp](arc137_e.cpp): 区間ごとの雇用と日ごとの販売量の最適化（LP）
