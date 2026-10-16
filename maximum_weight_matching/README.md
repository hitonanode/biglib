# 一般グラフの最大重みマッチング

**[maximum_weight_matching.hpp](maximum_weight_matching.hpp) だけで使える、外部依存なしの C++17 ソルバ。** 一般の無向グラフを primal-dual blossom 法で解く。整数重みでは多項式時間で厳密な最適解を返す。負重み、多重辺、任意の初期マッチング、時間・処理量による停止と再開に対応する。

GCC / Clang が対象。整数の総重みと証明書は処理系拡張 `__int128` で返す。内部演算は `32 max(1,n) max|w| <= INT64_MAX` なら64 bit、それ以外は128 bitへ自動切替し、64 bit全域の入力を厳密に扱う。`double` も使えるが、こちらは許容誤差つきの数値計算となる。

## 最初の例

```cpp
#include "maximum_weight_matching.hpp"
#include <iostream>

using Solver = maximum_weight_matching::Solver<>; // 重み: int64_t
int main() {
    std::vector<Solver::Edge> edges = {{0,1,6}, {0,2,10}, {1,3,1}};
    Solver solver(4, edges, {0}); // 任意。辺ID 0を初期マッチングにする
    auto r = solver.solve();
    std::cout << maximum_weight_matching::to_string(r.weight) << '\n'; // 11
    // r.edges = 選択した元の辺ID。r.mate[v] = 相手頂点、未マッチなら -1。
}
```

頂点は `0..n-1`、辺 ID は入力 `vector` の添字。`Solver<long long>`、`Solver<int>` も使える。自己ループ、不正な端点、重複頂点を含む初期マッチングは `std::invalid_argument`。空グラフ・孤立点は許す。目的関数は採用辺の重みの総和。重み0の辺を採用するかは保証しない。

## 任意の初期解を渡す

コンストラクタの第3引数に、頂点を共有しない元の辺 ID の配列を渡す。部分マッチングでよく、負重み辺や最良ではない多重辺も指定できる。初期解は固定制約ではなく、辺を取り替えて改善するための開始情報である。

```cpp
Solver first(n, edges, initial_edge_ids);
auto a = first.solve();

// 重みを変更した別の問題。辺の順序・端点を維持していれば同じIDを使える。
edges[7].weight += 20;
Solver changed(n, edges, a.edges);
auto b = changed.solve();

// 同じグラフで初期化からやり直す場合。
first.reset(other_edge_ids);
```

初期解は暫定解として保持し、未使用頂点を貪欲に補完して探索を始める。返す解の重みは初期解以上で、同じ Solver の再開間でも非減少となる。

`stats.retained_seed_edges` は探索開始時に保持できた指定辺数、`stats.initial_search_edges` は貪欲補完を含む開始辺数。

## 停止と再開

```cpp
maximum_weight_matching::Options options;
options.time_limit = 0.01; // 今回の solve() 内の10 ms
auto r = first.solve(options);
if (r.status == maximum_weight_matching::Status::Limit) {
    use_matching(r.edges);          // 常に実行可能、初期解以上
    r = first.solve(options);       // 双対・blossom・交互森・走査位置から再開
}
auto snapshot = first;              // 通常のコピーで独立した状態を保存できる
```

`step_limit` は1回の呼び出しで許す内部ステップ数。辺走査1本、双対更新・森操作等が1ステップに相当し、固定計算量の単位ではない。既定は時間・ステップとも無制限。`time_limit=0` または `step_limit=0` では初期化済みの暫定解を返す。`reset()` は状態を破棄する操作で、再開には不要。

時間制限は soft limit。256ステップごとの時計確認に加え、1回の縮約・展開・増大・候補修復は完了させる。**コンストラクタ／`reset()` の入力検査・貪欲初期化・双対初期化と、結果のコピーは `time_limit` の外**。初期化と返却の計算量は下表を参照。短い全体予算を使う場合はこれらの時間も見込む。

## 結果と最適性の証明書

| 値 | 内容 |
|---|---|
| `status` | `Optimal`、`Limit`、`NumericalError` |
| `edges` | 元の辺 ID の配列 |
| `mate` / `mate_edge` | 各頂点の相手／採用辺 ID。未マッチは `-1` |
| `weight` | 総重み。整数版は `__int128`、小数版は `long double` |
| `dual_bound_twice` | 双対目的値の2倍。整数版では最適値の厳密な上界の2倍 |
| `tolerance` | 小数版の倍化 slack 判定許容誤差。整数版では0 |
| `stats` | 累積の各操作数、初期辺の保持数、`solve()` 内の経過秒、実際に使用した `algorithm` |

整数版の `Optimal` では `dual_bound_twice == 2*weight`。`Limit` でも `2*weight <= dual_bound_twice` なので上下界が分かる。`certificate()` は倍化した頂点双対と、正の双対を持つ奇数集合を返す。全辺の双対制約、奇数集合のサイズ、primal/dual の一致を別プログラムで検査できる。証明書の奇数集合を書き出す量は最悪 `O(n²)`。

`NumericalError` では暫定マッチングは利用できるが、双対値を保証として使わない。整数は負の最小値・正の最大値を含む64 bit全域に対応し、総和が64 bitを越えても正確に扱う。

## 小数重み

```cpp
using FloatingSolver = maximum_weight_matching::Solver<double>;
FloatingSolver f(3, {{0,1,0.1}, {1,2,0.25}, {2,0,0.2}});
auto r = f.solve(); // 重みは long double。Optimal は数値許容誤差の範囲内。
```

内部は `long double`。macOS arm64 などでは `long double` が `double` と同じ精度になる。許容誤差は最大正重み・頂点数・内部型の epsilon に応じて決める。極小値では非正規化数の最小刻みも考慮する。NaN / 無限大は拒否し、中間計算を安全な範囲に収められない極大値には `overflow_error` を投げる。

差が許容誤差以下の候補の優劣や、整数版と同じ意味での厳密最適性・厳密上界は保証しない。`dual_bound_twice` も小数版では数値的な評価値である。小数桁数が決まっていて厳密性が必要なら、整数にスケールして渡す。

## アルゴリズムと計算量

Edmonds の primal-dual blossom 法で、奇数閉路を縮約しながら増大路を探索する。頂点数を `n`、入力辺数を `m` とし、密グラフ向けの `Dense` と疎グラフ向けの `Sparse` を使い分ける。

| 方式 | 求解の最悪時間 | 内部空間 | 整数の初期化 | 小数の初期化 |
|---|---|---|---|---|
| `Dense` | `O(n³+nm)`、単純グラフでは `O(n³)` | `O(n²+m)` | `O(n²+m)` | `O(n²+m log(m+1))` |
| `Sparse` | `O(n(n+m) log(n+m+2))` | `O(n+m)` | `O(n+m)` | `O(n+m log(m+1))` |

整数初期化の計算量は、固定64 bit幅の radix sort を使う場合の上界。

既定の `Automatic` は `m < n² / (⌊log₂ n⌋+1)` で疎版、それ以外で密版を選ぶ（空グラフも処理する）。理論上界に基づく目安であり、実測最速の選択を保証するものではない。実際の速度は頂点数や重み分布にも依存する。第4引数で指定できる。

```cpp
Solver dense(n, edges, initial_edge_ids, maximum_weight_matching::Algorithm::Dense);
Solver sparse(n, edges, initial_edge_ids, maximum_weight_matching::Algorithm::Sparse);
auto r = sparse.solve();
// r.stats.algorithm == maximum_weight_matching::Algorithm::Sparse
```

`result()` と `solve()` の結果返却は密版 `O(n)`、疎版 `O(n log(n+1))`。`certificate()` は集合の明示出力を含めて最悪 `O(n²)` 時間・出力空間を使う。上表の内部空間に、この出力先やコピーした solver の空間は含めない。細かく停止して多数回結果を返す場合は、その都度の返却費用が加わる。

## 使用例

[basic.cpp](basic.cpp): 初期マッチングを指定して解き、総重みと採用辺を出力する。
