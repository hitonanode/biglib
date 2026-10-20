# 最小シュタイナー木ソルバ

[`steiner_tree.hpp`](steiner_tree.hpp) だけで使用できる、無向・非負辺重みグラフの最小シュタイナー木ソルバです。C++17、標準ライブラリのみ。数十〜数百頂点、10 ms〜2 秒程度の予算を想定しています。整数の内部費用に `__int128` を使用するため、GCC / Clang 向けです。

最短路に基づく初期構築、局所探索、摂動と解の再結合を組み合わせます。端子が少ない場合は部分集合 DP による厳密解も求めます。

## 使用例

```cpp
#include "steiner_tree.hpp"

using Solver = steiner_tree::Solver<long long>;
std::vector<Solver::Edge> edges{
    {0, 3, 2}, {1, 3, 2}, {2, 3, 2}, {0, 1, 5}, {1, 2, 5}
};
std::vector<int> terminals{0, 1, 2};

Solver solver(4, edges, terminals);
auto r = solver.solve(0.010);       // 今回の呼び出しに10 ms
if (r.feasible) {
    // r.edges は元の edges 配列の添字。r.cost は合計費用。
}
r = solver.solve(0.100);            // 状態を保持して追加の100 ms
```

入力の頂点・辺 ID は0始まりです。端子の重複、平行辺、自己ループ、零費用を受け付けます。自己ループは解に使いません。端子が異なる連結成分にあれば `feasible=false`。空の端子集合と端子1個では費用0・辺集合が空の最適解です。

整数は符号付き・符号なし64 bitまで、負の費用は拒否します。`Solver<double>` / `Solver<long double>` にも対応します。範囲外の頂点・辺 ID、辺費用のNaN・無限大・負値、無効なオプションは `invalid_argument`。内部和・摂動を安全に計算できない極大の小数費用は `overflow_error` です。

## 初期解・不完全解と再開

```cpp
std::vector<int> fragments{3};      // 一部の辺だけでよい
Solver solver(4, edges, terminals, fragments, 42); // 最後は乱数seed
auto r = solver.solve(0.010);

auto checkpoint = solver;          // 探索状態の独立したコピー
r = checkpoint.solve(0.050);

Solver another(4, edges, terminals, r.edges); // 辺集合だけを引き継ぐ方法
auto completed = another.solve(0.010);
```

初期辺は、複数の離れた断片、閉路、重複、不要な枝を含んでも構いません。コンストラクタで、初期辺を優先して DSU で森林を作り、残りの入力辺で接続を補い、不要な葉を除きます。端子が接続可能なら、この時点で実行可能な木があるため `solve(0)` でも取得できます。初期辺集合が既に実行可能なら、その辺費用の総和を超える結果にはなりません。初期辺は変更可能なヒントであり、必ず採用する制約ではありません。

同じオブジェクトへの `solve()` は最良解、現在の探索木、乱数、elite pool、最短路の前処理、途中の近傍、厳密 DP を引き継ぎます。最良費用は悪化しません。Dijkstra の heap と辺走査位置、Floyd の pivot と行、DP の mask・分割・頂点、MST と経路復元の途中位置も保存します。完成していない候補が外へ返ることはありません。コピーで探索状態を独立に保持できます。

## 時間指定

```cpp
steiner_tree::Options o;
o.time_limit = 0.2;
o.deadline = steiner_tree::Clock::now() + std::chrono::milliseconds(150);
auto r = solver.solve(o);          // 相対予算と絶対締切の早い方
```

時間は `steady_clock` の壁時計です。32 work ごとに時計を確認し、終了処理の余裕を確保します。既定の余裕は残り時間の1%（最大2 ms）と0.2 msの大きい方で、残り時間の10%を上限とします。`safety_margin` で0.2 msの部分を変更できます。

`solve()` 内の探索・前処理・結果コピーを計測します。**コンストラクタの入力コピー・検査・初期解補完は、メンバー関数 `solve()` の予算外**です。準備も含めるには、次の補助関数を使います。

```cpp
steiner_tree::Options o;
o.time_limit = 0.01;
auto r = steiner_tree::solve<long long>(4, edges, terminals, o, fragments, 42);
```

補助関数は準備開始時に絶対締切を設定し、残り時間で探索します。ただし入力検査・コピー・最初の木の生成自体は中断しません。予算0、非常に短い予算、対象規模を大きく超えた入力でも、その準備費用は必要です。厳しい予算を繰り返し使う場合はオブジェクトを先に構築してください。

Floyd の1行、木の葉刈り、局所探索候補の準備、DP の最終復元整理などは途中で分断しない処理です。OSによる停止やメモリ確保もあるためハードリアルタイム保証ではありません。`stats.elapsed` はメンバー関数の結果構築までを含み、返却後のコピー・破棄や one-shot 内のソルバ破棄は含みません。

## オプション・結果

オプションは呼び出しごとに指定します。

| `Options` | 既定値・意味 |
| --- | --- |
| `time_limit` | 1.0 秒。0では探索を進めない。正の無限大で時間制限なし |
| `deadline` | 絶対締切。既定は無制限 |
| `safety_margin` | 0.0002 秒。上記の余裕計算に使用 |
| `max_work` | 無制限。1回の呼び出しで進める仕事単位の上限 |
| `exact_terminal_limit` | 11。0で部分集合 DP を停止。自明な厳密判定は有効 |
| `exact_memory_limit` | 64 MiB。部分集合 DP のセルの予測使用量上限 |
| `local_search` | true。局所探索を使用 |
| `recombination` | true。elite の辺を利用した再結合を使用 |
| `cascaded_recombination` | true。生成木とeliteの合成を、合計3回失敗するまで継続。falseで単発の弱い辺バイアス方式 |
| `vertex_perturbation` | true。新規構築時の半分で頂点ごとの係数を辺両端で平均。falseで辺ごとの独立摂動のみ |
| `fast_key_paths` | true。木の全頂点対を走査し、最良のkey-path交換をまとめて評価 |
| `steiner_reconnection` | true。分岐点除去後の成分MSTと、新しい中心頂点を経由する星型の再接続を比較 |
| `vertex_moves` | false。元グラフの誘導MST、頂点挿入・削除も評価する。追加評価の負荷があるためopt-in |

work は辺緩和、heap操作、Floydの1行、MSTの1頂点、DPの1組合せなどであり、一定のCPU時間ではありません。`time_limit=+∞` と `max_work` を使うと、同じ seed・オプションで再現可能な探索になります。仕事数を分割した呼び出しと1回の呼び出しは同じ探索経路を通ります。進行中の DP は `exact_terminal_limit=0` で一時的に止まり、正の値で再開します。未起動なら次回の上限・メモリ設定で起動可否を再評価します。

| `Result` | 意味 |
| --- | --- |
| `edges` | 元の辺 ID。昇順、重複なし。実行可能なら閉路・非端子葉のない木 |
| `cost` | 合計費用。整数版は `__int128_t`、小数版は `long double` |
| `feasible` | 全端子を接続できたか。false の費用を解として利用しない |
| `optimal` | ソルバ自身が最適性を確認したか |
| `lower_bound` | 下界。通常は0から開始し、自明な厳密ケースでは構築時から最適値。整数版では厳密 |
| `reason` | `time_limit` / `work_limit` / `optimal` / `infeasible` |
| `stats` | 今回の経過時間、累積 work・構築回数・近傍評価数・改善数・再結合回数・APSP/DP進捗 |

`result()` は探索せず現在の最良解を返します。`__int128_t` は標準の `ostream` に直接渡せないので、必要に応じて変換・出力処理を用意してください。

小数版では比較・総和・最短路に丸めがあり、`optimal` と下界は浮動小数点演算上の結果です。局所探索では相対的に32 epsilon以下の差を改善とみなしません。64 bit整数を小数へ変換して計算することはありません。正確な証明が必要で小数桁数が決まっている場合は、整数へスケールして渡せます。

## 採用したアルゴリズム

1. 初期辺を補完した木を確保し、SPH で改善します。最初の SPH より良い初期木は保持します。
2. 全点対最短路を計算します。`m > n²/8` なら Floyd–Warshall、それ以外は Dijkstra を各始点から実行します。
3. 端子間 metric-MST を元の辺へ展開し、MST と葉刈りで正規化します。これを完了すれば通常の KMB 型近似解以上の解を保持します。
4. key-path の交換、非端子分岐点と隣接 key-path 群の除去・再接続、端子＋分岐点の metric-MST における補助頂点の追加・削除を巡回します。key-path交換は現木の頂点対を一巡して最良候補を探します。分岐点除去では成分MSTに加え、任意の中心頂点から各成分への最短路を結ぶ星型候補を評価します。改善ごとに実木を再構成します。
5. 局所最適の後、最大16個のeliteから異なる木を選び、共通辺・片方だけの辺・その他の辺を区別した強い費用摂動でSPHと局所探索を行います。合成が改善したら親を更新し、合計3回の失敗で新規構築に戻ります。poolが満杯なら候補以上の費用のうち最も似た木を置換します。新規SPHは辺または頂点の重み摂動で多様化します。評価は常に元の費用です。
6. `vertex_moves=true` の場合は、誘導部分グラフのMST化と、元辺を使う頂点の挿入・削除も評価します。挿入は費用順にキャッシュした現木の辺と新しい隣接辺を併合し、同じ木辺を毎回heapへ投入する負荷を省きます。

端子数を `k`、`q=k-1` とすると、部分集合 DP は根の端子を除いた集合について「同じ頂点で2部分集合を併合→多始点 Dijkstra」で計算します。`k<=exact_terminal_limit`、`q<=20`、セルメモリ上限、`n*3^q<=40,000,000` を満たす場合に起動します。探索仕事3回につき DP 仕事1回を割り当て、完成したときだけ最適性を確定します。この比率は実行時間の比率ではありません。

端子0〜1個、端子を含む連結成分が木、費用0の実行可能解は直ちに厳密です。端子の連結成分の全頂点が端子なら Kruskal MST を求めます。また、端子間距離の最大値と `k * MST_metric / (2(k-1))` の大きい方を下界とし、整数版は後者を切り上げ、最良費用と一致すれば最適とします。端子2個の SPH は最短路ですが、`optimal` はDPまたは下界等による確認後に立ちます。

計算量の目安は次のとおりです。`k` は重複除去後の端子数、`L` は重複を含む入力端子数、`r` は metric-MST の必須頂点数、`c` は切断後の成分数です。表は完了までの上界であり、途中で予算により停止できます。

| 処理 | 時間 | 主な追加空間 |
| --- | --- | --- |
| 準備 | `O(L + (n+m+seed数) α(n) + n log n)` | `O(n+m)` |
| SPH | `O(k(n+m) log(n+m+1))` | 最悪 `O(km+n)` の遅延heapを含む |
| 全点対最短路 | 疎: `O(n(n+m) log(n+m+1))`、密: `O(n³+m)` | `O(n²+m)` |
| metric-MST候補 | `O(r² + rn log(rn+1) + n log n)` | `O(rn+n)` |
| key除去候補 | `O(n² + cn log(cn+1) + n log n)` | `O(c²+cn+n)` |
| key-path交換一巡 | `O(s² h + n log n)`。`s`は現木の頂点数、`h`は根付き木の高さ。経路展開・正規化は別 | `O(n+m)` |
| 元頂点挿入候補（opt-in） | `O(s + d log(d+1) + n log n)`。`d`は候補の元グラフ次数 | `O(n+d)` |
| 元頂点削除候補（opt-in） | `O(m log(m+1) + n log n)` | `O(n+m)` |
| 部分集合DP | `O(3^q n + 2^q(n+m) log(n+m+1))` | `O(2^q n + m)` |

## 実行可能な使用例

[basic.cpp](basic.cpp): 初期辺の補完、探索の継続、状態コピーの例。

```sh
c++ -std=c++17 -O2 steiner_tree/basic.cpp -o /tmp/steiner_example
/tmp/steiner_example
```
