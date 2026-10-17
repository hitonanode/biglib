# Symmetric TSP ソルバ

[symmetric_tsp.hpp](symmetric_tsp.hpp) だけで使える、C++17・標準ライブラリのみの対称巡回セールスマン問題ソルバ。最大数千頂点、10 ms〜2 秒程度の利用を想定し、時間指定・完全解と部分解からの開始・内部状態を保持した再開に対応する。

有限な対称コストの完全グラフを扱う。負コスト、小数、三角不等式を満たさないコストも入力できる。最適性や近似率の保証はない。

## 使い方

```cpp
#include "symmetric_tsp.hpp"

std::vector<std::vector<int>> cost = {
    {0, 1, 9, 1}, {1, 0, 1, 9},
    {9, 1, 0, 1}, {1, 9, 1, 0}
};
symmetric_tsp::Options o;
o.time_limit = 0.05;  // 関数呼び出し全体で50 ms
o.seed = 42;
auto r = symmetric_tsp::solve(cost, o);
// r.tour: 頂点0..n-1の順列。最後から先頭への辺も含む。
// r.cost: 巡回路の費用。上の例では4。
```

入力行列をコピーしない。平坦な行列や高速な距離関数も渡せる。

```cpp
int n = 1000;
std::vector<double> flat(static_cast<std::size_t>(n) * n);
// flatを埋める。
auto distance = [&](int i, int j) {
    return flat[static_cast<std::size_t>(i) * n + j];
};
auto r = symmetric_tsp::solve(n, distance, o);
```

## 途中解・不完全な解から始める

`initial_tour` には、全頂点の順列でも、一部の頂点だけの訪問順でも渡せる。先頭頂点を末尾へ重複して書かない。

```cpp
o.initial_tour = r.tour;  // 前回の完全解。それより悪い解は返さない。
o.time_limit = 0.02;
r = symmetric_tsp::solve(cost, o);

o.initial_tour = {2, 0};  // 部分解。未出現頂点1,3を補い、完全巡回路にして改善。
r = symmetric_tsp::solve(cost, o);
```

複数の断片がある場合は `initial_paths` を使う。

```cpp
o.initial_tour.clear();
o.initial_paths = {{2, 0}, {3}};
r = symmetric_tsp::solve(cost, o);
```

両オプションは同時指定しない。断片内・断片間の重複や範囲外の頂点は例外。空の断片は無視する。未出現頂点を単頂点の断片として追加し、予算内で断片の結合順を改善する。時間がなくても入力順の連結と未出現頂点の追加で完全巡回路を確保する。

断片の内部順は初期構成で保つが、その後の探索では全頂点を自由に並べ替える。

## 探索状態を保持して再開する

何度も短い予算を追加する場合は `make_solver` を使う。

```cpp
auto solver = symmetric_tsp::make_solver(cost, o);
auto r1 = solver.run(0.01);
auto r2 = solver.run(0.01);  // 前処理も探索も続きから。r2.cost <= r1.cost。
auto r3 = solver.run(0.20);
```

現在解・最良解・候補辺・乱数・前処理と探索の進行状態を保持する。

Solver はコピー可能なので、`auto snapshot = solver;` でプロセス内の状態を複製できる。巡回路だけを別プロセスに引き継ぐ場合は、保存した `r.tour` を `initial_tour` に指定する。

距離関数は値として保持するが、行列用の overload は行列を参照する。**行列や callback が参照するデータは Solver より長く生存させ、探索中に変更しない。** 一時行列からの `make_solver` はコンパイルエラーにする。変更した問題では、新しい Solver を構築して巡回路だけを渡す。同じ Solver への並行した呼び出しは対応しない。

## 時間制限の契約

`time_limit` / `run(seconds)` は `steady_clock` による実経過時間。初期補完・初期配置・候補生成・1-tree ascent・探索を含む。行列の生成や入出力は呼び出し元の処理。`solve` は形状検査と Solver 構築も同じ締切に含む。`make_solver` の構築時にはオプションのコピーと検査を行い、O(n²) の前処理は行わない。

呼び出し元の締切も指定できる。

```cpp
const auto end = symmetric_tsp::Clock::now() + std::chrono::milliseconds(1900);
o.deadline = end;
auto r = symmetric_tsp::solve(cost, o);

// 再利用するSolverでは、呼び出しごとの絶対締切も渡せる。
auto solver = symmetric_tsp::make_solver(cost);
auto resumed = solver.run(0.05, end);
```

指定秒数、構築時の `Options::deadline`、`run` に渡した締切のうち最も早い時刻を使う。構築時の絶対締切は再開でも延長されない。既定で最後の2 ms、短い予算では最大10%を返却用に確保する。密な走査や探索節点では128回程度ごとに時計を確認する。

**ハードリアルタイムの保証ではない。** OS による停止、メモリ確保・解放、callback の1回、採用交換の O(n) 反映・費用検算には割り込めない。callback は同じ引数に同じ有限値を返し、短時間で終了すること。0秒や期限切れでも、最初の完全巡回路の検査・構築・評価の O(n) は必要。部分入力の空断片を大量に渡す場合は、その入力サイズ分の処理も必要になる。

## 探索の構成と計算量

1. 初期巡回路を構成し、部分解がある場合は未訪問頂点を補う。
2. 最近傍候補と、1-tree の下界計算から得られる候補辺を作る。
3. 2-opt、長さ1〜3の区間移動、候補を制限した逐次交換、経路回転で改善する。
4. 現在解と最良解の間で、同じ頂点集合を持つ部分経路を交換する。
5. 停滞時は double bridge による摂動で探索を再開する。

候補数や探索深さを制限するヒューリスティックであり、最適性や近似率は保証しない。

頂点数を `n`、候補数を `K` とすると、初期入力の保管分を除く追加メモリは `O(nK+n)`。最近傍・α候補の生成は `O(n² log K)`、疎グラフでの1-tree更新は1回 `O(nK log(nK))`、最後の完全グラフでの1-tree計算は `O(n²)`。前処理は呼び出し間で再開できる。

採用する交換の巡回路再構築は `O(n)`。部分経路の照合は1走査 `O(n min(n,512))` で、交換時には費用を再計算する。

## 設定と結果

構造に関するオプションは Solver の構築時に固定される。`run(seconds)` の追加予算と絶対締切は毎回指定できる。

| Options | 既定 | 意味 |
|---|---:|---|
| `time_limit` | 1.0 | 秒数。0と正の無限大も指定可 |
| `deadline` | 最大時刻 | 絶対締切 |
| `safety_margin` | 0.002 | 返却に残す秒数。予算の10%が上限 |
| `seed` | 1 | 乱数シード |
| `max_trials` | uint64最大 | 1回のrunで完了させる局所探索試行数。0で前処理も省略 |
| `candidates` | 24 | 最近傍・α 各行の上位数、2〜256。α行には最近傍を追加で一部保持 |
| `search_candidates` | 8 | 深い交換・経路回転に使う先頭候補数、2〜256 |
| `max_depth` | 10 | 逐次交換の削除辺数上限、2〜20 |
| `search_breadth` | 5 | 深い段で選ぶ分岐数。最初の追加辺はより広く探索 |
| `search_nodes` | 160 | 逐次k-optの1起点に対する展開数上限。patch・経路回転は別上限 |
| `ascent_iterations` | 80 | ペナルティ反復数の上限、0〜1000。疎版では時間配分でも打ち切る。0でα生成も無効 |
| `sparse_ascent` | true | 疎グラフで調整し、最後に完全グラフで確認 |
| `relocate` | true | 長さ1〜3の区間移動、反転も含む |
| `patching` | true | 2閉路を接合する限定4-opt |
| `chained_rotations` | true | 最大48段の経路回転 |
| `merge_tours` | true | IPTによる部分経路転写 |
| `kicks` | true | 停滞後の摂動。falseなら局所探索1試行で終了 |
| `guided_kicks` | false | 候補に沿う局所的なdouble bridgeも併用 |
| `active_search` | true | 段階別に深い探索の再訪起点を制限 |
| `partition_seed` | true | 初期解なし・512頂点以上の均衡二分構成 |
| `staged_search` | true | 深さ5までを先に探索 |
| `relative_tolerance` | 1e-12 | 利得判定の相対許容誤差 |
| `initial_tour` | 空 | 完全または部分訪問順 |
| `initial_paths` | 空 | 互いに頂点を共有しない経路断片 |

`Result::tour`, `cost` は最良の完全巡回路と費用。`lower_bound` は完全グラフの最小1-treeで確認した下界で、未計算なら負の無限大。疎グラフの値は含めない。浮動小数点計算の参考値であり、厳密な証明用ではない。

`reason` は `time_limit` / `trial_limit` / `trivial`（0〜3頂点）。`stats.elapsed` は今回の内部計測時間で、その他の統計は再開をまたぐ累積値。`trials` は完了した試行、`moves` は現在解に採用した改善（初期構成の改善も含む）、`two_opt_moves` / `relocate_moves` / `exchange_moves` は各近傍の改善数。`patching_moves` / `chained_moves` / `merging_moves` は追加した各処理の改善数で、`merging_moves` には最良解だけを変えた転写も含む。`improvements` は最良解更新、`kicks` は摂動回数、`search_nodes` は逐次k-optの分岐展開数。`nearest_rows` / `alpha_rows` と `ascent_iterations` で前処理の進み具合を確認できる。疎版の最後のdense確認は `ascent_iterations` に追加計上しない。

時間によらない再現実験には、`time_limit=+infinity` と有限の `max_trials`、同じ入力・seedを使う。同じビルドと実行環境での再現を想定する。有限時間による停止や異なる浮動小数点環境では一致を要求しない。

コストと総和は `long double`。プラットフォームによって `double` と同じ精度の場合があり、巨大な整数の細かい差や、極端に異なる桁を混ぜた厳密整数計算は対象外。自己辺は読まず、0/1頂点の費用は0、2頂点では往復を数える。非有限値・内部演算がオーバーフローしうる値はアクセス時に例外とする。全辺を先に検証する API ではなく、**対称性は呼び出し元の前提条件**。非対称入力には [asymmetric_tsp](../asymmetric_tsp/README.md) を使う。

## 使用例

[solve_matrix.cpp](solve_matrix.cpp): `n 秒数` と `n × n` の費用行列を標準入力から読み、費用と巡回順を出力する。
