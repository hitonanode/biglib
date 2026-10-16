# Asymmetric TSP ソルバ

[asymmetric_tsp.hpp](asymmetric_tsp.hpp) だけで使える、C++17・標準ライブラリのみのヘッダオンリー実装。数百〜数千頂点の、有限コストの完全有向グラフを対象とする。頂点を1回ずつ通り、出発点に戻る巡回路の費用を最小化する。非対称、非距離的、負のコストに対応する。

AHC などで、短い時間予算を指定して何度も使う用途を想定している。ソースをそのまま提出コードの先頭に貼り付けられる。

## 使い方

```cpp
#include "asymmetric_tsp.hpp"

std::vector<std::vector<long long>> cost = {
    {0, 1, 9, 9},
    {9, 0, 1, 9},
    {9, 9, 0, 1},
    {1, 9, 9, 0},
};
asymmetric_tsp::Options o;
o.time_limit = 1.8;  // solve() 全体に1.8秒
o.seed = 42;
auto r = asymmetric_tsp::solve(cost, o);
// r.tour: 頂点0..n-1の順列。末尾から先頭へ戻る辺も含む。
// r.cost: その費用（この例では4）。先頭頂点の重複は含まない。
```

行列のコピーは行わない。コストをその場で高速計算できる場合や、平坦な配列を使いたい場合は、コールバックでも渡せる。

```cpp
int n = /* 頂点数 */;
std::vector<int> flat(static_cast<std::size_t>(n) * n);
// flat を埋める
auto r = asymmetric_tsp::solve(n, [&](int i, int j) {
    return flat[static_cast<std::size_t>(i) * n + j];
}, o);
```

コールバックは同じ引数に同じ有限値を返し、短時間で終了すること。I/O やロック待ちには向かない。探索中も元のデータを変更しない。自己ループの値は読まない。頂点数0/1の費用は0、2なら往復の費用を返す。

## 時間制限と再利用

`time_limit` は `steady_clock` による実経過時間。初期解生成、行・列の減算、割当問題、閉路結合、1-tree ascent、候補集合の生成も同じ予算に含む。入力行列の作成・読み込みは呼び出し元の処理。

大きな走査は約256要素/探索節点ごとに時計を確認する。最良の巡回路を保持し、前処理の途中でも打ち切って返す。既定では最後の2 msを返却のために残す。短い予算の場合は余裕を予算の10%までに制限する。

プロセス全体の締切に合わせる場合:

```cpp
const auto deadline = asymmetric_tsp::Clock::now() + std::chrono::milliseconds(1900);
// 入力や別の処理……
o.deadline = deadline;
auto r = asymmetric_tsp::solve(cost, o);
// time_limit と deadline の、早い方までに探索を切り上げる。
```

**ハードリアルタイムの保証ではない。** OS による停止、メモリ確保・解放、コールバック1回、採用した交換の反映・巡回路評価・起点の再活性化には割り込めない。`time_limit=0` や期限切れでも、初期巡回路の構築・検査・費用評価の O(n) は実行する。呼び出し元の入出力時間を含めた余裕は別途確保する。

前回の解を渡すと、それより悪い巡回路は返さない。

```cpp
o.initial_tour = r.tour;
o.time_limit = 0.05;
r = asymmetric_tsp::solve(cost, o);
```

保存・再利用するのは巡回路。割当の双対変数、1-treeのポテンシャル、候補集合は呼び出しごとに作り直す。初期解なしでも、時間切れ時には `0,1,...,n-1` またはそれを改善した巡回路がある。初期巡回路の重複・範囲外の頂点・サイズ違いは例外とする。

## 探索の構成と計算量

1. 最近傍法と割当問題から初期巡回路を作る。割当で得られた複数の閉路は、費用増分の小さい接続で結合する。
2. 固定辺付き1-treeの下界計算を使って、局所探索で優先する候補辺を選ぶ。前処理が時間内に完了しなければ、利用可能な情報で探索を続ける。
3. 1〜3頂点のブロック移動、非対称2-opt、候補を制限した逐次交換で巡回路を改善する。
4. 改善が止まったら、二閉路への分割と再結合、最良解との部分経路の合成、より深い逐次交換を試す。
5. 巡回路の一部を組み替える摂動や初期解の再構築で探索を再開する。最良の巡回路は別に保持する。

候補・深さ・分岐数を制限するヒューリスティックであり、最適性や完全な k-opt 局所最適性は保証しない。品質は入力構造と時間予算に依存する。

以下では頂点数を `n`、頂点あたりの候補辺数を `K` とする。

距離行列を内部に複製しない。追加メモリは主に候補辺の O(nK) と巡回路・双対・全域木の O(n)。前処理には O(n²) の走査、中断可能な割当計算（最悪 O(n³)）、最大80回の O(n²) の1-tree構築がある。深い交換の作業領域は O(31K)。部分路の合成は O(n min(n,256)) の走査で、途中打切りできる。交換ごとに O(n) で費用・位置・反転差分を再計算し、変更した辺の本数を m として O(mK) で候補近傍を再活性化する。これらの更新は中断しないため、最大 O(nK) の反映時間を返却余裕に含める。

## オプションと結果

| 項目 | 既定値 | 意味 |
|---|---:|---|
| `time_limit` | 1.0 | 呼び出し全体の秒数。0で初期巡回路を返す。正の無限大も指定可 |
| `deadline` | 最大時刻 | 呼び出し元の絶対締切 |
| `safety_margin` | 0.002 | 返却に残す秒数 |
| `seed` | 1 | 乱数シード |
| `max_trials` | 最大値 | 摂動と局所探索の試行数。0で前処理も省略 |
| `candidates` | 24 | 各頂点から出る候補辺数 |
| `max_exchange` | 31 | 逐次交換の深さ上限、3〜31。完成する交換は奇数深さ |
| `search_breadth` | 6 | 2本目以降の追加辺で試す分岐数 |
| `search_nodes` | 256 | 浅い降下の1起点からの探索節点数上限 |
| `deep_search_nodes` | 1024 | 深い交換の1起点からの探索節点数上限 |
| `staged_search` | true | まず深さmin(5,max_exchange)で降下し、停滞後に深い探索を行う |
| `adaptive_depth` | true | coverの閉路数qが正で16q≤nのとき深い段階を使う |
| `lookahead` | true | 次に削除できる辺の費用を使って分岐順を決める |
| `merge_tours` | true | 最良解へ、安い同集合・同方向端点の部分路を移す |
| `ascent_candidates` | true | 固定辺付き1-tree ascentと感度による候補を使う |
| `assignment` | true | 割当を解いて初期解を強化する |
| `assignment_fraction` | 0.25 | 割当に使える時刻を全予算のこの割合に置く。ここまでの初期化時間も含む |
| `reversal` | true | 非対称2-optを使う |
| `structural_candidates` | false | 最小全域木の感度で候補を順位づける任意設定 |
| `active_nodes` | true | 変更のない起点の再探索を省く |
| `patching` | true | 局所探索後に二閉路を修復する非逐次4-optを試す |
| `guided_kicks` | true | 偶数試行の最初の摂動で、高い reduced cost の辺を切る |
| `restarts` | true | 12試行ごとにランダム化した貪欲巡回路から再開する |
| `relative_tolerance` | 1e-12 | 交換の改善判定の相対許容誤差 |
| `initial_tour` | 空 | 最初から保持する巡回路 |

小さい時間予算なら既定値を推奨する。深さや候補数を増やすと1試行が重くなり、同じ時間で品質が上がるとは限らない。時間によらない再現実験には `time_limit=+infinity`、有限の `max_trials`、同じ seed を使う。同じビルドと実行環境で再現できる。標準ライブラリの `shuffle` や浮動小数点形式が異なる環境間では一致を要求しない。

| 結果 | 意味 |
|---|---|
| `tour`, `cost` | 最良の巡回路とその費用 |
| `lower_bound` | 割当双対と完成した1-tree反復による参考下界の最大値。有効な下界をまだ得ていなければ負の無限大 |
| `reason` | `time_limit` / `trial_limit` / `trivial`（0〜3頂点を直接解いた） |
| `stats.elapsed` | 内部で計測した秒数。返却後の内部メモリ解放は含まない |
| `stats.trials` | 局所探索の試行数。最後の途中打切りを含む |
| `stats.moves`, `stats.improvements` | 採用した局所交換数、最良巡回路の更新数 |
| `stats.search_nodes` | 逐次交換の探索節点数 |
| `stats.patching_moves` | 採用した非逐次4-optの数（`moves`にも含む） |
| `stats.merging_moves` | 採用した部分経路合成の数（`moves`にも含む） |
| `stats.assignment_cycles` | 結合前のcoverの閉路数。取得できなければ0 |
| `stats.ascent_iterations` | 完了した1-tree反復数 |
| `stats.one_tree_candidates` | 固定辺付き1-treeを候補に利用できたか |
| `stats.tree_candidates` | 1-treeまたは任意MSTを候補に利用できたか。候補全体の完成は `candidates_complete` |
| `stats.assignment_complete`, `stats.candidates_complete` | 各前処理が完了したか |
| `stats.assignment_repaired` | 打ち切った部分割当を補完し、cycle cover を作れたか。単一閉路への結合完了は意味しない |

コストは内部で `long double` に変換する。環境によって `double` と同じ精度の場合がある。総和が大きく整数が正確に表せない場合や、巨大な値と微小な差を混ぜる場合の厳密整数最適化は対象外。採用前に巡回路費用を再計算し、丸め誤差による偽の改善で巡回し続けることを防ぐ。`lower_bound` も浮動小数点計算なので厳密な証明用の値ではない。

`NaN`・無限大・内部計算がオーバーフローしうる値は、読んだ時点で `invalid_argument` とする。時間切れのときに未読の全辺の妥当性まで検査する API ではない。欠辺を `+infinity` で渡す形式には対応しない。有限の大きな罰金を使う場合、その辺を通らない保証はなく、返された巡回路を呼び出し元で確認すること。

## 使用例

- [solve_matrix.cpp](solve_matrix.cpp): 標準入力から費用行列を読み、費用と巡回路を出力
