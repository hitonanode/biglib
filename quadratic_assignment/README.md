# 二次割当問題（QAP）ソルバ

[quadratic_assignment.hpp](quadratic_assignment.hpp) だけで使える C++17・標準ライブラリのみのソルバ。最大数百要素、10 ms〜2 秒程度を想定し、時間制限、完全解・部分割当からの開始、内部状態を保持した再開に対応する。

施設 i の配置先を `p[i]` として、次の Koopmans–Beckmann 型 QAP を最小化する。

```text
sum(i=0..n-1, j=0..n-1) A[i][j] * B[p[i]][p[j]]
```

`p` は 0〜n−1 の順列。二重和には両方向と対角を含む。非対称、負の値、非ゼロ対角、小数に対応。一般の四添字費用、容量制約、禁止割当、固定割当、追加の線形費用は扱わない。大きい問題では最適性・近似率の保証はない。

## 使用例

```cpp
#include "quadratic_assignment.hpp"

std::vector<std::vector<int>> flow = {{0, 8, 2}, {8, 0, 3}, {2, 3, 0}};
std::vector<std::vector<int>> distance = {{0, 1, 4}, {1, 0, 2}, {4, 2, 0}};
quadratic_assignment::Options o;
o.time_limit = 0.01;  // 呼出し全体で10ms。小規模なら全列挙を終えて早く返る。
o.seed = 42;
auto r = quadratic_assignment::solve(flow, distance, o);
// r.assignment[i]: 施設iの配置先。r.cost: 目的値（この例では44）。
```

入力行列は所有コピーする。整数・小数の異なる型を混ぜてもよい。平坦な配列や関数でも渡せる。関数は構築時に各 n² 要素を一度ずつ読み、探索中には呼ばない。

```cpp
auto solver = quadratic_assignment::make_solver(n,
    [&](int i, int j) { return flat_flow[i * n + j]; },
    [&](int i, int j) { return flat_distance[i * n + j]; }, o);
auto r = solver.run(0.1);
```

[入出力例](solve_matrix.cpp) は `n 秒数`、A の n² 個の値、B の n² 個の値を読み、目的値と順列を出力する。

```sh
clang++ -std=c++17 -O3 quadratic_assignment/solve_matrix.cpp -o /tmp/qap
```

## 完全解・不完全な解から再開

`initial_assignment` は空、または n 要素の配列。`-1` は未割当を表す。既に割り当てた配置先の重複、不正な番号、サイズ違いは `std::invalid_argument`。

```cpp
o.initial_assignment = r.assignment;  // 完全解: これより悪い解を返さない。
o.time_limit = 0.02;
r = quadratic_assignment::solve(flow, distance, o);

o.initial_assignment = {2, -1, -1};   // 施設0→配置先2、残りを補完。
r = quadratic_assignment::solve(flow, distance, o);
```

未使用の配置先をシャッフルし、未割当施設へ一対一に補完する。**与えた割当は初期値であり、固定制約ではない。** 補完時は指定を守り、その後の探索では全施設を動かせる。部分解自体には完全な QAP の比較目的値がない。

追加時間で探索を続ける場合は、同じソルバを再利用する。

```cpp
o.initial_assignment.clear();
auto solver = quadratic_assignment::make_solver(flow, distance, o);
auto first = solver.run(0.01);
auto second = solver.run(0.05);  // 追加50ms。返却最良値は悪化しない。

auto checkpoint = solver;      // 行列を含む全探索状態のインメモリ複製。
auto best = solver.best();
auto current = solver.current_assignment();  // 探索中の配置。最良解とは異なりうる。
solver.set_start({-1, 0, -1});  // 新しい開始点。履歴・最良解・統計はリセット。
```

現在解・最良解、差分表と各行の最小値、配置順の行列キャッシュ、各表の初期化・更新位置、近傍走査位置、再評価中の総和、tabu、乱数、摂動回数、探索モード、再始動時の履歴消去位置を保持する。差分表を作り終える前の時間切れでも、次回は続きから進む。`set_start()` は乱数状態を維持する。完全な初期化には新しいソルバを作る。

`Result` は所有コピーなので次回の実行後も保持できる。`current_assignment()` はソルバ内部の `const vector<int>&`。`auto current = ...` ならコピーになる。配置配列を保存して再投入できる。小規模全列挙中の現在解は次の評価候補を指すことがある。

## 時間制限

単調時計 `std::chrono::steady_clock` を使う。O(n³) の初期差分計算、O(n²) の更新・候補走査・最良候補再評価を途中で中断できる。既定では期限の手前に `min(1ms, 残り予算×10%)` の返却余裕を取る。

```cpp
const auto deadline = quadratic_assignment::Clock::now() + std::chrono::milliseconds(1900);
o.deadline = deadline;
auto r = quadratic_assignment::solve(flow, distance, o);
// time_limit と deadline の早い方を使う。

auto more = solver.run_until(deadline);
// run(seconds, deadline) でも指定可。
```

- `solve(...)` は入力コピー・検証・初期解評価から計時する。
- `make_solver(...)` / コンストラクタと `set_start(...)` は、探索とは別に O(n²) の準備を行う。`run(...)` の秒数にはこれらを含まない。
- `time_limit=0` や期限切れでも、一括 `solve` の入力検証・メモリ確保・完全解構築・評価は必要。既存ソルバの `run(0)` は探索せず保存済み最良解を返す。
- コンストラクタへ渡した `Options.deadline` は継続実行にも適用する。期限を後へ延ばす用途では、ここを既定値にして各 `run` へ期限を渡す。
- 一般OSのスケジューラ停止、メモリ確保・解放、入力関数1回、O(n) の配列操作等への割込みはできず、ハードリアルタイムの上限保証ではない。

時間切れ時は再評価済みの最良解を返し、未完了の更新や候補検証は次回へ引き継ぐ。呼び出し元の入出力時間は別途確保する。

## 探索の構成と計算量

2施設の交換差分を初期化し、最良改善降下と適応摂動を組み合わせる `breakout`、tabu 探索、両者を切り替える `hybrid` を選べる。停滞時は最良解を保持してランダムな順列から探索を再開する。

交換差分の初期化は `O(n³)`、交換後の更新と候補走査は `O(n²)`。行列・差分・履歴・キャッシュの追加メモリは `O(n²)`。初期化や更新は途中で中断して再開できる。目的値と差分を定期的に再計算し、丸め誤差の蓄積を抑える。

既定では `n <= 8` を全列挙し、完了すれば評価した全順列の最小値を返す。全列挙も途中で停止・再開できる。

## オプションと結果

| オプション | 既定値 | 意味 |
| --- | ---: | --- |
| `time_limit` | 1.0 | 一回の秒数。0または正の無限大も可 |
| `deadline` | 最大時刻 | 絶対期限。各呼出しで適用 |
| `safety_margin` | 0.001 | 返却に残す秒数。実際は残り予算の10%以下 |
| `seed` | 1 | 乱数初期値 |
| `max_iterations` | 最大値 | 一回の交換数上限。全列挙中は評価置換数の上限 |
| `search` | `Search::hybrid` | `hybrid` / `breakout` / `tabu` |
| `exact_threshold` | 8 | 全列挙する最大サイズ。0〜10、0ならn≥2はヒューリスティック |
| `perturbation_fraction` | 0.15 | 基本摂動回数を n 倍した値、最低2交換。0より大きく1以下 |
| `stagnation_optima` | 100 | 強い摂動を行う非改善局所解数 |
| `hybrid_epoch` | 12 | 探索方式を切り替える交換数の n 倍係数 |
| `restart_stagnation` | 100 | 最良値更新・前回再始動から再始動までの交換数の n 倍係数。0で無効 |
| `relative_tolerance` | 1e-12 | 降下・aspiration の相対判定。目的項の絶対値総和でスケール |
| `initial_assignment` | 空 | 初期配置。n要素、-1で部分解 |

`run(seconds)` は `time_limit` だけをその呼出しで上書きする。反復上限は呼出しごとに適用し、統計は `set_start()` 以降の累積。前回の交換の表更新や最良値検証が未完なら、それを処理してから反復上限を判定する。

`Result.assignment`, `cost` が保存最良解。`reason` は `time_limit`、`iteration_limit`、`exhaustive`。`exhaustive` は全列挙の完了または n≤1 を表し、浮動小数入力の厳密算術による証明を意味しない。

`stats.elapsed` は今回の内部経過秒数（返却後のメモリ解放を除く）。ほかは累積で、`iterations` は交換/全列挙の評価数、`permutations` は全列挙だけの評価数、`improvements` は最良値更新数、`local_optima` はBLSの降下終了数、`perturbations` / `perturbation_moves` は摂動開始数/交換数、`mode_switches` は探索方式の切替数。`delta_evaluations` は初期化した差分と直接再評価した差分の個数、`delta_updates` は有効な差分の定数時間更新数（後で上書きするマスク要素への一時演算は含めない）、`rebuilds` は定期再構築開始数、`restarts` は停滞による再始動数。`symmetric` は対称経路、`batched_initialization` は安全な整数用初期化の選択、`delta_ready` は差分表の完成を表す（全列挙では表を作らない）。

時間によらない再現実験では `time_limit=+infinity` と有限の `max_iterations` を使う。同じビルド・環境・seedなら同じ探索になる。無限時間かつ無制限反復では、n≥2のヒューリスティック探索は停止しない。

## 数値と入力の範囲

内部の `Cost` は `double`。目的値の保存前評価は、安全な整数範囲では厳密な加算、それ以外では `long double` と補償和で行う。環境によって `long double` も `double` と同精度である。巨大整数の厳密計算、巨大な定数項に埋もれる微小な差の識別、最適性証明は対象外。整数入力も格納時に `double` へ変換する。

全入力要素を検査し、NaN・無限大と `sqrt(DBL_MAX)/(32*max(n,1))` を超える絶対値を拒否する。欠辺を無限大で表す形式は受け付けない。`-ffast-math` は使用しない。最大化する場合は片方の行列を安全な型で符号反転し、返却目的値の符号を戻す。
