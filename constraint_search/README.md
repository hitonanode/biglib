# 短時間の整数・大域制約ヒューリスティック

**[constraint_search.hpp](constraint_search.hpp) だけで動く**、外部依存のない C++17 ヘッダです。有限整数・Boolean・線形制約に加え、AllDifferent、NoOverlap、Circuit を直接扱います。10 ms〜2秒で使う近似探索を想定し、部分解・制約違反のある割当からの修復と、探索状態を保持した再開を実装しています。GCC / Clang の `__int128` を使用します。

## 使い方

```cpp
#include "constraint_search.hpp"
#include <iostream>

int main() {
    using namespace constraint_search;
    Model m;
    int a = m.add_integer(0, 30);
    int b = m.add_integer(0, 30);
    int c = m.add_integer(0, 30);
    m.add_no_overlap({Interval(a, 4), Interval(b, 7), Interval(c, 3)});
    m.add_ge({{c, 1}, {a, -1}}, 4); // c >= a + 4
    int end = m.add_maximum({Affine(a, 1, 4), Affine(b, 1, 7), Affine(c, 1, 3)});
    m.minimize({{end, 1}});

    Options o;
    o.time_limit = 0.010;
    o.seed = 42;
    o.initial_solution = {0, unknown, 4}; // b は未指定、end は自動計算
    auto r = solve(m, o);
    if (r.feasible) std::cout << r.objective << '\n';
}
```

```sh
c++ -std=c++17 -O3 constraint_search/solve.cpp -o /tmp/cp-example
/tmp/cp-example
```

## モデルの契約

変数 ID は追加順の0始まりです。`Terms` は `{変数ID, 整数係数}` の配列、`Affine(v,a,b)` は `a*x[v]+b`、定数は `Affine::constant(b)`。制約値・目的値は整数で計算し、入力・式の安全な範囲を保守的に `±2^60` に制限します。探索用の重み・スコアと公開する違反量の合計は浮動小数点です。実行可能性判定は浮動小数点の許容誤差に依存しません。線形式では `Σ |a_i| max(|lo_i|,|hi_i|)` も制限するため、相殺で収まる一部の大きな式は拒否します。

| API | 意味 |
|---|---|
| `add_integer(lo,hi)` / `add_bool()` | 閉区間整数 / 0,1 の決定変数 |
| `add_linear(terms,lo,hi,guards={})` | 両側線形制約。重複変数を集約 |
| `add_le` / `add_ge` / `add_eq` | 片側線形制約・等式 |
| `pos(v)` / `neg(v)` | Boolean literal |
| `add_bool_or(literals,guards={})` | 少なくとも1つ真 |
| `add_implication(a,b)` | `a ⇒ b` |
| `add_bool_xor(literals)` | 真の個数が奇数 |
| `add_all_different(vars)` | 変数の値が相異なる |
| `add_all_different(affines,guards={})` | affine 式の値が相異なる |
| `add_no_overlap(intervals,guards={})` | 存在する区間が順序を付けて重ならない |
| `add_circuit(successors,allow_self_loops=false)` | 後続変数による巡回路 |
| `add_element(index,constant_table)` | 配列参照を表す**従属変数**を作り、その ID を返す |
| `add_maximum(affines)` | 最大値を表す**従属変数**を作り、その ID を返す |
| `minimize(terms,offset=0)` / `maximize(...)` | 整数線形目的。再指定は置き換え |
| `evaluate(full_assignment)` | 全制約・従属値の再計算による検証 |

`guards` は全 literal が真のときだけ制約を課す含意です。同値条件ではありません。cardinality は Boolean の和の線形制約で書けます。

`Interval(start_var,fixed_duration,presence={})` と `Interval(start_affine,duration_affine,presence={})` を用意しています。duration の定義域は非負にしてください。省略した presence は真。長さ0の点も順序制約に参加し、例えば `[5,5]` は `[0,10]` と衝突しますが、`[0,0]` と `[0,10]` は衝突しません。

Circuit の各後続値は、渡した配列内の頂点番号 `0..n-1` です。既定は全頂点を通る1巡回（n=1のみ自己ループ可）。`allow_self_loops=true` は自己ループを不参加とする subcircuit で、非自明な閉路は高々1つ、全自己ループも実行可能です。

element / maximum は既存変数から一方向に計算します。従属変数への hint は再計算され、探索の自由度にはなりません。index の全定義域が配列内に収まる必要があります。


## warm start と途中再開

`Options::initial_solution` は変数順の prefix または全長配列です。未指定箇所を `unknown` にできます。決定変数の hint は定義域内であれば、制約違反を含んでも受け付けます。hint は固定ではなく、修復・改善の出発点です。定義域外は `invalid_argument`、長すぎる配列もエラーです。

```cpp
Options options;
options.initial_solution = {0, unknown, 4};
Solver search(model, options);
search.run(0.010);
const Result& r = search.run(0.010); // キャッシュ・重み・乱数状態も継続
// r は次の run / take_result / search の破棄まで有効
```

`Model` は `Solver` より長く生存させ、探索中は変更しないでください。モデル変更は次の `run` で `logic_error` になります。異なる `Solver` は独立ですが、同じオブジェクトを複数スレッドから同時に呼べません。

`run` が返すのは探索中に保存した最良結果です。実行可能解を得た後は実行可能性を失わず、目的値も悪化しません。実行可能解がなければ、生の違反量が最も小さい検証済み割当を返します。終了時の探索地点と最良結果は異なる場合があります。

`take_result()` は結果を移動して取り出す最終操作です。その後の同じ Solver の再利用はしないでください。通常は `result()` または `run()` の参照を使います。

## 時間と結果

`solve(model,options)` は呼び出し開始から初期化・探索・キャッシュ解放まで、`Solver::run(seconds)` はその呼び出し内の処理を時間予算に含めます。モデル構築、呼び出し元による `Result` のコピー、`Model::evaluate` は別です。`Options::deadline` に `steady_clock` の絶対締切を渡すと、相対制限との早い方を使います。絶対締切は再開後も維持されます。

内部では128回の短い処理ごとと主要ループで時計を確認し、中断した候補評価を巻き戻します。初期化と採用 move の一部は次回へ持ち越します。返却のための全制約スキャンは行いません。`run` は150 µs、one-shot `solve` はさらにモデル規模に応じた解放時間を先取りします。SAT符号化の対象になり得るモデルには、学習節・watch listの解放分も予約します。この予約は絶対締切にも適用します。

**1 ms以内の超過を目標とした協調的な停止で、ハードリアルタイム保証ではありません。** OS のプリエンプション、メモリ確保・解放、極端なモデルでは超過し得ます。実行可能解が必要な呼び出しでは、あらかじめ `Solver` を作って初期化した上で短い `run` を繰り返す使い方が適しています。巨大な1制約の初期化より短い slice だけを繰り返すと、その初期化が完了しない場合があります。

| 結果 | 解釈 |
|---|---|
| `evaluated && feasible` | 全制約・定義域・従属値を満たす保存解 |
| `evaluated && !feasible` | 違反のある保存割当。充足不能の証明ではない |
| `!evaluated` | 初期評価まで完了していない。`x` は空であり得る |
| `objective` | 元の整数目的値。違反罰則を混ぜない |
| `violation` / `violated` | 生の違反量の合計 / 違反制約数 |
| `stats.elapsed` / `total_elapsed` | 今回 / 累積の内部求解時間 |
| `stats.first_feasible` | 累積求解時間内の初実行可能解時刻、未発見は -1 |

冷開始の `run(0)` / `solve` 時間0は hint の実行可能性を未検証のまま返しません。既存 Solver の `run(0)` は既に検証済みの最良解をそのまま返します。`StopReason::ObjectiveBound` は変数の定義域だけで得られる目的下界／上界に達した場合です。一般的な最適性・充足不能証明機能は持ちません。

時間指定の結果は実行負荷に依存します。`seed` と `iteration_limit` を指定すれば、十分長い時間予算下で反復数を固定した検証もできます。`iteration_limit` は各 `run` の上限です。

## 探索と設定

線形行、AllDifferent、NoOverlap、Circuit と従属値を増分評価し、境界へのジャンプ、交換、巡回路の移動、適応重み、tabu、摂動を組み合わせる。モデルの構造に応じて、順列の初期割当と領域削減、SAT による残余 exact cover の探索、packing の小規模な分枝限定、MaxSAT の1変数・2変数反転、スケジュールの機械順序の改善などを使う。

| `Options` | 既定値 | 用途 |
|---|---:|---|
| `time_limit` | 1.0 | `solve` の時間予算（秒） |
| `seed` | 1 | 乱数シード |
| `iteration_limit` | uint64 の最大値 | 各 `run` の反復上限 |
| `deadline` | 時刻の最大値 | 複数回の再開に共通する絶対締切 |
| `initial_solution` | 空 | 部分割当・初期解 |
| `candidates` | 24 | 一般近傍の候補数（1〜1024） |
| `clause_pair_samples` | 8 | 2変数反転の候補収集に使う不満足節の標本数（1〜64） |

個別の探索機構を切り替えるフラグはすべて既定で有効。詳細は [Options](constraint_search.hpp) を参照。

## 使用例

- [solve.cpp](solve.cpp): 作業の重なりと先行関係を制約にした、完了時刻の最小化。部分割当から探索を始め、途中再開する。
