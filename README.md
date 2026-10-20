# biglib

C++17 以降に対応する、標準ライブラリだけで利用可能な C++ ヘッダーオンリーライブラリ。
最大重みマッチングとシュタイナー木は整数の総重み・費用に `__int128` を使うため、GCC / Clang が対象。

## ライブラリ一覧

| 分野 | ヘッダー | ドキュメント |
| --- | --- | --- |
| 線形計画・混合整数計画問題 | [mip.hpp](mip/mip.hpp) | [使い方](mip/README.md) |
| 一般化割当問題 | [generalized_assignment.hpp](generalized_assignment/generalized_assignment.hpp) | [使い方](generalized_assignment/README.md) |
| 対称巡回セールスマン問題 | [symmetric_tsp.hpp](symmetric_tsp/symmetric_tsp.hpp) | [使い方](symmetric_tsp/README.md) |
| 非対称巡回セールスマン問題 | [asymmetric_tsp.hpp](asymmetric_tsp/asymmetric_tsp.hpp) | [使い方](asymmetric_tsp/README.md) |
| 一般グラフの最大重みマッチング | [maximum_weight_matching.hpp](maximum_weight_matching/maximum_weight_matching.hpp) | [使い方](maximum_weight_matching/README.md) |
| 最小シュタイナー木 | [steiner_tree.hpp](steiner_tree/steiner_tree.hpp) | [使い方](steiner_tree/README.md) |
| 重み付き集合被覆・集合分割 | [set_cover_partition.hpp](set_cover_partition/set_cover_partition.hpp) | [使い方](set_cover_partition/README.md) |
| 二次割当問題 | [quadratic_assignment.hpp](quadratic_assignment/quadratic_assignment.hpp) | [使い方](quadratic_assignment/README.md) |

## 使い方

このリポジトリのルートをインクルード検索パスに追加する。

```sh
c++ -std=c++17 -O2 -I /path/to/biglib main.cpp -o main
```

各ヘッダーは単独で使え、必要なファイルの内容を提出コードへ貼り付けることもできる。
ヘッダーをソースと同じディレクトリに置く場合、`-I` の指定は不要。

## 構成

ヘッダーと同名のディレクトリに、`.hpp`、`README.md`、使用例の `.cpp` を配置する。

## コード整形

[.clang-format](.clang-format) の設定で clang-format を使用する。

```sh
make format        # .hpp / .cpp を整形
make format-check  # 整形済みか確認
```
