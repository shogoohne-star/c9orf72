# ドキュメント

鋳型指向型 IEDDA ライゲーション設計の計算結果をまとめた 2 種類の HTML。

| ファイル | 位置づけ |
|---|---|
| `summary.html` | ラボミーティング発表用。手法の要点・結果・そこから分かったこと |
| `detail.html`  | 詳解版。使用ツールと出典、鋳型モデルの検証、判定式、EM の定義、全数値、感度解析 |

## 生成方法

```
cd docs
python3 mkcharts.py     # グラフ SVG（assets/chart_*.svg, legend_*.html）
python3 mkfigs.py       # 概念図 SVG（assets/fig_*.svg）
python3 build.py        # *.tpl.html + style.css + assets -> summary.html / detail.html
```

`assets/{A,B,C}_{tz,cp}.svg` は RDKit で描いたリンカーの構造式。`build.py` が
XML 宣言と固定サイズを除去し、線の色を `currentColor` に置換してテーマ追従にする。

数値の出どころは `../iedda_template_sim_v18.py` の Colab 実行結果。
