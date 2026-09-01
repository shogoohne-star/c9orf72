# ドキュメント

鋳型指向型 IEDDA ライゲーション設計の計算結果をまとめた 2 種類の HTML。

| ファイル | 位置づけ |
|---|---|
| `summary.html` / `鋳型連結ASO_設計比較.docx` | ラボミーティング発表用。手法の要点・結果・そこから分かったこと |
| `detail.html` / `鋳型連結ASO_計算詳解.docx`  | 詳解版。使用ツールと出典、鋳型モデルの検証、判定式、EM の定義、全数値、感度解析 |

HTML 版と Word 版は同じ内容・同じ図表で、Word 版は A4 縦・余白 20/18 mm、
本文 10 pt（游ゴシック）。Word の見出しスタイルを使っているのでナビゲーション
ウィンドウで目次として使えます。

## 生成方法

```
cd docs
python3 mkcharts.py     # グラフ SVG（assets/chart_*.svg, legend_*.html）
python3 mkfigs.py       # 概念図 SVG（assets/fig_*.svg）
python3 build.py        # *.tpl.html + style.css + assets -> summary.html / detail.html

# Word 版
npm install             # docx (docx-js)
node render.js          # SVG -> png/*.png（Chromium でライトテーマ固定・3 倍解像度）
node mkdocx_summary.js  # -> 鋳型連結ASO_設計比較.docx
node mkdocx_detail.js   # -> 鋳型連結ASO_計算詳解.docx
python3 fixdocx.py *.docx   # Normal を既定スタイルにするなどの仕上げ
```

`docxlib.js` が見出し・表・図・用語集などの共通部品。本文は `mkdocx_*.js` に
`**太字**` `` `等幅` `` `^{上付き}` `_{下付き}` の簡易記法で書く。

`assets/{A,B,C}_{tz,cp}.svg` は RDKit で描いたリンカーの構造式。`build.py` が
XML 宣言と固定サイズを除去し、線の色を `currentColor` に置換してテーマ追従にする。

数値の出どころは `../iedda_template_sim_v18.py` の Colab 実行結果。
