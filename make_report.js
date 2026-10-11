const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  LevelFormat, PageBreak, convertInchesToTwip,
} = require('docx');

// A4 (11906 DXA) - 左右マージン 1134 x2 => 本文幅 9638
const W = 9638;
const MONO = "Consolas";
const JP = "Yu Gothic";

const P = (text, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: opts.line ?? 276 },
  alignment: opts.align,
  indent: opts.indent,
  children: [new TextRun({ text, font: opts.font ?? JP, size: opts.size ?? 20,
                           bold: opts.bold, italics: opts.italics, color: opts.color })],
});

const Rich = (runs, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: 276 },
  indent: opts.indent,
  children: runs.map(r => new TextRun({
    text: r.t, font: r.mono ? MONO : JP, size: r.mono ? 17 : (r.size ?? 20),
    bold: r.b, italics: r.i, color: r.c })),
});

const H1 = (t) => new Paragraph({ text: t, heading: HeadingLevel.HEADING_1,
  spacing: { before: 320, after: 160 } });
const H2 = (t) => new Paragraph({ text: t, heading: HeadingLevel.HEADING_2,
  spacing: { before: 260, after: 120 } });
const H3 = (t) => new Paragraph({ text: t, heading: HeadingLevel.HEADING_3,
  spacing: { before: 200, after: 100 } });

const Code = (lines) => lines.map((l, i) => new Paragraph({
  spacing: { after: i === lines.length - 1 ? 140 : 0, line: 240 },
  shading: { type: ShadingType.CLEAR, fill: "F4F4F4" },
  indent: { left: 170 },
  children: [new TextRun({ text: l || " ", font: MONO, size: 16 })],
}));

const Bullet = (t, lvl = 0) => new Paragraph({
  numbering: { reference: "bl", level: lvl },
  spacing: { after: 70, line: 270 },
  children: [new TextRun({ text: t, font: JP, size: 20 })],
});

let NUMI = 0;
const NewList = () => { NUMI += 1; };
const Num = (t) => new Paragraph({
  numbering: { reference: "nm", level: 0, instance: NUMI },
  spacing: { after: 70, line: 270 },
  children: [new TextRun({ text: t, font: JP, size: 20 })],
});

function table(headers, rows, widths, opts = {}) {
  const total = widths.reduce((a, b) => a + b, 0);
  const cw = widths.map(w => Math.round(w / total * W));
  const diff = W - cw.reduce((a, b) => a + b, 0);
  cw[cw.length - 1] += diff;

  const cell = (txt, i, head, hi) => new TableCell({
    width: { size: cw[i], type: WidthType.DXA },
    shading: { type: ShadingType.CLEAR,
               fill: head ? "2F5496" : (hi ? "FFF2CC" : "FFFFFF") },
    margins: { top: 40, bottom: 40, left: 55, right: 55 },
    children: [new Paragraph({
      spacing: { after: 0, line: 240 },
      alignment: i === 0 ? AlignmentType.LEFT : (opts.center ? AlignmentType.CENTER : AlignmentType.LEFT),
      children: [new TextRun({
        text: String(txt), bold: head, color: head ? "FFFFFF" : "000000",
        font: /^[ -~]*$/.test(String(txt)) && !head ? MONO : JP,
        size: head ? 16 : 16 })],
    })],
  });

  return new Table({
    columnWidths: cw,
    width: { size: W, type: WidthType.DXA },
    rows: [
      new TableRow({ tableHeader: true,
        children: headers.map((h, i) => cell(h, i, true, false)) }),
      ...rows.map(r => new TableRow({
        children: r.cells.map((c, i) => cell(c, i, false, r.hi)) })),
    ],
  });
}

const Spacer = () => new Paragraph({ spacing: { after: 100 }, children: [] });
const Rule = () => new Paragraph({
  spacing: { before: 120, after: 160 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "BBBBBB" } },
  children: [],
});

// =====================================================================
const body = [];

// ---------- 表紙 ----------
body.push(new Paragraph({
  spacing: { before: 1200, after: 200 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: "二価自壊性マスク鎖による mRNA 翻訳スイッチ",
                           font: JP, size: 34, bold: true })] }));
body.push(new Paragraph({
  spacing: { after: 700 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: "設計計算のまとめと、先行研究を加味した考察",
                           font: JP, size: 26, color: "444444" })] }));
body.push(P("2026年10月11日", { align: AlignmentType.CENTER, size: 20 }));
body.push(P("標的: 83 nt（5'UTR 1–53 / CDS 54–83）", { align: AlignmentType.CENTER, size: 20 }));
body.push(P("計算コード: zone_screen.py（ViennaRNA 2.7.2）", { align: AlignmentType.CENTER, size: 20 }));
body.push(Spacer());
body.push(new Paragraph({
  spacing: { before: 500, after: 120 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({ text: "— 本書の構成 —", font: JP, size: 20, color: "666666" })] }));
body.push(P("第1部　計算結果のまとめ（パラメータ、モデル、候補）", { align: AlignmentType.CENTER, size: 20 }));
body.push(P("第2部　先行研究を加味した考察（文献レビューと設計への含意）", { align: AlignmentType.CENTER, size: 20 }));
body.push(new Paragraph({ children: [new PageBreak()] }));

// =====================================================================
// 第1部
// =====================================================================
body.push(new Paragraph({
  spacing: { after: 240 },
  children: [new TextRun({ text: "第1部　計算結果のまとめ", font: JP, size: 30, bold: true, color: "2F5496" })] }));

body.push(H1("1. 設計と目的"));
body.push(P("5'UTR 上の2箇所を 2'-OMe マスク鎖で同時に塞ぎ（二価＝高親和性）、NQO1 による Q3PA（トリメチルロック型キノンプロピオン酸）の還元を引き金に自壊性コアが分解して一価化し、アビディティを失ってマスクが外れる。"));
body.push(P("無細胞翻訳系で cap近傍 / 中央 / Kozak近傍 の3ゾーンを比較し、最も阻害効率の高い位置を決めたうえで、その周辺を再スクリーニングする。"));

body.push(H2("1.1 カスケード機構（提示構造から読み取り）"));
body.push(NewList() || Num("NQO1 がキノンを還元 → ヒドロキノン"));
body.push(Num("トリメチルロックのラクトン化（速い）→ アニリン遊離"));
body.push(Num("ベンジルエーテルの 1,6-脱離 → オリゴ1（DBCO-トリアゾール側）放出"));
body.push(Num("生じたベンジルアルコールが ortho のアミドを攻撃しフタリド化 → オリゴ2 放出"));
body.push(Rich([{ t: "ON への切り替わりはステップ3で完了する", b: true },
                { t: "（ステップ4は後始末。本解析では無視）。" }]));

body.push(H2("1.2 標的配列と二次構造"));
body.push(...Code([
  "         1         11        21        31        41        51        61        71        81",
  "         |         |         |         |         |         |         |         |         |",
  "    AGAAUAAACU AGUAUUCUUC UGGUCCCCAC AGACUCAGAG AGAACCCGCC ACCAUGUUCG UGUUCCUGGU GCUGCUGCCU CUG",
  "    |------------------ 5'UTR (1-53) ------------------||--------- CDS (54-83) ---------|",
  "",
  "MFE 構造 (37 C, 83 nt 単独):",
  "   1 AGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACCAUGUUCG",
  "     ((((((......)))))).((((((......))).)))((((....((.(((((.(....",
  "  61 UGUUCCUGGUGCUGCUGCCUCUG",
  "     ....).))))).))....)))).",
  "",
  "MFE = -17.50 kcal/mol   アンサンブル dG = -19.40 kcal/mol",
]));

body.push(H1("2. パラメータ設定"));
body.push(Rich([{ t: "確度で4段階に分けて記載する。どれがどれか分からないまま数値だけ使うと誤用されるため。", b: true }]));

body.push(H2("2.1 確度の高い値（文献の標準パラメータ）"));
body.push(table(["パラメータ", "値", "出典・根拠"], [
  { cells: ["RNA 二次構造", "ViennaRNA 2.7.2 標準", "Turner 2004"] },
  { cells: ["最近接塩基対 ΔH, ΔS", "Xia et al. 1998 (RNA/RNA)", "duplex init ΔH 3.61 / ΔS −1.5、末端 A-U ペナルティ ΔH 3.72 / ΔS 10.5"] },
  { cells: ["温度", "37.0 °C", "翻訳反応温度"] },
  { cells: ["RT", "0.6163 kcal/mol", "—"] },
  { cells: ["A 型 rise / twist", "2.8 Å/bp / 32.7 °/bp", "11 bp/turn"] },
  { cells: ["リン酸バックボーン半径", "8.8 Å", "—"] },
], [20, 24, 56]));

body.push(H2("2.2 実験条件（実験に合わせて変更すること）"));
body.push(table(["パラメータ", "値", "備考"], [
  { cells: ["C_MRNA", "10 nM", "IVT 系で 5–50 nM 程度"] },
  { cells: ["MASK_EXCESS", "10 倍（= 100 nM）", "30倍まで可。gap 0–1 では上げすぎ注意（§7.1）"] },
  { cells: ["NA_CONC / MG_CONC", "0.100 M / 0.002 M", "→ Na+ 等価 0.248 M、Debye 長 6.1 Å"] },
  { cells: ["CLEAVAGE_FRAC", "1.0", "完全切断を仮定"] },
], [24, 24, 52]));

body.push(H2("2.3 仮定値（実測で置き換えるべき）"));
body.push(table(["パラメータ", "値", "不確実性と影響"], [
  { hi: true, cells: ["DDG_OME_PER_NT", "−0.14 kcal/mol/残基", "最大の不確実性。文献 ΔTm +0.5〜+1.5 °C/mod と幅。全予測の精度を支配"] },
  { cells: ["RIGID_SEGMENT_A", "14.1 Å", "m-フェニレン 4.9 + トリアゾール 2.2 + DBCO 7.0。gap の許容範囲に効く"] },
  { cells: ["SI_CONTOUR_A", "4.8 Å", "CH2-O 2.4 + アミド 2.4（柔軟部のみ）"] },
  { cells: ["TERM_PER_END_A", "9.0 Å", "ホスホジエステル 4 + C6 アミノ/アジド 5"] },
  { cells: ["DG_COAX", "−1.3 kcal/mol", "nick での coaxial stacking。gap 0–1 のみ"] },
  { cells: ["COAX_FRAC_BY_GAP", "{0: 1.00, 1: 0.35}", "1-nt ギャップでの減衰"] },
  { cells: ["COAX_LOSS_PER_SUBST", "0.5", "nick 末端の置換基1本あたりの coax 損失"] },
  { cells: ["ELEC_Q_EFF / ELEC_EPS", "0.5 e / 60", "全電荷 1.0 と Manning 凝縮 0.25 の中間。影響小（§7.4）"] },
  { hi: true, cells: ["eps（阻害効率）", "cap 0.90 / mid 0.60 / kozak 0.90", "完全な仮定。ON/OFF の絶対値を支配。これを測るのが本実験"] },
], [24, 22, 54]));

body.push(H2("2.4 未定量（モデルに入っていない）"));
body.push(table(["項目", "状態"], [
  { cells: ["接合部の立体障害", "STERIC_DG_BY_GAP = {}（= 0）。第一原理では出せない。gap 0–1 で約 1 kcal/mol を超えると推奨が gap 2 に動くことは確認済み"] },
  { cells: ["全長 mRNA の文脈", "83 nt 単独で折り畳んでいる。FULL_MRNA で差し替え可能"] },
  { cells: ["5' cap 構造", "未反映"] },
  { cells: ["速度論", "平衡のみ。切断直後の geminate recombination は未モデル"] },
], [26, 74]));

body.push(H2("2.5 探索条件とフィルタ"));
body.push(table(["項目", "値", "項目", "値"], [
  { cells: ["腕長", "6–13 nt（各腕独立）", "G 連続 / 同一塩基連続", "3 / 4 まで"] },
  { cells: ["gap", "0–45 nt（固定しない）", "マスク鎖の自己構造", "MFE > −1.5 kcal/mol"] },
  { cells: ["係留末端", "inner（M1-5' / M2-3'）", "鎖間二量体", "分子間 ΔG > −4.0 kcal/mol"] },
  { cells: ["最小被覆", "機能要素に 6 nt 以上", "オフターゲット余裕", "3.0 kcal/mol 以上"] },
], [20, 30, 24, 26]));

body.push(H2("2.6 ゾーン定義"));
body.push(table(["ゾーン", "region", "must_cover", "eps", "阻害機構（仮定）"], [
  { cells: ["cap", "1–28", "1–12", "0.90", "43S ロード / eIF4F 結合の阻害"] },
  { cells: ["mid", "16–45", "22–38", "0.60", "スキャニング阻害。43S のヘリカーゼ活性で解かれうる"] },
  { cells: ["kozak", "38–66", "48–56", "0.90", "Kozak + AUG。開始コドン認識の阻害"] },
], [12, 13, 15, 9, 51], { center: true }));

body.push(new Paragraph({ children: [new PageBreak()] }));

body.push(H1("3. 計算モデル"));
body.push(H2("3.1 自由エネルギー"));
body.push(...Code([
  "dG_mono,k = dG_hyb,k + dG_open,k - dG_self,k + dG_elec,k",
  "",
  "dG_biv    = dG_mono,1 + dG_mono,2 + ddG_coop - dG_cross",
  "            - RT*ln(C_eff) + dG_steric + dG_elec,biv",
  "",
  "ddG_coop  = (dG_open,both - dG_open,1 - dG_open,2) + dG_coax",
]));
body.push(Bullet("dG_hyb: Xia 1998 NN + 塩濃度補正 + 2'-OMe 補正"));
body.push(Bullet("dG_open: ViennaRNA hard constraint による制約付き分配関数差（RNAup 相当）"));
body.push(Bullet("ddG_coop: ΔG_open の非加算性。負なら「片方を開くともう片方も開く」= 二価が得"));
body.push(Bullet("C_eff: 半屈曲性鎖（Thirumalai–Ha 型）+ 剛直ロッドの球面平均畳み込み"));

body.push(H2("3.2 占有率（結合平衡を厳密に解く）"));
body.push(P("無細胞系は濃度が既知なので「マスク大過剰」近似を使わない。"));
body.push(Bullet("OFF 状態: T + M ⇌ TM の 1:1 平衡（二次方程式）"));
body.push(Bullet("ON 状態（切断後）: 2サイト4化学種（T, TA, TB, TAB）の連立平衡を交互求解。協同性 w = exp(−ΔΔG_coop/RT) を含む"));
body.push(Rich([{ t: "検証済み: ", b: true },
  { t: "質量保存（化学種の和 = 1.000000000）、1:1 極限との一致（差 2×10⁻¹³）、独立サイト極限との一致（差 8×10⁻⁶）、協同性因子の単調性。" }]));

body.push(H2("3.3 翻訳活性と評価指標"));
body.push(...Code([
  "活性     = SUM_species [分率] x (1 - 阻害効率)",
  "阻害効率 = eps_zone x min(1, 被覆長 / REQ_COVER_NT)",
  "ON/OFF   = 活性_ON / 活性_OFF",
]));
body.push(Rich([{ t: "eps は熱力学から予測できないため、ON/OFF は eps の関数として出力する。", b: true },
  { t: " eps_from_onoff() で測定値から eps を逆算できる。" }]));
body.push(Spacer());
body.push(table(["指標", "定義", "意味"], [
  { cells: ["occ_off", "二価状態での標的占有率", "マスク効率。高いほど良い"] },
  { cells: ["occ_on", "切断後の残存占有率", "リーク。低いほど良い"] },
  { hi: true, cells: ["S", "occ_off × (1 − occ_on)", "スイッチ品質。eps 非依存の熱力学指標"] },
  { cells: ["ON/OFF", "活性比", "eps の仮定に依存"] },
  { cells: ["ok_frac", "45条件での合格率", "頑健性（§6）"] },
], [14, 28, 58]));

body.push(new Paragraph({ children: [new PageBreak()] }));

body.push(H1("4. 計算結果 — ゾーン別候補"));

body.push(H2("4.1 cap（region 1–28、cap近傍 1–12 を被覆）"));
body.push(table(["標的位置", "gap", "MASK1", "MASK2", "Tm1/Tm2", "occOFF", "occON", "S", "ON/OFF", "頑健"], [
  { hi: true, cells: ["4-13 / 15-24", "1", "ACUAGUUUAU", "ACCAGAAGAA", "22.4/22.3", "0.9974", "0.0171", "0.9803", "9.66", "0.71"] },
  { cells: ["4-13 / 14-23", "0", "ACUAGUUUAU", "CCAGAAGAAU", "22.4/20.4", "0.9931", "0.0137", "0.9795", "9.31", "0.67"] },
  { hi: true, cells: ["3-12 / 15-24", "2", "CUAGUUUAUU", "ACCAGAAGAA", "18.6/22.3", "0.9861", "0.0072", "0.9790", "8.86", "0.64"] },
  { cells: ["2-11 / 15-24", "3", "UAGUUUAUUC", "ACCAGAAGAA", "21.2/22.3", "0.9876", "0.0112", "0.9765", "8.94", "0.60"] },
  { cells: ["4-13 / 16-24", "2", "ACUAGUUUAU", "ACCAGAAGA", "22.4/17.8", "0.9850", "0.0132", "0.9720", "8.71", "0.58"] },
  { cells: ["2-11 / 16-24", "4", "UAGUUUAUUC", "ACCAGAAGA", "21.2/17.8", "0.9269", "0.0078", "0.9197", "6.00", "0.40"] },
], [15, 5, 15, 15, 11, 9, 8, 8, 7, 7], { center: true }));
body.push(Rich([{ t: "gap 0–3 がほぼ等価（S 0.977–0.980）、gap 4 から劣化。3ゾーンで最も gap に寛容。MASK2 " },
  { t: "ACCAGAAGAA", mono: true }, { t: " が gap 1/2/3 で共通 → MASK1 を替えるだけで gap を振れる。" }]));

body.push(H2("4.2 mid（region 16–45、22–38 を被覆）"));
body.push(table(["標的位置", "gap", "MASK1", "MASK2", "Tm1/Tm2", "occOFF", "occON", "S", "ON/OFF", "頑健"], [
  { hi: true, cells: ["27-33 / 34-42", "0", "UCUGUGG", "CUCUCUGAG", "15.8/24.7", "0.9795", "0.0138", "0.9660", "2.41", "0.60"] },
  { cells: ["27-33 / 35-44", "1", "UCUGUGG", "UUCUCUCUGA", "15.8/25.8", "0.9838", "0.0199", "0.9642", "2.42", "0.53"] },
  { cells: ["28-34 / 35-44", "0", "GUCUGUG", "UUCUCUCUGA", "14.1/25.8", "0.9665", "0.0184", "0.9487", "2.36", "0.47"] },
  { cells: ["28-35 / 37-45", "1", "AGUCUGUG", "GUUCUCUCU", "22.1/19.0", "0.9530", "0.0134", "0.9402", "2.32", "0.51"] },
  { hi: true, cells: ["27-33 / 36-45", "2", "UCUGUGG", "GUUCUCUCUG", "15.8/29.3", "0.9941", "0.0720", "0.9225", "2.42", "0.40"] },
  { cells: ["27-34 / 38-45", "3", "GUCUGUGG", "GUUCUCUC", "29.8/11.9", "0.9211", "0.1603", "0.7734", "2.02", "0.00"] },
], [15, 5, 15, 15, 11, 9, 8, 8, 7, 7], { center: true }));
body.push(P("gap 0–2 に限られる。gap 3 で頑健性 0.00。ON/OFF が 2.4 止まりなのは eps = 0.6 の仮定による。"));

body.push(H2("4.3 kozak（region 38–66、Kozak+AUG 48–56 を被覆）"));
body.push(P("全候補とも AUG 3塩基を完全被覆。"));
body.push(table(["標的位置", "gap", "MASK1", "MASK2", "Tm1/Tm2", "occOFF", "occON", "S", "ON/OFF", "頑健"], [
  { cells: ["51-59 / 60-66", "0", "GAACAUGGU", "GGAACAC", "24.2/22.6", "0.9985", "0.0420", "0.9566", "9.73", "0.73"] },
  { hi: true, cells: ["42-48 / 50-57", "1", "CGGGUUC", "ACAUGGUG", "19.7/19.8", "0.9581", "0.0121", "0.9465", "7.22", "0.53"] },
  { hi: true, cells: ["42-48 / 49-56", "0", "CGGGUUC", "CAUGGUGG", "19.7/26.3", "0.9964", "0.0513", "0.9453", "9.29", "0.60"] },
  { cells: ["45-51 / 53-60", "1", "UGGCGGG", "CGAACAUG", "25.2/14.3", "0.9783", "0.0438", "0.9355", "8.14", "0.51"] },
  { hi: true, cells: ["38-46 / 49-56", "2", "GGUUCUCUC", "CAUGGUGG", "25.8/26.3", "0.9908", "0.0611", "0.9303", "8.88", "0.33"] },
  { cells: ["38-46 / 50-57", "3", "GGUUCUCUC", "ACAUGGUG", "25.8/19.8", "0.9122", "0.0230", "0.8912", "5.56", "0.27"] },
  { cells: ["48-55 / 60-66", "4", "AUGGUGGC", "GGAACAC", "27.4/22.6", "0.8965", "0.0923", "0.8138", "4.86", "0.00"] },
], [15, 5, 15, 15, 11, 9, 8, 8, 7, 7], { center: true }));
body.push(Rich([{ t: "gap 0–3。gap 4 で頑健性 0.00。3ゾーンで最も gap に敏感。" },
  { t: "なお最上段 51-59/60-66 は第2部 §9.3 で撤回する。", b: true, c: "C00000" }]));

body.push(new Paragraph({ children: [new PageBreak()] }));

body.push(H1("5. 試すべき位置と gap（計算のみに基づく推奨）"));
body.push(H2("5.1 第一選択 — gap 2 に揃える"));
body.push(P("立体障害が未定量であること、gap 0–1 では切断後の協同再結合リークが濃度依存で効くこと（§7.1）を踏まえ、gap 2 を標準とする。"));
body.push(table(["ゾーン", "MASK1（5'→コア）", "MASK2（3'→コア）", "標的位置", "gap", "スパン", "S", "ON/OFF"], [
  { cells: ["cap", "CUAGUUUAUU (10nt)", "ACCAGAAGAA (10nt)", "3-12 / 15-24", "2", "15.7 Å", "0.979", "8.86"] },
  { cells: ["mid", "UCUGUGG (7nt)", "GUUCUCUCUG (10nt)", "27-33 / 36-45", "2", "15.7 Å", "0.923", "2.43"] },
  { cells: ["kozak", "GGUUCUCUC (9nt)", "CAUGGUGG (8nt)", "38-46 / 49-56", "2", "15.7 Å", "0.930", "8.88"] },
], [10, 22, 22, 15, 7, 10, 7, 7], { center: true }));
body.push(P("全残基 2'-OMe。リンカーは MASK2 の 3' 末端と MASK1 の 5' 末端を繋ぐ（逆平行配向）。合成分子は 5'-[MASK2]-リンカー-コア-リンカー-[MASK1]-3'。"));

body.push(H2("5.2 第二選択 — S を最大化する（gap を揃えない）"));
body.push(table(["ゾーン", "MASK1", "MASK2", "標的位置", "gap", "S", "ON/OFF", "頑健性"], [
  { cells: ["cap", "ACUAGUUUAU", "ACCAGAAGAA", "4-13 / 15-24", "1", "0.980", "9.66", "0.71"] },
  { cells: ["mid", "UCUGUGG", "CUCUCUGAG", "27-33 / 34-42", "0", "0.966", "2.41", "0.60"] },
  { cells: ["kozak", "GAACAUGGU", "GGAACAC", "51-59 / 60-66", "0", "0.957", "9.73", "0.73"] },
], [10, 20, 20, 16, 7, 9, 9, 9], { center: true }));
body.push(Rich([{ t: "※ kozak の 51-59/60-66 は第2部 §9.3 で撤回する（MASK2 が CDS 内）。", b: true, c: "C00000" }]));

body.push(H2("5.3 gap の許容範囲"));
body.push(table(["ゾーン", "使える gap", "推奨", "破綻する gap"], [
  { cells: ["cap", "0–3", "1–2", "4 以上で劣化"] },
  { cells: ["mid", "0–2", "0–1", "3 で頑健性 0.00"] },
  { cells: ["kozak", "0–3", "0–2", "4 で頑健性 0.00"] },
], [20, 24, 20, 36], { center: true }));
body.push(P("gap 2 なら3ゾーンとも成立する。共通化するなら gap 2。"));

body.push(H2("5.4 リンカー"));
body.push(P("現構造のテザー contour は PEG なしで 36.9 Å（うち剛直部 14.1 Å）。gap 0–5 の全候補で PEG スペーサー追加は不要。"));
body.push(P("剛直な DBCO 側の腕（中心環から約 16 Å）がコアを二重鎖表面から離す。コアはアミド側オリゴの近くに座るため、NQO1 のアクセスを確保したい場合はアミド側オリゴの末端修飾を長くする（C6 → C12）のが有効。余裕が 21 Å あるので熱力学的コストはほぼゼロ。"));

body.push(H2("5.5 必要な対照"));
body.push(NewList() || Num("マスクなし mRNA — ON の上限、活性 1.0 の基準"));
body.push(Num("切断済み相当品 — MASK1 と MASK2 を別々の鎖として等モル添加。自壊性リンカーの切断反応を介さずに ON 状態の熱力学だけを検証でき、ON が出ない場合に「切断が進んでいない」のか「切断しても外れない」のかを切り分けられる"));
body.push(Num("スクランブル配列の二価構築体 — 非特異的な翻訳阻害の対照"));

body.push(H1("6. 頑健性"));
body.push(P("2'-OMe 補正（5水準）× マスク過剰率（5/10/30倍）× 末端修飾長（6/9/12 Å）の45条件で occ_OFF ≥ 0.95 かつ occ_ON ≤ 0.10 を満たす割合。"));
body.push(P("最良候補でも 0.60–0.73 に留まる。主因は 2'-OMe 補正を −0.05 まで振ると二価が弱くなりすぎること。この補正値の実測が最も効率的な精度改善。"));

body.push(new Paragraph({ children: [new PageBreak()] }));

body.push(H1("7. 設計上の注意点（計算で定量した5件）"));

body.push(H2("7.1 gap 0–1 ではマスク過剰率を上げられない"));
body.push(P("gap 0–1 では、切断後の2断片が隣接結合して coaxial stacking で協同的に再結合する三元複合体（TAB）が形成される。TAB は [mask]² に比例するため濃度で急増する。"));
body.push(table(["マスク過剰率", "×10", "×30", "×50", "×100"], [
  { cells: ["occ_ON", "0.018", "0.073", "0.144", "0.343"] },
  { cells: ["うち TAB 由来", "20%", "43%", "56%", "72%"] },
], [24, 19, 19, 19, 19], { center: true }));
body.push(P("gap 2 以上なら TAB 項が小さく、濃度を上げる自由度が残る。"));

body.push(H2("7.2 inner 係留では gap 0–1 の coax ボーナスは得られない"));
body.push(P("gap ≤ 1 では nick が Mask1-5' と Mask2-3' の間にできる。inner 係留ではその2末端とも置換基が付くため coaxial stacking はほぼ消失する。gap 0–1 を選ぶ積極的理由がない。"));

body.push(H2("7.3 鎖間二量体は繋がれた分子では桁違いに効く"));
body.push(P("2本が自由鎖なら ΔG −4.9 kcal/mol（Kd 360 µM）は 100 nM で無視できる。しかしリンカーで繋がれていると相手との実効濃度は約 40 mM になり、"));
body.push(...Code(["自己折りたたみ率 = 40 mM / (40 mM + 0.36 mM) = 約 99%"]));
body.push(P("構築体のほぼ全てが標的ではなく自分自身と対合する。CROSS_DIMER_MAX = −4.0 のフィルタはこのために必須。"));

body.push(H2("7.4 接合部の電荷反発はこの塩濃度では小さい"));
body.push(table(["成分", "大きさ (kcal/mol)", "Kd への影響"], [
  { cells: ["自分のマスク鎖との反発", "0.82", "なし（結合前後で不変）"] },
  { cells: ["標的鎖との反発（各腕）", "0.18", "一価・二価の両方に等しく効く"] },
  { cells: ["相手腕との反発（二価のみ）", "0.40 (inner) / 0.15 (3'-3')", "スイッチ窓を狭める"] },
], [30, 26, 44]));
body.push(P("Mg²⁺ 込みのイオン強度 0.248 M で Debye 長 6.1 Å と遮蔽が強いため小さい。低塩条件の無細胞系を使う場合は再計算が必要。"));

body.push(H2("7.5 ゾーン間の交絡要因（誘起構造変化）"));
body.push(P("マスクが結合すると残りの領域が再フォールドし、他の機能要素のアクセシビリティが変わる。"));
body.push(table(["ゾーン", "cap15 への影響", "Kozak/AUG への影響"], [
  { cells: ["cap", "—", "−1.17（AUG が開く）"] },
  { cells: ["mid", "+1.19（新構造で閉じる）", "−1.00（AUG が開く）"] },
  { hi: true, cells: ["kozak", "+0.35（影響小）", "—"] },
], [20, 40, 40], { center: true }));
body.push(P("cap と mid は自分が覆っていない領域の構造を約 1 kcal/mol 動かす。測定された ON/OFF を eps として読む際は、この副作用が混入していることに注意。kozak が最も解釈がクリーン。"));

body.push(H1("8. 既知の限界"));
body.push(NewList() || Num("83 nt 単独で折り畳んでいる（全長 mRNA の文脈は未反映）。ΔΔG_coop はこの設計の駆動力であり、全長では変わりうる"));
body.push(Num("10-mer 断片ベースの duplex ΔG は、全長 RNA 上での二本鎖形成を直接解くものではない"));
body.push(Num("ViennaRNA は 2'-OMe をモデル化していない（ΔΔG で補正）"));
body.push(Num("平衡のみ。切断直後の geminate recombination、リボソームによる置換は未モデル"));
body.push(Num("eps（阻害効率）は仮定値。ON/OFF の絶対値は予測ではない"));
body.push(Num("接合部の立体障害は未定量"));

body.push(new Paragraph({ children: [new PageBreak()] }));

// =====================================================================
// 第2部
// =====================================================================
body.push(new Paragraph({
  spacing: { after: 240 },
  children: [new TextRun({ text: "第2部　先行研究を加味した考察", font: JP, size: 30, bold: true, color: "2F5496" })] }));

body.push(P("第1部の計算は熱力学的なスイッチ品質（occ_OFF と occ_ON）を与えるが、「マスクされた mRNA の翻訳がどれだけ落ちるか」（阻害効率 eps）は熱力学の外にある。本部ではその eps について先行研究から読み取れることを整理し、設計への含意を導く。"));

body.push(H1("9. 位置依存性に関する先行研究"));

body.push(H2("9.1 キャップ近傍の構造は 43S ロードを阻害する"));
body.push(Rich([{ t: "Kozak (1989)", b: true },
  { t: "　−30 kcal/mol のヘアピンは、キャップから 12 塩基の位置では開始を阻害するが、52 塩基に離すと阻害しない。細胞内実験（Kozak 1986）では −30 kcal/mol は無効、−50 kcal/mol で 85–95% の低下。つまり安定性だけでなく位置が効く。" }]));
body.push(Rich([{ t: "Babendure, Babendure, Ding & Tsien (2006) RNA 12:851–861", b: true },
  { t: "　二色蛍光による単一細胞アッセイ。ヘアピン安定性 −10〜−50 kcal/mol、キャップ–ヘアピン距離 1–46 塩基を系統的に走査した。" }]));
body.push(Bullet("翻訳効率は ΔG が −25 → −35 kcal/mol に変わる狭い窓で急激に低下する"));
body.push(Bullet("ヘアピンをキャップに対してわずか 9 塩基ずらすだけで翻訳が 50 倍以上変動する"));
body.push(Bullet("安定性と距離を固定しても、GC 含量が独立に効く"));
body.push(Rich([{ t: "Baker et al. (1997) J Biol Chem 272:11994–12000", b: true },
  { t: "　2'-O-(2-methoxyethyl) 修飾オリゴを ICAM-1 のキャップ領域（nt −20）に標的化。RNase H 非依存の機構で IC50 2–6 nM。同じ標的の RNase H 依存型化合物より「exceptionally effective」と報告されている。ステリックブロックがキャップ近傍で実際に強く効くことの直接的な実例。" }]));

body.push(H2("9.2 中央部はスキャニングに融かされる"));
body.push(Rich([{ t: "Wang, Shin, Alvarado, Kim, Bohlen, Dever & Puglisi (2022) Cell", b: true },
  { t: "　単分子計測で 43S 走査をリアルタイム追跡。43S は約 100 nt/s で一方向に走査し、" },
  { t: "5'UTR 中の −52.1 kcal/mol までのヘアピンを初回遭遇で融かして開始コドンに到達する", b: true },
  { t: "。これは「強い構造＝絶対的な障壁」という見方を否定する。" }]));
body.push(P("本設計の二価マスクは 20 bp 相当の二重鎖で、ヘアピン換算では約 −26 kcal/mol。Wang らが融かせると示した −52 kcal/mol の半分でしかなく、スキャニングに対する裕度がない。"));

body.push(H2("9.3 開始コドン近傍は別機構で効く"));
body.push(P("同じ Wang (2022) は、ヘアピンが開始コドンに隣接している場合は挙動が異なることも示している。43S PIC が開始部位に留まる時間が約 3 秒と短くなり、PIC の 3'→5' 方向への後退が誘起されて再走査が必要になる。つまり AUG 近傍の構造は「融かされにくい」のではなく「開始複合体の安定化を妨げる」ことで効く。"));
body.push(P("EMCV を用いた古典的研究でも、5' 非翻訳領域と AUG コドンに対するオリゴは RNase H なしで有意な翻訳阻害を示し、RNase H 処理を加えても阻害はほとんど増えなかった（= ステリックブロックで完結している）。"));

body.push(H2("9.4 CDS 内のステリックブロックは伸長リボソームに剥がされる"));
body.push(Rich([{ t: "これが本設計に直接影響する最も重要な知見である。", b: true }]));
body.push(P("Nucleic Acids Research (1991) および上記 EMCV 研究は、コーディング領域を標的とするオリゴは RNase H によって mRNA が切断されない限り翻訳を阻害しないと報告している。RNase H 活性を持たない修飾オリゴをコーディング領域に標的化しても、タンパク質合成は阻害されない。伸長中の 80S リボソームがオリゴを mRNA から排除するためと説明されている。"));

body.push(H1("10. 本設計への含意"));

body.push(H2("10.1 eps の仮定（cap 0.90 / mid 0.60 / kozak 0.90）は文献と整合する"));
body.push(P("ただし3ゾーンの阻害機構はそれぞれ異なる。"));
body.push(table(["ゾーン", "機構", "支持する文献"], [
  { cells: ["cap", "43S / eIF4F のロードを阻害", "Kozak 1989、Babendure 2006、Baker 1997"] },
  { cells: ["mid", "スキャニングの物理的阻害", "Wang 2022（否定的）"] },
  { cells: ["kozak", "開始複合体の安定化を妨げ、PIC を後退させる", "Wang 2022、EMCV 研究"] },
], [12, 36, 52]));

body.push(H2("10.2 mid は仮定（0.60）より悪い可能性が高い"));
body.push(P("Wang 2022 の「−52 kcal/mol まで融かす」という数字に対し、本設計の二価マスクは −26 kcal/mol 相当にすぎない。ヘアピンと分子間二重鎖は厳密には異なり（ヘアピンは融けても即座に再形成するのに対し、剥がされたオリゴは溶液に戻るため再結合に時間がかかる）、2'-OMe/RNA 二重鎖は A 型で剥離に抵抗する可能性もあるが、裕度がないことは確かである。"));
body.push(Rich([{ t: "ただしこれは mid を落とす理由にはならない。", b: true },
  { t: " むしろ mid は「スキャニング融解モデルが本系でも成り立つか」を確かめる陰性対照として価値がある。3ゾーン比較で mid だけが効かなければ、Wang らのモデルが本系でも成立していることの確認になり、以降の設計方針が固まる。" }]));

body.push(H2("10.3 【設計修正】kozak の S 最大候補 51-59/60-66 は撤回する"));
body.push(Rich([{ t: "第1部 §4.3 および §5.2 に挙げた kozak 候補 51-59 / 60-66 を撤回する。", b: true, c: "C00000" }]));
body.push(P("この候補は MASK2 が 60–66、すなわち完全に CDS 内にある。§9.4 の知見によれば、CDS に結合したステリックブロックは伸長中のリボソームに剥がされるため、MASK2 は阻害にほとんど寄与しない。"));
body.push(Rich([{ t: "さらに悪いことに、暴走リークの経路になる。", b: true },
  { t: " 漏れ初期化が1回起きると、伸長するリボソームが MASK2 を剥がす。すると構築体は切断前でも実質的に一価となり、MASK1 のアビディティも失われて外れる。1回の漏れが OFF 状態の崩壊を引き起こす正のフィードバックが成立してしまう。" }]));
body.push(P("代替として、AUG より下流の CDS 侵入を最小にした候補を計算した。"));
body.push(table(["AUG 下流の CDS 侵入", "標的位置", "gap", "MASK1", "MASK2", "S", "ON/OFF"], [
  { hi: true, cells: ["0 nt", "42-48 / 49-56", "0", "CGGGUUC", "CAUGGUGG", "0.9453", "9.29"] },
  { hi: true, cells: ["0 nt", "38-46 / 49-56", "2", "GGUUCUCUC", "CAUGGUGG", "0.9303", "8.88"] },
  { cells: ["1 nt", "42-48 / 50-57", "1", "CGGGUUC", "ACAUGGUG", "0.9465", "7.22"] },
  { cells: ["4 nt", "45-51 / 53-60", "1", "UGGCGGG", "CGAACAUG", "0.9355", "8.14"] },
  { cells: ["10 nt（撤回）", "51-59 / 60-66", "0", "GAACAUGGU", "GGAACAC", "0.9566", "9.73"] },
], [20, 18, 7, 16, 16, 11, 12], { center: true }));
body.push(Rich([{ t: "第1部 §5.1 の gap 2 推奨 38-46 / 49-56 は、MASK2 が AUG の最終塩基 56 でちょうど終わるため CDS 侵入ゼロであり、修正は不要である。", b: true }]));

body.push(H2("10.4 cap 候補の位置は文献の窓に合っている"));
body.push(P("Kozak・Babendure の知見は「キャップから約 12 塩基以内の構造が 43S ロードを阻害する」ことを示す。本設計の cap 候補は MASK1 が 3–12 または 4–13、すなわち 5' 末端から 2–3 塩基の位置から始まり、この窓の中に完全に収まっている。"));
body.push(P("一方 MASK2（14–24 / 15–24）はこの窓の外にある。したがって cap ゾーンの阻害は実質的に MASK1 が担い、MASK2 はアビディティを稼ぐ役割に専念していると解釈すべきである。これは設計として合理的で、切断後に外れるべき腕（MASK1、Kd 8–32 µM）が阻害を担っている点も整合的である。"));
body.push(P("なお 5' 末端ちょうどから始める候補（1-11 / 15-24）は熱力学的に劣る（S 0.862）。2–4 塩基目から始めるのが最適で、これは文献の要求とも両立する。"));

body.push(H2("10.5 Babendure の「9 塩基で 50 倍」は位置の微調整に大きな価値があることを意味する"));
body.push(P("Babendure らは、ヘアピンをキャップに対して 9 塩基ずらすだけで翻訳効率が 50 倍以上変動することを示した。本設計の cap 候補は 1–3 塩基違いで並んでおり（2-11 / 3-12 / 4-13）、熱力学的には S 0.977–0.980 とほぼ同等だが、翻訳阻害効率では大きく異なる可能性がある。"));
body.push(Rich([{ t: "実験計画への提言: ", b: true },
  { t: "第一段で1ゾーン1候補しか作らないと、位置の最適値を外すリスクが高い。少なくとも cap ゾーンについては、3-12 / 15-24 と 4-13 / 15-24 の2点（MASK2 は共通なので MASK1 を2本合成するだけ）を同時に試すことを推奨する。合成コストの増分に対して得られる情報が大きい。" }]));

body.push(H2("10.6 二価→一価の遷移が Babendure の急峻域に重なる"));
body.push(P("Babendure らは翻訳効率が ΔG −25〜−35 kcal/mol の狭い窓で急激に低下することを見出した。本設計では、"));
body.push(Bullet("二価状態 = 20 bp の二重鎖 ≈ −26 kcal/mol 相当"));
body.push(Bullet("一価状態 = 10 bp の二重鎖 ≈ −13 kcal/mol 相当"));
body.push(P("であり、設計のスイッチ遷移がちょうど応答曲線の最も急峻な部分をまたぐ。これはスイッチ振幅にとって有利な偶然である。逆に言えば、腕長をこれ以上短くすると二価状態でも閾値を下回り、マスクが効かなくなる恐れがある。腕長 6–7 nt の候補（mid の UCUGUGG など）はこの観点では余裕が小さい。"));

body.push(H2("10.7 2'-OMe の ΔTm は文献値の中〜下寄りに設定してある"));
body.push(P("2'-OMe による二重鎖安定化は一般に +1.3 °C/修飾残基とされ、配列・骨格文脈により 0.5〜3.0 °C の幅がある。本モデルの −0.14 kcal/mol/残基 は約 +1 °C/mod に相当し、文献範囲の中〜下寄りで保守的である。"));
body.push(P("ただし Freier & Altmann (1997) が指摘するように Tm は熱力学量ではなく、親和性の指標としては ΔG°37 を用いるべきである。本モデルは ΔG ベースで補正しているため方針としては正しい。それでも実測を推奨する立場は変わらない（頑健性 0.60–0.73 の主因がこのパラメータであるため）。"));

body.push(H2("10.8 5'UTR ステリックブロックには再現性の問題の前例がある"));
body.push(P("uORF を標的とするステリックブロック ASO が RNASEH1 の発現を再現性よく上昇させなかったという 2024 年の報告がある。同報告は、それ以前の複数の研究が 5'UTR 標的 ASO による発現上昇を uORF 非依存的に報告していたことにも言及している。"));
body.push(P("方向（抑制 vs 上昇）は本設計と逆だが、5'UTR を標的とするステリックブロックの効果が文献上必ずしも安定していないことは留意すべきである。第1部 §5.5 に挙げた3種の対照群は省略できない。"));

body.push(new Paragraph({ children: [new PageBreak()] }));

body.push(H1("11. 文献を加味した最終推奨"));

body.push(H2("11.1 ゾーンの優先順位を変更する"));
body.push(P("第1部では熱力学に優劣がないことから「解釈のクリーンさで kozak > cap > mid」としたが、先行研究を加味すると順位が変わる。"));
body.push(table(["順位", "ゾーン", "根拠", "懸念"], [
  { hi: true, cells: ["1", "cap", "3つの独立した文献ライン（Kozak の位置依存性、Babendure の定量、Baker の実例 IC50 2–6 nM）が支持。機構が最も確立している", "AUG を 1.17 kcal/mol 開く副作用（§7.5）"] },
  { cells: ["2", "kozak", "Wang 2022 の後退・再走査機構が支持。他要素への誘起構造変化が最小で解釈がクリーン", "CDS 侵入を避ける制約（§10.3）"] },
  { cells: ["3", "mid", "熱力学的には成立する", "Wang 2022 が否定的（−52 kcal/mol まで融かす）。陰性対照としての価値はある"] },
], [8, 12, 46, 34]));

body.push(H2("11.2 文献を加味した最終候補"));
body.push(table(["ゾーン", "MASK1（5'→コア）", "MASK2（3'→コア）", "標的位置", "gap", "S", "ON/OFF", "備考"], [
  { hi: true, cells: ["cap", "CUAGUUUAUU", "ACCAGAAGAA", "3-12 / 15-24", "2", "0.979", "8.86", "第一選択"] },
  { hi: true, cells: ["cap", "ACUAGUUUAU", "ACCAGAAGAA", "4-13 / 15-24", "1", "0.980", "9.66", "位置振り用。MASK2 共通"] },
  { hi: true, cells: ["kozak", "GGUUCUCUC", "CAUGGUGG", "38-46 / 49-56", "2", "0.930", "8.88", "CDS 侵入ゼロ"] },
  { cells: ["mid", "UCUGUGG", "GUUCUCUCUG", "27-33 / 36-45", "2", "0.923", "2.43", "陰性対照的な位置づけ"] },
], [9, 19, 19, 15, 6, 8, 8, 16], { center: true }));
body.push(P("全残基 2'-OMe。合成分子は 5'-[MASK2]-リンカー-コア-リンカー-[MASK1]-3'。gap 2 に揃えた場合の必要スパンはいずれも 15.7 Å（伸長時 18.0 Å）で、現構造のテザー（contour 36.9 Å）なら PEG 追加不要。"));
body.push(Rich([{ t: "cap を2点作ることを推奨する。", b: true },
  { t: " MASK2 が共通なので追加合成は MASK1 1本のみであり、§10.5 の「9 塩基で 50 倍」を踏まえると費用対効果が高い。" }]));

body.push(H2("11.3 予想される結果と、その解釈"));
body.push(table(["観測される結果", "示唆"], [
  { cells: ["cap > kozak ≫ mid", "文献の予想どおり。キャップ近傍の 43S ロード阻害が最も強い。以降は cap 周辺を細かく振る"] },
  { cells: ["kozak > cap", "本系では開始コドン認識の阻害が優位。Wang 2022 の後退機構が強く働いている可能性。kozak 周辺を振る"] },
  { cells: ["mid も同程度に効く", "スキャニング融解モデルが本系（無細胞系、2'-OMe 二重鎖）では当てはまらない。分子間二重鎖がヘアピンより剥がれにくいことの証拠になり、設計自由度が大きく広がる"] },
  { cells: ["3ゾーンとも ON/OFF < 2", "マスク自体は掛かっている（occ_OFF 0.92–0.98 は熱力学的に確か）ので、eps が小さい。腕長を伸ばして二価の安定性を上げるか、3価化を検討する"] },
  { cells: ["OFF が掛からない", "2'-OMe 補正の実測値が想定より小さい。§6 の頑健性解析のとおり、このパラメータが最大の不確実性"] },
], [26, 74]));

body.push(H1("12. 次にやるべきこと（優先順）"));
body.push(NewList() || Num("2'-OMe マスク鎖の Tm 実測（UV melting、10-mer 1本で可）→ DDG_OME_PER_NT を確定。現状で最大の不確実性であり、全予測の精度が一段上がる"));
body.push(Num("無細胞翻訳で3ゾーン（cap は2点）を比較し、eps_from_onoff() で各ゾーンの eps を確定"));
body.push(Num("勝ったゾーンの周辺で再スクリーニング。Babendure の 9 塩基 / 50 倍を踏まえ、1–2 塩基刻みで振る"));
body.push(Num("自壊性コアの実寸が分かれば RIGID_SEGMENT_A / SI_CONTOUR_A を更新"));
body.push(Num("立体障害の見積りが得られれば STERIC_DG_BY_GAP に入れて再計算（gap 0–1 で約 1 kcal/mol を超えると推奨が gap 2 に動く）"));

body.push(H1("13. 参照文献"));
const refs = [
  ["Kozak M. (1986, 1989)", "キャップからの距離とヘアピン安定性の位置依存性。−30 kcal/mol のヘアピンはキャップから 12 塩基では阻害、52 塩基では阻害しない。"],
  ["Babendure JR, Babendure JL, Ding JH, Tsien RY. (2006) Control of mammalian translation by mRNA structure near caps. RNA 12(5):851–861.", "二色蛍光単一細胞アッセイ。ΔG −25〜−35 kcal/mol で急激な低下。9 塩基のシフトで 50 倍の変動。"],
  ["Baker BF et al. (1997) J Biol Chem 272:11994–12000.", "2'-MOE オリゴを ICAM-1 のキャップ領域に標的化。非 RNase H 機構で IC50 2–6 nM。"],
  ["Wang J, Shin BS, Alvarado C, Kim JR, Bohlen J, Dever TE, Puglisi JD. (2022) Rapid 40S scanning and its regulation by mRNA structure during eukaryotic translation initiation. Cell.", "単分子計測。43S は約 100 nt/s で走査し −52.1 kcal/mol までのヘアピンを融かす。開始コドン隣接のヘアピンは PIC の後退を誘起する。"],
  ["Nucleic Acids Research 19(5):1113 (1991) および EMCV 研究 (PMID 2476310)", "コーディング領域のステリックブロックは RNase H なしでは翻訳を阻害しない。5'UTR と AUG コドン標的は RNase H なしで阻害する。"],
  ["Freier SM, Altmann KH. (1997) Nucleic Acids Res 25:4429.", "2'-修飾による二重鎖安定化。Tm は熱力学量ではなく ΔG°37 を用いるべきとの指摘。"],
  ["Molecular Therapy Nucleic Acids (2024)", "uORF 標的ステリックブロック ASO は RNASEH1 発現を再現性よく上昇させなかった。5'UTR 標的 ASO の再現性に関する注意喚起。"],
  ["Xia T et al. (1998) Biochemistry 37:14719.", "RNA/RNA 二重鎖の最近接塩基対パラメータ。本計算の ΔH/ΔS の出典。"],
];
body.push(table(["文献", "本書での用途"], refs.map(r => ({ cells: r })), [44, 56]));

body.push(Spacer());
body.push(Rule());
body.push(Rich([{ t: "注記: ", b: true },
  { t: "第2部の文献情報は公開データベースの検索により収集した。Kozak 1989 および Baker 1997 については一次文献の全文を確認しておらず、二次文献の記述に基づく。具体的な数値を論文として引用する際は原著を確認されたい。Babendure 2006 と Wang 2022 については要旨および本文の記述を確認している。", i: true }]));

// =====================================================================
const doc = new Document({
  creator: "mask design pipeline",
  title: "二価自壊性マスク鎖 設計まとめと考察",
  numbering: {
    config: [
      { reference: "bl", levels: [
        { level: 0, format: LevelFormat.BULLET, text: "・", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 400, hanging: 200 } } } },
        { level: 1, format: LevelFormat.BULLET, text: "–", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 760, hanging: 200 } } } },
      ]},
      { reference: "nm", levels: [
        { level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 400, hanging: 220 } } } },
      ]},
    ],
  },
  styles: {
    default: {
      document: { run: { font: JP, size: 20 } },
      heading1: { run: { font: JP, size: 26, bold: true, color: "2F5496" } },
      heading2: { run: { font: JP, size: 22, bold: true, color: "2F5496" } },
      heading3: { run: { font: JP, size: 20, bold: true, color: "444444" } },
    },
  },
  sections: [{
    properties: { page: {
      size: { width: 11906, height: 16838 },
      margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 },
    }},
    children: body,
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync("マスク鎖設計_計算まとめと考察.docx", b);
  console.log("written");
});
