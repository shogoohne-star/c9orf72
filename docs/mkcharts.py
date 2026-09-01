# -*- coding: utf-8 -*-
"""ドキュメント用のインライン SVG グラフを生成する。
色は CSS 変数を参照するのでライト/ダーク両テーマで成立する。"""
import math, os, html

W_PAD_L, W_PAD_R, PAD_T, PAD_B = 118, 18, 26, 46

def esc(s): return html.escape(str(s), quote=True)

def frame(w, h, title=None):
    return [f'<svg viewBox="0 0 {w} {h}" role="img" class="chart" preserveAspectRatio="xMidYMid meet">']

def grouped_bars(name, cats, series, ymax, unit="", h=None, colors=("s1","s2","s3"),
                 fmt=lambda v: f"{v:.1f}", ticks=5, note=None):
    """cats: [ラベル], series: [(系列名, [値])]"""
    n_c, n_s = len(cats), len(series)
    rowh = 30 if n_s == 1 else 22 * n_s + 16
    h = h or PAD_T + n_c * rowh + PAD_B
    w = 660
    plot_w = w - W_PAD_L - W_PAD_R
    o = frame(w, h)
    # 目盛
    for i in range(ticks + 1):
        v = ymax * i / ticks
        x = W_PAD_L + plot_w * i / ticks
        o.append(f'<line x1="{x:.1f}" y1="{PAD_T-6}" x2="{x:.1f}" y2="{h-PAD_B+4}" class="grid"/>')
        o.append(f'<text x="{x:.1f}" y="{h-PAD_B+20}" class="tick" text-anchor="middle">{fmt(v)}</text>')
    y = PAD_T
    for ci, c in enumerate(cats):
        o.append(f'<text x="{W_PAD_L-10}" y="{y+rowh/2+4}" class="catlab" text-anchor="end">{esc(c)}</text>')
        bh = 15 if n_s > 1 else 17
        by = y + (rowh - (bh * n_s + 4 * (n_s - 1))) / 2
        for si, (sname, vals) in enumerate(series):
            v = vals[ci]
            bw = max(plot_w * (v / ymax), 0.8) if v > 0 else 0
            yy = by + si * (bh + 4)
            o.append(f'<rect x="{W_PAD_L}" y="{yy:.1f}" width="{bw:.1f}" height="{bh}" rx="2" class="bar {colors[si]}"/>')
            lab = fmt(v) + unit if v > 0 else "0"
            inside = bw > 62
            lx = W_PAD_L + bw - 7 if inside else W_PAD_L + bw + 7
            o.append(f'<text x="{lx:.1f}" y="{yy+bh-3.5:.1f}" class="vlab{" inv" if inside else ""}" '
                     f'text-anchor="{"end" if inside else "start"}">{esc(lab)}</text>')
        y += rowh
    if note:
        o.append(f'<text x="{W_PAD_L}" y="{h-6}" class="note">{esc(note)}</text>')
    o.append('</svg>')
    return "\n".join(o)

def log_bars(name, cats, series, lo, hi, colors=("s1","s2","s3"), h=None, note=None,
             fmtv=lambda v: f"{v:.2e}"):
    n_c, n_s = len(cats), len(series)
    rowh = 22 * n_s + 18
    h = h or PAD_T + n_c * rowh + PAD_B + (14 if note else 0)
    w = 660
    plot_w = w - W_PAD_L - W_PAD_R
    L0, L1 = math.log10(lo), math.log10(hi)
    def px(v):
        if v <= 0: return 0.0
        return plot_w * (math.log10(v) - L0) / (L1 - L0)
    o = frame(w, h)
    d = int(round(L1 - L0))
    for i in range(d + 1):
        e = int(L0) + i
        x = W_PAD_L + plot_w * i / d
        o.append(f'<line x1="{x:.1f}" y1="{PAD_T-6}" x2="{x:.1f}" y2="{h-PAD_B+4}" class="grid"/>')
        o.append(f'<text x="{x:.1f}" y="{h-PAD_B+20}" class="tick" text-anchor="middle">10<tspan dy="-5" font-size="9">{e}</tspan></text>')
    y = PAD_T
    for ci, c in enumerate(cats):
        o.append(f'<text x="{W_PAD_L-10}" y="{y+rowh/2+4}" class="catlab" text-anchor="end">{esc(c)}</text>')
        bh = 15
        by = y + (rowh - (bh * n_s + 4 * (n_s - 1))) / 2
        for si, (sname, vals) in enumerate(series):
            v = vals[ci]
            yy = by + si * (bh + 4)
            if v <= 0:
                o.append(f'<text x="{W_PAD_L+3}" y="{yy+bh-3.5:.1f}" class="vlab zero" text-anchor="start">到達なし</text>')
                continue
            bw = max(px(v), 1.0)
            o.append(f'<rect x="{W_PAD_L}" y="{yy:.1f}" width="{bw:.1f}" height="{bh}" rx="2" class="bar {colors[si]}"/>')
            inside = bw > 76
            lx = W_PAD_L + bw - 7 if inside else W_PAD_L + bw + 7
            o.append(f'<text x="{lx:.1f}" y="{yy+bh-3.5:.1f}" class="vlab{" inv" if inside else ""}" '
                     f'text-anchor="{"end" if inside else "start"}">{esc(fmtv(v))}</text>')
        y += rowh
    if note:
        o.append(f'<text x="{W_PAD_L}" y="{h-6}" class="note">{esc(note)}</text>')
    o.append('</svg>')
    return "\n".join(o)

def legend(items):
    sp = []
    for cls, lab in items:
        sp.append(f'<span class="lg"><i class="sw {cls}"></i>{esc(lab)}</span>')
    return '<div class="legend">' + "".join(sp) + '</div>'

OUT = {}

# --- A: 立体的に許される配置の割合 -----------------------------------------
OUT["chart_allowed"] = grouped_bars(
    "allowed",
    ["末端リン酸 (gap 0)", "塩基 C7 (2,5)", "ホスホトリエステル"],
    [("テトラジン側", [0.75, 0.60, 23.48]), ("シクロプロペン側", [0.91, 7.27, 29.69])],
    ymax=32, unit=" %", fmt=lambda v: f"{v:.0f}",
    note="RNA との立体反発が 4 kcal/mol 以下に収まる配置の割合（24 本の中央値）")
OUT["legend_allowed"] = legend([("s1","テトラジン側アーム"), ("s2","シクロプロペン側アーム")])

# --- B: 必要なひずみ ---------------------------------------------------------
OUT["chart_strain"] = grouped_bars(
    "strain",
    ["末端リン酸 (gap 0)", "塩基 C7 (2,5)", "ホスホトリエステル"],
    [("ひずみ", [15.7, 4.9, 2.5])],
    ymax=18, unit=" kcal/mol", fmt=lambda v: f"{v:.0f}", colors=("s3",),
    note="小さいほど良い。反応できる配置に届くまでに払うエネルギー（24 本の中央値）")

# --- C: EM 3 指標 ------------------------------------------------------------
OUT["chart_em"] = log_bars(
    "em",
    ["末端リン酸 (gap 0)", "塩基 C7 (2,5)", "ホスホトリエステル"],
    [("届きやすさ", [0.0, 3.88e-2, 4.93e-2]),
     ("中間",       [6.96e-9, 1.05e-2, 4.85e-2]),
     ("ひずみ込み", [7.21e-10, 8.08e-4, 2.58e-2])],
    lo=1e-10, hi=1e-1,
    note="実効モル濃度 EM (M)。3 本が揃うほど推定が信用できる")
OUT["legend_em"] = legend([("s1","届きやすさ重視 (EM_steric)"), ("s2","中間 (EM_unif)"), ("s3","ひずみ込み (EM_boltz)")])

# --- D: 到達できた鋳型の割合 -------------------------------------------------
OUT["chart_hit"] = grouped_bars(
    "hit",
    ["末端リン酸 (gap 0)", "塩基 C7 (2,5)", "ホスホトリエステル"],
    [("到達率", [0, 46, 100])],
    ymax=100, unit=" %", fmt=lambda v: f"{v:.0f}", colors=("s1",),
    note="24 本の二重鎖のうち、無理なく反応配座に届いたものの割合")

# --- E: ON/OFF 半減期 --------------------------------------------------------
OUT["chart_onoff"] = log_bars(
    "onoff",
    ["0.1 µM", "1 µM", "10 µM"],
    [("鋳型あり", [67.2, 43.6, 41.6]), ("鋳型なし", [1.066e7, 1.066e6, 1.066e5])],
    lo=1e1, hi=1e8,
    fmtv=lambda v: (f"{v:.0f} 秒" if v < 120 else (f"{v/3600:.0f} 時間" if v < 2*86400 else f"{v/86400:.0f} 日")),
    note="ホスホトリエステル型、k2 = 0.65 M⁻¹s⁻¹。反応が半分進むまでの時間")
OUT["legend_onoff"] = legend([("s1","鋳型あり（隣接 2 本が結合）"), ("s2","鋳型なし（溶液中の 2 分子衝突）")])

# --- F: C7 修飾位置の走査 ----------------------------------------------------
OUT["chart_pos"] = grouped_bars(
    "pos",
    ["(2,5) Δ=3", "(2,4) Δ=4", "(3,5) Δ=4", "(2,3) Δ=5", "(3,4) Δ=5", "(4,5) Δ=5"],
    [("到達率", [50, 69, 81, 75, 62, 62])],
    ymax=100, unit=" %", fmt=lambda v: f"{v:.0f}", colors=("s1",),
    note="塩基 C7 型。16 本、傍観アームなし。Δ = 接合部をまたぐ 2 ハンドルの鎖内間隔")

# --- G: 末端型のギャップ走査 -------------------------------------------------
OUT["chart_gap"] = grouped_bars(
    "gap",
    ["gap 0", "gap 1", "gap 2", "gap 3", "gap 4"],
    [("到達率", [0, 100, 100, 100, 100])],
    ymax=100, unit=" %", fmt=lambda v: f"{v:.0f}", colors=("s2",),
    note="末端リン酸型。12 本。gap = ASO 2 本の間に残る未占有の鋳型塩基数")

# --- H: 頑健性（不確かさの幅、倍） -------------------------------------------
OUT["chart_robust"] = log_bars(
    "robust",
    ["3 つの EM 指標の開き", "鋳型間の 95% 信頼区間", "反応判定を厳↔緩", "経験パラメータを振る"],
    [("塩基 C7 (2,5)", [48, 20400, 104, 1e6]), ("ホスホトリエステル", [1.9, 3.2, 1.6, 1.3])],
    lo=1, hi=1e6,
    fmtv=lambda v: ("到達 0 が出る" if v >= 1e6 else f"{v:,.1f} 倍" if v < 10 else f"{v:,.0f} 倍"),
    note="小さいほど良い。答えがどれだけブレるか")
OUT["legend_robust"] = legend([("s1","塩基 C7 (2,5)"), ("s2","ホスホトリエステル")])

os.makedirs("assets", exist_ok=True)
for k, v in OUT.items():
    open(f"assets/{k}.svg" if k.startswith("chart") else f"assets/{k}.html", "w").write(v)
    print(f"{k}: {len(v)} bytes")
