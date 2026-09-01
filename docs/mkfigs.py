# -*- coding: utf-8 -*-
import math

def sine_path(x0, x1, cy, amp, lam, phase, step=4):
    pts = []
    x = x0
    while x <= x1:
        pts.append((x, cy + amp * math.sin(2 * math.pi * (x - x0) / lam + phase)))
        x += step
    return "M " + " L ".join(f"{px:.1f} {py:.1f}" for px, py in pts)

def sine_y(x, x0, cy, amp, lam, phase):
    return cy + amp * math.sin(2 * math.pi * (x - x0) / lam + phase)

# =========== 図 1: 3 つの接続位置 ===========================================
W, H = 720, 340
X0, X1, CY, AMP, LAM = 56, 664, 168, 46, 202
pA = sine_path(X0, X1, CY, -AMP, LAM, 0)          # ASO 鎖
pB = sine_path(X0, X1, CY, AMP, LAM, 0)           # 鋳型鎖
o = [f'<svg viewBox="0 0 {W} {H}" role="img" class="figure" preserveAspectRatio="xMidYMid meet">']
# 溝のハイライト（strand が離れている区間 = 溝）
for k in range(3):
    xs = X0 + LAM * (k + 0.25)
    xe = X0 + LAM * (k + 0.75)
    if xe > X1: xe = X1
    if xs < X1:
        o.append(f'<rect x="{xs:.0f}" y="{CY-AMP-4:.0f}" width="{xe-xs:.0f}" height="{2*AMP+8:.0f}" class="groove"/>')
# 塩基対のはしご
x = X0
while x <= X1:
    ya = sine_y(x, X0, CY, -AMP, LAM, 0)
    yb = sine_y(x, X0, CY, AMP, LAM, 0)
    o.append(f'<line x1="{x:.1f}" y1="{ya:.1f}" x2="{x:.1f}" y2="{yb:.1f}" class="rung"/>')
    x += 13.2
o.append(f'<path d="{pB}" class="bb tmpl"/>')
o.append(f'<path d="{pA}" class="bb aso"/>')
# ニック（ASO 鎖の切れ目）
XN = X0 + LAM * 1.5
yn = sine_y(XN, X0, CY, -AMP, LAM, 0)
o.append(f'<circle cx="{XN:.1f}" cy="{yn:.1f}" r="7" class="nickdot"/>')
o.append(f'<text x="{XN:.0f}" y="{yn-14:.0f}" class="figlab" text-anchor="middle">ニック</text>')

def callout(x, y, tx, ty, num, text, cls):
    o.append(f'<path d="M {x:.0f} {y:.0f} L {tx:.0f} {ty:.0f}" class="lead {cls}"/>')
    o.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="5.5" class="dot {cls}"/>')
    anc = "start" if tx > x else "end"
    o.append(f'<circle cx="{tx:.0f}" cy="{ty:.0f}" r="10" class="badge {cls}"/>')
    o.append(f'<text x="{tx:.0f}" y="{ty+4:.0f}" class="badgenum" text-anchor="middle">{num}</text>')
    dx = 16 if anc == "start" else -16
    o.append(f'<text x="{tx+dx:.0f}" y="{ty+4:.0f}" class="figlab" text-anchor="{anc}">{text}</text>')

# ① 塩基 C7 = 主溝の底（はしごの中点、内側）
xc = X0 + LAM * 0.5
callout(xc, CY, xc - 6, 44, "1", "塩基 C7：主溝の底（狭い）", "c-steel")
# ② 末端リン酸 = ニック
callout(XN, yn, XN + 10, 44, "2", "末端リン酸：ニックの隙間", "c-amber")
# ③ 骨格リン酸 = 外側
xp = X0 + LAM * 2.25
yp = sine_y(xp, X0, CY, -AMP, LAM, 0)
callout(xp, yp, xp + 14, H - 46, "3", "ホスホトリエステル：骨格の外側", "c-accent")
o.append(f'<text x="{X0}" y="{H-10}" class="fignote">ASO 鎖（上）と鋳型鎖（下）。網掛けが溝。</text>')
o.append('</svg>')
open("assets/fig_anchor.svg", "w").write("\n".join(o))
print("fig_anchor", sum(len(s) for s in o))

# =========== 図 2: 反応できる配置（NAC）の判定 ==============================
W, H = 660, 290
o = [f'<svg viewBox="0 0 {W} {H}" role="img" class="figure" preserveAspectRatio="xMidYMid meet">']
# テトラジン環（下）を楕円で、シクロプロペン（上）を短い線で
cx, cyT, cyC = 250, 200, 118
o.append(f'<ellipse cx="{cx}" cy="{cyT}" rx="74" ry="17" class="ringT"/>')
o.append(f'<text x="{cx}" y="{cyT+5}" class="figlab dim" text-anchor="middle">テトラジン環</text>')
o.append(f'<line x1="{cx-26}" y1="{cyC}" x2="{cx+26}" y2="{cyC}" class="alkene"/>')
o.append(f'<text x="{cx}" y="{cyC-14}" class="figlab dim" text-anchor="middle">シクロプロペン C=C</text>')
# 形成結合
for sx, ex in ((cx-64, cx-26), (cx+64, cx+26)):
    o.append(f'<line x1="{sx}" y1="{cyT-6}" x2="{ex}" y2="{cyC+6}" class="forming"/>')
o.append(f'<text x="{cx-96}" y="{(cyT+cyC)/2+4}" class="figlab" text-anchor="end">形成する 2 本の結合</text>')
o.append(f'<line x1="{cx-92}" y1="{(cyT+cyC)/2}" x2="{cx-62}" y2="{(cyT+cyC)/2}" class="lead c-accent"/>')
# 高さの矢印
o.append(f'<line x1="{cx+112}" y1="{cyC}" x2="{cx+112}" y2="{cyT}" class="dimline"/>')
o.append(f'<text x="{cx+120}" y="{(cyT+cyC)/2+4}" class="figlab" text-anchor="start">面からの高さ</text>')
# 判定条件
conds = ["形成する 2 本の結合が 2.2–3.8 Å",
         "2 つの環の面のなす角が 40° 以内",
         "C=C 軸と環の軸のなす角が 45° 以内",
         "アルケン中点の面内ずれが 1.6 Å 以内",
         "アルケンの 2 炭素が環面の同じ側"]
for i, c in enumerate(conds):
    y = 40 + i * 21
    o.append(f'<circle cx="{W-262}" cy="{y-4}" r="3" class="dot c-accent"/>')
    o.append(f'<text x="{W-252}" y="{y}" class="figlab" text-anchor="start">{c}</text>')
o.append(f'<text x="{W-262}" y="{H-14}" class="fignote">5 つすべてを満たす配置だけを「反応できる」と数える</text>')
o.append('</svg>')
open("assets/fig_nac.svg", "w").write("\n".join(o))
print("fig_nac ok")

# =========== 図 3: 実効モル濃度 EM の考え方 =================================
W, H = 660, 230
o = [f'<svg viewBox="0 0 {W} {H}" role="img" class="figure" preserveAspectRatio="xMidYMid meet">']
o.append(f'<text x="150" y="26" class="figlab hd" text-anchor="middle">鋳型につながれた 2 つ</text>')
o.append(f'<text x="500" y="26" class="figlab hd" text-anchor="middle">自由な溶液中の 2 つ</text>')
o.append(f'<rect x="34" y="42" width="232" height="120" rx="6" class="panel"/>')
o.append(f'<rect x="384" y="42" width="232" height="120" rx="6" class="panel"/>')
# 左: 二重鎖 + 2 本のアーム
o.append(f'<rect x="52" y="128" width="196" height="16" rx="4" class="duplex"/>')
o.append(f'<text x="150" y="157" class="fignote" text-anchor="middle">鋳型 RNA</text>')
o.append(f'<path d="M 96 128 C 88 96 118 84 132 92" class="arm"/>')
o.append(f'<path d="M 204 128 C 212 96 182 84 168 92" class="arm"/>')
o.append(f'<circle cx="134" cy="92" r="7" class="dot c-accent"/>')
o.append(f'<circle cx="166" cy="92" r="7" class="dot c-steel"/>')
# 右: 自由な 2 分子
for (px, py, cls) in ((430, 76, "c-accent"), (566, 132, "c-steel")):
    o.append(f'<circle cx="{px}" cy="{py}" r="7" class="dot {cls}"/>')
o.append(f'<path d="M 442 84 Q 500 100 556 126" class="lead dashed"/>')
o.append(f'<text x="500" y="152" class="fignote" text-anchor="middle">出会うには拡散が必要</text>')
o.append(f'<text x="330" y="106" class="figlab hd" text-anchor="middle">÷</text>')
o.append(f'<text x="{W/2}" y="196" class="figlab" text-anchor="middle">'
         f'EM = （鋳型上で反応できる形になる確率）÷（自由な溶液で同じ形になる確率／体積）</text>')
o.append(f'<text x="{W/2}" y="216" class="fignote" text-anchor="middle">'
         f'単位は M。「相手が実際に何 M の濃度でそばにいるのと同じか」を表す</text>')
o.append('</svg>')
open("assets/fig_em.svg", "w").write("\n".join(o))
print("fig_em ok")
