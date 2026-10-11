# -*- coding: utf-8 -*-
# =============================================================================
#  Zone-wise design of bivalent self-immolative mask strands for a cell-free
#  translation ON/OFF assay
#  ---------------------------------------------------------------------------
#  目的:
#    5'UTR を 3 つのゾーン (cap近傍 / 中央 / Kozak近傍) に分け、それぞれについて
#    無細胞翻訳系で試すべき MASK1 / MASK2 の配列・腕長・リンカー長を確定させる。
#
#  前スクリプト (mask_screen.py) からの変更点:
#    1. 濃度を実験条件そのものにした。無細胞系では [mRNA] と [mask] が既知なので
#       「マスク鎖大過剰」近似をやめ、結合平衡を厳密に解く。
#       切断後は 2 サイト 4 化学種 (T, TA, TB, TAB) の連立平衡を解く。
#    2. ΔG_open の非加算性 (協同性) を ddG_coop として独立に出力し、
#       平衡計算にも正しい結合定数として反映する。
#    3. Tm 窓はハードフィルタから外した。恣意的な閾値が順位を歪めたため、
#       ランキングは物理量ベース、Tm 窓は情報として併記する。
#    4. ロバストネス解析 (2'-OMe 補正 / 自壊性ユニット長 / mask:mRNA 比) を標準装備。
#    5. 腕長の探索範囲を 6-13 nt に拡大。GC 含量の高いゾーンでは 10 nt では
#       一価が強すぎるため、ゾーンごとに最適長が変わる。
#    6. 全長 mRNA 配列を任意で受け付け、局所折り畳みで accessibility を計算できる。
#
#  ON/OFF について:
#    占有率 (どれだけマスクされるか) は熱力学で予測できるが、
#    「マスクされた mRNA の翻訳がどれだけ落ちるか」= 阻害効率 eps は予測できない。
#    それが今回の実験で測る量そのものなので、ON/OFF は eps の関数として出力する。
#
#  Colab: このセルをそのまま貼って実行
# =============================================================================

try:
    import RNA
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "ViennaRNA"], check=True)
    import RNA

import math, os, itertools
import numpy as np
import pandas as pd

try:
    from scipy.optimize import brentq
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
try:
    import matplotlib.pyplot as plt
    HAVE_PLT = True
except ImportError:
    HAVE_PLT = False

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 120)


# =============================================================================
# 1. CONFIG
# =============================================================================
class CFG:
    # ---- 標的 ------------------------------------------------------------
    TARGET = ("AGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACC"
              "AUGUUCGUGUUCCUGGUGCUGCUGCCUCUG")
    UTR_RANGE = (1, 53)
    CDS_RANGE = (54, 83)

    # 全長 mRNA 配列 (任意)。与えると accessibility を全長の文脈で計算する。
    # 83 nt だけで折り畳むと実際と違う構造になりうるため、分かるなら入れること。
    FULL_MRNA = ""          # 例: "GGG...（5'UTR+CDS+3'UTR 全長）..."
    MAX_BP_SPAN = 120       # 全長使用時の局所折り畳み窓 [nt]

    # ---- ゾーン定義 ------------------------------------------------------
    #   region     : 腕を置いてよい範囲 (1-based inclusive)
    #   must_cover : この機能要素に腕が掛かることを要求する
    #   eps        : マスクされた時の翻訳阻害効率の「仮定値」。実験で測る対象。
    ZONES = [
        dict(name="cap",   region=(1, 28),  must_cover=(1, 12),
             eps=0.90, note="43S ロード / eIF4F 結合の阻害"),
        dict(name="mid",   region=(16, 45), must_cover=(22, 38),
             eps=0.60, note="スキャニング阻害。43S のヘリカーゼ活性で解かれうる"),
        dict(name="kozak", region=(38, 62), must_cover=(48, 56),
             eps=0.90, note="Kozak + AUG。開始コドン認識の阻害"),
    ]
    EPS_SWEEP = [0.5, 0.7, 0.9, 0.99]   # ON/OFF を eps の関数として出す
    # 機能要素にこの長さ以上掛かれば阻害は飽和するとみなす [nt]。
    # must_cover の窓幅で割る方式だと窓幅という恣意的な値が予測 ON/OFF を決めて
    # しまいゾーン間が比較できなくなるため、絶対長で頭打ちにする。
    REQ_COVER_NT = 6
    # マスクによる誘起構造変化を監視する機能要素 (name, start, end)。
    # 自分が覆っていない要素についてのみ計算する。
    WATCH_ELEMENTS = [("cap15", 1, 15), ("kozakAUG", 46, 56)]

    # ---- マスク鎖 --------------------------------------------------------
    ARM_LENGTHS = range(6, 14)      # 各腕の長さ探索範囲 [nt]
    ALL_2OME = True
    DDG_OME_PER_NT = -0.14          # kcal/mol / 修飾残基 (ΔTm ~ +1 degC/mod 相当)

    # ---- 無細胞翻訳系の条件 (ここが前回からの主な改善点) -------------------
    C_MRNA  = 10.0e-9       # [M] 反応中の mRNA 濃度 (IVT 系で 5-50 nM 程度)
    MASK_EXCESS = 10.0      # マスク鎖 / mRNA のモル比
    ANNEAL_TEMP_C = 37.0    # 評価温度 (翻訳反応温度に合わせる)
    CLEAVAGE_FRAC = 1.0     # 自壊性リンカーの切断率 (1.0 = 完全切断)

    # ---- 溶液条件 --------------------------------------------------------
    NA_CONC = 0.100         # 一価カチオン [M] (無細胞系の K+/ NH4+ 相当)
    MG_CONC = 0.0020        # Mg2+ [M]
    USE_MG  = True

    # ---- リンカー (PEG 系) -----------------------------------------------
    EO_CONTOUR_A   = 3.5
    EO_PERSIST_A   = 3.8
    SI_CONTOUR_A   = 12.0   # 自壊性ユニットの contour length [Å]。要実測/要見積り
    TERM_LINK_A    = 5.0
    EO_SCAN        = range(0, 61)
    CEFF_FRAC      = 0.50
    REQUIRE_EXTENDED_REACH = True   # ギャップ伸長時にも届く長さのみ推奨

    NICK_OFFSET_A  = 6.0
    RISE_SS_A      = 3.4
    RISE_SS_MAX_A  = 6.0

    # ---- 構造・配列フィルタ ----------------------------------------------
    # gap の下限は 0 にして物理に判定させる。gap<=COAX_GAP_MAX では
    # 切断後に 2 断片が隣接結合し coaxial stacking (DG_COAX) で協同再結合する
    # リークが生じるが、それは ddG_coop 経由で切断後の 2 サイト平衡にも
    # 正しく入っているので、手でフィルタせずスコアに淘汰させる。
    MIN_GAP_NT      = 0
    MAX_GAP_NT      = 45     # 上限は C_eff (アビディティ) が物理的に決めるので緩く取る
    MAX_G_RUN       = 3
    MAX_ANY_RUN     = 4
    SELF_FOLD_MAX   = -1.5
    CROSS_DIMER_MAX = -4.0
    OFFTARGET_MARGIN = 3.0
    DG_COAX         = -1.3
    COAX_GAP_MAX    = 1

    # ---- 目標値 (ランキングには使わず、情報として併記) ---------------------
    TM_MONO_TARGET = (18.0, 33.0)
    OCC_OFF_MIN    = 0.95    # OFF 状態で目標としたい占有率
    OCC_ON_MAX     = 0.10    # ON 状態 (切断後) で許容する残存占有率

    # ---- ロバストネス解析 -------------------------------------------------
    RUN_ROBUSTNESS   = True
    ROB_TOP_K        = 5             # 各ゾーンの上位何件を解析するか
    ROB_OME          = [-0.05, -0.10, -0.14, -0.20, -0.25]
    ROB_SI           = [6.0, 9.0, 12.0, 15.0]
    ROB_EXCESS       = [2.0, 5.0, 10.0, 50.0]

    TOP_N   = 12
    OUTDIR  = "/content" if os.path.isdir("/content") else "."


R_GAS = 1.9872e-3
T0K   = 273.15
TK    = CFG.ANNEAL_TEMP_C + T0K
RT    = R_GAS * TK
AVOGADRO = 6.02214076e23


# =============================================================================
# 2. 配列ユーティリティ
# =============================================================================
COMP = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"}

def to_rna(s):
    return "".join(c for c in s.upper().replace("T", "U") if c in "AUGC")

def revcomp(s):
    return "".join(COMP[c] for c in reversed(to_rna(s)))

def gc_frac(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0

def max_run(s):
    best = cur = 1
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best

def max_base_run(s, base):
    best = cur = 0
    for c in s:
        cur = cur + 1 if c == base else 0
        best = max(best, cur)
    return best

def overlap_len(a, b, c, d):
    return max(0, min(b, d) - max(a, c) + 1)


TARGET = to_rna(CFG.TARGET)
N_T = len(TARGET)
ZONE_BY_NAME = {z["name"]: z for z in CFG.ZONES}

# 折り畳みに使う配列と、標的座標 -> 折り畳み座標のオフセット
if CFG.FULL_MRNA:
    FOLD_SEQ = to_rna(CFG.FULL_MRNA)
    _off = FOLD_SEQ.find(TARGET)
    if _off < 0:
        raise ValueError("FULL_MRNA の中に TARGET 配列が見つかりません")
    OFFSET = _off
    USE_FULL = True
else:
    FOLD_SEQ, OFFSET, USE_FULL = TARGET, 0, False


# =============================================================================
# 3. 最近接塩基対熱力学 (Xia et al. 1998, RNA/RNA)
# =============================================================================
NN_XIA = {
    "AA": (-6.82, -19.0), "UU": (-6.82, -19.0),
    "AU": (-9.38, -26.7), "UA": (-7.69, -20.5),
    "CU": (-10.48, -27.1), "AG": (-10.48, -27.1),
    "CA": (-10.44, -26.9), "UG": (-10.44, -26.9),
    "GU": (-11.40, -29.5), "AC": (-11.40, -29.5),
    "GA": (-12.44, -32.5), "UC": (-12.44, -32.5),
    "CG": (-10.64, -26.7),
    "GG": (-13.39, -32.7), "CC": (-13.39, -32.7),
    "GC": (-14.88, -36.9),
}
DH_INIT, DS_INIT = 3.61, -1.5
DH_TERM_AU, DS_TERM_AU = 3.72, 10.5


def na_equivalent():
    na = CFG.NA_CONC
    if CFG.USE_MG and CFG.MG_CONC > 0:
        na += 3.3 * math.sqrt(CFG.MG_CONC)
    return min(max(na, 0.01), 1.2)


def nn_duplex(site):
    """完全相補二重鎖の ΔH [kcal/mol], ΔG37 [kcal/mol] を返す。"""
    s = to_rna(site)
    n = len(s)
    dH, dS = DH_INIT, DS_INIT
    for k in range(n - 1):
        h, ds = NN_XIA[s[k:k + 2]]
        dH += h
        dS += ds
    for end in (s[0], s[-1]):
        if end in "AU":
            dH += DH_TERM_AU
            dS += DS_TERM_AU
    dS += 0.368 * (n - 1) * math.log(na_equivalent())
    dG = dH - TK * dS / 1000.0
    if CFG.ALL_2OME:
        dG += n * CFG.DDG_OME_PER_NT
    return dH, dG


def tm_apparent(dH, dG, ct):
    """ΔH と ΔG(T) から見かけ Tm [degC]。ct = 全鎖濃度。"""
    dS = (dH - dG) / TK * 1000.0
    denom = dS + R_GAS * math.log(ct / 4.0) * 1000.0
    return dH * 1000.0 / denom - T0K if denom < 0 else float("nan")


def kd_from_dg(dG):
    return math.exp(dG / RT)


# =============================================================================
# 4. ViennaRNA ラッパ
# =============================================================================
RNA.cvar.temperature = CFG.ANNEAL_TEMP_C

def _md():
    md = RNA.md()
    md.temperature = CFG.ANNEAL_TEMP_C
    if USE_FULL:
        md.max_bp_span = CFG.MAX_BP_SPAN
    return md

_CACHE = {}

def ens_dg(seq):
    key = ("pf", seq)
    if key not in _CACHE:
        _CACHE[key] = RNA.fold_compound(seq, _md()).pf()[1]
    return _CACHE[key]

def mfe_fold(seq):
    return RNA.fold_compound(seq, _md()).mfe()

def dg_open(regions):
    """
    標的座標 (1-based inclusive) の領域を強制的に unpaired にするコスト。
    FULL_MRNA が与えられていれば全長の文脈で計算する。返り値 >= 0。
    """
    key = ("open", tuple(sorted(regions)))
    if key in _CACHE:
        return _CACHE[key]
    g_free = ens_dg(FOLD_SEQ)
    fc = RNA.fold_compound(FOLD_SEQ, _md())
    for (a, b) in regions:
        for k in range(a + OFFSET, b + OFFSET + 1):
            try:
                fc.hc_add_up(k)
            except Exception:
                fc.hc_add_up(k, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    val = fc.pf()[1] - g_free
    _CACHE[key] = val
    return val

def self_fold_ome(seq):
    ss, e = mfe_fold(seq)
    if CFG.ALL_2OME:
        e += 2 * ss.count("(") * CFG.DDG_OME_PER_NT
    return min(e, 0.0), ss

def intermolecular_dg(a, b):
    fc = RNA.fold_compound(a + "&" + b, _md())
    try:
        ss, e_ab = fc.mfe_dimer()
    except Exception:
        ss, e_ab = RNA.cofold(a + "&" + b)
    inter = e_ab - mfe_fold(a)[1] - mfe_fold(b)[1]
    if CFG.ALL_2OME:
        inter += 2 * ss.count("(") * CFG.DDG_OME_PER_NT
    return inter

_OT_CACHE = {}

def best_offtarget_dg(masker, exclude):
    """
    意図外サイトとの最良ハイブリ ΔG。
    FULL_MRNA が与えられていれば全長に対して走査する。無細胞系の反応液には
    全長 mRNA が入るので、83 nt だけを見ていたのでは off-target を見落とす。
    exclude は標的座標なので、折り畳み座標へ OFFSET を足して比較する。
    """
    key = (masker, exclude, CFG.DDG_OME_PER_NT)
    if key in _OT_CACHE:
        return _OT_CACHE[key]
    L = len(masker)
    ex_a, ex_b = exclude[0] + OFFSET, exclude[1] + OFFSET
    best = 0.0
    for i in range(0, len(FOLD_SEQ) - L + 1):
        a, b = i + 1, i + L
        if overlap_len(a, b, ex_a, ex_b) > 0:
            continue
        e = RNA.duplexfold(masker, FOLD_SEQ[i:i + L]).energy
        if CFG.ALL_2OME:
            e += L * CFG.DDG_OME_PER_NT
        best = min(best, e)
    _OT_CACHE[key] = best
    return best


# =============================================================================
# 5. リンカー / アビディティ (半屈曲性鎖モデル)
# =============================================================================
def required_span(gap_nt, rise):
    return CFG.NICK_OFFSET_A + gap_nt * rise

def _wlc_shape(x, L_over_lp):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    t = 1.0 - x * x
    expo = -0.75 * L_over_lp / t
    return 0.0 if expo < -700 else 4 * math.pi * x * x * t ** (-4.5) * math.exp(expo)

_WLCN = {}
def _wlc_norm(Llp, n=400):
    key = round(Llp, 4)
    if key not in _WLCN:
        h = 1.0 / n
        acc = sum((1 if i in (0, n) else (4 if i % 2 else 2)) * _wlc_shape(i * h, Llp)
                  for i in range(n + 1))
        _WLCN[key] = acc * h / 3.0
    return _WLCN[key]

def ceff_wlc(r_A, contour_A):
    """テザー末端間の実効濃度 [M]。理想 Gauss 鎖と違い r -> L_c で正しく 0 に落ちる。"""
    if contour_A <= 0 or r_A <= 0 or r_A >= contour_A:
        return 0.0
    Llp = contour_A / CFG.EO_PERSIST_A
    Z = _wlc_norm(Llp)
    if Z <= 0:
        return 0.0
    P_r = _wlc_shape(r_A / contour_A, Llp) / (Z * contour_A)
    return P_r / (4 * math.pi * r_A * r_A) * 1e27 / AVOGADRO

def optimize_linker(gap_nt):
    """
    C_eff(L_c) は極端に浅い最大しか持たないので「最大化」では無意味に長い
    リンカーが選ばれる。最大値の CEFF_FRAC に達する最短長を推奨値とし、
    さらにギャップ伸長時にも届くことを要求する (余裕のコストはほぼゼロ)。
    """
    r_req = required_span(gap_nt, CFG.RISE_SS_A)
    r_max = required_span(gap_nt, CFG.RISE_SS_MAX_A)
    scan = []
    for n in CFG.EO_SCAN:
        c = n * CFG.EO_CONTOUR_A + CFG.SI_CONTOUR_A + CFG.TERM_LINK_A
        if c > r_req:
            scan.append((n, c, ceff_wlc(r_req, c)))
    if not scan:
        return None
    ce_max = max(z[2] for z in scan)
    if ce_max <= 0:
        return None
    pool = scan
    if CFG.REQUIRE_EXTENDED_REACH:
        reach = [z for z in scan if z[1] >= r_max]
        if reach:
            pool = reach
    cand = [z for z in pool if z[2] >= CFG.CEFF_FRAC * ce_max] or \
           [max(pool, key=lambda z: z[2])]
    n_rec, c_rec, ce_rec = cand[0]
    return dict(n_eo=n_rec, contour=c_rec, ceff=ce_rec,
                r_req=r_req, r_req_max=r_max,
                reaches_extended=int(c_rec >= r_max),
                dG_link=-RT * math.log(ce_rec))


# =============================================================================
# 6. 結合平衡 (無細胞系の実濃度で厳密に解く)
# =============================================================================
def solve_1to1(Tt, Mt, Kd):
    """T + M <-> TM。[TM] を返す。"""
    b = Tt + Mt + Kd
    disc = max(b * b - 4.0 * Tt * Mt, 0.0)
    return (b - math.sqrt(disc)) / 2.0


def solve_two_site(Tt, At, Bt, K1, K2, w):
    """
    2 サイト 4 化学種の連立平衡:
        T + A <-> TA   (K1)
        T + B <-> TB   (K2)
        TA + B <-> TAB (K2/w),  TB + A <-> TAB (K1/w)
      w = exp(-ddG_coup/RT) : 協同性因子 (w > 1 で正の協同性)
    自由濃度 a, b を 2 次元で数値求解し、各化学種の分率を返す。
    """
    def species(a, b):
        ta = a / K1
        tb = b / K2
        tab = ta * tb * w
        denom = 1.0 + ta + tb + tab
        T0 = Tt / denom
        return T0, T0 * ta, T0 * tb, T0 * tab

    def resid(a, b):
        T0, TA, TB, TAB = species(a, b)
        return (a + TA + TAB - At, b + TB + TAB - Bt)

    # a, b それぞれについて単調なので、交互に brentq / 二分法で詰める
    a, b = At, Bt
    for _ in range(60):
        a_old, b_old = a, b

        def fa(x):
            return resid(x, b)[0]
        def fb(x):
            return resid(a, x)[1]

        for f, lo_hi in ((fa, "a"), (fb, "b")):
            lo, hi = 1e-24, max(At, Bt) * 1.0
            hi = max(hi, 1e-12)
            if f(hi) < 0:
                hi = max(At, Bt) * 10.0
            try:
                if HAVE_SCIPY:
                    root = brentq(f, lo, hi, xtol=1e-26, rtol=1e-10, maxiter=100)
                else:
                    for _ in range(120):
                        mid = math.sqrt(lo * hi)
                        if f(mid) > 0:
                            hi = mid
                        else:
                            lo = mid
                    root = math.sqrt(lo * hi)
            except ValueError:
                root = hi
            if lo_hi == "a":
                a = root
            else:
                b = root

        # 1e-14 では break がほぼ発火せず常に最大回数回って遅かった。
        # 占有率の精度としては 1e-10 で十分 (検証済み)。
        if abs(a - a_old) <= 1e-10 * max(a, 1e-30) and \
           abs(b - b_old) <= 1e-10 * max(b, 1e-30):
            break

    T0, TA, TB, TAB = species(a, b)
    return dict(T=T0 / Tt, TA=TA / Tt, TB=TB / Tt, TAB=TAB / Tt)


def induced_access_change(sites, element):
    """
    マスクが結合したときに、別の機能要素のアクセシビリティがどう変わるか。

    返り値 < 0 : マスクによってその要素が「開く」(構造が壊れる)
    返り値 > 0 : マスクによってその要素に新しい構造ができ、より閉じる

    どちらも実験の交絡要因になる。例えば cap を塞いだ結果 AUG が開いて
    翻訳が上がれば、「cap のマスクが効かなかった」のか「副作用で AUG が
    開いた」のか区別できない。
    マスク鎖に占有された領域は分子内対合に使えないので、強制 unpaired 制約が
    結合状態の正しい表現になっている。
    """
    free = dg_open([element])
    cond = dg_open(list(sites) + [element]) - dg_open(list(sites))
    return cond - free


def residence_time_s(dG, k_on=1.0e6):
    """
    k_off = k_on * Kd から OFF 状態の滞留半減期を出す [s]。
    翻訳反応時間 (1-2 h) と比べて十分長いかを見るための目安。
    k_on は 2'-OMe/RNA の典型値 1e6 M^-1 s^-1 を仮定 (桁の目安)。
    """
    k_off = k_on * kd_from_dg(dG)
    return math.log(2) / k_off if k_off > 0 else float("inf")


# =============================================================================
# 7. 翻訳活性モデル
# =============================================================================
def block_eff(zone, covered, eps_zone=None):
    """
    結合している腕の「合算フットプリント」が機能要素をどれだけ覆うかから阻害効率を出す。

    must_cover の幅で割る形にすると、ゾーン定義の窓幅という恣意的な値がそのまま
    予測 ON/OFF に効いてゾーン間が比較できなくなる。実際には機能要素に一定長
    掛かれば阻害は飽和するはずなので、REQ_NT で頭打ちにする。
    また 43S のロードやスキャニングは「どちらの腕が覆ったか」ではなく占拠された
    領域全体を見るので、腕ごとに独立な阻害を掛けるのではなく合算して評価する。

    covered: [(a, b), ...] その化学種で二重鎖になっている標的領域
    """
    eps = zone["eps"] if eps_zone is None else eps_zone
    ma, mb = zone["must_cover"]
    req = zone.get("req_nt", CFG.REQ_COVER_NT)
    pos = set()
    for (a, b) in covered:
        pos |= set(range(max(a, ma), min(b, mb) + 1))
    return eps * min(1.0, len(pos) / req)


def activity_from_species(zone, frac, up, down, eps_zone=None):
    """各化学種の分率から相対翻訳活性 (マスクなし = 1.0) を出す。"""
    return (frac.get("T", 0.0) * 1.0
            + frac.get("TA", 0.0) * (1 - block_eff(zone, [up], eps_zone))
            + frac.get("TB", 0.0) * (1 - block_eff(zone, [down], eps_zone))
            + frac.get("TAB", 0.0) * (1 - block_eff(zone, [up, down], eps_zone)))


def predict_on_off(rec, eps_scale=None):
    """
    OFF 状態 (二価、切断前) と ON 状態 (切断後) の翻訳活性と ON/OFF 比。
    eps_scale を与えるとゾーンの eps をその値で置き換える (eps 感度解析用)。
    """
    zone = ZONE_BY_NAME[rec["zone"]]
    up = (rec["up_start"], rec["up_end"])
    down = (rec["down_start"], rec["down_end"])

    Tt = CFG.C_MRNA
    Mt = CFG.C_MRNA * CFG.MASK_EXCESS

    # --- OFF: 二価構築体が 1:1 で結合 (両腕が同時に乗る = TAB 相当)
    bound = solve_1to1(Tt, Mt, kd_from_dg(rec["dG_biv"]))
    f_off = {"T": 1 - bound / Tt, "TAB": bound / Tt}
    act_off = activity_from_species(zone, f_off, up, down, eps_scale)

    # --- ON: 切断後。各腕が独立に、かつ協同性 w を伴って結合
    w = math.exp(-rec["ddG_coop"] / RT)
    f_cut = solve_two_site(Tt, Mt, Mt,
                           kd_from_dg(rec["dG_mono_up"]),
                           kd_from_dg(rec["dG_mono_down"]), w)
    act_cut = activity_from_species(zone, f_cut, up, down, eps_scale)

    # 部分切断を線形に混合
    c = CFG.CLEAVAGE_FRAC
    act_on = c * act_cut + (1 - c) * act_off

    blk_full = block_eff(zone, [up, down], eps_scale)
    return dict(
        occ_off=round(bound / Tt, 4),
        occ_on=round(1 - f_cut["T"], 4),
        f_TA=round(f_cut["TA"], 4), f_TB=round(f_cut["TB"], 4),
        f_TAB=round(f_cut["TAB"], 4),
        block_full=round(blk_full, 3),
        act_off=round(act_off, 5), act_on=round(act_on, 5),
        on_off=round(act_on / act_off, 2) if act_off > 1e-12 else np.inf,
        # 完全マスク・完全解離が達成できた場合の上限
        on_off_ceiling=round(1.0 / (1 - blk_full), 2) if blk_full < 1 else np.inf,
    )


# =============================================================================
# 8. 候補の評価
# =============================================================================
def evaluate(zone, a1, b1, a2, b2):
    lo, hi = zone["region"]
    if a1 < lo or b2 > hi:
        return None, "region"
    gap = a2 - b1 - 1
    if not (CFG.MIN_GAP_NT <= gap <= CFG.MAX_GAP_NT):
        return None, "gap"

    site_up, site_dn = TARGET[a1 - 1:b1], TARGET[a2 - 1:b2]
    m_up, m_dn = revcomp(site_up), revcomp(site_dn)

    # 機能要素を実際に覆っているか (両腕合算のフットプリントで判定)
    ma, mb = zone["must_cover"]
    cov_nt = overlap_len(a1, b1, ma, mb) + overlap_len(a2, b2, ma, mb)
    if cov_nt <= 0:
        return None, "no-coverage"

    for m in (m_up, m_dn):
        if max_base_run(m, "G") > CFG.MAX_G_RUN:
            return None, "G-run"
        if max_run(m) > CFG.MAX_ANY_RUN:
            return None, "homopolymer"

    dH_up, dG_up = nn_duplex(site_up)
    dH_dn, dG_dn = nn_duplex(site_dn)

    # --- 安価なフィルタを先に通す。2 サイト同時の dg_open はペアごとに一意で
    #     キャッシュが効かず (83 nt の制約付き分配関数) ここが律速になるため、
    #     落ちる候補にはそれを計算させない。
    s_up, _ = self_fold_ome(m_up)
    s_dn, _ = self_fold_ome(m_dn)
    if s_up < CFG.SELF_FOLD_MAX or s_dn < CFG.SELF_FOLD_MAX:
        return None, "self-structure"

    cross = intermolecular_dg(m_up, m_dn)
    if cross < CFG.CROSS_DIMER_MAX:
        return None, "cross-dimer"

    if (best_offtarget_dg(m_up, (a1, b1)) - dG_up) < CFG.OFFTARGET_MARGIN:
        return None, "offtarget"
    if (best_offtarget_dg(m_dn, (a2, b2)) - dG_dn) < CFG.OFFTARGET_MARGIN:
        return None, "offtarget"

    lk = optimize_linker(gap)
    if lk is None:
        return None, "linker"

    # 単一領域の dg_open は領域ごとにキャッシュされるので安価。
    # 2 サイト同時のものだけが一意で高価。
    o_up, o_dn = dg_open([(a1, b1)]), dg_open([(a2, b2)])
    o_both = dg_open([(a1, b1), (a2, b2)])

    coax = CFG.DG_COAX if gap <= CFG.COAX_GAP_MAX else 0.0

    # 一価 (切断後に各腕が単独で結合するときの ΔG)
    dG_mono_up = dG_up + o_up - s_up
    dG_mono_dn = dG_dn + o_dn - s_dn

    # 2 サイト同時結合の協同性。ΔG_open の非加算性がここに入る。
    #   ddG_coop < 0 なら正の協同性 (片方を開くともう片方も開く)
    ddG_coop = (o_both - o_up - o_dn) + coax

    # 二価 (切断前): 一価 2 つ + 協同性 + アビディティ項
    dG_biv = dG_mono_up + dG_mono_dn + ddG_coop - cross + lk["dG_link"]

    ct = CFG.C_MRNA * (1 + CFG.MASK_EXCESS)
    rec = dict(
        zone=zone["name"],
        up_start=a1, up_end=b1, L_up=b1 - a1 + 1,
        down_start=a2, down_end=b2, L_down=b2 - a2 + 1,
        gap_nt=gap,
        site_up=site_up, site_down=site_dn,
        MASK1=m_up, MASK2=m_dn,
        GC_up=round(gc_frac(site_up), 2), GC_down=round(gc_frac(site_dn), 2),
        dG_open_up=round(o_up, 2), dG_open_down=round(o_dn, 2),
        dG_open_both=round(o_both, 2),
        ddG_coop=round(ddG_coop, 2),
        dG_self_up=round(s_up, 2), dG_self_down=round(s_dn, 2),
        dG_cross=round(cross, 2),
        dG_mono_up=dG_mono_up, dG_mono_down=dG_mono_dn, dG_biv=dG_biv,
        dG37_mono_up=round(dG_mono_up, 2), dG37_mono_down=round(dG_mono_dn, 2),
        dG37_biv=round(dG_biv, 2),
        Tm_mono_up=round(tm_apparent(dH_up, dG_mono_up, ct), 1),
        Tm_mono_down=round(tm_apparent(dH_dn, dG_mono_dn, ct), 1),
        Tm_biv=round(tm_apparent(dH_up + dH_dn, dG_biv, ct), 1),
        Kd_mono_up_uM=round(kd_from_dg(dG_mono_up) * 1e6, 3),
        Kd_mono_down_uM=round(kd_from_dg(dG_mono_dn) * 1e6, 3),
        Kd_biv_nM=round(kd_from_dg(dG_biv) * 1e9, 4),
        n_EO=lk["n_eo"], linker_contour_A=round(lk["contour"], 1),
        span_req_A=round(lk["r_req"], 1), span_req_max_A=round(lk["r_req_max"], 1),
        reaches_extended=lk["reaches_extended"],
        C_eff_mM=round(lk["ceff"] * 1e3, 3), dG_link=round(lk["dG_link"], 2),
        cov_element_nt=cov_nt,
        # OFF 状態の滞留半減期 [min]。翻訳反応時間より十分長い必要がある
        t_half_off_min=round(residence_time_s(dG_biv) / 60.0, 1),
        t_half_mono_s=round(min(residence_time_s(dG_mono_up),
                                residence_time_s(dG_mono_dn)), 3),
        # 他の機能要素への誘起構造変化 (交絡要因の検出)
        **{f"induced_{nm}": round(
               induced_access_change([(a1, b1), (a2, b2)], (ea, eb)), 2)
           for nm, ea, eb in CFG.WATCH_ELEMENTS
           if overlap_len(a1, b1, ea, eb) == 0 and overlap_len(a2, b2, ea, eb) == 0},
    )
    rec.update(predict_on_off(rec))
    rec["construct"] = f"5'-{m_dn}-[PEG{lk['n_eo']}+SI]-{m_up}-3'"
    rec["Tm_window_ok"] = int(
        CFG.TM_MONO_TARGET[0] <= rec["Tm_mono_up"] <= CFG.TM_MONO_TARGET[1] and
        CFG.TM_MONO_TARGET[0] <= rec["Tm_mono_down"] <= CFG.TM_MONO_TARGET[1])
    return rec, ""


def screen_zone(zone):
    lo, hi = zone["region"]
    rows, rej = [], {}
    for L1 in CFG.ARM_LENGTHS:
        for a1 in range(lo, hi - L1 + 2):
            b1 = a1 + L1 - 1
            for L2 in CFG.ARM_LENGTHS:
                for a2 in range(b1 + 1 + CFG.MIN_GAP_NT, hi - L2 + 2):
                    b2 = a2 + L2 - 1
                    if b2 > hi:
                        break
                    rec, why = evaluate(zone, a1, b1, a2, b2)
                    if rec is None:
                        rej[why] = rej.get(why, 0) + 1
                    else:
                        rows.append(rec)
    return pd.DataFrame(rows), rej


def rank(df):
    """
    ランキングは物理量ベース。Tm 窓のような恣意的な閾値は順位に使わない
    (前版ではそれが順位を歪めた)。
    """
    if df.empty:
        return df
    d = df.copy()
    def nz(x):
        r = x.max() - x.min()
        return (x - x.min()) / r if r > 1e-12 else pd.Series(0.5, index=x.index)
    d["_onoff"] = nz(np.log10(d["on_off"].replace(np.inf, np.nan)).fillna(0))
    d["_off"]   = nz(d["occ_off"])
    d["_on"]    = 1 - nz(d["occ_on"])
    d["_cov"]   = nz(d["cov_element_nt"].clip(upper=CFG.REQ_COVER_NT))
    d["score"]  = (0.40 * d["_onoff"] + 0.25 * d["_off"]
                   + 0.20 * d["_on"] + 0.15 * d["_cov"]).round(4)
    return d.drop(columns=[c for c in d.columns if c.startswith("_")])


def positional_scan(region=None, gaps=None, arms=None):
    """
    ゾーン境界を使わずに、マスクのフットプリント位置を 5'UTR 全域で動かし、
    各位置で到達可能な「スイッチ品質」の上限を出す。

    ゾーンの region は私が手で引いた恣意的な値であり、その中だけを探索して
    いては境界の妥当性を検証できない。この関数はゾーン定義から独立に
    「どの位置なら良いスイッチが組めるか」を出すので、ゾーンの切り方が
    妥当だったか、取りこぼした領域がないかを確認できる。

    スイッチ品質 S = occ_off * (1 - occ_on)  (0-1、1 が理想)
      eps に依存しない = 翻訳阻害効率の仮定を含まない純粋な熱力学指標。
      「塞げるか」と「外れるか」は別の軸なので、eps の仮定と分離して見る。
    """
    lo, hi = region or (1, max(z["region"][1] for z in CFG.ZONES))
    # gap も腕長も決め打ちしない。ここで絞ると「ゾーン境界に依存しない検証」
    # という目的そのものが崩れる (別の隠れた制約が入る)。
    # 長い gap はアビディティ項 -RT ln(C_eff) が物理的に淘汰する。
    gaps = range(CFG.MIN_GAP_NT, CFG.MAX_GAP_NT + 1) if gaps is None else gaps
    arms = CFG.ARM_LENGTHS if arms is None else arms
    probe = dict(name="_scan", region=(lo, hi), must_cover=(lo, hi), eps=0.9)
    ZONE_BY_NAME[probe["name"]] = probe   # predict_on_off が名前で引くため
    best = {}
    for L1 in arms:
        for a1 in range(lo, hi - L1 + 2):
            b1 = a1 + L1 - 1
            for g in gaps:
                a2 = b1 + 1 + g
                for L2 in arms:
                    b2 = a2 + L2 - 1
                    if b2 > hi:
                        continue
                    rec, _ = evaluate(probe, a1, b1, a2, b2)
                    if rec is None:
                        continue
                    S = rec["occ_off"] * (1 - rec["occ_on"])
                    ctr = (a1 + b2) // 2
                    if ctr not in best or S > best[ctr]["S"]:
                        best[ctr] = dict(center=ctr, S=round(S, 4),
                                         occ_off=rec["occ_off"], occ_on=rec["occ_on"],
                                         span=f"{a1}-{b1}/{a2}-{b2}",
                                         L_up=L1, L_down=L2, gap=g,
                                         ddG_coop=rec["ddG_coop"])
    return pd.DataFrame(sorted(best.values(), key=lambda r: r["center"]))


def uncovered_peaks(scan_df, zone_best_S=None, s_min=0.90):
    """
    位置スキャンで見つかった良い位置のうち、どのゾーンにも入っていないものを
    列挙する。ゾーンの region は手で引いた恣意的な値なので、境界をまたぐ
    良い位置を取りこぼしていないかをここで検証する。

    zone_best_S を渡すと、ゾーン内で実際に見つかった最良 S を上回るものだけに
    絞る。単に「ゾーン外にも S>=0.90 の場所がある」だけでは行動に繋がらない
    (ゾーン内の方が良いなら取りこぼしではない)。
    """
    thr = s_min if zone_best_S is None else max(s_min, zone_best_S)
    out = []
    for _, r in scan_df[scan_df.S >= thr].iterrows():
        a1 = int(r["span"].split("/")[0].split("-")[0])
        b2 = int(r["span"].split("/")[1].split("-")[1])
        inside = [z["name"] for z in CFG.ZONES
                  if a1 >= z["region"][0] and b2 <= z["region"][1]]
        if not inside:
            out.append(dict(span=r["span"], center=r["center"], S=r["S"],
                            occ_off=r["occ_off"], occ_on=r["occ_on"],
                            ddG_coop=r["ddG_coop"],
                            straddles="/".join(
                                z["name"] for z in CFG.ZONES
                                if overlap_len(a1, b2, *z["region"]) > 0) or "-"))
    return pd.DataFrame(out)


# =============================================================================
# 9. ロバストネス解析
# =============================================================================
def robustness(zone, top):
    """2'-OMe 補正 / 自壊性ユニット長 / mask:mRNA 比を振って候補の頑健性を見る。"""
    keep = (CFG.DDG_OME_PER_NT, CFG.SI_CONTOUR_A, CFG.MASK_EXCESS)
    rows = []
    grid = list(itertools.product(CFG.ROB_OME, CFG.ROB_SI, CFG.ROB_EXCESS))
    for _, r in top.iterrows():
        tag = f"{int(r['up_start'])}-{int(r['up_end'])}/{int(r['down_start'])}-{int(r['down_end'])}"
        for ome, si, exc in grid:
            CFG.DDG_OME_PER_NT, CFG.SI_CONTOUR_A, CFG.MASK_EXCESS = ome, si, exc
            rec, _ = evaluate(zone, int(r["up_start"]), int(r["up_end"]),
                              int(r["down_start"]), int(r["down_end"]))
            if rec is None:
                rows.append(dict(cand=tag, ome=ome, si=si, exc=exc,
                                 on_off=np.nan, occ_off=np.nan, occ_on=np.nan,
                                 reach=0, ok=False))
                continue
            ok = (rec["occ_off"] >= CFG.OCC_OFF_MIN and rec["occ_on"] <= CFG.OCC_ON_MAX)
            rows.append(dict(cand=tag, ome=ome, si=si, exc=exc,
                             on_off=rec["on_off"], occ_off=rec["occ_off"],
                             occ_on=rec["occ_on"],
                             reach=rec["reaches_extended"], ok=ok))
    CFG.DDG_OME_PER_NT, CFG.SI_CONTOUR_A, CFG.MASK_EXCESS = keep
    d = pd.DataFrame(rows)
    g = d.groupby("cand").agg(
        ok_frac=("ok", "mean"), onoff_min=("on_off", "min"),
        onoff_med=("on_off", "median"), occ_off_min=("occ_off", "min"),
        occ_on_max=("occ_on", "max"), reach_frac=("reach", "mean")
    # ok_frac だけだと同点が出て 1 位が恣意的に決まるので、
    # 最悪条件での ON/OFF で同点を割る
    ).round(3).sort_values(["ok_frac", "onoff_min"], ascending=[False, False])
    return g, d


# =============================================================================
# 10. レポート
# =============================================================================
def show_structure():
    ss, mfe = mfe_fold(FOLD_SEQ)
    print("=" * 104)
    print(f"TARGET {N_T} nt | folding context: "
          f"{'FULL mRNA %d nt (max_bp_span=%d)' % (len(FOLD_SEQ), CFG.MAX_BP_SPAN) if USE_FULL else '83-mer のみ'}"
          f" | T={CFG.ANNEAL_TEMP_C} C | Na+eq={na_equivalent():.3f} M")
    print("=" * 104)
    sub = ss[OFFSET:OFFSET + N_T]
    for o in range(0, N_T, 60):
        sl = slice(o, min(o + 60, N_T))
        print(f"{o+1:>4} {TARGET[sl]}")
        print(f"     {sub[sl]}")
    print(f"\nMFE (folding seq) = {mfe:.2f} kcal/mol")
    if not USE_FULL:
        print("[warn] 83 nt だけで折り畳んでいます。実物の mRNA では両端の文脈で")
        print("       構造が変わりえます。CFG.FULL_MRNA に全長を入れると局所折り畳みで")
        print("       計算し直します。ΔG_open はこの設計の支配因子なので影響は大きいです。")
    print("\nzones:")
    for z in CFG.ZONES:
        a, b = z["region"]; ma, mb = z["must_cover"]
        print(f"  {z['name']:<6} region {a:>3}-{b:<3}  must_cover {ma:>3}-{mb:<3} "
              f"({TARGET[ma-1:mb]})  eps={z['eps']}  {z['note']}")
    print(f"\nassay: [mRNA]={CFG.C_MRNA*1e9:.1f} nM, mask excess x{CFG.MASK_EXCESS:.0f} "
          f"(= {CFG.C_MRNA*CFG.MASK_EXCESS*1e9:.0f} nM), cleavage {CFG.CLEAVAGE_FRAC*100:.0f}%")


SHOW = ["zone", "up_start", "up_end", "L_up", "gap_nt", "down_start", "down_end", "L_down",
        "MASK1", "MASK2", "GC_up", "GC_down",
        "dG_open_up", "dG_open_down", "ddG_coop", "dG_cross",
        "dG37_mono_up", "dG37_mono_down", "dG37_biv",
        "Tm_mono_up", "Tm_mono_down", "Tm_biv", "Tm_window_ok",
        "n_EO", "linker_contour_A", "reaches_extended", "C_eff_mM", "dG_link",
        "Kd_mono_up_uM", "Kd_biv_nM",
        "occ_off", "occ_on", "block_full", "act_off", "act_on", "on_off",
        "on_off_ceiling", "cov_element_nt", "score"]


def eps_from_onoff(rec, R, lo=1e-4, hi=0.9999):
    """
    測定された ON/OFF (= R) から阻害効率 eps を逆算する。
    実験後にこれを使って各ゾーンの eps を確定させ、次段の最適化の入力にする。
    ON/OFF は eps に対して単調増加なので二分法でよい。
    """
    f = lambda e: predict_on_off(rec, eps_scale=e)["on_off"] - R
    if f(lo) > 0:
        return None            # R が小さすぎる (eps ~ 0 でも届かない)
    if f(hi) < 0:
        return None            # R が大きすぎる (占有率の方が律速)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def inverse_table(rec, measured=(1.5, 2, 3, 5, 10, 20, 50)):
    rows = []
    for R in measured:
        e = eps_from_onoff(rec, R)
        rows.append(dict(measured_on_off=R,
                         implied_eps=round(e, 4) if e is not None else np.nan,
                         note="" if e is not None else "この構築体では到達不能"))
    return pd.DataFrame(rows)


def eps_table(rec):
    out = []
    for e in CFG.EPS_SWEEP:
        p = predict_on_off(rec, eps_scale=e)
        out.append(dict(eps=e, act_off=p["act_off"], act_on=p["act_on"],
                        on_off=p["on_off"], ceiling=p["on_off_ceiling"]))
    return pd.DataFrame(out)


def detail(rec):
    print("-" * 96)
    print(f"[{rec['zone']}]  MASK1 (上流腕) 標的 {rec['up_start']}-{rec['up_end']} "
          f"({rec['L_up']} nt)  |  gap {rec['gap_nt']} nt  |  "
          f"MASK2 (下流腕) 標的 {rec['down_start']}-{rec['down_end']} ({rec['L_down']} nt)")
    print(f"  MASK1 : 5'-{rec['MASK1']}-3'   GC={rec['GC_up']}  (all 2'-OMe)")
    print(f"  MASK2 : 5'-{rec['MASK2']}-3'   GC={rec['GC_down']}  (all 2'-OMe)")
    print(f"  構築体: {rec['construct']}")
    print(f"          テザー contour {rec['linker_contour_A']} A / 必要スパン "
          f"{rec['span_req_A']} A (伸長時 {rec['span_req_max_A']} A, "
          f"{'届く' if rec['reaches_extended'] else '届かない'})")
    print(f"  構造  : dG_open {rec['dG_open_up']} / {rec['dG_open_down']} "
          f"(both {rec['dG_open_both']}),  協同性 ddG_coop {rec['ddG_coop']} kcal/mol")
    print(f"  一価  : dG37 {rec['dG37_mono_up']} / {rec['dG37_mono_down']}, "
          f"Tm {rec['Tm_mono_up']} / {rec['Tm_mono_down']} C, "
          f"Kd {rec['Kd_mono_up_uM']} / {rec['Kd_mono_down_uM']} uM")
    print(f"  二価  : dG37 {rec['dG37_biv']}, Tm {rec['Tm_biv']} C, "
          f"Kd {rec['Kd_biv_nM']} nM  (avidity {rec['dG_link']}, C_eff {rec['C_eff_mM']} mM)")
    print(f"  占有率: OFF {rec['occ_off']}  ->  ON {rec['occ_on']} "
          f"(切断後の内訳 TA {rec['f_TA']} / TB {rec['f_TB']} / TAB {rec['f_TAB']})")
    print(f"  翻訳  : act_OFF {rec['act_off']}  act_ON {rec['act_on']}  "
          f"ON/OFF {rec['on_off']} (上限 {rec['on_off_ceiling']})")
    print("  阻害効率 eps 依存性 (eps は実験で測る量。熱力学からは予測できない):")
    for _, r in eps_table(rec).iterrows():
        print(f"      eps={r['eps']:<5} act_OFF={r['act_off']:<9.5f} "
              f"act_ON={r['act_on']:<9.5f} ON/OFF={r['on_off']:<8} (上限 {r['ceiling']})")
    print("  逆引き (実験後: 測定した ON/OFF から eps を読む):")
    for _, r in inverse_table(rec).iterrows():
        e = r["implied_eps"]
        print(f"      測定 ON/OFF = {r['measured_on_off']:<5} -> eps = "
              f"{('%.3f' % e) if e == e else '---'}  {r['note']}")


def order_sheet(picks):
    print("\n" + "=" * 104)
    print("発注シート: 無細胞翻訳 ON/OFF スクリーニング用 (3 ゾーン)")
    print("=" * 104)
    print("全残基 2'-O-methyl RNA。リンカーは MASK2 の 3' 末端と MASK1 の 5' 末端を繋ぐ")
    print("(逆平行配向のため)。合成分子は 5'-MASK2-linker-SI-linker-MASK1-3' の 1 本鎖。\n")
    rows = []
    for r in picks:
        rows.append(dict(
            zone=r["zone"], strand="MASK1 (上流腕)", seq_5to3=r["MASK1"],
            length=r["L_up"], GC=r["GC_up"], target_pos=f"{r['up_start']}-{r['up_end']}",
            Tm_mono=r["Tm_mono_up"], attach="5' 末端をリンカーへ"))
        rows.append(dict(
            zone=r["zone"], strand="MASK2 (下流腕)", seq_5to3=r["MASK2"],
            length=r["L_down"], GC=r["GC_down"], target_pos=f"{r['down_start']}-{r['down_end']}",
            Tm_mono=r["Tm_mono_down"], attach="3' 末端をリンカーへ"))
    print(pd.DataFrame(rows).to_string(index=False))
    print("\nリンカー仕様 (自壊性ユニット長が確定したら PEG 単位数で差分を埋める):")
    for r in picks:
        print(f"  {r['zone']:<6} ギャップ {r['gap_nt']} nt -> テザー contour 目標 "
              f"{max(r['span_req_max_A'], r['span_req_A']*1.3):.0f}-"
              f"{max(r['span_req_max_A'], r['span_req_A']*1.3)+8:.0f} A "
              f"(SI={CFG.SI_CONTOUR_A:.0f} A 仮定なら PEG{r['n_EO']})")
    print("\n対照として必要なもの:")
    print("  - マスクなし mRNA (ON の上限 / 活性 1.0 の基準)")
    print("  - 切断済み相当品: MASK1 と MASK2 を別々の鎖として等モル添加")
    print("    (自壊性リンカーの切断反応を介さずに ON 状態の熱力学だけを検証できる)")
    print("  - スクランブル配列の二価構築体 (非特異的な翻訳阻害の対照)")


def plots(all_df):
    """
    図のラベルは英語で書く。Colab には日本語フォントが入っておらず、
    日本語を使うと字形が落ちて豆腐になるため。
    """
    if not HAVE_PLT or all_df.empty:
        return
    zones = [z["name"] for z in CFG.ZONES if (all_df.zone == z["name"]).any()]

    fig, axes = plt.subplots(1, len(zones), figsize=(5.2 * len(zones), 4.4), squeeze=False)
    for ax, zn in zip(axes[0], zones):
        d = all_df[all_df.zone == zn]
        sc = ax.scatter(d["occ_on"], d["occ_off"],
                        c=np.log10(d["on_off"].clip(1e-3, 1e6)), cmap="viridis", s=16)
        ax.axhline(CFG.OCC_OFF_MIN, color="r", ls="--", lw=1)
        ax.axvline(CFG.OCC_ON_MAX, color="r", ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_xlabel("occupancy ON, post-cleavage (want low)")
        ax.set_ylabel("occupancy OFF (want high)")
        ax.set_title(f"zone: {zn}  (n={len(d)})")
        fig.colorbar(sc, ax=ax, label="log10 ON/OFF")
    plt.tight_layout(); plt.show()

    # --- 逆引き図: 測定 ON/OFF から eps を読む ---
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    es = np.linspace(0.05, 0.995, 120)
    for zn in zones:
        d = all_df[all_df.zone == zn].nlargest(1, "score")
        if d.empty:
            continue
        rec = d.iloc[0].to_dict()
        ys = [predict_on_off(rec, eps_scale=e)["on_off"] for e in es]
        ax.plot(es, ys, lw=2,
                label=f"{zn}  {rec['up_start']}-{rec['up_end']}/"
                      f"{rec['down_start']}-{rec['down_end']}")
    for R in (2, 5, 10, 20):
        ax.axhline(R, color="grey", ls=":", lw=0.8)
        ax.annotate(f"ON/OFF={R}", (0.06, R), fontsize=7, color="grey",
                    va="bottom")
    ax.set_xlabel("blocking efficiency  eps   (measured by THIS experiment)")
    ax.set_ylabel("ON/OFF")
    ax.set_yscale("log")
    ax.grid(alpha=.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Inverse lookup: measured ON/OFF  ->  eps of that zone\n"
                 "(curves overlap = the three constructs are thermodynamically "
                 "equivalent,\n so any measured difference is eps, not binding)")
    plt.tight_layout(); plt.show()


def plot_positional(scan, picks=None):
    """位置スキャンの landscape。ゾーン境界の妥当性を目で確認するための図。"""
    if not HAVE_PLT or scan is None or scan.empty:
        return
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    colors = {"cap": "tab:blue", "mid": "tab:orange", "kozak": "tab:green"}
    for z in CFG.ZONES:
        a, b = z["region"]
        ax.axvspan(a, b, alpha=.10, color=colors.get(z["name"], "grey"))
        ax.text((a + b) / 2, 1.04, z["name"], ha="center", fontsize=9,
                color=colors.get(z["name"], "grey"))
    ax.plot(scan["center"], scan["S"], "o-", color="k", ms=4, lw=1.2,
            label="best achievable switch quality S")
    ax.axhline(0.90, color="r", ls="--", lw=1, label="S = 0.90")
    for _, r in scan.nlargest(3, "S").iterrows():
        ax.annotate(r["span"], (r["center"], r["S"]), textcoords="offset points",
                    xytext=(0, -16), ha="center", fontsize=7.5)
    for nm, a, b in CFG.WATCH_ELEMENTS:
        ax.axvspan(a, b, ymin=0, ymax=0.045, color="purple", alpha=.55)
        ax.text((a + b) / 2, -0.055, nm, ha="center", fontsize=7, color="purple")
    ax.set_xlabel("center of the masked footprint (nt, 1-based)")
    ax.set_ylabel("S = occ_OFF x (1 - occ_ON)")
    ax.set_ylim(-0.08, 1.12)
    ax.grid(alpha=.3)
    ax.legend(fontsize=8, loc="lower center")
    ax.set_title("Switch quality landscape, independent of zone boundaries and of eps\n"
                 "(gaps = no viable pair exists at that position)")
    plt.tight_layout(); plt.show()


# =============================================================================
# 11. MAIN
# =============================================================================
def main():
    show_structure()

    all_df, picks = [], []
    for zone in CFG.ZONES:
        print("\n" + "=" * 104)
        print(f"ZONE: {zone['name']}  region {zone['region']}  "
              f"must_cover {zone['must_cover']}")
        print("=" * 104)
        df, rej = screen_zone(zone)
        print(f"candidates kept: {len(df)}")
        if rej:
            print("rejected:", ", ".join(f"{k}={v}" for k, v in sorted(rej.items())))
        if df.empty:
            print("このゾーンには候補がありません。region を広げるか "
                  "ARM_LENGTHS / フィルタを緩めてください。")
            continue
        df = rank(df).sort_values("score", ascending=False)
        df.to_csv(os.path.join(CFG.OUTDIR, f"zone_{zone['name']}.csv"), index=False)
        all_df.append(df)

        print(f"\n--- TOP {CFG.TOP_N} ---")
        print(df[SHOW].head(CFG.TOP_N).to_string(index=False))
        print(f"\n--- 推奨候補の詳細 ---")
        detail(df.iloc[0].to_dict())
        picks.append(df.iloc[0].to_dict())

        if CFG.RUN_ROBUSTNESS:
            print(f"\n--- ロバストネス解析 (OMe x SI長 x mask比 = "
                  f"{len(CFG.ROB_OME)*len(CFG.ROB_SI)*len(CFG.ROB_EXCESS)} 条件) ---")
            g, _ = robustness(zone, df.head(CFG.ROB_TOP_K))
            print(g.to_string())
            best = g.index[0]
            cur = f"{picks[-1]['up_start']}-{picks[-1]['up_end']}/" \
                  f"{picks[-1]['down_start']}-{picks[-1]['down_end']}"
            if best != cur:
                print(f"\n[note] score 1 位 ({cur}) とロバストネス 1 位 ({best}) が異なります。")
                print("       score は 1 点のパラメータ設定での値です。パラメータに不確実性が")
                print("       ある以上、ok_frac と onoff_min が高い方を選ぶのが安全です。")

    if not all_df:
        print("\n全ゾーンで候補ゼロ。CFG を見直してください。")
        return

    # --- ゾーン定義から独立した位置スキャン ---
    print("\n" + "=" * 104)
    print("位置スキャン (ゾーン境界を使わない検証)")
    print("=" * 104)
    print("region は手で引いた恣意的な値なので、その中だけを探索していては")
    print("境界の妥当性を検証できない。フットプリント位置を 5'UTR 全域で動かし、")
    print("各位置で到達可能なスイッチ品質 S = occ_off * (1 - occ_on) の上限を出す。")
    print("S は eps に依存しない = 翻訳阻害効率の仮定を含まない純粋な熱力学指標。\n")
    scan = positional_scan()
    scan.to_csv(os.path.join(CFG.OUTDIR, "positional_scan.csv"), index=False)
    print(scan.to_string(index=False))

    zone_best_S = max(r["occ_off"] * (1 - r["occ_on"]) for r in picks)
    print(f"\nゾーン内で実際に得られた最良 S = {zone_best_S:.4f}")
    unc = uncovered_peaks(scan, zone_best_S=zone_best_S)
    if not unc.empty:
        print("[重要] ゾーン内の最良を上回る位置が、どのゾーンにも収まっていません。")
        print("       ゾーン境界の引き方で取りこぼしています:")
        print(unc.to_string(index=False))
    else:
        print("ゾーン外に、ゾーン内の最良を上回る位置はありません。")
        print("= 現在のゾーン境界で取りこぼしは生じていない。")
    allcat = pd.concat(all_df, ignore_index=True)
    allcat.to_csv(os.path.join(CFG.OUTDIR, "zones_all.csv"), index=False)

    print("\n" + "=" * 104)
    print("ゾーン間比較 (各ゾーンの推奨候補)")
    print("=" * 104)
    cmp_cols = ["zone", "MASK1", "MASK2", "L_up", "L_down", "gap_nt", "n_EO",
                "ddG_coop", "dG37_biv", "Kd_biv_nM", "Tm_mono_up", "Tm_mono_down",
                "occ_off", "occ_on", "t_half_off_min", "t_half_mono_s",
                "block_full", "on_off", "on_off_ceiling", "cov_element_nt"]
    cmp_df = pd.DataFrame(picks)
    ind_cols = [c for c in cmp_df.columns if c.startswith("induced_")]
    print(cmp_df[cmp_cols + ind_cols].to_string(index=False))
    print("\n  t_half_off_min : OFF 複合体の解離半減期。翻訳反応時間 (60-120 min) と")
    print("                   同程度以下だと、平衡が保たれない状況で OFF が漏れる。")
    print("  t_half_mono_s  : 切断後の解離半減期。短いほど ON への切替が速い。")
    print("  induced_*      : マスクが覆っていない機能要素への誘起構造変化。")
    print("                   負 = その要素が開く / 正 = 新たな構造ができ閉じる。")
    print("                   |値| が大きいほど実験解釈の交絡要因になる。")

    order_sheet(picks)
    plots(allcat)
    plot_positional(scan, picks)

    print("\n" + "=" * 104)
    print("出力ファイル:")
    for f in [f"zone_{z['name']}.csv" for z in CFG.ZONES] + ["zones_all.csv"]:
        p = os.path.join(CFG.OUTDIR, f)
        if os.path.exists(p):
            print("  " + p)
    print("=" * 104)
    print("""
読み方
------
* occ_off / occ_on : 熱力学から予測できる量。OFF で高く、ON で低いほどよい。
                     無細胞系は濃度が既知なので、ここは信頼してよい。
* on_off           : eps (マスクされた mRNA の翻訳がどれだけ落ちるか) の仮定に
                     依存する。eps は熱力学からは出ないので、これが実験で測る量。
                     eps 依存性の表を見て、どの範囲なら差が見えるか確認すること。
* ddG_coop         : ΔG_open の非加算性。負なら「片方のサイトを開くともう片方も
                     開く」= 二価が得をする。この設計の主な駆動力であり、
                     アビディティ項より大きいことが多い。
* reaches_extended : ギャップ ssRNA が伸長してもリンカーが届くか。0 なら要再設計。
* ロバストネス     : score 1 位とロバストネス 1 位が食い違ったら後者を採ること。
""")


if __name__ == "__main__":
    main()
