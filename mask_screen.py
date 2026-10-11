# -*- coding: utf-8 -*-
# =============================================================================
#  Bivalent self-immolative mask strand screening for 5'UTR masking
#  ---------------------------------------------------------------------------
#  設計:   5'-[Arm_down]-(linker)-[SI unit]-(linker)-[Arm_up]-3'
#          2本の 10 nt 2'-OMe マスク鎖が 5'UTR 上の2箇所に同時結合（二価 = 高親和性）
#          バイオマーカーで自壊性リンカーが切断 → 一価化 → アビディティ消失 → 解離
#
#  本スクリプトがやること:
#   STAGE 0: 全 10-nt ウィンドウの単独特性（accessibility / Tm / GC / 自己構造 / off-target）
#   STAGE 1: 2サイトの全組合せを列挙し、リンカー長を各組合せごとに最適化
#   STAGE 2: 一価/二価の ΔG37, Kd, 見かけ Tm, 占有率, スイッチ窓を計算してランキング
#   STAGE 3: 上位候補について各腕の鎖長 (6-14 nt) を再スキャンし、目標 Tm に合わせる
#
#  Colab: このセルをそのまま貼って実行（ViennaRNA は自動 pip install）
# =============================================================================

# ----------------------------------------------------------------------------
# 0. 依存関係
# ----------------------------------------------------------------------------
try:
    import RNA
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "ViennaRNA"], check=True)
    import RNA

import math, itertools, os, sys
import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAVE_PLT = True
except ImportError:
    HAVE_PLT = False

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 100)

# =============================================================================
# 1. CONFIG  --- ここだけ編集すれば挙動が変わります
# =============================================================================
class CFG:
    # ---- 標的 ------------------------------------------------------------
    TARGET = ("AGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACC"
              "AUGUUCGUGUUCCUGGUGCUGCUGCCUCUG")
    UTR_RANGE   = (1, 53)    # 1-based, inclusive
    CDS_RANGE   = (54, 83)

    # マスク鎖を置いてよい領域（1-based, inclusive）
    # 5'UTR 内に限る場合は (1, 53)。AUG を覆わせたい場合は (1, 60) 等に広げる
    MASK_REGION = (1, 53)

    # ---- マスク鎖 --------------------------------------------------------
    L1 = 10                  # 上流側腕の長さ (STAGE1 ではこれで固定)
    L2 = 10                  # 下流側腕の長さ
    ALL_2OME = True          # 全残基 2'-OMe

    # 2'-OMe 補正: 2'-OMe RNA : RNA は RNA:RNA より安定
    #   文献値 ΔTm ≒ +0.5 ~ +1.5 degC / modification (Freier & Altmann 1997 等)
    #   ここでは ΔΔG37 として与える（ΔH は不変、ΔS を再計算 = エントロピー寄与扱い）
    DDG_OME_PER_NT = -0.14   # kcal/mol / 修飾残基（対合している残基のみカウント）

    # ---- 溶液条件 --------------------------------------------------------
    TEMP_C      = 37.0       # ViennaRNA / ΔG37 の基準温度
    NA_CONC     = 0.150      # 一価カチオン等価濃度 [M]（細胞内 K+ 相当）
    MG_CONC     = 0.0010     # 遊離 Mg2+ [M]  (0 にすると無効)
    USE_MG      = True       # True なら Na+ 等価に繰り込む（粗い近似）

    # 濃度は 2 つに分ける。ここが設計の成否を大きく左右する:
    #   C_MASK_LOAD : in vitro でマスクを装着するときの濃度（二価の占有率 = マスク効率）
    #   C_MASK_CELL : 送達後・切断後に遊離腕が実効的に感じる濃度（一価のリーク）
    #                 細胞内では希釈され拡散で逃げるので LOAD より 1-2 桁低いのが現実的
    C_MASK_LOAD = 1.0e-6     # [M]
    C_MASK_CELL = 5.0e-8     # [M]
    C_TOTAL_TM  = 1.0e-6     # Tm 計算の全鎖濃度 C_T [M]（UV melting 条件を想定）

    # ---- リンカー（PEG 系を想定）-----------------------------------------
    # 全テザー = spacer(PEG n1) + 自壊性ユニット + spacer(PEG n2)
    EO_CONTOUR_A   = 3.5     # -CH2CH2O- 1 単位あたりの contour length [Å]
    EO_PERSIST_A   = 3.8     # PEG の持続長 l_p [Å]（Kuhn 長 ~7.6 Å の半分）
    SI_CONTOUR_A   = 12.0    # 自壊性ユニット自体の contour length [Å]（設計依存、要調整）
    TERM_LINK_A    = 5.0     # 末端修飾（3'-amino / 5'-hexynyl 等）の寄与 [Å]
    EO_SCAN        = range(0, 81)   # 探索する PEG 単位数（片側合計）
    # C_eff は L_c に対して非常に浅い最大を持つため「最大化」で選ぶと
    # 無意味に長いリンカーが選ばれる。最大値の CEFF_FRAC に達する最短長を推奨値とする。
    CEFF_FRAC      = 0.50

    # 末端間距離の見積り
    NICK_OFFSET_A  = 6.0     # ギャップ 0（nick）時の末端間距離 [Å]
    RISE_SS_A      = 3.4     # ギャップ ssRNA の 1 nt あたり伸長 [Å]（stacked 寄り）
    RISE_SS_MAX_A  = 6.0     # 完全伸長時 [Å]（上限報告用）

    # ---- 構造・配列フィルタ ----------------------------------------------
    MIN_GAP_NT      = 2      # 2サイト間の最小ギャップ。小さすぎると切断後に
                             # 2断片が隣接結合 + coaxial stack して協同的に再結合する
    MAX_GAP_NT      = 40     # 最大ギャップ
    # GC は Tm 窓フィルタで実質的に制約されるので、ここは極端な配列の除去のみ
    GC_MIN, GC_MAX  = 0.15, 0.75
    MAX_G_RUN       = 3      # GGGG 以上を禁止（G4 / 合成性）
    MAX_ANY_RUN     = 4      # 同一塩基の連続上限
    SELF_FOLD_MAX   = -1.5   # マスク鎖単独の MFE がこれより安定なら不可 [kcal/mol]
    CROSS_DIMER_MAX = -4.0   # 2本のマスク鎖間の分子間 ΔG がこれより安定なら不可
    OFFTARGET_MARGIN = 3.0   # 意図サイト以外の最良ハイブリ ΔG との差 [kcal/mol]

    # coaxial stacking（nick 部）のボーナス（二価状態 / 切断後協同状態の両方に効く）
    DG_COAX         = -1.3   # kcal/mol、ギャップ <= COAX_GAP_MAX のとき適用
    COAX_GAP_MAX    = 1

    # ---- 目標熱力学値（フィルタ & ランキング）-----------------------------
    TM_MONO_TARGET  = (18.0, 33.0)   # 一価の見かけ Tm 目標窓 [degC]
    # 二価の融解は ΔH が 2 本鎖分（~ -140 kcal/mol）で極端に急峻なため、
    # Kd が sub-nM でも見かけ Tm は 45 degC 程度にしかならない。
    # 従って Tm_biv は補助指標とし、判定は occ_biv / Kd_biv を主にする。
    TM_BIV_MIN      = 42.0           # 二価の見かけ Tm 下限 [degC]
    OCC_BIV_MIN     = 0.90           # 二価での標的占有率（= マスク効率）下限
    OCC_MONO_MAX    = 0.10           # 一価での残存占有率（= リーク）上限
    SWITCH_MIN_LOG  = 2.0            # log10(Kd_mono/Kd_biv) 下限

    # ---- 機能要素（翻訳阻害に効く位置）------------------------------------
    # (name, start, end, weight)   1-based inclusive
    FUNC_ELEMENTS = [
        ("cap_proximal", 1,  15, 1.0),   # 43S ロード / eIF4F
        ("mid_utr",      16, 40, 0.3),
        ("kozak",        46, 53, 1.0),   # GCCACC
        ("AUG",          54, 56, 1.2),
    ]

    # ---- 実行制御 --------------------------------------------------------
    TOP_N          = 25
    RUN_STAGE3     = True
    STAGE3_TOP_K   = 6
    STAGE3_LENGTHS = range(6, 15)   # 10 nt + 2'-OMe は一価でも強すぎる傾向 → 短い側も探索
    OUTDIR         = "/content" if os.path.isdir("/content") else "."
    VERBOSE        = True


R_GAS  = 1.9872e-3            # kcal/(mol K)
T0     = 273.15
TK     = CFG.TEMP_C + T0
RT     = R_GAS * TK           # ~0.6163 kcal/mol at 37 C
AVOGADRO = 6.02214076e23


# =============================================================================
# 2. 配列ユーティリティ
# =============================================================================
COMP = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}

def to_rna(s):
    return s.upper().replace("T", "U").replace(" ", "").replace("\n", "")

def revcomp(s):
    return "".join(COMP[c] for c in reversed(to_rna(s)))

def gc_frac(s):
    s = to_rna(s)
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

TARGET = to_rna(CFG.TARGET)
N_T = len(TARGET)


# =============================================================================
# 3. 最近接塩基対 (NN) 熱力学  --- Xia et al. 1998 (RNA/RNA, 1 M NaCl)
#    ΔH [kcal/mol], ΔS [cal/(mol K)]
#    完全相補二重鎖なので、標的側配列の dinucleotide step から計算して等価
# =============================================================================
NN_XIA = {
    "AA": (-6.82, -19.0), "UU": (-6.82, -19.0),
    "AU": (-9.38, -26.7),
    "UA": (-7.69, -20.5),
    "CU": (-10.48, -27.1), "AG": (-10.48, -27.1),
    "CA": (-10.44, -26.9), "UG": (-10.44, -26.9),
    "GU": (-11.40, -29.5), "AC": (-11.40, -29.5),
    "GA": (-12.44, -32.5), "UC": (-12.44, -32.5),
    "CG": (-10.64, -26.7),
    "GG": (-13.39, -32.7), "CC": (-13.39, -32.7),
    "GC": (-14.88, -36.9),
}
DH_INIT, DS_INIT = 3.61, -1.5          # duplex initiation
DH_TERM_AU, DS_TERM_AU = 3.72, 10.5    # terminal A-U penalty (ΔG37 = +0.45)


def na_equivalent():
    """Mg2+ を一価カチオン等価に粗く繰り込む。"""
    na = CFG.NA_CONC
    if CFG.USE_MG and CFG.MG_CONC > 0:
        # Ahsen et al. 風の経験式（粗い近似）。過大評価を避けるため上限をかける
        na = na + 3.3 * math.sqrt(max(CFG.MG_CONC, 0.0))
    return min(max(na, 0.01), 1.2)


def nn_duplex(site, n_ome_paired=None):
    """
    完全相補二重鎖の ΔH, ΔS, ΔG37, Tm を返す。
    site: 標的側の 5'->3' 配列（マスク鎖はその完全相補鎖）
    n_ome_paired: 2'-OMe 補正を掛ける修飾残基数（既定 = マスク鎖長）
    """
    s = to_rna(site)
    n = len(s)
    dH, dS = DH_INIT, DS_INIT
    for k in range(n - 1):
        step = s[k:k + 2]
        if step not in NN_XIA:
            raise ValueError(f"unknown NN step {step}")
        h, ds = NN_XIA[step]
        dH += h
        dS += ds
    for end in (s[0], s[-1]):
        if end in "AU":
            dH += DH_TERM_AU
            dS += DS_TERM_AU

    # 塩濃度補正 (phosphate 数 = N-1)
    dS += 0.368 * (n - 1) * math.log(na_equivalent())

    dG37 = dH - TK * dS / 1000.0

    # 2'-OMe 補正（ΔG ベース。ΔH は不変として ΔS を再計算）
    if CFG.ALL_2OME:
        nmod = n if n_ome_paired is None else n_ome_paired
        dG37 += nmod * CFG.DDG_OME_PER_NT
        dS = (dH - dG37) / TK * 1000.0

    return dH, dS, dG37


def tm_from_hs(dH, dS, ct=None, bimolecular=True):
    """非自己相補二重鎖の Tm [degC]。ct = 全鎖濃度 (equimolar)。"""
    ct = CFG.C_TOTAL_TM if ct is None else ct
    conc_term = R_GAS * math.log(ct / 4.0) * 1000.0 if bimolecular else 0.0
    denom = dS + conc_term
    if denom >= 0:
        return float("nan")
    return dH * 1000.0 / denom - T0


def tm_from_dg(dH, dG37, ct=None, bimolecular=True):
    """ΔH と（補正後の）ΔG37 から見かけ Tm を出す。"""
    dS_eff = (dH - dG37) / TK * 1000.0
    return tm_from_hs(dH, dS_eff, ct=ct, bimolecular=bimolecular)


def kd_from_dg(dG):
    """ΔG37 [kcal/mol] -> Kd [M]"""
    return math.exp(dG / RT)


def occupancy(kd, c_free):
    return c_free / (kd + c_free)


# =============================================================================
# 4. ViennaRNA ラッパ
# =============================================================================
RNA.cvar.temperature = CFG.TEMP_C

def _md():
    md = RNA.md()
    md.temperature = CFG.TEMP_C
    return md

_PF_CACHE = {}

def ens_fold_dg(seq):
    """単鎖のアンサンブル自由エネルギー（分配関数ベース）[kcal/mol]"""
    key = ("pf", seq)
    if key in _PF_CACHE:
        return _PF_CACHE[key]
    fc = RNA.fold_compound(seq, _md())
    _, g = fc.pf()
    _PF_CACHE[key] = g
    return g

def mfe_fold(seq):
    fc = RNA.fold_compound(seq, _md())
    ss, e = fc.mfe()
    return ss, e

def _force_unpaired(fc, regions):
    for (a, b) in regions:                  # 1-based inclusive
        for k in range(a, b + 1):
            try:
                fc.hc_add_up(k)
            except Exception:
                fc.hc_add_up(k, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

def dg_open(seq, regions):
    """
    指定領域を強制的に unpaired にするコスト（RNAup 的 opening energy）。
    regions: [(start1,end1), ...]  1-based inclusive
    返り値 >= 0
    """
    key = ("open", seq, tuple(sorted(regions)))
    if key in _PF_CACHE:
        return _PF_CACHE[key]
    g_free = ens_fold_dg(seq)
    fc = RNA.fold_compound(seq, _md())
    _force_unpaired(fc, regions)
    _, g_con = fc.pf()
    val = g_con - g_free
    _PF_CACHE[key] = val
    return val

def dimer_dg(a, b):
    """2本鎖の MFE（分子内も含む全体）[kcal/mol]"""
    fc = RNA.fold_compound(a + "&" + b, _md())
    try:
        ss, e = fc.mfe_dimer()
    except Exception:
        ss, e = RNA.cofold(a + "&" + b)
    return ss, e

def intermolecular_dg(a, b):
    """
    分子間成分のみ: ΔG(AB) - ΔG(A) - ΔG(B)
    2'-OMe 補正は、形成された bp 数 x 2 残基分を掛ける
    """
    ss, e_ab = dimer_dg(a, b)
    _, e_a = mfe_fold(a)
    _, e_b = mfe_fold(b)
    inter = e_ab - e_a - e_b
    nbp = ss.count("(")
    if CFG.ALL_2OME:
        inter += 2 * nbp * CFG.DDG_OME_PER_NT
    return inter, ss

def self_fold_dg_ome(seq):
    """マスク鎖単独の自己構造 MFE（2'-OMe 補正込み）[kcal/mol], <= 0"""
    ss, e = mfe_fold(seq)
    nbp = ss.count("(")
    if CFG.ALL_2OME:
        e += 2 * nbp * CFG.DDG_OME_PER_NT
    return min(e, 0.0), ss

def best_offtarget_dg(masker, target, exclude):
    """
    マスク鎖が標的上の意図外サイトと作るハイブリの最良 ΔG。
    exclude: (start, end) 1-based inclusive — 意図サイト（この範囲と重なる解は除外）
    duplexfold を窓スライドで評価
    """
    L = len(masker)
    best = 0.0
    best_pos = None
    for s in range(0, len(target) - L + 1):
        a, b = s + 1, s + L
        if not (b < exclude[0] or a > exclude[1]):
            continue
        d = RNA.duplexfold(masker, target[s:s + L])
        e = d.energy
        if CFG.ALL_2OME:
            e += L * CFG.DDG_OME_PER_NT        # 粗い補正
        if e < best:
            best, best_pos = e, (a, b)
    return best, best_pos


# =============================================================================
# 5. リンカー / アビディティモデル
# =============================================================================
def required_span(gap_nt, rise=None):
    """
    二価結合時に、上流腕の 5' 末端と下流腕の 3' 末端が取るべき距離 [Å]。
    （逆平行なので、ギャップに面する末端は up-arm の 5' と down-arm の 3'）
    """
    rise = CFG.RISE_SS_A if rise is None else rise
    return CFG.NICK_OFFSET_A + gap_nt * rise


def _wlc_shape(x, L_over_lp):
    """
    Thirumalai-Ha / Becker 型の半屈曲性鎖 末端間距離分布（規格化前の形）。
    x = r / L_c
      P(x) ~ 4 pi x^2 (1-x^2)^{-9/2} exp( -3 L/(4 l_p (1-x^2)) )
    理想 Gauss 鎖と違い x -> 1 で正しくゼロに落ちるため、
    「ぎりぎり届く長さ」のエントロピーコストを過小評価しない。
    """
    if x <= 0.0 or x >= 1.0:
        return 0.0
    t = 1.0 - x * x
    expo = -0.75 * L_over_lp / t
    if expo < -700.0:
        return 0.0
    return 4.0 * math.pi * x * x * t ** (-4.5) * math.exp(expo)


_WLC_NORM_CACHE = {}

def _wlc_norm(L_over_lp, n=400):
    """Simpson 法で ∫0^1 P(x) dx を求める（規格化定数）。"""
    key = round(L_over_lp, 4)
    if key in _WLC_NORM_CACHE:
        return _WLC_NORM_CACHE[key]
    h = 1.0 / n
    acc = 0.0
    for i in range(n + 1):
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        acc += w * _wlc_shape(i * h, L_over_lp)
    val = acc * h / 3.0
    _WLC_NORM_CACHE[key] = val
    return val


def ceff_wlc(r_A, contour_A, lp_A=None):
    """
    テザー末端が距離 r に存在する確率密度から実効濃度 C_eff [M] を出す。
      density[1/Å^3] = P(r) / (4 pi r^2)      (P は径方向分布, 1/Å に規格化)
      C_eff[M]       = density * 1e27 / N_A
    """
    lp_A = CFG.EO_PERSIST_A if lp_A is None else lp_A
    if contour_A <= 0 or r_A <= 0 or r_A >= contour_A:
        return 0.0
    Llp = contour_A / lp_A
    Z = _wlc_norm(Llp)
    if Z <= 0:
        return 0.0
    P_r = _wlc_shape(r_A / contour_A, Llp) / (Z * contour_A)   # [1/Å]
    dens = P_r / (4.0 * math.pi * r_A * r_A)                   # [1/Å^3]
    return dens * 1.0e27 / AVOGADRO


def optimize_linker(gap_nt):
    """
    ギャップに対して PEG 単位数を走査する。
    C_eff(L_c) は非常に浅い最大しか持たないので、
      n_eo_max : C_eff 最大の長さ（参考値）
      n_eo     : 最大値の CFG.CEFF_FRAC 以上を与える「最短」長 = 推奨値
    を返し、スコアリングには推奨値を使う。
    """
    r_req     = required_span(gap_nt, CFG.RISE_SS_A)
    r_req_max = required_span(gap_nt, CFG.RISE_SS_MAX_A)

    scan = []
    for n_eo in CFG.EO_SCAN:
        contour = n_eo * CFG.EO_CONTOUR_A + CFG.SI_CONTOUR_A + CFG.TERM_LINK_A
        if contour <= r_req:                 # 物理的に届かない
            continue
        scan.append((n_eo, contour, ceff_wlc(r_req, contour)))

    base = dict(n_eo=np.nan, n_eo_max=np.nan, contour=np.nan, ceff=0.0,
                r_req=r_req, r_req_max=r_req_max, dG_link=99.0, at_bound=0)
    if not scan:
        return base

    n_max, c_max, ce_max = max(scan, key=lambda z: z[2])
    if ce_max <= 0:
        return base
    thr = CFG.CEFF_FRAC * ce_max
    n_rec, c_rec, ce_rec = next(z for z in scan if z[2] >= thr)

    return dict(
        n_eo=n_rec, contour=c_rec, ceff=ce_rec,
        n_eo_max=n_max, ceff_max=ce_max,
        r_req=r_req, r_req_max=r_req_max,
        dG_link=-RT * math.log(ce_rec),
        # 最大値が走査範囲の端で得られている場合は EO_SCAN を広げる必要がある
        at_bound=int(n_rec >= max(CFG.EO_SCAN)),
    )


# =============================================================================
# 6. 機能要素カバー率
# =============================================================================
def functional_score(sites):
    """sites: [(a,b), ...] 1-based inclusive。重み付きカバー率 0..1"""
    covered = set()
    for (a, b) in sites:
        covered.update(range(a, b + 1))
    num = den = 0.0
    detail = {}
    for name, a, b, w in CFG.FUNC_ELEMENTS:
        tot = b - a + 1
        hit = len(covered & set(range(a, b + 1)))
        detail[name] = hit / tot
        num += w * (hit / tot)
        den += w
    return (num / den if den else 0.0), detail


# =============================================================================
# 7. STAGE 0 : 全単一ウィンドウの特性
# =============================================================================
def stage0(length):
    lo, hi = CFG.MASK_REGION
    rows = []
    for a in range(lo, hi - length + 2):
        b = a + length - 1
        site = TARGET[a - 1:b]
        masker = revcomp(site)

        dH, dS, dG37 = nn_duplex(site)
        tm_ideal = tm_from_hs(dH, dS)
        d_open = dg_open(TARGET, [(a, b)])
        d_self, ss_self = self_fold_dg_ome(masker)
        dG_eff = dG37 + d_open - d_self
        tm_eff = tm_from_dg(dH, dG_eff)

        ot, ot_pos = best_offtarget_dg(masker, TARGET, (a, b))

        rows.append(dict(
            start=a, end=b, site=site, masker=masker,
            GC=round(gc_frac(site), 3),
            dH=round(dH, 2), dG37_duplex=round(dG37, 2),
            Tm_ideal=round(tm_ideal, 1),
            dG_open=round(d_open, 2),
            dG_self=round(d_self, 2), self_struct=ss_self,
            dG37_eff=round(dG_eff, 2), Tm_eff=round(tm_eff, 1),
            Kd_eff_uM=round(kd_from_dg(dG_eff) * 1e6, 3),
            dG_offtarget=round(ot, 2),
            offtarget_margin=round(ot - dG37, 2),
            max_run=max_run(site), G_run=max_base_run(site, "G"),
        ))
    return pd.DataFrame(rows)


# =============================================================================
# 8. STAGE 1-2 : 2サイト組合せの列挙と評価
# =============================================================================
def seq_filters_ok(masker, site):
    g = gc_frac(site)
    if not (CFG.GC_MIN <= g <= CFG.GC_MAX):
        return False, "GC"
    if max_base_run(masker, "G") > CFG.MAX_G_RUN:
        return False, "G-run"
    if max_run(masker) > CFG.MAX_ANY_RUN:
        return False, "homopolymer"
    return True, ""


def evaluate_pair(a1, b1, a2, b2):
    """
    a1<b1 < a2<b2 (1-based inclusive)。上流腕 = arm_up (site1), 下流腕 = arm_down (site2)
    """
    site1, site2 = TARGET[a1 - 1:b1], TARGET[a2 - 1:b2]
    m_up, m_down = revcomp(site1), revcomp(site2)
    gap = a2 - b1 - 1

    # --- 配列フィルタ
    for m, s in ((m_up, site1), (m_down, site2)):
        ok, why = seq_filters_ok(m, s)
        if not ok:
            return None, f"seq:{why}"

    # --- 単腕の熱力学
    dH1, dS1, dG1 = nn_duplex(site1)
    dH2, dS2, dG2 = nn_duplex(site2)
    open1 = dg_open(TARGET, [(a1, b1)])
    open2 = dg_open(TARGET, [(a2, b2)])
    open12 = dg_open(TARGET, [(a1, b1), (a2, b2)])
    self1, ss1 = self_fold_dg_ome(m_up)
    self2, ss2 = self_fold_dg_ome(m_down)
    if self1 < CFG.SELF_FOLD_MAX or self2 < CFG.SELF_FOLD_MAX:
        return None, "self-structure"

    cross, ss_cross = intermolecular_dg(m_up, m_down)
    if cross < CFG.CROSS_DIMER_MAX:
        return None, "cross-dimer"

    # --- off-target
    ot1, _ = best_offtarget_dg(m_up, TARGET, (a1, b1))
    ot2, _ = best_offtarget_dg(m_down, TARGET, (a2, b2))
    if (ot1 - dG1) < CFG.OFFTARGET_MARGIN or (ot2 - dG2) < CFG.OFFTARGET_MARGIN:
        return None, "offtarget"

    # --- 一価（切断後）: 各腕が独立に結合
    dG_mono1 = dG1 + open1 - self1
    dG_mono2 = dG2 + open2 - self2
    tm_mono1 = tm_from_dg(dH1, dG_mono1)
    tm_mono2 = tm_from_dg(dH2, dG_mono2)
    kd_mono1, kd_mono2 = kd_from_dg(dG_mono1), kd_from_dg(dG_mono2)

    # 切断後に両断片が隣接結合 + coaxial stack して協同的に再結合するリーク
    coax = CFG.DG_COAX if gap <= CFG.COAX_GAP_MAX else 0.0
    dG_coop_cut = dG1 + dG2 + open12 - self1 - self2 + coax     # 2分子が同時結合（3体）
    # 2段階結合として扱い、実効 Kd（見かけ 1 分子目基準）を出す
    cc = CFG.C_MASK_CELL
    occ_cut_coop = (cc / (kd_from_dg(dG_mono1) + cc)) * \
                   (cc / (kd_from_dg(dG_coop_cut - dG_mono1) + cc))

    # --- リンカー最適化 & アビディティ
    if not (CFG.MIN_GAP_NT <= gap <= CFG.MAX_GAP_NT):
        return None, "gap"
    lk = optimize_linker(gap)
    if lk["ceff"] <= 0:
        return None, "linker"

    # --- 二価
    dG_biv = (dG1 + dG2 + open12 - self1 - self2 - cross
              + lk["dG_link"] + coax)
    dH_biv = dH1 + dH2
    tm_biv = tm_from_dg(dH_biv, dG_biv)
    kd_biv = kd_from_dg(dG_biv)

    occ_biv = occupancy(kd_biv, CFG.C_MASK_LOAD)
    occ_mono_worst = max(occupancy(kd_mono1, CFG.C_MASK_CELL),
                         occupancy(kd_mono2, CFG.C_MASK_CELL))
    occ_leak = max(occ_mono_worst, occ_cut_coop)

    kd_mono_eff = min(kd_mono1, kd_mono2)
    switch_log = math.log10(kd_mono_eff / kd_biv) if kd_biv > 0 else np.nan

    fs, fdet = functional_score([(a1, b1), (a2, b2)])

    rec = dict(
        # 位置
        up_start=a1, up_end=b1, down_start=a2, down_end=b2, gap_nt=gap,
        site_up=site1, site_down=site2,
        mask_up=m_up, mask_down=m_down,
        # 構築体（5'->3' で合成される分子）
        construct=f"5'-{m_down}-[PEG{lk['n_eo']}+SI]-{m_up}-3'",
        GC_up=round(gc_frac(site1), 3), GC_down=round(gc_frac(site2), 3),
        # 構造コスト
        dG_open_up=round(open1, 2), dG_open_down=round(open2, 2),
        dG_open_both=round(open12, 2),
        dG_self_up=round(self1, 2), dG_self_down=round(self2, 2),
        dG_cross=round(cross, 2),
        # 一価
        dG37_mono_up=round(dG_mono1, 2), dG37_mono_down=round(dG_mono2, 2),
        Tm_mono_up=round(tm_mono1, 1), Tm_mono_down=round(tm_mono2, 1),
        Kd_mono_up_uM=round(kd_mono1 * 1e6, 3), Kd_mono_down_uM=round(kd_mono2 * 1e6, 3),
        # リンカー
        n_EO=lk["n_eo"], n_EO_ceffmax=lk["n_eo_max"],
        linker_contour_A=round(lk["contour"], 1),
        span_req_A=round(lk["r_req"], 1), span_req_max_A=round(lk["r_req_max"], 1),
        C_eff_mM=round(lk["ceff"] * 1e3, 4), dG_link=round(lk["dG_link"], 2),
        eo_scan_at_bound=lk["at_bound"],
        # 二価
        dG37_biv=round(dG_biv, 2), Tm_biv=round(tm_biv, 1),
        Kd_biv_nM=round(kd_biv * 1e9, 4),
        # 性能
        occ_biv=round(occ_biv, 4),
        occ_leak=round(occ_leak, 4),
        occ_cut_coop=round(occ_cut_coop, 4),
        switch_log10=round(switch_log, 2),
        coax_flag=int(coax != 0.0),
        func_score=round(fs, 3),
        **{f"cov_{k}": round(v, 2) for k, v in fdet.items()},
    )
    return rec, ""


def passes_targets(r):
    lo, hi = CFG.TM_MONO_TARGET
    if not (lo <= r["Tm_mono_up"] <= hi):
        return False
    if not (lo <= r["Tm_mono_down"] <= hi):
        return False
    if r["Tm_biv"] < CFG.TM_BIV_MIN:
        return False
    if r["occ_biv"] < CFG.OCC_BIV_MIN:
        return False
    if r["occ_leak"] > CFG.OCC_MONO_MAX:
        return False
    if r["switch_log10"] < CFG.SWITCH_MIN_LOG:
        return False
    return True


def funnel(df, label=""):
    """各設計目標を1つずつ重ねたときの残存候補数を表示（どの条件が律速かを見る）。"""
    if df.empty:
        return
    lo, hi = CFG.TM_MONO_TARGET
    steps = [
        ("all candidates",            pd.Series(True, index=df.index)),
        (f"Tm_mono in [{lo},{hi}]",   df["Tm_mono_up"].between(lo, hi) &
                                      df["Tm_mono_down"].between(lo, hi)),
        (f"Tm_biv >= {CFG.TM_BIV_MIN}",   df["Tm_biv"] >= CFG.TM_BIV_MIN),
        (f"occ_biv >= {CFG.OCC_BIV_MIN}", df["occ_biv"] >= CFG.OCC_BIV_MIN),
        (f"leak <= {CFG.OCC_MONO_MAX}",   df["occ_leak"] <= CFG.OCC_MONO_MAX),
        (f"switch >= 10^{CFG.SWITCH_MIN_LOG}", df["switch_log10"] >= CFG.SWITCH_MIN_LOG),
    ]
    print(f"\n--- filter funnel {label} ---")
    cum = pd.Series(True, index=df.index)
    for name, m in steps:
        alone = int(m.sum())
        cum = cum & m
        print(f"  {name:<28} alone={alone:>5}   cumulative={int(cum.sum()):>5}")


def composite_score(df):
    """スイッチ窓・マスク効率・リーク・機能要素カバーの合成スコア"""
    if df.empty:
        return df
    d = df.copy()
    def nz(x):
        rng = x.max() - x.min()
        return (x - x.min()) / rng if rng > 1e-12 else pd.Series(0.5, index=x.index)
    d["_sw"]   = nz(d["switch_log10"])
    d["_mask"] = nz(d["occ_biv"])
    d["_leak"] = 1.0 - nz(d["occ_leak"])
    d["_func"] = nz(d["func_score"])
    d["_tmw"]  = 1.0 - nz((d[["Tm_mono_up", "Tm_mono_down"]].mean(axis=1)
                           - sum(CFG.TM_MONO_TARGET) / 2).abs())
    d["score"] = (0.35 * d["_sw"] + 0.20 * d["_mask"] + 0.20 * d["_leak"]
                  + 0.15 * d["_func"] + 0.10 * d["_tmw"]).round(4)
    return d.drop(columns=[c for c in d.columns if c.startswith("_")])


def stage12(L1=None, L2=None):
    L1 = CFG.L1 if L1 is None else L1
    L2 = CFG.L2 if L2 is None else L2
    lo, hi = CFG.MASK_REGION
    rows, rejects = [], {}
    starts1 = range(lo, hi - L1 + 2)
    for a1 in starts1:
        b1 = a1 + L1 - 1
        for a2 in range(b1 + 1 + CFG.MIN_GAP_NT, hi - L2 + 2):
            b2 = a2 + L2 - 1
            if b2 > hi:
                break
            rec, why = evaluate_pair(a1, b1, a2, b2)
            if rec is None:
                rejects[why] = rejects.get(why, 0) + 1
                continue
            rec["pass_targets"] = passes_targets(rec)
            rows.append(rec)
    return pd.DataFrame(rows), rejects


# =============================================================================
# 9. STAGE 3 : 上位候補について腕長を最適化
#    ギャップに面する内側末端を固定し、外側へ伸縮させる
# =============================================================================
def stage3(top_df):
    lo, hi = CFG.MASK_REGION
    out = []
    for _, r in top_df.head(CFG.STAGE3_TOP_K).iterrows():
        b1_fix = int(r["up_end"])      # 上流腕の 3'UTR 側端（ギャップに面する）
        a2_fix = int(r["down_start"])  # 下流腕の 5'UTR 側端（ギャップに面する）
        base_tag = f"{int(r['up_start'])}-{b1_fix}/{a2_fix}-{int(r['down_end'])}"
        for l1 in CFG.STAGE3_LENGTHS:
            a1 = b1_fix - l1 + 1
            if a1 < lo:
                continue
            for l2 in CFG.STAGE3_LENGTHS:
                b2 = a2_fix + l2 - 1
                if b2 > hi:
                    continue
                rec, why = evaluate_pair(a1, b1_fix, a2_fix, b2)
                if rec is None:
                    continue
                rec["parent"] = base_tag
                rec["L_up"], rec["L_down"] = l1, l2
                rec["pass_targets"] = passes_targets(rec)
                out.append(rec)
    return pd.DataFrame(out)


# =============================================================================
# 10. レポート
# =============================================================================
def show_target_structure():
    ss, mfe = mfe_fold(TARGET)
    fc = RNA.fold_compound(TARGET, _md())
    _, gens = fc.pf()
    print("=" * 100)
    print(f"TARGET ({N_T} nt)   T = {CFG.TEMP_C} C,  Na+eq = {na_equivalent():.3f} M")
    print("=" * 100)
    for off in range(0, N_T, 60):
        seg = slice(off, min(off + 60, N_T))
        print(f"{off+1:>4} {TARGET[seg]}")
        print(f"     {ss[seg]}")
    print(f"\nMFE          = {mfe:.2f} kcal/mol")
    print(f"ensemble dG  = {gens:.2f} kcal/mol")

    # 1 nt 単位の unpaired probability
    try:
        fc2 = RNA.fold_compound(TARGET, _md())
        fc2.pf()
        bppm = fc2.bpp()
        pu = []
        for i in range(1, N_T + 1):
            p = 0.0
            for j in range(1, N_T + 1):
                if i < j:
                    p += bppm[i][j]
                elif j < i:
                    p += bppm[j][i]
            pu.append(1.0 - p)
        print("\nunpaired probability (0-9 scale, 9 = most accessible):")
        for off in range(0, N_T, 60):
            seg = slice(off, min(off + 60, N_T))
            print(f"{off+1:>4} {TARGET[seg]}")
            print("     " + "".join(str(min(9, int(x * 10))) for x in pu[seg]))
    except Exception as e:
        print(f"(unpaired prob skipped: {e})")

    print("\nannotation:")
    for name, a, b, w in CFG.FUNC_ELEMENTS:
        print(f"  {name:<14} {a:>3}-{b:<3}  w={w}  {TARGET[a-1:b]}")
    print(f"  5'UTR          {CFG.UTR_RANGE[0]}-{CFG.UTR_RANGE[1]}")
    print(f"  CDS            {CFG.CDS_RANGE[0]}-{CFG.CDS_RANGE[1]}")


def print_design_targets():
    print("\n" + "=" * 100)
    print("DESIGN TARGETS (CFG で変更可)")
    print("=" * 100)
    print(f"  一価 Tm 目標窓        : {CFG.TM_MONO_TARGET[0]} - {CFG.TM_MONO_TARGET[1]} degC "
          f"(C_T = {CFG.C_TOTAL_TM*1e6:.1f} uM)")
    print(f"  二価 Tm 下限          : {CFG.TM_BIV_MIN} degC")
    print(f"  二価での標的占有率    : >= {CFG.OCC_BIV_MIN}   "
          f"(@ [mask]_load = {CFG.C_MASK_LOAD*1e6:.2f} uM)")
    print(f"  切断後の残存占有率    : <= {CFG.OCC_MONO_MAX}  "
          f"(@ [mask]_cell = {CFG.C_MASK_CELL*1e9:.1f} nM)")
    print(f"  スイッチ窓            : log10(Kd_mono/Kd_biv) >= {CFG.SWITCH_MIN_LOG}")
    print(f"  最小ギャップ          : {CFG.MIN_GAP_NT} nt  "
          f"(これ未満だと切断後に coaxial stack で協同再結合するリーク)")
    print(f"  2'-OMe 補正           : {CFG.DDG_OME_PER_NT} kcal/mol / 修飾残基")


MAIN_COLS = ["up_start", "up_end", "down_start", "down_end", "gap_nt",
             "mask_up", "mask_down", "GC_up", "GC_down",
             "dG_open_up", "dG_open_down", "dG_open_both",
             "dG_self_up", "dG_self_down", "dG_cross",
             "Tm_mono_up", "Tm_mono_down", "Kd_mono_up_uM", "Kd_mono_down_uM",
             "n_EO", "n_EO_ceffmax", "linker_contour_A", "span_req_A", "C_eff_mM", "dG_link",
             "dG37_biv", "Tm_biv", "Kd_biv_nM",
             "occ_biv", "occ_leak", "switch_log10", "func_score",
             "pass_targets", "score"]


def plots(pair_df, L1, L2):
    if not HAVE_PLT or pair_df.empty:
        return
    lo, hi = CFG.MASK_REGION
    s1 = sorted(pair_df["up_start"].unique())
    s2 = sorted(pair_df["down_start"].unique())
    idx1 = {v: i for i, v in enumerate(s1)}
    idx2 = {v: i for i, v in enumerate(s2)}

    for metric, title in [("score", "composite score"),
                          ("switch_log10", "log10(Kd_mono / Kd_biv)"),
                          ("Tm_biv", "bivalent apparent Tm (degC)"),
                          ("occ_leak", "post-cleavage leak (occupancy)")]:
        if metric not in pair_df.columns:
            continue
        M = np.full((len(s1), len(s2)), np.nan)
        for _, r in pair_df.iterrows():
            M[idx1[r["up_start"]], idx2[r["down_start"]]] = r[metric]
        fig, ax = plt.subplots(figsize=(9, 6))
        im = ax.imshow(M, aspect="auto", origin="lower", cmap="viridis",
                       extent=[min(s2) - .5, max(s2) + .5, min(s1) - .5, max(s1) + .5])
        ax.set_xlabel(f"downstream arm start (1-based, L={L2})")
        ax.set_ylabel(f"upstream arm start (1-based, L={L1})")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.show()

    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(pair_df["Tm_mono_up"].clip(-20, 90),
                    pair_df["Tm_biv"].clip(-20, 110),
                    c=pair_df["switch_log10"], cmap="plasma", s=18)
    ax.axvspan(*CFG.TM_MONO_TARGET, color="g", alpha=.12, label="target Tm_mono")
    ax.axhline(CFG.TM_BIV_MIN, color="r", ls="--", lw=1, label="Tm_biv min")
    ax.axhline(37, color="k", ls=":", lw=1, label="37 C")
    ax.axvline(37, color="k", ls=":", lw=1)
    ax.set_xlabel("Tm monovalent, upstream arm (degC)")
    ax.set_ylabel("Tm bivalent, apparent (degC)")
    ax.legend(fontsize=8)
    fig.colorbar(sc, ax=ax, label="log10 switch window")
    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(7, 5))
    g = pair_df.groupby("gap_nt").agg(sw=("switch_log10", "max"),
                                      ne=("n_EO", "median")).reset_index()
    ax.plot(g["gap_nt"], g["sw"], "o-", label="max log10 switch window")
    ax2 = ax.twinx()
    ax2.plot(g["gap_nt"], g["ne"], "s--", color="tab:orange", label="optimal PEG units")
    ax.set_xlabel("gap between the two sites (nt)")
    ax.set_ylabel("max log10(Kd_mono/Kd_biv)")
    ax2.set_ylabel("optimal n_EO (PEG units)")
    ax.set_title("gap / linker length trade-off")
    ax.grid(alpha=.3)
    plt.tight_layout()
    plt.show()


# =============================================================================
# 11. MAIN
# =============================================================================
def main():
    show_target_structure()
    print_design_targets()

    print("\n" + "=" * 100)
    print(f"STAGE 0 : single-window survey (L = {CFG.L1} nt)")
    print("=" * 100)
    s0 = stage0(CFG.L1)
    cols0 = ["start", "end", "masker", "GC", "Tm_ideal", "dG_open",
             "dG_self", "dG37_eff", "Tm_eff", "Kd_eff_uM",
             "offtarget_margin", "G_run"]
    print(s0[cols0].to_string(index=False))
    s0.to_csv(os.path.join(CFG.OUTDIR, "stage0_windows.csv"), index=False)

    print("\n" + "=" * 100)
    print(f"STAGE 1-2 : pairwise screen (L_up = {CFG.L1}, L_down = {CFG.L2})")
    print("=" * 100)
    pair, rej = stage12()
    print(f"evaluated pairs kept : {len(pair)}")
    if rej:
        print("rejected :", ", ".join(f"{k}={v}" for k, v in sorted(rej.items())))
    if pair.empty:
        print("\n候補がゼロです。CFG のフィルタ（GC_MIN/MAX, SELF_FOLD_MAX, "
              "OFFTARGET_MARGIN, MIN_GAP_NT）を緩めて再実行してください。")
        return

    pair = composite_score(pair)
    pair = pair.sort_values(["pass_targets", "score"], ascending=[False, False])
    pair.to_csv(os.path.join(CFG.OUTDIR, "stage12_pairs.csv"), index=False)

    funnel(pair, "(STAGE 1-2)")
    npass = int(pair["pass_targets"].sum())
    print(f"\nmeeting all design targets : {npass}")
    if int(pair["eo_scan_at_bound"].sum()) > 0:
        print(f"[warn] {int(pair['eo_scan_at_bound'].sum())} 件で推奨 PEG 長が EO_SCAN の端に "
              f"達しています -> CFG.EO_SCAN を広げてください")
    if npass == 0:
        print("\n[診断] 設計目標を満たす組合せがありません。上の funnel で律速条件を確認し、")
        print("       多くの場合は『一価が強すぎる』のが原因です → STAGE 3 で腕を短く、")
        print("       あるいは CFG.C_MASK_CELL を下げる / MIN_GAP_NT を小さくして")
        print("       C_eff を稼ぐ、のいずれかで窓が開きます。以下は目標未達でも上位の候補です。")
    print(f"\n--- TOP {CFG.TOP_N} ---")
    print(pair[MAIN_COLS].head(CFG.TOP_N).to_string(index=False))

    print("\n--- TOP 5 の詳細 ---")
    for _, r in pair.head(5).iterrows():
        print("-" * 90)
        print(f"sites      : up {r['up_start']}-{r['up_end']}  |  gap {r['gap_nt']} nt  |  "
              f"down {r['down_start']}-{r['down_end']}")
        print(f"target up  : 5'-{r['site_up']}-3'    GC={r['GC_up']}")
        print(f"target down: 5'-{r['site_down']}-3'  GC={r['GC_down']}")
        print(f"mask up    : 5'-{r['mask_up']}-3'   (all 2'-OMe)")
        print(f"mask down  : 5'-{r['mask_down']}-3' (all 2'-OMe)")
        print(f"construct  : {r['construct']}")
        print(f"             PEG {r['n_EO']} units (contour {r['linker_contour_A']} A), "
              f"required span {r['span_req_A']} A "
              f"(extended-gap upper bound {r['span_req_max_A']} A)")
        print(f"             C_eff = {r['C_eff_mM']} mM  ->  dG_avidity = {r['dG_link']} kcal/mol")
        print(f"monovalent : dG37 {r['dG37_mono_up']} / {r['dG37_mono_down']} kcal/mol, "
              f"Tm {r['Tm_mono_up']} / {r['Tm_mono_down']} C, "
              f"Kd {r['Kd_mono_up_uM']} / {r['Kd_mono_down_uM']} uM")
        print(f"bivalent   : dG37 {r['dG37_biv']} kcal/mol, Tm {r['Tm_biv']} C, "
              f"Kd {r['Kd_biv_nM']} nM")
        print(f"performance: occ_biv {r['occ_biv']}, leak {r['occ_leak']}, "
              f"switch 10^{r['switch_log10']}, func {r['func_score']}, "
              f"pass={bool(r['pass_targets'])}")

    plots(pair, CFG.L1, CFG.L2)

    if CFG.RUN_STAGE3:
        print("\n" + "=" * 100)
        print(f"STAGE 3 : arm-length optimization for top {CFG.STAGE3_TOP_K} geometries "
              f"(L = {min(CFG.STAGE3_LENGTHS)}-{max(CFG.STAGE3_LENGTHS)} nt)")
        print("=" * 100)
        s3 = stage3(pair)
        if s3.empty:
            print("STAGE 3 で候補なし。")
        else:
            s3 = composite_score(s3)
            s3 = s3.sort_values(["pass_targets", "score"], ascending=[False, False])
            s3.to_csv(os.path.join(CFG.OUTDIR, "stage3_lengths.csv"), index=False)
            c3 = ["parent", "L_up", "L_down", "up_start", "up_end",
                  "down_start", "down_end", "gap_nt", "mask_up", "mask_down",
                  "Tm_mono_up", "Tm_mono_down", "Tm_biv", "n_EO",
                  "Kd_mono_up_uM", "Kd_biv_nM", "occ_biv", "occ_leak",
                  "switch_log10", "func_score", "pass_targets", "score"]
            print(s3[c3].head(40).to_string(index=False))
            print(f"\nSTAGE3 passing design targets : {int(s3['pass_targets'].sum())}")

    print("\n" + "=" * 100)
    print("出力ファイル:")
    for f in ["stage0_windows.csv", "stage12_pairs.csv", "stage3_lengths.csv"]:
        p = os.path.join(CFG.OUTDIR, f)
        if os.path.exists(p):
            print("  " + p)
    print("=" * 100)
    print("""
解釈のポイント
--------------
* dG_open  : そのサイトを 1 本鎖にするコスト。大きいサイトは二次構造の幹に
             埋まっており、一価でも二価でも結合しにくい。
* dG_link  : アビディティ項 -RT ln(C_eff)。ギャップが小さいほど C_eff が大きく
             (= dG_link が小さく) 二価が有利 → スイッチ窓が広がる。
* occ_leak : 切断後の残存マスク。一価 2 本の独立結合と、隣接サイトでの
             coaxial stacking による協同再結合の大きい方を取っている。
* switch_log10 が窓の広さ。1 価/2 価の Tm 差ではなく Kd 比で見るのが本質。
* Tm は ΔG_open の ΔH を無視した「見かけ Tm」なので、設計判断は
  dG37 / Kd / occupancy を主、Tm を従にしてください。
""")


if __name__ == "__main__":
    main()
