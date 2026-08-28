# =============================================================================
#  IEDDA テンプレート連結シミュレーション v17
#  ―― C9orf72 (C4G2)n 標的・7-deaza-G C7 リンカー型 ASO の鋳型上ライゲーション
#
#  Google Colab にそのまま貼り付けて実行してください。
#  ランタイム: 標準 CPU で 20-35 分 (QUICK=True なら 4-7 分)
#
# -----------------------------------------------------------------------------
#  v16 からの主な変更点（なぜ変えたかは README セクション参照）
#
#  [1] 標的配列そのものの二重鎖を構築する
#      v16 は無関係な PDB 5 本を「A 型らしい」という理由だけで使っていた。
#      v17 は A 型 RNA らせんから (a) 塩基対ステップ変換 (b) G / C ヌクレオチド
#      の理想テンプレートを抽出し、標的配列
#         ASO 鎖   5'-CGGGGC CGGGGC-3'  (6mer ASO 2 本、間にニック)
#         鋳型鎖   3'-GCCCCG GCCCCG-5'  = (C4G2)n アンチセンス反復
#      の正則 A 型二重鎖を de novo 構築する。ニック（A6-B1 間のリン酸欠損）も再現。
#
#  [2] 反応判定を TS 類似の NAC (near-attack conformation) に厳格化
#      v16: 形成結合 4.0 Å 以内 + 環面から 1.2 Å 以上、だけ。角度拘束が皆無で
#      f_free = 6.08% という非物理的な値になっていた（実際の NAC 確率は 1e-3 以下）。
#      v17: 形成結合 2 本とも 2.2-3.8 Å、両アルケン炭素がテトラジン環面の同じ側で
#      面からの高さ 1.9 Å 以上、環法線とシクロプロペン環法線のなす角 <= 40°、
#      C3···C6 軸と C=C 軸のなす角 <= 45°、アルケン中点の環面内オフセット <= 1.6 Å。
#      さらに厳/中/緩の 3 段階で感度を出す。
#
#  [3] ハード足切りをやめ、ボルツマン重み付き統計に置き換え
#      v16: MMFF 配座を等重みで扱い、vdW 半径 0.75 倍の 0/1 判定で受理していた。
#      v17: 配座 MMFF エネルギー + RNA との相互作用エネルギー
#      (Lennard-Jones + Debye-Huckel 静電) から重みを作り、全量を重み付き平均。
#
#  [4] 静電相互作用を導入
#      テトラジン側リンカーの 2 級アミンは pH 7.4 でプロトン化して +1。
#      リン酸 (-1) との引力で主溝に張り付く。150 mM 塩の Debye-Huckel で評価。
#
#  [5] ペア列挙を厳密化・高速化
#      v16 は 500 万組の乱数サンプリング + Python ループ。v17 は KD-tree で
#      NAC になりうる組だけを列挙し、重み付き厳密平均を取る（統計誤差ゼロ）。
#
#  [6] 実際の分子設計を反映
#      6mer の 2 位と 5 位が両方修飾されている = 鋳型上で交互共重合する設計。
#      鎖内 (2->5) も接合部越え (5->2) も間隔 3 nt になる。反応する対のほかに
#      未反応ハンドル（傍観者アーム）が両隣にいるので、その立体障害も入れる。
#
#  [7] 熱ゆらぎアンサンブル
#      実測 A 型ステップ変換をブートストラップして 二重鎖を N_TEMPLATES 本生成。
#      結晶構造 5 本のばらつきではなく、標的配列のらせん熱ゆらぎで誤差を出す。
#
#  [8] 観測量まで接続
#      EM だけで半減期を出すのは誤り。6mer RNA 二重鎖は解離が速いので、
#      Turner 2004 最近接塩基対法で Kd を出し、三元複合体の寿命と化学反応速度の
#      競合（reactive-encounter モデル）から実際の擬一次速度と半減期を計算する。
# =============================================================================

import sys, subprocess
try:
    import rdkit
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "rdkit"], check=True)

import numpy as np, urllib.request, math, io, warnings
from collections import defaultdict
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")
import matplotlib
import matplotlib.pyplot as plt

# =============================================================================
#  設定
# =============================================================================
QUICK = False          # True にすると粗いサンプリングで数分で終わる

# --- 標的・設計 -------------------------------------------------------------
ASO_SEQ        = "CGGGGC"   # ASO (5'->3')。鋳型は (C4G2)n アンチセンス反復
MOD_POSITIONS  = (2, 5)     # 修飾する塩基の ASO 内 1-based 位置（どちらも G）
N_FLANK_BP     = 3          # 構築する二重鎖の両端に足す余分な塩基対
DELTA_SCAN     = [2, 3, 4, 5, 6]   # 反応する 2 ハンドル間の鎖内 nt 間隔の走査
TEMP_C         = 37.0       # 温度 (°C)

# --- 構造ソース -------------------------------------------------------------
# A 型 RNA らせんを含む高分解能構造。ここから G/C ヌクレオチド template と
# 塩基対ステップ変換を抽出する（配列は問わない。幾何だけ使う）。
HARVEST_PDB_IDS = ["1QC0", "1QCU", "157D", "397D", "280D", "1SDR", "353D", "1RNA",
                   "354D", "413D", "255D", "438D"]
PDB_LOCAL_FILES = []        # ネットが無い環境用: ローカル .pdb のパスを入れる

# --- アーム分子 -------------------------------------------------------------
# 三重結合クリック後の 1,2,3-トリアゾール経由。トリアゾール C4 上のメチル基が
# 7-deaza-G の C7 に置き換わる（= キャップ原子を C7 位置に重ねる）
SMI_TZ = "Cc1cn(C(C)CCNCc2ccc(C3=NN=C(C)N=N3)cc2)nn1"      # テトラジン側
SMI_CP = "Cc1cn(CCNC(=O)OCC2C=C2C)nn1"                      # シクロプロペン側

N_CONFS      = 2500 if not QUICK else 300
PRUNE_RMS    = 0.4
E_WINDOW     = 12.0      # kcal/mol。これを超える配座は捨てる
MAX_CONFS_USE= 600 if not QUICK else 120
N_TORSION    = 30 if not QUICK else 12   # C7-トリアゾール結合まわりの回転
N_TEMPLATES  = 16 if not QUICK else 4    # 熱ゆらぎ二重鎖の本数
N_TEMPL_SCAN = 6  if not QUICK else 2    # 修飾位置の走査に使う本数
N_FREE_MC    = 4_000_000 if not QUICK else 500_000

USE_BYSTANDERS = True     # 未反応の隣接ハンドルを立体障害として入れるか
N_BYS          = 3 if not QUICK else 2   # 傍観アームの配置をボルツマン分布から何点引くか

# --- 物理定数・力場 ---------------------------------------------------------
KCAL = 0.0019872041           # R (kcal/mol/K)
IONIC_M       = 0.15          # 塩濃度 (M)
DIELECTRIC    = 78.0
AMINE_CHARGE  = +1.0          # プロトン化 2 級アミン
PHOS_O_CHARGE = -0.5          # OP1 / OP2 それぞれ
LJ_CUTOFF     = 9.0
LJ_KMAX       = 48            # KD-tree で拾う近傍数
E_CLAMP       = 15.0          # 最良配置に対する相互作用エネルギーの上限 (kcal/mol)
LJ_REP_CAP    = 3.0           # 1 接触あたりの反発の飽和値 (kcal/mol)。RNA とアームの
                              # 局所緩和を剛体近似のまま取り込むソフトコア（0 で無効化）
ATTR_SCALE    = 0.20          # 引力項の減衰係数。真空 LJ の引力は水中では脱溶媒和で
                              # ほぼ相殺される。0.15-0.35 が妥当な範囲（感度を [11] で確認）
DIR_WOBBLE_DEG = 8.0          # C7 置換基の生える向きのゆらぎ（塩基のプロペラ・呼吸）

# united-atom 的に少し膨らませた Rmin/2 (Å) と ε (kcal/mol)
RMIN2 = {"C": 2.00, "N": 1.85, "O": 1.75, "P": 2.10, "S": 2.10}
EPS   = {"C": 0.10, "N": 0.16, "O": 0.19, "P": 0.20, "S": 0.20}

# --- NAC (反応可能配座) 判定 -----------------------------------------------
# (d_max, h_min, ang_normal_max, ang_axis_max, offset_max)
NAC_LEVELS = {
    "tight":  (3.4, 1.9, 30.0, 35.0, 1.2),
    "medium": (3.8, 1.9, 40.0, 45.0, 1.6),
    "loose":  (4.2, 1.9, 50.0, 55.0, 2.2),
}
NAC_DMIN   = 2.2
NAC_MAIN   = "medium"
#  EM の 3 通りの重み付け
#   steric : 立体的に許される配置を等重みで数える（幾何学的な到達可能性。統計が安定）
#   unif   : RNA との相互作用はボルツマン、配座エネルギーは等重み
#   boltz  : 配座も相互作用もボルツマン（ひずみを全部払った現実的な下限側の値）
WMODES     = ("steric", "unif", "boltz")
WMODE_MAIN = "boltz"
WMODE_GEOM = "steric"
E_ALLOW      = 5.0    # kcal/mol : 最良配置からこの範囲を「立体的に許される」とみなす
E_ALLOW_PAIR = 3.0    # kcal/mol : アーム同士の許容相互作用
CLASH_HARD   = 2.0    # Å : これより近い重原子接触があれば即座に棄却（高速化も兼ねる）
RUN_SENSITIVITY = True
R_FREE     = 6.0     # 非鋳型参照系の球半径 (Å)

SEED = 20260828
rng  = np.random.default_rng(SEED)
T_K  = TEMP_C + 273.15
RT   = KCAL * T_K
KAPPA = 0.329 * math.sqrt(IONIC_M)      # Debye 遮蔽 (1/Å) @25°C 近似
TO_M  = 1660.5391                        # 1 molecule/Å^3 -> M

print("=" * 92)
print("  IEDDA テンプレート連結シミュレーション v17")
print("=" * 92)
print(f"  RDKit {Chem.rdBase.rdkitVersion} / T = {TEMP_C:.0f} °C / QUICK = {QUICK}")
print(f"  標的: ASO 5'-{ASO_SEQ}-3' x2 (修飾 {MOD_POSITIONS}) on (C4G2)n")
print()

# =============================================================================
#  Part 1  A 型 RNA らせんの幾何を抽出し、標的配列の二重鎖を de novo 構築する
# =============================================================================

BB_ATOMS   = ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
BASE_ATOMS = {
    "G": ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
    "A": ["N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
    "C": ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    "U": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
}
PURINES = ("G", "A")
COMPL   = {"G": "C", "C": "G", "A": "U", "U": "A"}

_RESMAP = {"A": "A", "C": "C", "G": "G", "U": "U",
           "RA": "A", "RC": "C", "RG": "G", "RU": "U",
           "ADE": "A", "CYT": "C", "GUA": "G", "URA": "U", "URI": "U"}

def _norm_resname(rn):
    rn = rn.strip().upper()
    if rn in _RESMAP:
        return _RESMAP[rn]
    core = rn.rstrip("35").lstrip("DR")
    return _RESMAP.get(core, _RESMAP.get(core[-1:], None) if core else None)

def _norm_atom(name):
    n = name.strip().upper().replace("*", "'")
    return {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3"}.get(n, n)

def parse_pdb(text):
    """-> list of residues in file order: dict(chain, seq, base, atoms{name: xyz})"""
    order, store = [], {}
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        if len(ln) < 54 or ln[16] not in (" ", "A"):
            continue
        base = _norm_resname(ln[17:20])
        if base is None:
            continue
        name = _norm_atom(ln[12:16])
        elem = (ln[76:78].strip() or name[0]).upper()
        if elem == "H" or name.startswith("H") or name[:1].isdigit():
            continue
        try:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        key = (ln[21], ln[22:27])
        if key not in store:
            store[key] = dict(chain=ln[21], seq=key[1], base=base, atoms={})
            order.append(key)
        store[key]["atoms"].setdefault(name, xyz)
    return [store[k] for k in order]

# --- Watson-Crick 塩基対（水素結合を全本数チェック）--------------------------
WC_HB = {("G", "C"): [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
         ("A", "U"): [("N1", "N3"), ("N6", "O4")]}

def wc_paired(r1, r2, lo=2.5, hi=3.5):
    key = (r1["base"], r2["base"])
    if key in WC_HB:
        hb, a, b = WC_HB[key], r1, r2
    elif key[::-1] in WC_HB:
        hb, a, b = WC_HB[key[::-1]], r2, r1
    else:
        return False
    for x, y in hb:
        if x not in a["atoms"] or y not in b["atoms"]:
            return False
        d = np.linalg.norm(a["atoms"][x] - b["atoms"][y])
        if not (lo <= d <= hi):
            return False
    return True

def torsion(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))

def delta_torsion(r):
    a = r["atoms"]
    if not all(k in a for k in ("C5'", "C4'", "C3'", "O3'")):
        return np.nan
    return torsion(a["C5'"], a["C4'"], a["C3'"], a["O3'"]) % 360.0

def base_normal(r):
    idx = [n for n in BASE_ATOMS[r["base"]] if n in r["atoms"]]
    if len(idx) < 5:
        return None
    P = np.array([r["atoms"][n] for n in idx])
    return np.linalg.svd(P - P.mean(0))[2][2]

def pair_frame(r1, r2):
    """WC 塩基対の局所座標系。x: 溝方向 / y: C1'->C1' / z: らせん軸方向"""
    if "C1'" not in r1["atoms"] or "C1'" not in r2["atoms"]:
        return None
    o = 0.5 * (r1["atoms"]["C1'"] + r2["atoms"]["C1'"])
    y = r1["atoms"]["C1'"] - r2["atoms"]["C1'"]
    y /= np.linalg.norm(y)
    n1, n2 = base_normal(r1), base_normal(r2)
    if n1 is None or n2 is None:
        return None
    if np.dot(n1, n2) < 0:
        n2 = -n2
    z = n1 + n2
    z -= np.dot(z, y) * y
    nz = np.linalg.norm(z)
    if nz < 1e-6:
        return None
    z /= nz
    x = np.cross(y, z)
    return o, np.stack([x, y, z], axis=1)   # 列が基底ベクトル

def orient_run_frames(residues, run):
    """らせん軸方向 z を strand I の 5'->3' 向きにそろえる（SVD 法線の符号不定を除去）"""
    frames = []
    for (i, j) in run:
        f = pair_frame(residues[i], residues[j])
        if f is None:
            return None
        frames.append(f)
    out = []
    for k, (o, F) in enumerate(frames):
        ref = frames[k + 1][0] - o if k + 1 < len(frames) else o - frames[k - 1][0]
        if np.dot(F[:, 2], ref) < 0:
            F = np.stack([-F[:, 0], F[:, 1], -F[:, 2]], axis=1)   # y 軸まわり 180°
        out.append((o, F))
    return out

def screw_params(R, t):
    """剛体変換 (R, t) をらせん運動に分解: (らせん twist °, らせん rise Å, 軸までの距離 Å)"""
    ang = math.acos(max(-1.0, min(1.0, (np.trace(R) - 1.0) / 2.0)))
    if ang < 1e-6:
        return 0.0, float(np.linalg.norm(t)), 0.0
    ax = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    ax /= np.linalg.norm(ax)
    rise = float(np.dot(t, ax))
    # 回転成分の不動点（らせん軸の位置）
    tperp = t - rise * ax
    A = np.eye(3) - R
    c, *_ = np.linalg.lstsq(A + np.outer(ax, ax), tperp, rcond=None)
    c = c - np.dot(c, ax) * ax
    return math.degrees(ang) * np.sign(np.dot(ax, [0, 0, 1])), rise * np.sign(np.dot(ax, [0, 0, 1])), float(np.linalg.norm(c))

def quat_from_R(R):
    m = R
    t = np.trace(m)
    if t > 0:
        s = math.sqrt(t + 1.0) * 2
        q = np.array([0.25 * s, (m[2, 1] - m[1, 2]) / s, (m[0, 2] - m[2, 0]) / s, (m[1, 0] - m[0, 1]) / s])
    else:
        i = int(np.argmax(np.diag(m)))
        j, k = (i + 1) % 3, (i + 2) % 3
        s = math.sqrt(1.0 + m[i, i] - m[j, j] - m[k, k]) * 2
        q = np.zeros(4)
        q[0] = (m[k, j] - m[j, k]) / s
        q[i + 1] = 0.25 * s
        q[j + 1] = (m[j, i] + m[i, j]) / s
        q[k + 1] = (m[k, i] + m[i, k]) / s
    return q / np.linalg.norm(q)

def R_from_quat(q):
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])

RX180 = np.diag([1.0, -1.0, -1.0])   # 塩基対の擬 2 回軸 (strand I <-> II)

def fetch_pdb_text(pid):
    url = f"https://files.rcsb.org/download/{pid}.pdb"
    with urllib.request.urlopen(url, timeout=40) as r:
        return r.read().decode("utf-8", errors="ignore")

def find_helical_runs(residues, min_len=4):
    """連続する WC 塩基対の逆平行スタック（らせん）を探す"""
    n = len(residues)
    idx = {id(r): i for i, r in enumerate(residues)}
    # 骨格でつながっている隣接関係
    nxt = {}
    for i in range(n - 1):
        a, b = residues[i], residues[i + 1]
        if a["chain"] == b["chain"] and "O3'" in a["atoms"] and "P" in b["atoms"]:
            if np.linalg.norm(a["atoms"]["O3'"] - b["atoms"]["P"]) < 2.0:
                nxt[i] = i + 1
    prv = {v: k for k, v in nxt.items()}

    # WC ペア（近傍のみ調べる）
    c1 = []
    keep = []
    for i, r in enumerate(residues):
        if "C1'" in r["atoms"]:
            c1.append(r["atoms"]["C1'"]); keep.append(i)
    if len(keep) < 4:
        return []
    tree = cKDTree(np.array(c1))
    pairs = set()
    for a, b in tree.query_pairs(12.0):
        i, j = keep[a], keep[b]
        if wc_paired(residues[i], residues[j]):
            pairs.add((i, j)); pairs.add((j, i))

    runs, used, seen_runs = [], set(), set()
    for (i, j) in sorted(pairs):
        if (i, j) in used:
            continue
        # 逆平行な起点か: (prev(i), next(j)) が対でない
        pi, nj = prv.get(i), nxt.get(j)
        if pi is not None and nj is not None and (pi, nj) in pairs:
            continue
        run = [(i, j)]
        used.add((i, j))
        a, b = i, j
        while True:
            na, pb = nxt.get(a), prv.get(b)
            if na is None or pb is None or (na, pb) not in pairs:
                break
            run.append((na, pb)); used.add((na, pb))
            a, b = na, pb
        if len(run) >= min_len:
            key = frozenset(frozenset(p) for p in run)
            if key not in seen_runs:
                seen_runs.add(key)
                runs.append(run)
    return runs

def harvest_geometry(pdb_texts, verbose=True):
    """A 型 RNA らせんから (1) ステップ変換の集合 (2) 残基テンプレートを集める"""
    steps, templates, stats = [], defaultdict(list), []
    for tag, txt in pdb_texts:
        res = parse_pdb(txt)
        if len(res) < 8:
            if verbose: print(f"    {tag}: 残基不足"); continue
        runs = find_helical_runs(res)
        n_ok = 0
        for run in runs:
            members = [k for p in run for k in p]
            dels = [delta_torsion(res[k]) for k in members]
            dels = [d for d in dels if not np.isnan(d)]
            if not dels:
                continue
            med_delta = float(np.median(dels))
            if not (65.0 <= med_delta <= 105.0):     # C3'-endo (A 型) でなければ捨てる
                continue
            frames = orient_run_frames(res, run)
            if frames is None or len(frames) < 4:
                continue
            # ステップ変換
            local_steps, twists, rises, xdisp = [], [], [], []
            for k in range(len(frames) - 1):
                o0, F0 = frames[k]; o1, F1 = frames[k + 1]
                Rs = F0.T @ F1
                ts = F0.T @ (o1 - o0)
                tw, ri, xd = screw_params(Rs, ts)
                local_steps.append((Rs, ts)); twists.append(tw); rises.append(ri); xdisp.append(xd)
            if not (26.0 <= np.median(twists) <= 38.0 and 2.3 <= np.median(rises) <= 3.4):
                continue
            if np.median(xdisp) < 4.0:          # A 型は軸から大きくずれる (B 型は小さい)
                continue
            steps.extend(local_steps)
            n_ok += 1
            # 残基テンプレート（strand I 規約に統一）
            for si, (i, j) in enumerate(run):
                o, F = frames[si]
                for side, k in ((0, i), (1, j)):
                    r = res[k]
                    need = BB_ATOMS + BASE_ATOMS[r["base"]]
                    if not all(a in r["atoms"] for a in need):
                        continue
                    L = np.array([F.T @ (r["atoms"][a] - o) for a in need])
                    if side == 1:
                        L = L @ RX180.T
                    templates[r["base"]].append((need, L))
            stats.append((tag, len(run), med_delta, float(np.median(twists)),
                          float(np.median(rises)), float(np.median(xdisp))))
        if verbose:
            print(f"    {tag}: A 型らせん {n_ok} 本 / 全 WC らせん {len(runs)} 本")
    return steps, templates, stats

def medoid_template(entries):
    """同一塩基種の局所座標群から代表（medoid）を選ぶ"""
    names = entries[0][0]
    arrs = np.stack([L for nm, L in entries if nm == names])
    d = np.sqrt(((arrs[:, None, :, :] - arrs[None, :, :, :]) ** 2).sum(-1).mean(-1))
    i = int(np.argmin(d.sum(1)))
    return names, arrs[i], float(np.median(d[i]))

def symmetrize_step(R, t):
    """二重鎖の擬 2 回対称 (S M S = M^-1, S = RX180) を満たすようステップ変換を対称化"""
    S = RX180
    R2 = S @ R.T @ S
    t2 = -(S @ (R.T @ t))
    q = quat_from_R(R); q2 = quat_from_R(R2)
    if q @ q2 < 0:
        q2 = -q2
    return R_from_quat(q + q2), 0.5 * (t + t2)

def mean_step(steps):
    """ステップ変換の平均（回転はクォータニオン平均）"""
    qs = np.stack([quat_from_R(R) for R, t in steps])
    qs *= np.sign(qs @ qs[0])[:, None]
    q = qs.mean(0)
    return symmetrize_step(R_from_quat(q), np.stack([t for R, t in steps]).mean(0))

# --- 二重鎖の構築 -----------------------------------------------------------
class Duplex:
    """strand I (ASO 鎖) と strand II (鋳型鎖) からなる正則 A 型二重鎖"""
    __slots__ = ("res1", "res2", "atoms", "elems", "charges", "nbp")

def build_duplex(seq1, step_list, tmpl, nick_after=None, deaza_pos=(), rng=None):
    """seq1: strand I の 5'->3' 配列。step_list: 使用するステップ変換のリスト
       nick_after: strand I の 1-based 位置。この後ろでリン酸を落とす（ニック）
       deaza_pos: 7-deaza 化する strand I の 1-based 位置（N7 -> C7）"""
    nbp = len(seq1)
    F = [(np.zeros(3), np.eye(3))]
    for k in range(nbp - 1):
        Rs, ts = step_list[k]
        o0, F0 = F[-1]
        F.append((o0 + F0 @ ts, F0 @ Rs))

    dx = Duplex()
    dx.res1, dx.res2 = [], []
    dx.nbp = nbp
    for k in range(nbp):
        o, Fk = F[k]
        b1 = seq1[k]
        b2 = COMPL[b1]
        names1, L1 = tmpl[b1]
        names2, L2 = tmpl[b2]
        X1 = (L1 @ Fk.T) + o
        X2 = ((L2 @ RX180.T) @ Fk.T) + o
        dx.res1.append(dict(base=b1, atoms={n: x for n, x in zip(names1, X1)}, pos=k + 1))
        dx.res2.append(dict(base=b2, atoms={n: x for n, x in zip(names2, X2)}, pos=k + 1))

    # 末端 / ニックのリン酸を落とす
    for strand in (dx.res1, dx.res2):
        for a in ("P", "OP1", "OP2"):
            strand[0]["atoms"].pop(a, None)
    if nick_after is not None and 1 <= nick_after < nbp:
        for a in ("P", "OP1", "OP2"):
            dx.res1[nick_after]["atoms"].pop(a, None)

    # 7-deaza 化（N7 を C7 に置換）
    for p in deaza_pos:
        r = dx.res1[p - 1]
        if "N7" in r["atoms"]:
            r["atoms"]["C7"] = r["atoms"].pop("N7")
    return dx

def duplex_arrays(dx, exclude=()):
    """(座標, 元素, 電荷, タグ) を平坦化。exclude = {(strand, pos, atomname)}"""
    X, E, Q, TAG = [], [], [], []
    for si, strand in ((0, dx.res1), (1, dx.res2)):
        for r in strand:
            for n, x in r["atoms"].items():
                if (si, r["pos"], n) in exclude:
                    continue
                X.append(x)
                E.append("C" if n.startswith("C") else n[0])
                Q.append(PHOS_O_CHARGE if n in ("OP1", "OP2") else 0.0)
                TAG.append((si, r["pos"], n))
    return np.array(X), np.array(E), np.array(Q), TAG

def write_model_pdb(path, dx, arms=()):
    """構築した二重鎖（＋アーム配置）を PDB として保存。PyMOL 等で確認できる"""
    lines, n = [], 0
    for ch, strand in (("A", dx.res1), ("B", dx.res2)):
        for r in strand:
            for name, x in r["atoms"].items():
                n += 1
                nm = name if len(name) == 4 else f" {name:<3s}"
                lines.append(f"ATOM  {n:5d} {nm}{r['base']:>4s} {ch}{r['pos']:4d}    "
                             f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          "
                             f"{('C' if name.startswith('C') else name[0]):>2s}")
    for k, (X, E, tag) in enumerate(arms):
        cnt = {}
        for x, e in zip(X, E):
            n += 1
            cnt[e] = cnt.get(e, 0) + 1
            an = f"{e}{cnt[e]}"[:3]
            lines.append(f"HETATM{n:5d}  {an:<3s}{tag[:3]:>4s} {chr(ord('L')+k)}{k+1:4d}    "
                         f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {e:>2s}")
    open(path, "w").write("\n".join(lines) + "\nEND\n")

def exocyclic_C7_dir(r):
    """7-deaza-G の C7 から環外へ伸びる方向（置換基が生える向き）"""
    a = r["atoms"]
    C7 = a["C7"] if "C7" in a else a["N7"]
    d = 2 * C7 - a["C5"] - a["C8"]
    return C7, d / np.linalg.norm(d)

# --- 骨格を閉じるための剛体微調整 -------------------------------------------
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial.transform import Rotation as _Rot

_BOND_O3P, _ANG_C3P, _ANG_O3O5 = 1.607, 2.64, 2.52
_CLASH_D, _CLASH_TARGET = 2.70, 2.85

def refine_duplex(dx, w_pos=0.6, w_bond=6.0, w_ang=2.0, w_clash=4.0):
    """テンプレート medoid と平均ステップの不整合で生じる骨格のずれを
       剛体微動（残基ごとに 6 自由度）で解消する。塩基位置は弱い拘束で保持。"""
    residues = list(dx.res1) + list(dx.res2)
    names = [sorted(r["atoms"]) for r in residues]
    X0 = [np.array([r["atoms"][n] for n in nm]) for r, nm in zip(residues, names)]
    cen = [x.mean(0) for x in X0]
    n1 = len(dx.res1)
    idx_of = {}
    for i, r in enumerate(residues):
        idx_of[(0 if i < n1 else 1, r["pos"])] = i

    links = []   # (i_donor(O3'), i_acceptor(P))
    for k in range(n1 - 1):
        if "P" in dx.res1[k + 1]["atoms"]:
            links.append((idx_of[(0, k + 1)], idx_of[(0, k + 2)]))
    for k in range(len(dx.res2) - 1):
        if "P" in dx.res2[k]["atoms"]:
            links.append((idx_of[(1, k + 2)], idx_of[(1, k + 1)]))
    pos_of = [{n: t for t, n in enumerate(nm)} for nm in names]

    # 残基間の近すぎる非結合接触（テンプレート合成で生じる 2.3 Å 級の重なり）を解く
    flat, own = [], []
    for i, x in enumerate(X0):
        for t in range(len(x)):
            flat.append(x[t]); own.append((i, t))
    flat = np.array(flat)
    tr0 = cKDTree(flat)
    clashes = []
    for a, b in tr0.query_pairs(_CLASH_D):
        ia_, ta = own[a]; ib_, tb = own[b]
        if ia_ == ib_:
            continue
        if any((ia_, ib_) == (u, v) or (ia_, ib_) == (v, u) for u, v in links):
            continue
        clashes.append((ia_, ta, ib_, tb))

    def transform(par):
        out = []
        for i, x in enumerate(X0):
            p = par[6 * i:6 * i + 6]
            R = _Rot.from_rotvec(p[:3]).as_matrix()
            out.append((x - cen[i]) @ R.T + cen[i] + p[3:])
        return out

    def resid(par):
        Xs = transform(par)
        r = [w_pos * (Xs[i] - X0[i]).ravel() for i in range(len(X0))]
        lr = []
        for a, b in links:
            if "O3'" not in pos_of[a] or "P" not in pos_of[b]:
                continue
            o3 = Xs[a][pos_of[a]["O3'"]]; c3 = Xs[a][pos_of[a]["C3'"]]
            pp = Xs[b][pos_of[b]["P"]];   o5 = Xs[b][pos_of[b]["O5'"]]
            lr += [w_bond * (np.linalg.norm(o3 - pp) - _BOND_O3P),
                   w_ang * (np.linalg.norm(c3 - pp) - _ANG_C3P),
                   w_ang * (np.linalg.norm(o3 - o5) - _ANG_O3O5)]
        cr = [w_clash * min(0.0, np.linalg.norm(Xs[a][ta] - Xs[b][tb]) - _CLASH_TARGET)
              for a, ta, b, tb in clashes]
        return np.concatenate(r + [np.array(lr), np.array(cr)])

    n_par = 6 * len(residues)
    nres = resid(np.zeros(n_par)).size
    S = lil_matrix((nres, n_par), dtype=int)
    off = 0
    for i, x in enumerate(X0):
        S[off:off + x.size, 6 * i:6 * i + 6] = 1
        off += x.size
    for a, b in links:
        S[off:off + 3, 6 * a:6 * a + 6] = 1
        S[off:off + 3, 6 * b:6 * b + 6] = 1
        off += 3
    for a, ta, b, tb in clashes:
        S[off, 6 * a:6 * a + 6] = 1
        S[off, 6 * b:6 * b + 6] = 1
        off += 1
    sol = least_squares(resid, np.zeros(n_par), jac_sparsity=S.tocsr(),
                        method="trf", xtol=1e-8, max_nfev=120, verbose=0)
    Xs = transform(sol.x)
    for i, r in enumerate(residues):
        for n, x in zip(names[i], Xs[i]):
            r["atoms"][n] = x
    return dx

def duplex_quality(dx):
    ds = []
    for k in range(len(dx.res1) - 1):
        if "P" in dx.res1[k + 1]["atoms"]:
            ds.append(np.linalg.norm(dx.res1[k]["atoms"]["O3'"] - dx.res1[k + 1]["atoms"]["P"]))
    for k in range(len(dx.res2) - 1):
        if "P" in dx.res2[k]["atoms"]:
            ds.append(np.linalg.norm(dx.res2[k + 1]["atoms"]["O3'"] - dx.res2[k]["atoms"]["P"]))
    WCA = {"G": ("N1", "N3"), "A": ("N1", "N3"), "C": ("N3", "N1"), "U": ("N3", "N1")}
    wc = [np.linalg.norm(dx.res1[k]["atoms"][WCA[dx.res1[k]["base"]][0]] -
                         dx.res2[k]["atoms"][WCA[dx.res1[k]["base"]][1]]) for k in range(dx.nbp)]
    X, E, Q, TAG = duplex_arrays(dx)
    tr = cKDTree(X)
    bad = sum(1 for i, j in tr.query_pairs(2.4)
              if TAG[i][0] != TAG[j][0] or abs(TAG[i][1] - TAG[j][1]) > 1)
    return dict(o3p=(float(np.min(ds)), float(np.max(ds)), float(np.mean(ds))),
                wc=(float(np.min(wc)), float(np.max(wc))), clashes=bad)

# =============================================================================
#  Part 2  アーム分子の配座アンサンブルと、主溝内での配置・エネルギー評価
# =============================================================================

P_TZS  = Chem.MolFromSmarts("[#6]1~[#7]~[#7]~[#6]~[#7]~[#7]1")   # 1,2,4,5-テトラジン
P_CPS  = Chem.MolFromSmarts("[#6]1=[#6][#6]1")                    # シクロプロペン
P_TRI  = Chem.MolFromSmarts("[#6]1~[#6]~[#7]~[#7]~[#7]1")         # 1,2,3-トリアゾール
P_AMIN = Chem.MolFromSmarts("[NX3;H1;!$(N[C,S]=[O,S,N]);!$(N~[#7,#8])]")  # プロトン化 2 級アミン

V2_TORSION = 2.0   # kcal/mol : C7-トリアゾール結合の 2 回対称ねじれ障壁（共平面が有利）

class Arm:
    pass

def build_arm(smi, label):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    p = AllChem.ETKDGv3()
    p.useMacrocycleTorsions = True; p.useSmallRingTorsions = True
    p.randomSeed = SEED; p.numThreads = 0; p.pruneRmsThresh = PRUNE_RMS
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_CONFS, params=p))
    props = AllChem.MMFFGetMoleculeProperties(mol)
    try:
        props.SetMMFFDielectricModel(2); props.SetMMFFDielectricConstant(DIELECTRIC)
    except Exception:
        pass
    E = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        ff.Minimize(maxIts=3000)
        E.append(ff.CalcEnergy())
    E = np.array(E); E -= E.min()
    keep = np.where(E <= E_WINDOW)[0]
    if len(keep) > MAX_CONFS_USE:      # エネルギー順ではなく無作為に間引く
        keep = np.sort(rng.choice(keep, MAX_CONFS_USE, replace=False))   # 配座の多様性を保つ
    P = np.stack([np.array(mol.GetConformer(cids[i]).GetPositions()) for i in keep])
    Ek = E[keep]

    a = Arm()
    a.mol, a.label, a.pos, a.Econf = mol, label, P, Ek
    a.wconf = np.exp(-(Ek - Ek.min()) / RT)

    # C7 に置き換わるキャップ（トリアゾール上のメチル基）とその結合先
    tri = set(mol.GetSubstructMatch(P_TRI))
    cap = None
    for at in mol.GetAtoms():
        if at.GetSymbol() != "C" or at.GetIdx() in tri:
            continue
        hv = [n for n in at.GetNeighbors() if n.GetAtomicNum() > 1]
        if len(hv) == 1 and hv[0].GetIdx() in tri and hv[0].GetSymbol() == "C":
            cap = (at.GetIdx(), hv[0].GetIdx()); break
    if cap is None:
        raise RuntimeError(f"{label}: C7 結合点（トリアゾールのメチル）が同定できません")
    a.cap = cap
    a.tri_ref = [n.GetIdx() for n in mol.GetAtomWithIdx(cap[1]).GetNeighbors()
                 if n.GetIdx() in tri][0]

    # 反応部位
    if mol.HasSubstructMatch(P_TZS):
        ring = list(mol.GetSubstructMatch(P_TZS))
        a.kind = "tz"
        a.ring = ring
        a.rC = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "C"]
        core = set(ring)
    else:
        m3 = list(mol.GetSubstructMatch(P_CPS))
        a.kind = "cp"
        a.ring = m3
        a.rC = [i for i in m3 if mol.GetAtomWithIdx(i).GetHybridization() == Chem.HybridizationType.SP2]
        core = set(m3)
    a.core = sorted(core)

    # 立体判定に使う重原子（キャップ = C7 は RNA 側の原子なので除く）
    drop = {cap[0]} | {n.GetIdx() for n in mol.GetAtomWithIdx(cap[0]).GetNeighbors() if n.GetAtomicNum() == 1}
    a.heavy = [x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum() > 1 and x.GetIdx() not in drop]
    a.rmin = np.array([RMIN2.get(mol.GetAtomWithIdx(i).GetSymbol(), 2.0) for i in a.heavy])
    a.eps  = np.array([EPS.get(mol.GetAtomWithIdx(i).GetSymbol(), 0.12) for i in a.heavy])
    # アーム間の立体判定では反応部位とその隣接は除く（結合形成中なので）
    excl = set(a.core)
    for i in a.core:
        excl |= {n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() > 1}
    a.inter = [k for k, i in enumerate(a.heavy) if i not in excl]

    # 形式電荷（プロトン化 2 級アミン）
    q = np.zeros(len(a.heavy))
    hmap = {i: k for k, i in enumerate(a.heavy)}
    for (n,) in mol.GetSubstructMatches(P_AMIN):
        if n in hmap:
            q[hmap[n]] = AMINE_CHARGE
    a.q = q
    print(f"  [{label}] {Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi))}  "
          f"配座 {len(cids)} -> 採用 {len(keep)} (E <= {E_WINDOW:.0f} kcal/mol)  "
          f"重原子 {len(a.heavy)}  正電荷 {int((q > 0).sum())}")
    return a

def rot_between(a, b):
    v = np.cross(a, b); c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else -np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))

def rot_axis(ax, th):
    ax = ax / np.linalg.norm(ax)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * (K @ K)

def dihedral_batch(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    return np.arctan2((np.cross(b1n, v) * w).sum(-1), (v * w).sum(-1))

COULOMB = 332.0637

def softcore(E):
    """反発は 1 接触あたり LJ_REP_CAP kcal/mol で飽和（局所緩和の近似）、
       引力は ATTR_SCALE 倍に減衰（脱溶媒和で相殺される分）。"""
    E = np.where(E < 0, ATTR_SCALE * E, E)
    if LJ_REP_CAP <= 0:
        return E
    return np.where(E > 0, LJ_REP_CAP * E / (E + LJ_REP_CAP), E)

def _lj_energy(Xh, rmin_a, eps_a, tree, rmin_r, eps_r, block=400):
    """アーム配置 (M, na, 3) と RNA の Lennard-Jones 相互作用エネルギー (M,)
       まず最近接 1 点だけ調べて明白な衝突を落とし、残りだけ本計算する"""
    M, na, _ = Xh.shape
    out = np.zeros(M)
    n_rna = len(rmin_r)
    d1, _i1 = tree.query(Xh.reshape(-1, 3), k=1, workers=-1)
    hard = (d1.reshape(M, na) < CLASH_HARD).any(axis=1)
    out[hard] = 1e4
    keep = np.where(~hard)[0]
    Xh = Xh[keep]
    M = len(keep)
    sub = np.zeros(M)
    for s in range(0, M, block):
        e = min(s + block, M)
        flat = Xh[s:e].reshape(-1, 3)
        d, idx = tree.query(flat, k=LJ_KMAX, distance_upper_bound=LJ_CUTOFF, workers=-1)
        d = d.reshape(e - s, na, LJ_KMAX); idx = idx.reshape(e - s, na, LJ_KMAX)
        ok = np.isfinite(d) & (idx < n_rna)
        idc = np.where(ok, idx, 0)
        dd = np.where(ok, d, 1e3)
        rm = rmin_a[None, :, None] + rmin_r[idc]
        ep = np.sqrt(eps_a[None, :, None] * eps_r[idc])
        x = (rm / np.maximum(dd, 0.7)) ** 6
        ev = np.clip(ep * (x * x - 2.0 * x), -5.0, 1e6)
        ev = softcore(ev)
        sub[s:e] = np.where(ok, ev, 0.0).sum(axis=(1, 2))
    out[keep] = sub
    return out

def _dh_energy(Xh, q_a, Xq_r, q_r):
    """Debye-Huckel 静電（プロトン化アミン <-> リン酸酸素）"""
    ia = np.where(q_a != 0)[0]
    if len(ia) == 0 or len(q_r) == 0:
        return np.zeros(len(Xh))
    D = np.linalg.norm(Xh[:, ia, None, :] - Xq_r[None, None, :, :], axis=-1)
    D = np.maximum(D, 1.5)
    E = COULOMB * (q_a[ia][None, :, None] * q_r[None, None, :]) * np.exp(-KAPPA * D) / (DIELECTRIC * D)
    return E.sum(axis=(1, 2))

class Placement:
    pass

def place_arm(arm, anchor, direction, ref_atom, rna):
    """C7 (= anchor) に arm を生やし、全配座 x 全ねじれ角の配置を作って重み付けする。
       rna = dict(tree, rmin, eps, xq, q, xyz)"""
    cap0, cap1 = arm.cap
    hmap = {i: k for k, i in enumerate(arm.heavy)}
    P = arm.pos
    bases = []
    for p in P:
        v0 = p[cap1] - p[cap0]; v0 /= np.linalg.norm(v0)
        bases.append((p - p[cap0]) @ rot_between(v0, direction).T)
    B = np.stack(bases)                                   # (nc, na, 3)
    th = (np.arange(N_TORSION) + 0.5) * 2 * np.pi / N_TORSION
    Rs = np.stack([rot_axis(direction, t) for t in th])   # (nt, 3, 3)
    C = np.einsum("rij,caj->crai", Rs, B).reshape(-1, B.shape[1], 3)
    if DIR_WOBBLE_DEG > 0:      # 塩基のプロペラ・呼吸による結合方向のゆらぎ
        n = len(C)
        v = rng.normal(size=(n, 3))
        v -= (v * direction).sum(1)[:, None] * direction
        v /= np.linalg.norm(v, axis=1)[:, None]
        ang = np.abs(rng.normal(0.0, math.radians(DIR_WOBBLE_DEG), n))
        K = np.zeros((n, 3, 3))
        K[:, 0, 1] = -v[:, 2]; K[:, 0, 2] = v[:, 1]; K[:, 1, 0] = v[:, 2]
        K[:, 1, 2] = -v[:, 0]; K[:, 2, 0] = -v[:, 1]; K[:, 2, 1] = v[:, 0]
        Rw = np.eye(3) + np.sin(ang)[:, None, None] * K + (1 - np.cos(ang))[:, None, None] * (K @ K)
        C = np.einsum("nij,naj->nai", Rw, C)
    C = C + anchor
    Xh = C[:, arm.heavy, :]
    nc, nt = B.shape[0], N_TORSION
    wconf = np.repeat(arm.wconf, nt)

    E = _lj_energy(Xh, arm.rmin, arm.eps, rna["tree"], rna["rmin"], rna["eps"])
    E = E + _dh_energy(Xh, arm.q, rna["xq"], rna["q"])
    # C7-トリアゾール結合まわりのねじれポテンシャル（共平面が有利）
    phi = dihedral_batch(np.broadcast_to(ref_atom, (len(C), 3)),
                         np.broadcast_to(anchor, (len(C), 3)),
                         C[:, cap1, :], C[:, arm.tri_ref, :])
    E = E + 0.5 * V2_TORSION * (1 - np.cos(2 * phi))

    E = np.minimum(E, E.min() + E_CLAMP)
    dE = E - E.min()
    bf = np.exp(-dE / RT)
    pl = Placement()
    pl.arm, pl.X, pl.E = arm, Xh, E
    pl.W = {"boltz": wconf * bf, "unif": bf, "steric": (dE <= E_ALLOW).astype(float)}
    w = pl.W["boltz"]
    ring = [hmap[i] for i in arm.ring]
    R3 = Xh[:, ring, :]
    pl.cen = R3.mean(1)
    if arm.kind == "tz":
        n = np.cross(R3[:, 2] - R3[:, 0], R3[:, 4] - R3[:, 0])
    else:
        n = np.cross(R3[:, 1] - R3[:, 0], R3[:, 2] - R3[:, 0])
    pl.nrm = n / np.linalg.norm(n, axis=1, keepdims=True)
    rc = [hmap[i] for i in arm.rC]
    pl.C1, pl.C2 = Xh[:, rc[0], :], Xh[:, rc[1], :]
    pl.mid = 0.5 * (pl.C1 + pl.C2)
    # 参考量: 主溝内で残る配座自由度（Kish の有効サンプル数比）
    pl.freedom = float((w.sum() ** 2 / (w ** 2).sum()) / len(w))
    pl.occupancy = float(bf.mean())                 # 自由な状態に対する相対分配関数
    pl.allowed = float(pl.W["steric"].mean())       # 立体的に許される配置の割合
    pl.Emin = float(E.min())
    pl.n = len(E)
    return pl

def nac_mask(Q, N, Ct1, Ct2, Ca, Cb, Mn, level):
    d_max, h_min, ang_n, ang_ax, off_max = NAC_LEVELS[level]
    d11 = np.linalg.norm(Ct1 - Ca, axis=-1); d22 = np.linalg.norm(Ct2 - Cb, axis=-1)
    d12 = np.linalg.norm(Ct1 - Cb, axis=-1); d21 = np.linalg.norm(Ct2 - Ca, axis=-1)
    mx1, mx2 = np.maximum(d11, d22), np.maximum(d12, d21)
    mn1, mn2 = np.minimum(d11, d22), np.minimum(d12, d21)
    use1 = mx1 <= mx2
    dmx = np.where(use1, mx1, mx2); dmn = np.where(use1, mn1, mn2)
    ok = (dmx <= d_max) & (dmn >= NAC_DMIN)
    ha = ((Ca - Q) * N).sum(-1); hb = ((Cb - Q) * N).sum(-1)
    ok &= (np.sign(ha) == np.sign(hb)) & (np.minimum(np.abs(ha), np.abs(hb)) >= h_min)
    ok &= np.abs((N * Mn).sum(-1)) >= math.cos(math.radians(ang_n))
    ax1 = Ct2 - Ct1; ax2 = Cb - Ca
    cosax = np.abs((ax1 * ax2).sum(-1)) / (np.linalg.norm(ax1, axis=-1) * np.linalg.norm(ax2, axis=-1))
    ok &= cosax >= math.cos(math.radians(ang_ax))
    M = 0.5 * (Ca + Cb); dv = M - Q
    inp = dv - ((dv * N).sum(-1))[..., None] * N
    ok &= np.linalg.norm(inp, axis=-1) <= off_max
    return ok

def interarm_pars(armA, armB):
    rm = armA.rmin[armA.inter][:, None] + armB.rmin[armB.inter][None, :]
    ep = np.sqrt(armA.eps[armA.inter][:, None] * armB.eps[armB.inter][None, :])
    return rm, ep

def interarm_E(A, B, rm, ep, block=20000):
    """反応配座での 2 本のアーム同士の相互作用エネルギー。A,B は (n, k, 3)"""
    out = np.empty(len(A))
    for s in range(0, len(A), block):
        e = min(s + block, len(A))
        D = np.linalg.norm(A[s:e][:, :, None, :] - B[s:e][:, None, :, :], axis=-1)
        x = (rm[None] / np.maximum(D, 0.7)) ** 6
        out[s:e] = softcore(np.clip(ep[None] * (x * x - 2 * x), -5.0, 1e6)).sum(axis=(1, 2))
    return np.minimum(out, E_CLAMP)

# =============================================================================
#  Part 3  非鋳型（自由溶液）参照系と、実効モル濃度 EM
# =============================================================================

NAC_SEARCH_R = 4.9

def local_geom(arm):
    """配座ごとの反応部位の局所幾何（反応中心を原点に取る）"""
    hmap = {i: k for k, i in enumerate(arm.heavy)}
    ring = [hmap[i] for i in arm.ring]
    rc = [hmap[i] for i in arm.rC]
    Xh = arm.pos[:, arm.heavy, :]
    R3 = Xh[:, ring, :]
    if arm.kind == "tz":
        origin = R3.mean(1)
        n = np.cross(R3[:, 2] - R3[:, 0], R3[:, 4] - R3[:, 0])
    else:
        origin = 0.5 * (Xh[:, rc[0], :] + Xh[:, rc[1], :])
        n = np.cross(R3[:, 1] - R3[:, 0], R3[:, 2] - R3[:, 0])
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    return dict(c1=Xh[:, rc[0], :] - origin, c2=Xh[:, rc[1], :] - origin,
                nrm=n, cen=R3.mean(1) - origin, inter=Xh[:, arm.inter, :] - origin[:, None, :])

def random_rot(n):
    u1, u2, u3 = rng.random(n), rng.random(n), rng.random(n)
    x = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2); y = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    z = np.sqrt(u1) * np.sin(2 * np.pi * u3);     w = np.sqrt(u1) * np.cos(2 * np.pi * u3)
    R = np.empty((n, 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w); R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w); R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w); R[:, 2, 1] = 2 * (y * z + x * w); R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R

def free_reference(arm_tz, arm_cp, n_mc=None, chunk=400_000):
    """自由溶液中で NAC になる確率 f_free。EM の参照状態。
       戻り値 {wmode: {level: f}}"""
    n_mc = n_mc or N_FREE_MC
    gt, gc = local_geom(arm_tz), local_geom(arm_cp)
    pt = arm_tz.wconf / arm_tz.wconf.sum()
    pu = np.full(len(pt), 1.0 / len(pt))
    pc = arm_cp.wconf / arm_cp.wconf.sum()
    pcu = np.full(len(pc), 1.0 / len(pc))
    rm, ep = interarm_pars(arm_tz, arm_cp)
    acc = {wm: {lv: 0.0 for lv in NAC_LEVELS} for wm in WMODES}
    done = 0
    while done < n_mc:
        m = min(chunk, n_mc - done); done += m
        ic = rng.integers(0, len(pt), m); jc = rng.integers(0, len(pc), m)   # 一様に引いて重みで補正
        wt = pt[ic] * len(pt); wc = pc[jc] * len(pc)
        RA, RB = random_rot(m), random_rot(m)
        d = R_FREE * rng.random(m) ** (1 / 3)
        u = rng.normal(size=(m, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
        t = u * d[:, None]
        N   = np.einsum("mij,mj->mi", RA, gt["nrm"][ic])
        Ct1 = np.einsum("mij,mj->mi", RA, gt["c1"][ic])
        Ct2 = np.einsum("mij,mj->mi", RA, gt["c2"][ic])
        Q   = np.einsum("mij,mj->mi", RA, gt["cen"][ic])
        Ca  = np.einsum("mij,mj->mi", RB, gc["c1"][jc]) + t
        Cb  = np.einsum("mij,mj->mi", RB, gc["c2"][jc]) + t
        Mn  = np.einsum("mij,mj->mi", RB, gc["nrm"][jc])
        ok = nac_mask(Q, N, Ct1, Ct2, Ca, Cb, Mn, "loose")
        sel = np.where(ok)[0]
        if len(sel) == 0:
            continue
        A = np.einsum("mij,maj->mai", RA[sel], gt["inter"][ic[sel]])
        B = np.einsum("mij,maj->mai", RB[sel], gc["inter"][jc[sel]]) + t[sel][:, None, :]
        Ei = interarm_E(A, B, rm, ep)
        fac = {"boltz": wt[sel] * wc[sel] * np.exp(-Ei / RT),
               "unif": np.exp(-Ei / RT),
               "steric": (Ei <= E_ALLOW_PAIR).astype(float)}
        for lv in NAC_LEVELS:
            mk = nac_mask(Q[sel], N[sel], Ct1[sel], Ct2[sel], Ca[sel], Cb[sel], Mn[sel], lv)
            for wm in WMODES:
                acc[wm][lv] += float((fac[wm] * mk).sum())
    return {wm: {lv: acc[wm][lv] / n_mc for lv in NAC_LEVELS} for wm in WMODES}

V_FREE = 4 / 3 * np.pi * R_FREE ** 3

def pair_EM(plA, plB, f_free, block=1500):
    """テンプレート上で繋がれた 2 アームの NAC 確率 -> 実効モル濃度 EM (M)
       f_free = {wmode: {level: f}}。戻り値も {wmode: {level: EM}}"""
    treeB = cKDTree(plB.mid)
    rm, ep = interarm_pars(plA.arm, plB.arm)
    num = {wm: {lv: 0.0 for lv in NAC_LEVELS} for wm in WMODES}
    n_cand = 0; n_nac = 0; d_best = np.inf; e_best = np.inf; best_xyz = None
    for s0 in range(0, plA.n, block):
        s1 = min(s0 + block, plA.n)
        lists = treeB.query_ball_point(plA.cen[s0:s1], r=NAC_SEARCH_R)
        cnt = [len(x) for x in lists]
        tot = int(np.sum(cnt))
        if tot == 0:
            continue
        n_cand += tot
        ia = np.repeat(np.arange(s0, s1), cnt)
        ib = np.fromiter((j for x in lists for j in x), dtype=np.int64, count=tot)
        Q, N = plA.cen[ia], plA.nrm[ia]
        Ct1, Ct2 = plA.C1[ia], plA.C2[ia]
        Ca, Cb, Mn = plB.C1[ib], plB.C2[ib], plB.nrm[ib]
        ok = nac_mask(Q, N, Ct1, Ct2, Ca, Cb, Mn, "loose")
        sel = np.where(ok)[0]
        if len(sel) == 0:
            continue
        n_nac += len(sel)
        A = plA.X[ia[sel]][:, plA.arm.inter, :]
        B = plB.X[ib[sel]][:, plB.arm.inter, :]
        Ei = interarm_E(A, B, rm, ep)
        masks = {lv: nac_mask(Q[sel], N[sel], Ct1[sel], Ct2[sel], Ca[sel], Cb[sel], Mn[sel], lv)
                 for lv in NAC_LEVELS}
        good = masks[NAC_MAIN] & (Ei <= E_ALLOW_PAIR)
        if good.any():
            strain = (plA.E[ia[sel]] - plA.E.min()) + (plB.E[ib[sel]] - plB.E.min()) + np.maximum(Ei, 0)
            g = np.where(good)[0]
            k = g[int(np.argmin(strain[g]))]
            if float(strain[k]) < e_best:
                e_best = float(strain[k])
                best_xyz = (plA.X[ia[sel][k]].copy(), plB.X[ib[sel][k]].copy())
            dd = np.maximum(np.linalg.norm(Ct1[sel] - Ca[sel], axis=1),
                            np.linalg.norm(Ct2[sel] - Cb[sel], axis=1))
            d_best = min(d_best, float(dd[g].min()))
        fac = {"boltz": np.exp(-Ei / RT), "unif": np.exp(-Ei / RT),
               "steric": (Ei <= E_ALLOW_PAIR).astype(float)}
        for wm in WMODES:
            wp = plA.W[wm][ia[sel]] * plB.W[wm][ib[sel]] * fac[wm]
            for lv in NAC_LEVELS:
                num[wm][lv] += float((wp * masks[lv]).sum())
    out = {wm: {} for wm in WMODES}
    for wm in WMODES:
        Z = plA.W[wm].sum() * plB.W[wm].sum()
        for lv in NAC_LEVELS:
            ff = f_free[wm][lv]
            out[wm][lv] = ((num[wm][lv] / Z) / (ff * V_FREE)) * TO_M if (ff > 0 and Z > 0) else np.nan
    return out, dict(n_cand=n_cand, n_nac=n_nac, d_min=d_best, e_strain=e_best, xyz=best_xyz)

# =============================================================================
#  Part 4  系の組み立て（ASO 2 本 + ニック + 4 個のハンドル）と実行
# =============================================================================

L_ASO = len(ASO_SEQ)

def make_system(step_pool, tmpl, design, rl):
    """design = (p, q): ASO 内の修飾位置。反応対は ASO_A の q と ASO_B の p。
       接合部越えの間隔 Δ = L - q + p"""
    p, q = design
    F = N_FLANK_BP
    total = 2 * F + 2 * L_ASO
    seq = "".join(ASO_SEQ[(i - F) % L_ASO] for i in range(total))
    gA = lambda k: F + k
    gB = lambda k: F + L_ASO + k
    deaza = (gA(p), gA(q), gB(p), gB(q))
    idx = rl.integers(0, len(step_pool), total - 1)
    dx = build_duplex(seq, [step_pool[i] for i in idx], tmpl,
                      nick_after=F + L_ASO, deaza_pos=deaza)
    refine_duplex(dx)
    info = dict(seq=seq, tz_react=gA(q), cp_react=gB(p),
                tz_bys=gA(p), cp_bys=gB(q), delta=L_ASO - q + p)
    return dx, info

def anchor_of(dx, pos):
    r = dx.res1[pos - 1]
    C7, u = exocyclic_C7_dir(r)
    return C7, u, r["atoms"]["C8"]

def rna_env(dx, host_pos, extra=None):
    excl = {(0, host_pos, a) for a in ("N9", "C8", "C7", "N7", "C5", "C4")}
    X, E, Q, TAG = duplex_arrays(dx, exclude=excl)
    if extra:
        X = np.vstack([X] + [e[0] for e in extra])
        E = np.concatenate([E] + [e[1] for e in extra])
        Q = np.concatenate([Q] + [np.zeros(len(e[0])) for e in extra])
    rmin = np.array([RMIN2.get(e, 2.0) for e in E])
    eps = np.array([EPS.get(e, 0.12) for e in E])
    m = Q != 0
    return dict(tree=cKDTree(X), rmin=rmin, eps=eps, xq=X[m], q=Q[m], X=X)

def arm_elements(arm):
    return np.array([arm.mol.GetAtomWithIdx(i).GetSymbol() for i in arm.heavy])

def bystander_ensemble(dx, pos, arm):
    C7, u, C8 = anchor_of(dx, pos)
    return place_arm(arm, C7, u, C8, rna_env(dx, pos))

def _eval_once(dx, info, arm_tz, arm_cp, f_free, extra):
    aT = anchor_of(dx, info["tz_react"]); aC = anchor_of(dx, info["cp_react"])
    plT = place_arm(arm_tz, aT[0], aT[1], aT[2], rna_env(dx, info["tz_react"], extra))
    plC = place_arm(arm_cp, aC[0], aC[1], aC[2], rna_env(dx, info["cp_react"], extra))
    em, stat = pair_EM(plT, plC, f_free)
    return dict(EM=em, d77=float(np.linalg.norm(aT[0] - aC[0])),
                allowT=plT.allowed, allowC=plC.allowed,
                occT=plT.occupancy, occC=plC.occupancy,
                EminT=plT.Emin, EminC=plC.Emin,
                dmin=stat["d_min"], strain=stat["e_strain"], xyz=stat["xyz"],
                ncand=stat["n_cand"], nnac=stat["n_nac"], delta=info["delta"])

def evaluate_system(dx, info, arm_tz, arm_cp, f_free, bystanders=True):
    """傍観ハンドルは「最安定配置に固定」ではなく、自身のボルツマン分布から
       N_BYS 点引いて平均する（固定すると反応相手を人為的に塞いでしまうため）"""
    if not bystanders:
        return _eval_once(dx, info, arm_tz, arm_cp, f_free, [])
    pT = bystander_ensemble(dx, info["tz_bys"], arm_tz)
    pC = bystander_ensemble(dx, info["cp_bys"], arm_cp)
    eT, eC = arm_elements(arm_tz), arm_elements(arm_cp)
    def draw(pl, n):
        w = pl.W["boltz"]; w = w / w.sum()
        return rng.choice(len(w), n, p=w, replace=True)
    iT, iC = draw(pT, N_BYS), draw(pC, N_BYS)
    outs = [_eval_once(dx, info, arm_tz, arm_cp, f_free,
                       [(pT.X[a], eT), (pC.X[b], eC)]) for a, b in zip(iT, iC)]
    res = dict(min(outs, key=lambda o: o["strain"] if np.isfinite(o["strain"]) else np.inf))
    res["EM"] = {wm: {lv: float(np.mean([o["EM"][wm][lv] for o in outs]))
                      for lv in NAC_LEVELS} for wm in WMODES}
    for k in ("allowT", "allowC", "occT", "occC", "EminT", "EminC", "ncand", "nnac"):
        res[k] = float(np.mean([o[k] for o in outs]))
    for k in ("dmin", "strain"):
        v = [o[k] for o in outs if np.isfinite(o[k])]
        res[k] = float(np.median(v)) if v else np.inf
    return res

def boot_ci(v, n=4000, lo=2.5, hi=97.5):
    v = np.asarray([x for x in v if np.isfinite(x)])
    if len(v) < 3:
        return (np.nan, np.nan)
    b = np.median(rng.choice(v, (n, len(v))), axis=1)
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))

# =============================================================================
#  Part 5  ハイブリダイゼーション熱力学と、観測可能な反応速度
# =============================================================================
#  EM だけで半減期を語るのは誤り。6mer RNA 二重鎖は寿命が短いので、
#  「三元複合体ができている間に化学反応が終わるか」の競合を解く必要がある。
#      速度/接合部 = θ_AB * 2k_off * k_chem / (k_chem + 2k_off),  k_chem = k2 * EM
#  k_chem << k_off なら平衡近似 θ_AB*k2*EM に、k_chem >> k_off なら
#  三元複合体の生成律速 θ_AB*2k_off に収束する。

NN_RNA = {   # Turner 2004 (1 M NaCl): 5'XY3'/3'X'Y'5' -> (dH kcal/mol, dS cal/mol/K)
    "AA": (-6.82, -19.0), "AU": (-9.38, -26.7), "UA": (-7.69, -20.5),
    "CU": (-10.48, -27.1), "CA": (-10.44, -26.9), "GU": (-11.40, -29.5),
    "GA": (-12.44, -32.5), "CG": (-10.64, -26.7), "GG": (-13.39, -32.7),
    "GC": (-14.88, -36.9),
}
INIT_H, INIT_S = 3.61, -1.5
AU_END_H, AU_END_S = 3.72, 10.5
SALT_COEF = 0.10          # ΔG salt 補正係数（RNA は DNA の 0.114 よりやや小さい）

def _revcomp(s):
    return "".join({"A": "U", "U": "A", "G": "C", "C": "G"}[c] for c in reversed(s))

def duplex_thermo(seq, na_M=IONIC_M, T=T_K):
    dH, dS = INIT_H, INIT_S
    for i in range(len(seq) - 1):
        st = seq[i:i + 2]
        v = NN_RNA.get(st) or NN_RNA.get(_revcomp(st))
        if v is None:
            raise KeyError(st)
        dH += v[0]; dS += v[1]
    for e in (seq[0], seq[-1]):
        if e in "AU":
            dH += AU_END_H; dS += AU_END_S
    dG = dH - T * dS / 1000.0
    dG += SALT_COEF * (len(seq) - 1) * math.log(na_M)   # 低塩は不安定化 (ln<0)
    Kd = math.exp(dG / (KCAL * T))
    return dict(dH=dH, dS=dS, dG=dG, Kd=Kd)

def ligation_kinetics(EM, k2, aso_M, Kd, f_acc, k_on=1.0e6):
    Kd_eff = Kd / max(f_acc, 1e-9)
    th = aso_M / (aso_M + Kd_eff)
    th_AB = th * th
    k_off = k_on * Kd_eff
    k_chem = k2 * EM
    rate = th_AB * 2 * k_off * k_chem / (k_chem + 2 * k_off) if (k_chem + k_off) > 0 else 0.0
    k_bg = k2 * aso_M            # 非鋳型（2 分子）の擬 1 次速度
    return dict(theta=th, theta_AB=th_AB, k_off=k_off, k_chem=k_chem,
                rate=rate, t_half=(math.log(2) / rate if rate > 0 else np.inf),
                k_bg=k_bg, templating=(rate / k_bg if k_bg > 0 else np.inf),
                eq_limit=th_AB * k_chem, res_limit=th_AB * 2 * k_off)

def fmt_time(s):
    if not np.isfinite(s): return "     ---"
    if s < 60: return f"{s:7.1f} 秒"
    if s < 3600: return f"{s/60:7.1f} 分"
    if s < 86400 * 2: return f"{s/3600:7.1f} 時間"
    return f"{s/86400:7.1f} 日"

# =============================================================================
#  Part 6  実行
# =============================================================================

print("=" * 92); print("  [1] A 型 RNA らせんの幾何を収集"); print("=" * 92)
texts = []
for path in PDB_LOCAL_FILES:
    try:
        texts.append((path.split("/")[-1], open(path).read()))
    except Exception as e:
        print(f"    {path}: 読み込み失敗 {e}")
for pid in HARVEST_PDB_IDS:
    try:
        texts.append((pid, fetch_pdb_text(pid)))
    except Exception as e:
        print(f"    {pid}: 取得失敗 ({type(e).__name__})")
assert texts, "構造を 1 つも取得できませんでした（PDB_LOCAL_FILES を指定してください）"

STEPS, TMPL_RAW, STATS = harvest_geometry(texts)
assert len(STEPS) >= 4, "A 型らせんが見つかりませんでした"
need = set(ASO_SEQ) | {COMPL[c] for c in ASO_SEQ}
missing = [b for b in need if len(TMPL_RAW.get(b, [])) == 0]
assert not missing, f"塩基テンプレートが不足: {missing}"
TMPL = {}
print()
for b in sorted(need):
    nm, Lc, sp = medoid_template(TMPL_RAW[b])
    TMPL[b] = (nm, Lc)
    print(f"    {b} テンプレート: {len(TMPL_RAW[b])} 残基から medoid を選択（残基間 RMSD 中央値 {sp:.2f} Å）")
tw = [screw_params(R, t)[0] for R, t in STEPS]
ri = [screw_params(R, t)[1] for R, t in STEPS]
print(f"\n    採用したステップ変換 {len(STEPS)} 個: らせん twist {np.mean(tw):.1f} ± {np.std(tw):.1f}°, "
      f"rise {np.mean(ri):.2f} ± {np.std(ri):.2f} Å")
STEPS_SYM = [symmetrize_step(R, t) for R, t in STEPS]
print()

print("=" * 92); print("  [2] アーム分子の配座アンサンブル"); print("=" * 92)
ARM_TZ = build_arm(SMI_TZ, "テトラジンアーム")
ARM_CP = build_arm(SMI_CP, "シクロプロペンアーム")
print()

print("=" * 92); print("  [3] 非鋳型（自由溶液）参照系"); print("=" * 92)
F_FREE = free_reference(ARM_TZ, ARM_CP)
for lv in ("tight", "medium", "loose"):
    d, h, an, ax, off = NAC_LEVELS[lv]
    print(f"    {lv:6s} (結合 <= {d:.1f} Å, 面間角 <= {an:.0f}°, 軸角 <= {ax:.0f}°, ずれ <= {off:.1f} Å) : "
          f"f_free = " + " / ".join(f"{wm} {F_FREE[wm][lv]*100:.4f} %" for wm in WMODES))
print("    ※ v16 の判定 (結合 4.0 Å + 面から 1.2 Å のみ) は上の loose よりさらに緩く、")
print("       f_free = 6.1 % という非物理的な値になっていた。EM は比なので一部相殺されるが、")
print("       角度拘束がないぶん鋳型側で歪んだ配置まで「反応可能」と数えてしまう。")
assert F_FREE[WMODE_MAIN][NAC_MAIN] > 0, "参照系で NAC が 1 つも出ませんでした（サンプル数を増やしてください）"
print()
import time

PURINE_POS = [i + 1 for i, c in enumerate(ASO_SEQ) if c in PURINES]
DESIGNS = [(p, q) for p in PURINE_POS for q in PURINE_POS if p < q]
MAIN_DESIGN = tuple(MOD_POSITIONS)
assert MAIN_DESIGN in DESIGNS, f"修飾位置 {MOD_POSITIONS} はプリンではありません（使用可: {PURINE_POS}）"

print("=" * 92); print("  [4] 標的配列 A 型二重鎖の構築と検証"); print("=" * 92)
_rl = np.random.default_rng(SEED + 1)
_dx0, _info0 = make_system(STEPS_SYM, TMPL, MAIN_DESIGN, _rl)
_q = duplex_quality(_dx0)
print(f"    ASO 鎖 5'-{_info0['seq'][:N_FLANK_BP].lower()}"
      f"{_info0['seq'][N_FLANK_BP:N_FLANK_BP+L_ASO]}|{_info0['seq'][N_FLANK_BP+L_ASO:N_FLANK_BP+2*L_ASO]}"
      f"{_info0['seq'][N_FLANK_BP+2*L_ASO:].lower()}-3'   ( | = ニック)")
print(f"    鋳型鎖 3'-{''.join(COMPL[c] for c in _info0['seq'])}-5'   = (C4G2)n アンチセンス反復")
print(f"    骨格 O3'-P : {_q['o3p'][0]:.3f} - {_q['o3p'][1]:.3f} Å (理想 1.607)")
print(f"    WC 水素結合 N1···N3 : {_q['wc'][0]:.2f} - {_q['wc'][1]:.2f} Å")
print(f"    残基間の近接衝突 (<2.4 Å) : {_q['clashes']} 個")
print(f"    反応する 2 ハンドル: ASO_A 第{MAIN_DESIGN[1]}位 (pos {_info0['tz_react']}, Tz) と "
      f"ASO_B 第{MAIN_DESIGN[0]}位 (pos {_info0['cp_react']}, Cp) / 鎖内間隔 Δ = {_info0['delta']} nt")
print(f"    未反応の傍観ハンドル: pos {_info0['tz_bys']} (Tz), pos {_info0['cp_bys']} (Cp)")
print()

write_model_pdb("duplex_v17.pdb", _dx0)
print("    構築した二重鎖を duplex_v17.pdb に保存しました（PyMOL 等で確認できます）")
print()

print("=" * 92); print(f"  [5] 実行設計 (p,q) = {MAIN_DESIGN} の実効モル濃度"); print("=" * 92)
print(f"    熱ゆらぎ二重鎖 {N_TEMPLATES} 本 x 傍観アーム有無")
RES = {True: [], False: []}
t0 = time.time()
for it in range(N_TEMPLATES):
    rl = np.random.default_rng(SEED + 100 + it)
    dx, info = make_system(STEPS_SYM, TMPL, MAIN_DESIGN, rl)
    for byst in ([True, False] if USE_BYSTANDERS else [False]):
        RES[byst].append(evaluate_system(dx, info, ARM_TZ, ARM_CP, F_FREE, bystanders=byst))
    if (it + 1) % max(1, N_TEMPLATES // 6) == 0:
        print(f"      {it+1}/{N_TEMPLATES} 本 ({time.time()-t0:.0f} s)")

def summarize(rows, label):
    em = np.array([r["EM"][WMODE_MAIN][NAC_MAIN] for r in rows])
    lo, hi = boot_ci(em)
    print(f"\n    --- {label} (n = {len(rows)}) ---")
    print(f"    C7···C7 距離           : {np.median([r['d77'] for r in rows]):.2f} Å")
    print(f"    主溝で許される配置の割合: Tz {np.median([r['allowT'] for r in rows])*100:.2f} % / "
          f"Cp {np.median([r['allowC'] for r in rows])*100:.2f} %  (最良配置 +{E_ALLOW:.0f} kcal/mol 以内)")
    print(f"    NAC 候補対             : {np.median([r['nnac'] for r in rows]):.0f} / "
          f"{np.median([r['ncand'] for r in rows]):.0f} 対")
    _dm = [r['dmin'] for r in rows if np.isfinite(r['dmin'])]
    _st = [r['strain'] for r in rows if np.isfinite(r['strain'])]
    print(f"    到達できた最短形成結合 : {np.median(_dm):.2f} Å" if _dm else
          "    到達できた最短形成結合 : 到達不能")
    print(f"    反応配座に必要なひずみ : {np.median(_st):.1f} kcal/mol (最良配置基準)" if _st else
          "    反応配座に必要なひずみ : ---")
    print(f"    {'':22s}{'steric':>12s} {'unif':>12s} {'boltz':>12s}")
    for lv in ("tight", "medium", "loose"):
        vals = {wm: np.array([r["EM"][wm][lv] for r in rows]) for wm in WMODES}
        print(f"    EM ({lv:6s}) 中央値 " + " ".join(f"{np.median(vals[wm]):12.3e}" for wm in WMODES))
    print(f"    {WMODE_MAIN}/{NAC_MAIN} の 95%CI : {lo:.3e} – {hi:.3e}   "
          f"ゼロ {int((em == 0).sum())}/{len(em)}")
    return float(np.median(em)), (lo, hi), em

_MAIN_ROWS = RES[True] if USE_BYSTANDERS else RES[False]
EM_MAIN, CI_MAIN, EM_ARR = summarize(_MAIN_ROWS,
                                     "傍観ハンドルあり（実際の設計）" if USE_BYSTANDERS else "傍観ハンドルなし")
EM_GEOM = float(np.median([r["EM"][WMODE_GEOM][NAC_MAIN] for r in _MAIN_ROWS]))
_cands = [(i, r) for i, r in enumerate(_MAIN_ROWS)
          if r.get("xyz") is not None and np.isfinite(r["strain"])]
_bi, _bestrow = min(_cands, key=lambda t: t[1]["strain"]) if _cands else (None, None)
if _bestrow is not None:
    _dxb, _ = make_system(STEPS_SYM, TMPL, MAIN_DESIGN, np.random.default_rng(SEED + 100 + _bi))
    write_model_pdb("nac_model_v17.pdb", _dxb,
                    arms=[(_bestrow["xyz"][0], arm_elements(ARM_TZ), "TZ"),
                          (_bestrow["xyz"][1], arm_elements(ARM_CP), "CP")])
    print(f"\n    最もひずみの小さい反応配座 (NAC) を nac_model_v17.pdb に保存しました "
          f"(ひずみ {_bestrow['strain']:.1f} kcal/mol, 形成結合 {_bestrow['dmin']:.2f} Å)")
if USE_BYSTANDERS:
    EM_NB, CI_NB, EM_ARR_NB = summarize(RES[False], "傍観ハンドルなし（参考）")
print()
print("=" * 92); print("  [6] 修飾位置 (p, q) の走査"); print("=" * 92)
print(f"    接合部越えの間隔 Δ = {L_ASO} - q + p。ASO 内の 2 つのハンドルは同種なので")
print(f"    反応するのは常に接合部越えの対のみ。傍観アームは省いた比較（{N_TEMPL_SCAN} 本平均）\n")
SCAN = {}
for design in sorted(DESIGNS, key=lambda d: (L_ASO - d[1] + d[0], d)):
    rows = []
    for it in range(N_TEMPL_SCAN):
        rl = np.random.default_rng(SEED + 500 + it)
        dx, info = make_system(STEPS_SYM, TMPL, design, rl)
        rows.append(evaluate_system(dx, info, ARM_TZ, ARM_CP, F_FREE, bystanders=False))
    v = np.array([r["EM"][WMODE_MAIN][NAC_MAIN] for r in rows])
    vg = np.array([r["EM"][WMODE_GEOM][NAC_MAIN] for r in rows])
    lo, hi = boot_ci(v)
    SCAN[design] = dict(em=v, emg=vg, med=float(np.median(v)), medg=float(np.median(vg)), ci=(lo, hi),
                        d77=float(np.median([r["d77"] for r in rows])), delta=rows[0]["delta"],
                        strain=float(np.median([r["strain"] for r in rows if np.isfinite(r["strain"])]))
                        if any(np.isfinite(r["strain"]) for r in rows) else np.nan)
    mark = "  <-- 今回の設計" if design == MAIN_DESIGN else ""
    print(f"    (p,q) = {design}  Δ = {rows[0]['delta']} nt  C7···C7 {SCAN[design]['d77']:5.2f} Å  "
          f"EM(steric) {SCAN[design]['medg']:.2e} M  EM(boltz) {SCAN[design]['med']:.2e} M  "
          f"ひずみ {SCAN[design]['strain']:5.1f} kcal/mol{mark}")
print()

print("=" * 92); print("  [7] ハイブリダイゼーション熱力学"); print("=" * 92)
TH = duplex_thermo(ASO_SEQ)
Tm_1uM = None
try:
    c_t = 1e-6
    Tm_1uM = (TH["dH"] * 1000.0) / (TH["dS"] + 1.987 * math.log(c_t / 4.0)) - 273.15
except Exception:
    pass
print(f"    ASO 5'-{ASO_SEQ}-3' : ΔH° = {TH['dH']:.1f} kcal/mol, ΔS° = {TH['dS']:.1f} cal/mol/K")
print(f"    {TEMP_C:.0f} °C, {IONIC_M*1000:.0f} mM Na+ : ΔG° = {TH['dG']:.2f} kcal/mol -> Kd = {TH['Kd']*1e9:.1f} nM")
if Tm_1uM: print(f"    Tm (全鎖 1 µM, 1 M NaCl 基準) ≈ {Tm_1uM:.0f} °C")
print("    ※ 修飾（C7 置換基）と標的の高次構造（G4 / ヘアピン）による不安定化は f_acc で表現する")
print()

print("=" * 92); print("  [8] 予測される反応速度"); print("=" * 92)
K2_LIST   = [1.0, 5.0, 20.0]      # M^-1 s^-1 : 3-メチル-6-アリールテトラジン x 1-メチルシクロプロペン
ASO_LIST  = [1e-7, 1e-6, 1e-5]
F_ACC     = 0.1                    # 標的の高次構造などで実効的に結合できる割合
K_ON      = 1.0e6
print(f"    k_on = {K_ON:.0e} M⁻¹s⁻¹, 標的接近可能割合 f_acc = {F_ACC}")
print(f"    EM は上下 2 通りで評価: EM_steric = {EM_GEOM:.3e} M（幾何学的到達性のみ・楽観側） /")
print(f"                            EM_boltz  = {EM_MAIN:.3e} M（反応配座のひずみを全部払う・悲観側）")
print(f"    以下の表は EM_boltz。EM_steric なら速度は約 {EM_GEOM/EM_MAIN if EM_MAIN>0 else float('nan'):.0f} 倍（ただし k_off 律速で頭打ち）")
print(f"\n    {'k2':>6s} {'[ASO]':>8s} {'θ(1鎖)':>8s} {'θ_AB':>8s} {'k_off':>10s} {'k_chem':>10s} "
      f"{'速度':>11s} {'半減期':>12s} {'鋳型効果':>10s}")
KIN = {}
for k2 in K2_LIST:
    for c in ASO_LIST:
        r = ligation_kinetics(EM_MAIN, k2, c, TH["Kd"], F_ACC, K_ON)
        KIN[(k2, c)] = r
        print(f"    {k2:6.1f} {c*1e6:7.2f}µM {r['theta']:8.3f} {r['theta_AB']:8.3f} "
              f"{r['k_off']:10.3f} {r['k_chem']:10.4f} {r['rate']:11.3e} {fmt_time(r['t_half'])} "
              f"{r['templating']:9.1f}x")
print("\n    速度 = θ_AB · 2k_off · k_chem/(k_chem+2k_off)  [s⁻¹ / 接合部]")
print("    鋳型効果 = (鋳型上の擬 1 次速度) / (非鋳型 2 分子反応の擬 1 次速度 k2[ASO])")
_r = KIN[(K2_LIST[1], ASO_LIST[1])]
print(f"    平衡近似の上限 {_r['eq_limit']:.3e} s⁻¹ / 三元複合体寿命による上限 {_r['res_limit']:.3e} s⁻¹")
print(f"    -> 律速は {'化学反応' if _r['k_chem'] < 2*_r['k_off'] else '二重鎖の解離（滞在時間）'}")
print()
if RUN_SENSITIVITY:
    print("=" * 92); print("  [9] 経験的パラメータへの感度"); print("=" * 92)
    print("    剛体近似を補うためのノブ（引力の減衰・反発の飽和・結合方向のゆらぎ）を振る。")
    print("    絶対値はこれらに依存するので、設計間の相対比較を主に見るべき理由を示す。\n")
    _base = dict(ATTR_SCALE=ATTR_SCALE, LJ_REP_CAP=LJ_REP_CAP, DIR_WOBBLE_DEG=DIR_WOBBLE_DEG)
    _cases = [("基準", {}),
              ("引力 x0.10", dict(ATTR_SCALE=0.10)),
              ("引力 x0.35", dict(ATTR_SCALE=0.35)),
              ("反発を硬く (cap 10)", dict(LJ_REP_CAP=10.0)),
              ("方向ゆらぎ 0°", dict(DIR_WOBBLE_DEG=0.0)),
              ("方向ゆらぎ 15°", dict(DIR_WOBBLE_DEG=15.0))]
    N_SENS = 3
    print(f"    {'条件':<22s}{'許容配置 Tz/Cp':>18s}{'ひずみ':>10s}{'EM_steric':>13s}{'EM_boltz':>13s}")
    _ff_cache = {}
    for name, over in _cases:
        for k, v in _base.items():
            globals()[k] = v
        for k, v in over.items():
            globals()[k] = v
        key = (ATTR_SCALE, LJ_REP_CAP)
        if key not in _ff_cache:
            _ff_cache[key] = free_reference(ARM_TZ, ARM_CP, n_mc=min(N_FREE_MC, 1_000_000))
        ff = _ff_cache[key]
        rr = []
        for it in range(N_SENS):
            rl = np.random.default_rng(SEED + 100 + it)
            dxs, infos = make_system(STEPS_SYM, TMPL, MAIN_DESIGN, rl)
            rr.append(evaluate_system(dxs, infos, ARM_TZ, ARM_CP, ff, bystanders=False))
        st = [r["strain"] for r in rr if np.isfinite(r["strain"])]
        print(f"    {name:<22s}{np.median([r['allowT'] for r in rr])*100:8.1f} %/"
              f"{np.median([r['allowC'] for r in rr])*100:5.1f} %"
              f"{(np.median(st) if st else np.nan):10.1f}"
              f"{np.median([r['EM'][WMODE_GEOM][NAC_MAIN] for r in rr]):13.2e}"
              f"{np.median([r['EM'][WMODE_MAIN][NAC_MAIN] for r in rr]):13.2e}")
    for k, v in _base.items():
        globals()[k] = v
    print()

print("=" * 92); print("  [10] 図"); print("=" * 92)
fig = plt.figure(figsize=(15, 9))

ax = fig.add_subplot(2, 3, 1)
ds = sorted(SCAN, key=lambda d: (SCAN[d]["delta"], d))
xs = np.arange(len(ds))
_pos = [SCAN[d]["medg"] for d in ds if SCAN[d]["medg"] > 0]
_floor = (min(_pos) / 30.0) if _pos else 1e-9
med = [max(SCAN[d]["medg"], _floor) for d in ds]
ax.bar(xs, med, bottom=0, color=["#c0392b" if d == MAIN_DESIGN else "#5b8db8" for d in ds])
for i, d in enumerate(ds):
    if SCAN[d]["medg"] <= 0:
        ax.text(i, _floor * 1.3, "0", ha="center", fontsize=8)
ax.set_ylim(_floor * 0.6, max(med) * 3)
ax.set_xticks(xs); ax.set_xticklabels([f"{d}\nΔ={SCAN[d]['delta']}" for d in ds], fontsize=7)
ax.set_yscale("log"); ax.set_ylabel("EM steric (M)"); ax.set_title("Modification positions (p,q)\nred = current design")
ax.grid(alpha=.3, axis="y")

ax = fig.add_subplot(2, 3, 2)
data = [EM_ARR] + ([EM_ARR_NB] if USE_BYSTANDERS else [])
ax.boxplot(data, tick_labels=(["with\nbystanders", "without"] if USE_BYSTANDERS else ["EM"]))
ax.set_yscale("log"); ax.set_ylabel("EM (M)")
ax.set_title(f"EM over {N_TEMPLATES} thermal duplexes"); ax.grid(alpha=.3, axis="y")

ax = fig.add_subplot(2, 3, 3)
rows = RES[True] if USE_BYSTANDERS else RES[False]
lv_names = ["tight", "medium", "loose"]
ax.boxplot([[r["EM"][WMODE_MAIN][lv] for r in rows] for lv in lv_names], tick_labels=lv_names)
ax.set_yscale("log"); ax.set_ylabel("EM (M)")
ax.set_title("Sensitivity to NAC stringency"); ax.grid(alpha=.3, axis="y")

ax = fig.add_subplot(2, 3, 4)
cc = np.logspace(-8, -4, 60)
for k2 in K2_LIST:
    ln, = ax.loglog(cc * 1e6, [ligation_kinetics(EM_MAIN, k2, c, TH["Kd"], F_ACC, K_ON)["rate"] for c in cc],
                    label=f"k2 = {k2:g} M⁻¹s⁻¹")
    ax.loglog(cc * 1e6, [k2 * c for c in cc], "--", lw=.8, alpha=.6, color=ln.get_color())
    ax.loglog(cc * 1e6, [ligation_kinetics(EM_GEOM, k2, c, TH["Kd"], F_ACC, K_ON)["rate"] for c in cc],
              ":", lw=1.0, alpha=.8, color=ln.get_color())
ax.set_xlabel("[ASO] (µM)"); ax.set_ylabel("pseudo-1st-order rate (s⁻¹)")
ax.set_title("templated EM_boltz (solid) / EM_steric (dotted)\nvs untemplated (dashed)"); ax.legend(fontsize=7); ax.grid(alpha=.3, which="both")

ax = fig.add_subplot(2, 3, 5)
for k2 in K2_LIST:
    ax.semilogx(cc * 1e6, [ligation_kinetics(EM_MAIN, k2, c, TH["Kd"], F_ACC, K_ON)["templating"] for c in cc],
                label=f"k2 = {k2:g}")
ax.set_xlabel("[ASO] (µM)"); ax.set_ylabel("templating factor"); ax.set_yscale("log")
ax.set_title("Rate enhancement by the template"); ax.legend(fontsize=7); ax.grid(alpha=.3, which="both")

ax = fig.add_subplot(2, 3, 6)
tt = np.logspace(0, 6, 200)
for k2 in K2_LIST:
    r = ligation_kinetics(EM_MAIN, k2, ASO_LIST[1], TH["Kd"], F_ACC, K_ON)["rate"]
    ax.semilogx(tt, 100 * (1 - np.exp(-r * tt)), label=f"k2 = {k2:g}")
ax.set_xlabel("time (s)"); ax.set_ylabel("ligated junctions (%)")
ax.set_title(f"Yield at [ASO] = {ASO_LIST[1]*1e6:g} µM"); ax.legend(fontsize=7); ax.grid(alpha=.3, which="both")

plt.tight_layout(); plt.savefig("iedda_v17.png", dpi=150); plt.show()

# =============================================================================
#  Part 7  まとめと限界
# =============================================================================
print("=" * 92); print("  [11] まとめ"); print("=" * 92)
best = max(SCAN, key=lambda d: SCAN[d]["medg"])
print(f"""
  ■ 幾何
    - 標的 (C4G2)n 上で ASO 2 本が隣接すると、設計 (p,q)={MAIN_DESIGN} の反応対は
      鎖内 Δ = {L_ASO - MAIN_DESIGN[1] + MAIN_DESIGN[0]} nt 離れ、C7···C7 は {np.median([r['d77'] for r in rows]):.1f} Å。
    - A 型 RNA の主溝は深く狭く、立体的に許される配置は全配置の
      {np.median([r['allowT'] for r in rows])*100:.1f} % (Tz) / {np.median([r['allowC'] for r in rows])*100:.1f} % (Cp) しかない。
    - 反応配座 (NAC) に届くには最良配置から中央値 {np.median([r['strain'] for r in rows if np.isfinite(r['strain'])]) if any(np.isfinite(r['strain']) for r in rows) else float('nan'):.1f} kcal/mol のひずみが要る。
    - 修飾位置の走査では (p,q) = {best} が最良 (EM_steric {SCAN[best]['medg']:.2e} M)。
      現行設計との比 = {SCAN[MAIN_DESIGN]['medg']/SCAN[best]['medg'] if SCAN[best]['medg']>0 else float('nan'):.2f}

  ■ 実効モル濃度（2 通りの見方）
    EM_steric = {EM_GEOM:.3e} M  幾何学的に到達できるかだけを見た値（統計が安定・楽観側）
    EM_boltz  = {EM_MAIN:.3e} M  反応配座のひずみを全部払った値（悲観側）
                95%CI {CI_MAIN[0]:.2e} – {CI_MAIN[1]:.2e}, NAC = {NAC_MAIN}
    真の値はこの間にある。v16 の 2.3e-2 M は判定が緩く EM_steric 寄りの評価だった。
    NAC 基準を厳/緩に振ると {np.median([r['EM'][WMODE_MAIN]['tight'] for r in rows]):.2e} – {np.median([r['EM'][WMODE_MAIN]['loose'] for r in rows]):.2e} M。
    ここは基準依存性が最大の不確かさなので、EM の絶対値より
    「設計間の相対比較」を信用してください。

  ■ 速度論（ここが v16 で最も欠けていた部分）
    6mer RNA:RNA の Kd = {TH['Kd']*1e9:.0f} nM、f_acc = {F_ACC} を仮定すると
    実効 Kd = {TH['Kd']/F_ACC*1e9:.0f} nM。1 µM では 1 サイトの占有率 θ = {KIN[(K2_LIST[1],ASO_LIST[1])]['theta']:.2f}、
    隣接 2 サイトが同時に埋まる確率 θ_AB = {KIN[(K2_LIST[1],ASO_LIST[1])]['theta_AB']:.2f}。
    k_off = {KIN[(K2_LIST[1],ASO_LIST[1])]['k_off']:.2f} s⁻¹（三元複合体の寿命 {1/max(KIN[(K2_LIST[1],ASO_LIST[1])]['k_off'],1e-12)/2:.1f} 秒）に対し
    k_chem = k2·EM = {KIN[(K2_LIST[1],ASO_LIST[1])]['k_chem']:.4f} s⁻¹。
    -> 一度の結合イベント中に反応が完了する確率は
       {KIN[(K2_LIST[1],ASO_LIST[1])]['k_chem']/(KIN[(K2_LIST[1],ASO_LIST[1])]['k_chem']+2*KIN[(K2_LIST[1],ASO_LIST[1])]['k_off'])*100:.2f} %。
    「EM が高い」ことと「速く進む」ことは別問題であることに注意。

  ■ この計算が答えていないこと（設計上の重要因子）
    1. 標的 r(C4G2)n 自身のヘアピン / G-四重鎖形成。f_acc に押し込んであるだけで、
       実際にはこれが律速になりうる。
    2. ASO 自身 (CGGGGC は G リッチ) の自己会合・四重鎖化。
    3. 二重鎖は剛体として扱っており、アームを収めるための局所的な変形
       （主溝の拡がり、糖パッカーの変化）を許していない -> EM は過小評価側。
    4. MMFF 配座エネルギーは気相ベース。水中では折れ畳んだ配座が過大評価され、
       伸びた配座が過小評価されている可能性がある。
    5. 連結後は 12mer になって二重鎖が一気に安定化する（協同的な鎖状伸長）。
       ここでは 1 接合部の初回反応のみを扱っている。
    6. テトラジンアームの 2 級アミンはテトラジンを求核攻撃・還元しうる。
       速度以前に安定性の検討が必要。
    7. 1-メチルシクロプロペンと 3-メチル-6-アリールテトラジンの k2 は
       文献値の幅が大きい。ここでは 1-20 M⁻¹s⁻¹ を仮定した。
""")
print("  図を iedda_v17.png に保存しました。")
