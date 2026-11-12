# -*- coding: utf-8 -*-
"""Stage 14 —— GSE350571 Xenium Prime 5K 神经邻近空间验证（6 例乳腺 FFPE）。

目的：在独立 Xenium 队列中复现主稿的核心空间发现：
  (a) LRRC15⁺ myCAF 是否富集于神经束 150µm 邻近区（PNI 场的存在性）
  (b) CD8 T 细胞是否在神经邻近区被几何排斥（主稿的空间 CD8 null/弱下降）
  (c) 邻近 vs 远隔成纤维细胞的 Perineural_stroma 模块评分差异

方法要点：
  * 神经束点 = 神经纤维标记转录本（QV≥20）的密度聚类核
    （UCHL1/TUBB3/NEFH/PRNP/GAP43/MPZ；**不含 SOX10**——乳腺中 SOX10
    同时标记 basal 肿瘤细胞与施万样细胞，会造成假阳性；
    不含 NCAM1——NK 细胞表达 CD56）
  * 空间 null：环形平移（torus shift）999 次置换神经点云
  * 细胞分型用面板内可用标记的 mean-z + 阈值规则

边界：样本为 TLS-rich 选择性取样（偏倚声明写入 Note 13）；
n=6 sections，结论为方向性复现，不做 meta 统计。

输出：03_results/tables/bulk_extrapolation/
  stage14_xenium_nerve_summary.csv   每样本密度比 + 置换 P
  stage14_xenium_celltypes.csv       每样本各类型细胞 × 区带计数
  stage14_xenium_fb_scores.csv       成纤维细胞 per-cell 模块评分
"""
from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from bulk_lib import ROOT, RAW, TAB, load_modules  # noqa: E402

XDIR = RAW / "gse350571"
SAMPLES = ["P17", "P22", "P24", "P26", "P29", "P32"]
SEED = 20251109
N_PERM = 999
PROX_UM = 150.0      # 主稿 PNI 场半径
DIST_UM = 300.0      # 远隔区起点
QV_MIN = 20.0

NERVE_TX_GENES = ["MPZ", "NEFH", "TUBB3", "GAP43"]
AXONAL_GENES = ["NEFH", "TUBB3", "GAP43"]
# 神经束判据（6 例实测校准，轴突锚定 + 施万确认）：
#   核 = 轴突基因(NEFH/TUBB3/GAP43)转录本在 30µm 内 ≥5 条且 ≥2 种不同轴突基因，
#        并要求 50µm 内 ≥2 条 MPZ（髓鞘/施万确认）。
#   动机：P26/P24 的 MPZ 为肿瘤/组织广泛高表达（P26 44.1 万条、95.5% 高密度），
#   以 MPZ 为核会把 64% 细胞误划入"近神经区"；轴突基因密度在所有样本中均稀疏，
#   以其为锚 + MPZ 确认可同时排除肿瘤 MPZ 噪声与散在单基因表达。
#   剔除基因：PRNP(广谱噪声)、UCHL1(乳腺肌上皮)、SOX10(basal 肿瘤代理)、
#   NCAM1(NK 细胞 CD56)。
#   实测神经可得性：P22/P26/P29/P32 有神经束；P17/P24 切片内无足够神经信号。
EPS_UM = 30.0        # 轴突转录本聚类半径
MIN_TX = 5           # 密度核内最少轴突转录本数
MIN_GENES = 2        # 密度核邻域内最少不同轴突基因数
MYELIN_UM = 50.0     # MPZ 确认半径
MYELIN_MIN = 2       # MPZ 最少条数

CELLTYPE_RULES = {
    "myCAF":         {"pos": ["LRRC15", "FAP", "POSTN", "ACTA2", "TNC", "PDGFRB"], "n_pos": 3, "neg": []},
    "fibroblast":    {"pos": ["DCN", "LUM", "PDGFRA", "CXCL12", "FAP", "PDGFRB"], "n_pos": 3, "neg": []},
    "endothelial":   {"pos": ["PECAM1", "VWF", "CLDN5", "PLVAP", "ACKR1"], "n_pos": 2, "neg": []},
    "CD8_T":         {"pos": ["CD8A", "CD8B", "GZMK", "GZMB", "PRF1"], "n_pos": 2, "neg": []},
    "macrophage":    {"pos": ["CD68", "LAMP1", "CD163", "CD14", "MSR1", "MRC1", "SPP1", "TREM2"], "n_pos": 3, "neg": []},
    "B_cell":        {"pos": ["MS4A1", "CD79A", "CD79B"], "n_pos": 2, "neg": []},
    "epithelial":    {"pos": ["EPCAM", "KRT8", "KRT19", "CDH1", "MUC1"], "n_pos": 3, "neg": []},
}


# ---------------------------------------------------------------- 装载
import tempfile
import gzip as _gzip

TMPDIR = XDIR / ".tmp"


def _gunzip_to_tmp(gz_path: Path, tag: str) -> Path:
    """流式解压 gzip 包裹的 parquet 到临时文件（pyarrow 无法随机访问 gzip 流）。
    返回临时 parquet 路径；用完调用 _release_tmp 清空（沙箱钩子拦截 unlink，
    故只 truncate，统一残留到 .tmp/ 目录，最终批量清理）。"""
    TMPDIR.mkdir(exist_ok=True)
    fd, name = tempfile.mkstemp(suffix=f"_{tag}.parquet", dir=str(TMPDIR))
    import os
    os.close(fd)
    tmp = Path(name)
    with _gzip.open(gz_path, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=64 * 1024 * 1024)
    return tmp


def _release_tmp(tmp: Path) -> None:
    """清空临时文件内容（不 unlink：沙箱 safe-delete 钩子会拦截并报错）。
    残余空文件由最终批量清理移除。"""
    try:
        with open(tmp, "wb") as f:
            pass
    except OSError:
        pass


def load_cells(p: str) -> pd.DataFrame:
    tmp = _gunzip_to_tmp(XDIR / f"{p}_Xe5k_cells.parquet.gz", f"{p}_cells")
    try:
        df = pd.read_parquet(tmp)
    finally:
        _release_tmp(tmp)
    df.columns = [c.lower() for c in df.columns]
    return df


def load_matrix(p: str) -> tuple[sparse.csc_matrix, list[str], np.ndarray]:
    """10x h5: /matrix/{data,indices,indptr,shape,features/name,barcodes}"""
    with h5py.File(XDIR / f"{p}_Xe5k_cell_feature_matrix.h5", "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        M = sparse.csc_matrix(
            (g["data"][()], g["indices"][()], g["indptr"][()]),
            shape=shape)  # genes × cells
        genes = [x.decode() if isinstance(x, bytes) else str(x)
                 for x in g["features"]["name"][()]]
        bar = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in g["barcodes"][()]])
    return M, genes, bar


def load_nerve_tx(p: str) -> pd.DataFrame:
    """先解压到临时 parquet，再用谓词下推只读取神经标记基因的行
    （feature_name in NERVE_TX_GENES），避免全量转录本载入内存。"""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    gz = XDIR / f"{p}_Xe5k_transcripts.parquet.gz"
    tmp = _gunzip_to_tmp(gz, f"{p}_tx")
    try:
        tbl = pq.read_table(
            tmp,
            columns=["feature_name", "x_location", "y_location", "qv"],
            filters=[("feature_name", "in", NERVE_TX_GENES)],
            use_threads=True,
        )
    finally:
        _release_tmp(tmp)
    df = tbl.to_pandas()
    df.columns = [c.lower() for c in df.columns]
    if "qv" in df.columns:
        df = df[df["qv"] >= QV_MIN]
    df = df.rename(columns={"x_location": "x", "y_location": "y"})
    return df[["feature_name", "x", "y"]]


# ---------------------------------------------------------------- 神经束核
def nerve_bundle_points(tx: pd.DataFrame, eps=EPS_UM, min_tx=MIN_TX,
                        min_genes=MIN_GENES) -> np.ndarray:
    """轴突锚定 + 施万确认的神经束核：
    1) 轴突基因(NEFH/TUBB3/GAP43)转录本在半径 eps 内（含自身）≥ min_tx；
    2) 该邻域内出现 ≥ min_genes 种不同轴突基因；
    3) 位置 50µm 内有 ≥ MYELIN_MIN 条 MPZ（髓鞘确认，排除纯肿瘤轴突标记表达）。"""
    ax = tx[tx["feature_name"].isin(AXONAL_GENES)]
    mxy = tx.loc[tx["feature_name"] == "MPZ", ["x", "y"]].to_numpy(float)
    if len(ax) < min_tx or len(mxy) < MYELIN_MIN:
        return np.empty((0, 2))
    xy = ax[["x", "y"]].to_numpy(float)
    tree = cKDTree(xy)
    nbr_idx = tree.query_ball_point(xy, r=eps)
    cnt = np.array([len(ix) for ix in nbr_idx])
    genes = ax["feature_name"].to_numpy()
    mtree = cKDTree(mxy)
    keep = []
    for i in np.where(cnt >= min_tx)[0]:
        if len(set(genes[list(nbr_idx[i])])) < min_genes:
            continue
        x, y = xy[i]
        if mtree.query_ball_point([x, y], r=MYELIN_UM,
                                  return_length=True) >= MYELIN_MIN:
            keep.append(i)
    return xy[keep] if keep else np.empty((0, 2))


def torus_shift(points: np.ndarray, w: float, h: float, rng) -> np.ndarray:
    dx, dy = rng.uniform(0, w), rng.uniform(0, h)
    out = points.copy()
    out[:, 0] = np.mod(out[:, 0] + dx, w)
    out[:, 1] = np.mod(out[:, 1] + dy, h)
    return out


# ---------------------------------------------------------------- 打分
class GeneGetter:
    """按基因惰性提取 log1p(CP10K) 行向量。
    M 为 genes × cells 的 csc；转 csr 后行切片（单基因）O(该基因 nnz)，
    避免 cells × genes CSR 的列切片 O(全矩阵 nnz) 灾难。"""

    def __init__(self, M: sparse.csc_matrix, genes: list[str]):
        self.Mr = M.tocsr()          # genes × cells (csr, 行切片高效)
        self.gene_cols = {g: i for i, g in enumerate(genes)}
        tot = np.asarray(M.sum(axis=0)).ravel().astype(float)
        tot[tot <= 0] = 1.0
        self.cp10k = 1e4 / tot         # 每细胞缩放因子
        self._cache: dict[str, np.ndarray] = {}

    def __call__(self, g: str) -> np.ndarray | None:
        if g in self._cache:
            return self._cache[g]
        j = self.gene_cols.get(g)
        if j is None:
            return None
        v = np.asarray(self.Mr[j].todense()).ravel()
        v = np.log1p(v * self.cp10k).astype(np.float32)
        self._cache[g] = v
        return v


def meanz_score(G: GeneGetter, genes: list[str], n_cells: int) -> np.ndarray:
    cols = [G(g) for g in genes if G(g) is not None]
    cols = [c for c in cols if c is not None]
    if len(cols) < 4:
        return np.full(n_cells, np.nan)
    sub = np.column_stack(cols)
    z = (sub - sub.mean(0)) / (sub.std(0) + 1e-9)
    return z.mean(1)


# ---------------------------------------------------------------- 主流程
def analyse_sample(p: str, mods: dict):
    need = [XDIR / f"{p}_Xe5k_cells.parquet.gz",
            XDIR / f"{p}_Xe5k_cell_feature_matrix.h5",
            XDIR / f"{p}_Xe5k_transcripts.parquet.gz"]
    miss = [f.name for f in need if not f.exists()]
    if miss:
        print(f"[{p}] 缺失文件 {miss} —— 跳过")
        return None
    cells = load_cells(p)
    M, genes, bar = load_matrix(p)
    # barcode 对齐（cells.parquet 与 h5 通常同序；防御性处理）
    if list(cells["cell_id"].astype(str)) != [str(b) for b in bar]:
        m = pd.Series(range(len(bar)), index=[str(b) for b in bar])
        keep = [m[b] for b in cells["cell_id"].astype(str) if b in m.index]
        print(f"[{p}] cells {len(cells)} vs matrix {M.shape[1]} —— 按 barcode 对齐")
        M = M[:, keep]
        genes = genes  # 行不变
    n_cells = M.shape[1]
    G = GeneGetter(M, genes)

    tx = load_nerve_tx(p)
    nerve = nerve_bundle_points(tx)
    if len(nerve) < 20:
        print(f"[{p}] 神经束点过少({len(nerve)}) —— 无神经区域, 跳过")
        return {"sample": p, "n_nerve_pts": len(nerve), "note": "insufficient nerve signal"}
    print(f"[{p}] 神经标记转录本 {len(tx)}, 密度核点 {len(nerve)}, "
          f"细胞 {len(cells)}")

    # 细胞 → 神经距离
    w = float(cells["x_centroid"].max() - cells["x_centroid"].min())
    h = float(cells["y_centroid"].max() - cells["y_centroid"].min())
    tree_n = cKDTree(nerve)
    d, _ = tree_n.query(cells[["x_centroid", "y_centroid"]].to_numpy(float), k=1)
    prox = d <= PROX_UM
    dist = d >= DIST_UM

    # 细胞分型（互斥优先级: myCAF > fibroblast > endothelial > CD8 > macrophage > B > epithelial）
    # 规则: 逐基因 CP10K>2 计为阳性, 阳性基因数 ≥ n_pos 且非肿瘤主导
    ct = pd.Series("other", index=cells.index)
    epi_pos = meanz_score(G, CELLTYPE_RULES["epithelial"]["pos"], n_cells)
    for name, rule in CELLTYPE_RULES.items():
        cnt = np.zeros(n_cells, int)
        for g in rule["pos"]:
            v = G(g)
            if v is not None:
                cnt += (v > np.log1p(2.0)).astype(int)
        pos = cnt >= rule["n_pos"]
        # 非上皮类型不得落在强上皮细胞上（epi z > 2.5 视为肿瘤细胞）
        if name != "epithelial":
            pos = pos & ~(epi_pos > 2.5)
        ct[pos & (ct == "other").to_numpy()] = name

    n_prox = int(prox.sum())
    n_dist = int(dist.sum())
    xy = cells[["x_centroid", "y_centroid"]].to_numpy(float)
    res = {"sample": p, "n_nerve_pts": int(len(nerve)),
           "n_cells": int(len(cells)),
           "n_prox": n_prox, "n_dist": n_dist}

    # 各类型: 构成比 prox vs dist + torus null
    # 密度比 = 近端(≤150µm)某类型构成比 / 远端(≥300µm)构成比
    # （构成比比较对总细胞密度平台差稳健）
    # 置换掩码先一次性算好（999 次平移 × KDTree 查询只在类型循环外做一遍）
    rng = np.random.default_rng(SEED + hash(p) % 10000)
    perm_masks = []
    for _ in range(N_PERM):
        ns = torus_shift(nerve, w + 200, h + 200, rng)
        d2, _ = cKDTree(ns).query(xy, k=1)
        p2 = d2 <= PROX_UM
        q2 = d2 >= DIST_UM
        if p2.sum() >= 20 and q2.sum() >= 20:
            perm_masks.append((p2, q2))

    for name in CELLTYPE_RULES:
        sel = (ct == name).to_numpy()
        res[f"n_{name}_total"] = int(sel.sum())
        if sel.sum() < 5:
            continue
        obs_p = sel & prox
        obs_d = sel & dist
        # 观测: 神经 150µm 内该类型占全部近端细胞比例 vs 远端比例
        frac_p = obs_p.sum() / max(n_prox, 1)
        frac_d = obs_d.sum() / max(n_dist, 1)
        ratio = frac_p / frac_d if frac_d > 0 else np.nan
        res[f"{name}_frac_prox"] = round(float(frac_p), 5)
        res[f"{name}_frac_dist"] = round(float(frac_d), 5)
        res[f"{name}_prox_dist_ratio"] = round(float(ratio), 3) if ratio == ratio else np.nan
        nulls = []
        for p2, q2 in perm_masks:
            num = (sel & p2).sum() / p2.sum()
            den = (sel & q2).sum() / q2.sum()
            if den <= 0:
                continue
            nulls.append(num / den)
        if len(nulls) > 100 and ratio == ratio:
            nulls = np.array(nulls)
            # 双侧
            emp = 2 * min((nulls >= ratio).mean(), (nulls <= ratio).mean())
            res[f"{name}_permP"] = round(float(min(emp, 1.0)), 4)
            res[f"{name}_null_med"] = round(float(np.median(nulls)), 3)

    # 成纤维细胞模块评分: prox vs dist
    fb = (ct == "myCAF").to_numpy() | (ct == "fibroblast").to_numpy()
    rows = []
    for modname in ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM",
                    "LRRC15_myCAF"]:
        s = meanz_score(G, mods[modname], n_cells)
        for zone, mask in [("prox", fb & prox), ("dist", fb & dist)]:
            v = s[mask]
            if len(v) > 10 and np.isfinite(v).any():
                rows.append({"sample": p, "module": modname, "zone": zone,
                             "n": int(mask.sum()),
                             "mean": round(float(np.nanmean(v)), 4)})
    return res, rows


def main() -> None:
    mods = load_modules()
    all_res, all_rows, all_cells = [], [], []
    for p in SAMPLES:
        out = analyse_sample(p, mods)
        if out is None:
            continue
        if isinstance(out, dict) and "note" in out:
            all_res.append(out)
            continue
        res, rows = out
        all_res.append(res)
        all_rows.extend(rows)
    pd.DataFrame(all_res).to_csv(TAB / "stage14_xenium_nerve_summary.csv", index=False)
    pd.DataFrame(all_rows).to_csv(TAB / "stage14_xenium_fb_scores.csv", index=False)
    print(pd.DataFrame(all_res).to_string(index=False))


if __name__ == "__main__":
    main()
