# -*- coding: utf-8 -*-
"""公共 bulk 队列外部外推分析库 —— Stage 1/2/4/6 共用。

数据契约
--------
expr : DataFrame, index = sample_id, columns = HGNC symbol, float (log 空间或 TPM)
clin : DataFrame, index = sample_id

模块评分一律走 `score_diff` / `score_mean`（mean-z），跨队列可比。
"""
from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"D:/嗜神经浸润")
RAW = ROOT / "01_data" / "raw"
PROC = ROOT / "01_data" / "processed"
TAB = ROOT / "03_results" / "tables" / "bulk_extrapolation"
TAB.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------- 模块定义
def load_modules() -> dict[str, list[str]]:
    """从稿件权威表读取模块基因，避免代码与文稿不一致。"""
    d = pd.read_csv(ROOT / "03_results" / "tables" / "supplementary"
                    / "TableS8_module_definitions.csv")
    mod = {r["module"]: r["genes"].split(";") for _, r in d.iterrows()}
    with open(PROC / "npp_signature.json", encoding="utf-8") as f:
        sig = json.load(f)
    out = {
        "NPP_up": sig["up"],
        "NPP_down": sig["down"],
        "Neural_intrinsic": mod["Neural_intrinsic"],
        "Perineural_stroma": mod["Perineural_stroma"],
        # 稿件 Methods 2.x 的 curated CAF/ECM（14 基因）
        "CAF_ECM": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "POSTN", "ACTA2",
                    "TAGLN", "FAP", "PDGFRB", "SPARC", "THBS2", "CTHRC1", "LOX"],
        "CD8_cytotox": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "IFNG"],
        # 免疫抑制读数（Stage 2）
        "Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
        "M2_macro": ["CD163", "MRC1", "MSR1", "IL10", "TGFB1", "CCL18"],
        "SPP1_macro": ["SPP1", "MMP9", "VEGFA", "TREM2"],
        "TGFb_core": ["TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SERPINE1", "JUNB"],
        "CXCL12_axis": ["CXCL12", "CXCR4", "ACKR3"],
        "LGALS1_axis": ["LGALS1", "LGALS9", "HAVCR2"],
        "IFN_gamma": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "IDO1"],
        "TLS": ["CCL21", "CCL19", "CXCL13", "LAMP3", "MS4A1", "CD79A"],
        "MHC_II": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "CIITA"],
        "LRRC15_myCAF": ["LRRC15", "SPARC", "POSTN", "ACTA2", "TNC", "COL1A1", "FAP", "PDGFRB"],
    }
    return out


# ------------------------------------------------- 纤维化程序（Stage 8 定量）
# 与 CAF_ECM（14 基因"总量"）正交设计：这里按"纤维化的质"拆分程序，
# 用来回答 desmoplasia 是"胶原多"还是"胶原成熟/交联/僵硬"。
PROGRAMS: dict[str, list[str]] = {
    "COLLAGEN_CORE": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2"],
    "COLLAGEN_MATURATION": ["LOX", "LOXL1", "LOXL2", "PLOD1", "PLOD2", "P4HA1", "P4HA2",
                            "SERPINH1", "ADAMTS2", "BMP1", "COL11A1", "COL12A1"],
    "SPARC_POSTN_CORE": ["SPARC", "POSTN", "THBS1", "THBS2", "TNC", "FN1", "CTHRC1", "SPP1"],
    "MECHANOTRANSDUCTION": ["YAP1", "WWTR1", "TAGLN", "ACTA2", "MYL9", "LIMS1", "ILK", "FLNA"],
    "EMT_CORE": ["VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "MMP2", "MMP14"],
    "TGFb_PROGRAM": ["TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD3", "SKIL", "SERPINE1"],
}

# 肿瘤细胞含量代理（纯度校正用；不含增殖基因，避免与侵袭性混淆）
PURITY_PROXY = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "MUC1"]


# ------------------------------------------------- 反卷积签名矩阵（Stage 10）
# 为什么不是 CIBERSORTx / xCell: 前者需 token + 网页上传，后者需 R 包；
# 本机两者都没有，且 raw.githubusercontent 在本网络被阻断（curl exit 7），
# LM22 原始矩阵无法取得。这里用 12 区室、每区室 4–7 个 canonical marker 的
# 自建矩阵做 NNLS —— 与 CIBERSORTx 同为「线性混合 + 最小二乘」框架，差别在
# 签名来源与 batch 校正（S-mode）。所有输出一律标注 deconv_proxy，不冒充 CIBERSORTx。
DECONV_PANEL: dict[str, list[str]] = {
    "CD8_T":       ["CD8A", "CD8B", "GZMK", "GZMA", "PRF1", "NKG7", "CD3D", "CD3E"],
    "CD4_T":       ["CD4", "IL7R", "CCR7", "LTB", "CD40LG", "ICOS", "CD2"],
    "Treg":        ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
    "NK":          ["KLRD1", "KLRF1", "NCAM1", "GNLY", "SH2D1B", "NCR1"],
    "B_plasma":    ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "IGKC"],
    "M1_macro":    ["NOS2", "CXCL9", "CXCL10", "IDO1", "SOCS3", "IL1B"],
    "M2_macro":    ["CD163", "MRC1", "MSR1", "CCL18", "IL10", "MRC2"],
    "DC":          ["CD1C", "LAMP3", "CLEC9A", "ITGAX", "HLA-DRA", "HLA-DPA1"],
    "Neutrophil":  ["FCGR3B", "CSF3R", "S100A8", "S100A9", "ELANE", "MPO"],
    "Mast":        ["TPSAB1", "CPA3", "MSA2", "KIT"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "PLVAP"],
    "Fibroblast":  ["COL1A1", "DCN", "LUM", "PDGFRB", "FAP", "COL3A1"],
    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
}
DECONV_PANEL["Mast"] = ["TPSAB1", "CPA3", "MS4A2", "KIT"]

# celltype_minor（GSE176078 注释）→ 反卷积区室
# 注意: 该参考集没有独立 Treg / M1 / M2 / Mast / Neutrophil 标签，
# 这五类不进签名矩阵（Treg 与 M2 另用 marker 打分读数，见 load_modules）。
MINOR_TO_COMPARTMENT: dict[str, str] = {
    "T cells CD8+": "CD8_T", "Cycling T-cells": "CD8_T",
    "T cells CD4+": "CD4_T",
    "NK cells": "NK", "NKT cells": "NK",
    "B cells Memory": "B_plasma", "B cells Naive": "B_plasma",
    "Plasmablasts": "B_plasma",
    "Macrophage": "Macrophage", "Cycling_Myeloid": "Monocyte",
    "Monocyte": "Monocyte", "DCs": "DC",
    "Endothelial ACKR1": "Endothelial", "Endothelial CXCL12": "Endothelial",
    "Endothelial RGS5": "Endothelial", "Endothelial Lymphatic LYVE1": "Endothelial",
    "CAFs myCAF-like": "Fibroblast", "CAFs MSC iCAF-like": "Fibroblast",
    "PVL Differentiated": "Fibroblast", "PVL Immature": "Fibroblast",
    "Cycling PVL": "Fibroblast",
    "Cancer LumA SC": "Epithelial", "Cancer LumB SC": "Epithelial",
    "Cancer Basal SC": "Epithelial", "Cancer Her2 SC": "Epithelial",
    "Cancer Cycling": "Epithelial", "Luminal Progenitors": "Epithelial",
    "Mature Luminal": "Epithelial", "Myoepithelial": "Epithelial",
}


def build_signature(scrna_root: Path, genes: list[str] | None = None,
                    out_csv: Path | None = None) -> pd.DataFrame:
    """从 GSE176078 单细胞参考构建区室 × 基因签名矩阵（CP10K 均值）。

    与 CIBERSORTx 的差别只在签名来源（自建 vs LM22）与是否做 batch 校正；
    混合模型与求解器（NNLS）同构。缺失区室不补零，直接缺列。
    """
    from scipy import sparse
    from scipy.io import mmread

    genes = genes or sorted({g for v in DECONV_PANEL.values() for g in v})
    acc: dict[str, dict[str, list[float]]] = {}
    for d in sorted([p for p in scrna_root.iterdir()
                     if p.is_dir() and (p / "count_matrix_sparse.mtx").exists()]):
        with open(d / "count_matrix_sparse.mtx", "rb") as fh:
            X = sparse.csc_matrix(mmread(fh))                       # genes × cells
        gl = pd.read_csv(d / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
        gi = {g: i for i, g in enumerate(gl)}
        keep = [(g, gi[g]) for g in genes if g in gi]
        meta = pd.read_csv(d / "metadata.csv", index_col=0)
        meta.index = meta.index.astype(str)
        bar = pd.read_csv(d / "count_matrix_barcodes.tsv", header=None)[0].astype(str).tolist()
        if len(bar) != X.shape[1]:
            continue
        meta = meta.reindex(bar)
        tot = np.asarray(X.sum(axis=0)).ravel().astype(float)
        tot[tot <= 0] = np.nan
        idx = [i for _, i in keep]
        sub = np.asarray(X[idx, :].todense(), float) / tot * 1e4      # 基因 × 细胞
        ct = meta["celltype_minor"].astype(str).map(MINOR_TO_COMPARTMENT)
        for comp in pd.unique(ct.dropna()):
            mask = (ct == comp).to_numpy()
            if mask.sum() < 10:
                continue
            prof = np.nanmean(sub[:, mask], axis=1)
            acc.setdefault(comp, {})
            for (g, _), v in zip(keep, prof):
                acc[comp].setdefault(g, []).append(float(v))
    S = pd.DataFrame({c: {g: np.mean(v) for g, v in d.items()} for c, d in acc.items()}).T
    S.index.name = "compartment"
    if out_csv is not None:
        S.to_csv(out_csv)
    return S


def specific_markers(S: pd.DataFrame, min_val: float = 0.5,
                     min_ratio: float = 2.0) -> dict[str, list[str]]:
    """从签名矩阵挑「区室特异」基因：本位室 CP10K ≥ min_val 且 ≥ 次高区室的 min_ratio 倍。

    这是 xCell / Danaher 式富集打分的前置步骤，不是反卷积。
    """
    arr = S.to_numpy(float)
    out: dict[str, list[str]] = {}
    for j, g in enumerate(S.columns):
        col = arr[:, j]
        i = int(np.nanargmax(col))
        v = col[i]
        sec = np.sort(np.delete(col, i))[-1]
        if v >= min_val and v / (sec + 1e-6) >= min_ratio:
            out.setdefault(S.index[i], []).append(g)
    return out


def compartment_scores(expr: pd.DataFrame, S: pd.DataFrame) -> pd.DataFrame:
    """区室富集分（marker 均值 z），附每个区室用到的基因数。

    为什么不做 NNLS 比例: 256 基因的靶向面板（ECM/神经为主）对 10 区室线性
    混合不可识别——实测三种标度化方案残差 R² < 0.02，权重退化为均匀分布。
    富集分只要求「标记基因相对上调」，不受混合可识别性约束，
    其与真实细胞占比的对应关系由配对单细胞 pseudobulk 直接验证（见 stage10）。
    """
    mk = specific_markers(S)
    out = {}
    for comp, genes in mk.items():
        s = score_mean(expr, genes, min_genes=2)
        if s.notna().any():
            out[comp] = s
    df = pd.DataFrame(out)
    df.attrs["markers"] = mk
    return df


# ----------------------------------------------------------------- 评分
def _z(df: pd.DataFrame) -> pd.DataFrame:
    return (df - df.mean()) / df.std(ddof=0).replace(0, np.nan)


def _sub(expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    keep = [g for g in genes if g in expr.columns]
    return expr[keep]


def score_mean(expr: pd.DataFrame, genes: list[str], min_genes: int = 4) -> pd.Series:
    sub = _sub(expr, genes)
    if sub.shape[1] < min_genes:
        return pd.Series(np.nan, index=expr.index, name="nan")
    return _z(sub).mean(axis=1)


def score_diff(expr: pd.DataFrame, up: list[str], down: list[str] | None = None,
               min_genes: int = 5) -> pd.Series:
    zu = _sub(expr, up)
    if zu.shape[1] < min_genes:
        raise ValueError(f"up 仅剩 {zu.shape[1]} 基因")
    s = _z(zu).mean(axis=1)
    zd = _sub(expr, down or [])
    if zd.shape[1]:
        s = s - _z(zd).mean(axis=1)
    return s


def score_all(expr: pd.DataFrame, mods: dict[str, list[str]]) -> pd.DataFrame:
    out = {}
    out["NPP58"] = score_diff(expr, mods["NPP_up"], mods["NPP_down"])
    for k in mods:
        if k.startswith("NPP"):
            continue
        s = score_mean(expr, mods[k])
        if s.notna().any():
            out[k] = s
    df = pd.DataFrame(out)
    return df


# ----------------------------------------------------------------- 统计工具
def partial_spearman(x, y, z, method="spearman") -> tuple[float, float, int]:
    """x~y 控制 z（Series 或 DataFrame）后的偏相关。纯 numpy 实现，避开 pandas 索引陷阱。"""
    df = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float)})
    if isinstance(z, pd.DataFrame):
        zcols = []
        for i, c in enumerate(z.columns):
            df[f"z{i}"] = np.asarray(z[c], float)
            zcols.append(f"z{i}")
    else:
        df["z0"] = np.asarray(z, float)
        zcols = ["z0"]
    df = df.dropna()
    if len(df) < 10:
        return np.nan, np.nan, len(df)
    if method == "spearman":
        df = df.rank()
    xv = df["x"].to_numpy(float)
    yv = df["y"].to_numpy(float)
    Z = np.column_stack([np.ones(len(df))] + [df[c].to_numpy(float) for c in zcols])

    def resid(v):
        b, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ b

    r, p = stats.pearsonr(resid(xv), resid(yv))
    return float(r), float(p), len(df)


def boot_ci(x: pd.Series, y: pd.Series, n=1000, seed=20251108, alpha=0.05):
    rng = np.random.default_rng(seed)
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    n_ = len(d)
    rs = np.empty(n)
    xa, ya = d["x"].to_numpy(), d["y"].to_numpy()
    for i in range(n):
        idx = rng.integers(0, n_, n_)
        rs[i] = stats.spearmanr(xa[idx], ya[idx])[0]
    return float(np.nanpercentile(rs, 100 * alpha / 2)), float(np.nanpercentile(rs, 100 * (1 - alpha / 2)))


# ------------------------------------------------------- 正交化（Stage 7 核心）
def rank_resid(y: pd.Series, x: pd.Series | pd.DataFrame) -> tuple[pd.Series, float]:
    """把 y 对 x 做秩空间 OLS，返回标准化残差与 R²。

    秩变换后做 OLS，等价于"扣除 x 的单调可解释部分"，残差即
    「给定纤维化水平后，PNI 模块仍剩余的那部分」。
    """
    yy = y.rank(pct=True)
    if isinstance(x, pd.DataFrame):
        X = x.rank(pct=True)
    else:
        X = x.rank(pct=True).to_frame()
    d = pd.concat([yy.rename("_y"), X], axis=1).dropna()
    if len(d) < 20:
        return pd.Series(np.nan, index=y.index), np.nan
    A = np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in X.columns])
    b, *_ = np.linalg.lstsq(A, d["_y"].to_numpy(float), rcond=None)
    fit = A @ b
    resid = d["_y"].to_numpy(float) - fit
    ss_res = float(resid @ resid)
    ss_tot = float(((d["_y"] - d["_y"].mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    resid = (resid - resid.mean()) / (resid.std(ddof=0) or np.nan)
    return pd.Series(resid, index=d.index).reindex(y.index), float(r2)


def cox_breslow(t, d, X, max_iter=80, tol=1e-10):
    """Breslow ties 的 Cox 比例风险，Newton-Raphson 自实现。

    为什么不用 statsmodels.PHReg: 当效应量接近 0（OR≈1.00）时其 Hessian
    数值退化，bse 返回 0/NaN，会产出 "HR=1.001, p=0" 这类伪显著。
    X 需已标准化（HR = per SD）。
    """
    t = np.asarray(t, float)
    d = np.asarray(d, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    order = np.argsort(-t, kind="mergesort")      # 时间降序 → 风险集为前缀累积
    Xo, to, do = X[order], t[order], d[order]
    n, p = Xo.shape
    b = np.zeros(p)
    H = None
    for _ in range(max_iter):
        eta = np.clip(Xo @ b, -40, 40)
        w = np.exp(eta)
        S0 = np.cumsum(w)
        SW = np.cumsum(w[:, None] * Xo, axis=0)
        xbar = SW / S0[:, None]
        # d_loglik = sum d_i (x_i - xbar_i)
        g = (do[:, None] * (Xo - xbar)).sum(axis=0)
        WXX = np.cumsum(w[:, None, None] * Xo[:, :, None] * Xo[:, None, :], axis=0)
        A = WXX / S0[:, None, None]
        H = -np.einsum("i,ijk->jk", do, A - xbar[:, :, None] * xbar[:, None, :])
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, g, rcond=None)[0]
        step = np.clip(step, -5, 5)
        b = b - step
        if np.max(np.abs(step)) < tol:
            break
    try:
        cov = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(-H)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    z = np.where(se > 0, b / np.where(se > 0, se, 1), np.nan)
    pv = 2 * (1 - stats.norm.cdf(np.abs(z)))
    return b, se, pv


def cox_hr(t, d, x, cov=None):
    """返回 (HR, lo, hi, p) per SD；x 已在调用方标准化。"""
    X = np.asarray(x, float).reshape(-1, 1)
    if cov is not None and np.size(cov):
        X = np.column_stack([X, np.asarray(cov, float)])
    b, se, pv = cox_breslow(t, d, X)
    return float(np.exp(b[0])), float(np.exp(b[0] - 1.96 * se[0])), \
        float(np.exp(b[0] + 1.96 * se[0])), float(pv[0])


def logrank(t, e, g) -> tuple[float, float]:
    """两组 log-rank（Mantel-Cox）。g 为 0/1。"""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    g = np.asarray(g, int)
    times = np.unique(t[e == 1])
    O1 = E1 = V = 0.0
    for tt in times:
        at = t >= tt
        n = at.sum()
        n1 = (at & (g == 1)).sum()
        dd = ((t == tt) & (e == 1)).sum()
        d1 = ((t == tt) & (e == 1) & (g == 1)).sum()
        if n <= 1:
            continue
        O1 += d1
        E1 += dd * n1 / n
        V += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi = (O1 - E1) ** 2 / V
    return float(chi), float(1 - stats.chi2.cdf(chi, 1))


def fdr_bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    if not len(pv):
        return q
    order = np.argsort(pv)
    m = len(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


# ----------------------------------------------------------------- 队列加载
def read_gpl_annot(path: Path) -> pd.DataFrame:
    """GEO *.annot.gz -> DataFrame(probe, symbol)。

    GEO annot 文件头是若干 `#列名 = 说明` 行（#ID / #Gene title / #Gene symbol ...），
    数据行无表头，因此必须先收集全部 # 行构造列名。
    """
    header, rows = [], []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                head = line[1:].split("=")[0].strip()
                if head:
                    header.append(head)
                continue
            if line.startswith(("!", "^")) or not header:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rows.append(parts[: len(header)] + [""] * max(0, len(header) - len(parts)))
    df = pd.DataFrame(rows, columns=header)
    df = df.rename(columns={df.columns[0]: "probe"})
    sym = [c for c in df.columns if "symbol" in c.lower()]
    if not sym:
        raise ValueError(f"{path} 无 Gene symbol 列: {list(df.columns)[:6]}")
    df = df.rename(columns={sym[0]: "symbol"})
    df["symbol"] = df["symbol"].astype(str).str.split("///").str[0].str.strip()
    df = df[df["symbol"].notna() & ~df["symbol"].isin(["", "nan", "---", "None"])]
    df = df[df["probe"].astype(str).str.upper() != "ID"]
    return df[["probe", "symbol"]]


def series_matrix_expr(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """GEO series matrix -> (expr probe×sample, sample meta 宽表)。

    注意: !Sample_characteristics_ch1 会出现多行(每行一个变量)，必须逐行解析成
    "key: value" 再展开，否则只保留第一行会丢掉 pCR 标签。
    """
    meta_rows, table_lines = [], []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                meta_rows.append(line.rstrip("\n"))
    expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    expr.index.name = "probe"
    meta_rows = [r for r in meta_rows if r.strip()]
    meta_rows = [r for r in meta_rows if "table" not in r.split("\t")[0]]
    n = len(expr.columns)
    sm = {}
    for line in meta_rows:
        parts = line.split("\t")
        key = parts[0][len("!Sample_"):].strip('"')
        vals = [p.strip().strip('"') for p in parts[1:]]
        if len(vals) != n:
            continue
        if key == "characteristics_ch1":
            if not any(":" in v for v in vals):
                continue
            for i, v in enumerate(vals):
                if ":" not in v:
                    continue
                kk, vv = v.split(":", 1)
                sm.setdefault(kk.strip(), [None] * n)[i] = vv.strip()
        else:
            sm[key] = vals
    meta = pd.DataFrame(sm, index=expr.columns)
    return expr, meta


def parse_gsm_chars(path: Path) -> pd.DataFrame:
    """GEO acc.cgi?targ=gsm&form=text&view=brief 输出 -> 样本特性表。"""
    txt = path.read_text(encoding="utf-8", errors="replace")
    rows = []
    for b in txt.split("^SAMPLE = ")[1:]:
        d = {}
        for line in b.split("\n"):
            m = re.match(r"!Sample_(title|characteristics_ch1|geo_accession)\s*=\s*(.*)", line)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k == "characteristics_ch1" and ":" in v:
                kk, vv = v.split(":", 1)
                d[kk.strip()] = vv.strip()
            elif k in ("title", "geo_accession"):
                d[k] = v
        rows.append(d)
    return pd.DataFrame(rows)
