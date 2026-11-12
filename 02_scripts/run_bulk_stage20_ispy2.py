"""Stage 20 (bulk): I-SPY2 chemo vs IO+chemo pCR 分层 (GSE194040, n=987).

用户第四轮方案第 5 条: "stroma-high 在 chemo-only 里 pCR 低; 在 IO+chemo 里差异收窄".
I-SPY2-990 数据资源: 987 例/10 臂 (Control 210, pembro 69, 其余 8 个非 IO 研究臂),
pre-treatment 表达 (Agilent 44k, ComBat 合并基因级), pCR + HR/HER2/MP 标签。

分析:
  1. TNBC (HR-/HER2-) 全体: stroma 模块 pCR logistic (age/race/HR 调整按可得性)
  2. 臂分层: Control (chemo-only) vs pembrolizumab (IO+chemo) — 模块 × 臂交互
  3. 对照: CD8 细胞毒 (预期阳性) + neural 模块
  4. 随机集基线 (510 组) 在 TNBC 全体层
口径: hypothesis-generating, cross-regimen (I-SPY2 各臂为不同试验方案, 非随机头对头),
      pembro 臂 n=69 中 TNBC 子集更小 — 全部标注探索性。

输出: 03_results/tables/bulk_extrapolation/stage20_ispy2_pcr.csv
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import bulk_lib as BL                                    # noqa: E402

RAW = Path(r"D:/嗜神经浸润/01_data/raw/gse194040")
TAB = Path(r"D:/嗜神经浸润/03_results/tables/bulk_extrapolation")
SEED = 20261120
N_RANDOM = 510


def parse_gsm(path: Path) -> pd.DataFrame:
    """family.soft.gz 全量样本 → 临床表 (patient, arm, hr, her2, pcr, mp)。"""
    import gzip
    t = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows = []
    for b in t.split("^SAMPLE = ")[1:]:
        gsm = b.split("\n")[0].strip()
        d = {}
        for line in b.splitlines():
            m = re.match(r"!Sample_characteristics_ch1 = (.+?): (.*)", line.strip())
            if m:
                d.setdefault(m.group(1).strip(), []).append(m.group(2).strip())
        rows.append(dict(
            gsm=gsm,
            patient=(d.get("patient id") or [None])[0],
            arm=(d.get("arm") or [None])[0],
            hr=(d.get("hr") or [None])[0],
            her2=(d.get("her2") or [None])[0],
            pcr=(d.get("pcr") or [None])[0],
            mp=(d.get("mp") or [None])[0],
        ))
    return pd.DataFrame(rows)


def load_expr() -> pd.DataFrame:
    """基因级矩阵: 行=基因, 列=ResID → 样本×基因。"""
    df = pd.read_csv(RAW / "GSE194040_gene_n988.txt.gz", sep="\t", index_col=0)
    # 探测结构: 第一列可能为基因名
    if df.shape[0] > df.shape[1]:      # 行=基因
        idx_name = df.index.name or "gene"
        df = df.T
    return df.astype(float)


def score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    avail = [g for g in genes if g in expr.columns]
    if not avail:
        return pd.Series(np.nan, index=expr.index)
    sub = expr[avail]
    z = (sub - sub.mean()) / (sub.std(ddof=0) + 1e-12)
    return z.mean(axis=1)


def logistic_or(x: np.ndarray, y: np.ndarray) -> tuple:
    """简单 logistic: 手写 IRLS, 返回 OR, CI, p, AUC(leave-one-out 近似用 in-sample)。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    X = np.column_stack([np.ones(len(x)), x])
    X = X[np.isfinite(X).all(1)]
    m = np.isfinite(y)
    Xf, yf = X[m], y[m]
    if yf.sum() < 3 or (1 - yf).sum() < 3:
        return (np.nan,) * 5
    try:
        lr = LogisticRegression(C=1e6, max_iter=500).fit(Xf[:, 1:2], yf)
        b = np.r_[lr.intercept_, lr.coef_[0]]
        # 标准误 via Hessian
        p = 1 / (1 + np.exp(-(Xf @ b)))
        W = p * (1 - p)
        H = Xf.T @ (Xf * W[:, None])
        cov = np.linalg.pinv(H)
        se = np.sqrt(np.diag(cov))[1]
        z = b[1] / se
        pval = 2 * stats.norm.sf(abs(z))
        orv = float(np.exp(b[1]))
        lo, hi = float(np.exp(b[1] - 1.96 * se)), float(np.exp(b[1] + 1.96 * se))
        auc = float(roc_auc_score(yf, Xf[:, 1]))
        return orv, lo, hi, pval, auc
    except Exception:
        return (np.nan,) * 5


def main() -> None:
    mods = BL.load_modules()
    clin = parse_gsm(RAW / "GSE194040_family.soft.gz")
    expr = load_expr()
    print(f"[data] GSM {len(clin)}, expr {expr.shape}")
    print(f"[arms] {clin.arm.value_counts().to_dict()}")
    print(f"[pcr] {clin.pcr.value_counts().to_dict()}")

    # 样本对齐: 表达矩阵 index = ISPY2 ResID (患者号)
    idx = set(expr.index.astype(str))
    clin["patient"] = clin["patient"].astype(str)
    hit = clin["patient"].isin(idx).mean()
    print(f"[align] patient->expr.index hit {hit:.1%}")
    clin = clin.drop_duplicates("patient").set_index("patient").reindex(expr.index.astype(str))
    expr.index = expr.index.astype(str)
    print(f"[aligned] n with arm = {clin.arm.notna().sum()}, pcr = {clin.pcr.notna().sum()}")

    # pCR 二值化
    pcr_raw = clin.pcr.astype(str).str.strip().str.lower()
    y = pcr_raw.map({"yes": 1, "1": 1, "no": 0, "0": 0, "pcr": 1, "non-pcr": 0, "rd": 0})
    hr = clin.hr.astype(str).str.strip()
    her2 = clin.her2.astype(str).str.strip()
    tnbc = (hr.isin(["0", "negative", "neg"])) & (her2.isin(["0", "negative", "neg"]))
    arm = clin.arm.astype(str).str.strip().str.lower()

    stroma = score(expr, mods["Perineural_stroma"])
    neural = score(expr, mods["Neural_intrinsic"])
    cd8 = score(expr, mods["CD8_cytotox"])
    caf = score(expr, mods["CAF_ECM"])

    def z(s):
        return (s - s.mean()) / (s.std(ddof=0) + 1e-12)

    rows = []
    # --- 层 1: TNBC 全体 (所有臂) ---
    for label, mask in [("TNBC_all", tnbc), ("TNBC_chemo_ctl", tnbc & (arm == "paclitaxel")),
                        ("TNBC_pembro", tnbc & arm.str.contains("pembrolizumab"))]:
        d = mask & y.notna() & stroma.notna()
        n, k = int(d.sum()), int(y[d].sum())
        print(f"\n[{label}] n={n} pCR={k}")
        if n < 20 or k < 3:
            continue
        for rname, s in [("Perineural_stroma", stroma), ("Neural_intrinsic", neural),
                         ("CAF_ECM", caf), ("CD8_cytotox", cd8)]:
            orv, lo, hi, p, auc = logistic_or(z(s[d].to_numpy()), y[d].to_numpy().astype(float))
            rows.append(dict(stratum=label, readout=rname, n=n, n_pcr=k,
                             OR=orv, ci_lo=lo, ci_hi=hi, p=p, AUC=auc))
            if orv == orv:
                print(f"  {rname:20s} OR={orv:.2f} ({lo:.2f}-{hi:.2f}) P={p:.3g} AUC={auc:.2f}")

    # --- 层 2: 交互 (TNBC, ctl vs pembro) ---
    d = tnbc & y.notna() & stroma.notna() & (arm.str.contains("pembrolizumab") | (arm == "paclitaxel"))
    if d.sum() > 60:
        io = arm[d].str.contains("pembrolizumab").astype(float).to_numpy()
        yv = y[d].to_numpy().astype(float)
        zs = z(stroma[d]).to_numpy()
        # 交互 logistic: pcr ~ stroma * pembro
        from sklearn.linear_model import LogisticRegression
        X = np.column_stack([np.ones(d.sum()), zs, io, zs * io])
        W = np.diag(np.ones(d.sum()))
        # IRLS 手写
        b = np.zeros(4)
        for _ in range(100):
            eta = X @ b
            p = 1 / (1 + np.exp(-eta))
            Wm = p * (1 - p)
            g = X.T @ (yv - p)
            H = X.T @ (X * Wm[:, None])
            step = np.linalg.pinv(H) @ g
            b = b + step
            if np.max(np.abs(step)) < 1e-8:
                break
        se = np.sqrt(np.diag(np.linalg.pinv(X.T @ (X * (p * (1 - p))[:, None]))))
        zint = b[3] / se[3]
        pint = 2 * stats.norm.sf(abs(zint))
        rows.append(dict(stratum="TNBC_interaction", readout="stroma_x_pembro",
                         n=int(d.sum()), n_pcr=int(yv.sum()),
                         OR=float(np.exp(b[3])), ci_lo=float(np.exp(b[3] - 1.96 * se[3])),
                         ci_hi=float(np.exp(b[3] + 1.96 * se[3])), p=float(pint), AUC=np.nan))
        print(f"\n[interaction] stroma x pembro (TNBC): OR={np.exp(b[3]):.2f} "
              f"({np.exp(b[3]-1.96*se[3]):.2f}-{np.exp(b[3]+1.96*se[3]):.2f}) P={pint:.3g}")
        # 主效应 (臂内 OR 差)
        for a, an in [(0, "chemo_ctl"), (1, "pembro")]:
            m = (io == a)
            if m.sum() > 15 and yv[m].sum() > 2:
                orv, lo, hi, p, auc = logistic_or(zs[m], yv[m])
                rows.append(dict(stratum=f"TNBC_{an}_only", readout="Perineural_stroma",
                                 n=int(m.sum()), n_pcr=int(yv[m].sum()),
                                 OR=orv, ci_lo=lo, ci_hi=hi, p=p, AUC=auc))
                print(f"  [{an}] stroma OR={orv:.2f} ({lo:.2f}-{hi:.2f}) P={p:.3g}")

    # --- 层 3: 随机基线 (TNBC 全体) ---
    d = tnbc & y.notna()
    if d.sum() > 60:
        rng = np.random.default_rng(SEED)
        exclude = set()
        for v in mods.values():
            exclude.update(v)
        pool = np.array([g for g in expr.columns if g not in exclude])
        yv = y[d].to_numpy().astype(float)
        sizes = [len([g for g in mods['Perineural_stroma'] if g in expr.columns])]
        aucs = []
        for _ in range(N_RANDOM):
            pick = rng.choice(pool, sizes[0], replace=False)
            s = score(expr, list(pick))
            if np.isnan(s[d]).all():
                continue
            _, _, _, _, auc = logistic_or(z(s[d].to_numpy()), yv)
            if auc == auc:
                aucs.append(auc)
        obs_rows = [r for r in rows if r.get("stratum") == "TNBC_all" and r.get("readout") == "Perineural_stroma"]
        if aucs and obs_rows:
            obs_auc = obs_rows[0]["AUC"]
            emp = float((1 + int((np.abs(np.array(aucs) - 0.5) >= abs(obs_auc - 0.5)).sum())) / (1 + len(aucs)))
            print(f"\n[random] {len(aucs)} sets, median AUC {np.median(aucs):.3f}; "
                  f"observed stroma AUC {obs_auc:.3f}, empirical P={emp:.3g}")
            rows.append(dict(stratum="TNBC_all_random", readout="Perineural_stroma",
                             n=int(d.sum()), n_pcr=int(yv.sum()), OR=np.nan, ci_lo=np.nan,
                             ci_hi=np.nan, p=emp, AUC=float(np.median(aucs))))

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "stage20_ispy2_pcr.csv", index=False)
    print("\nsaved:", TAB / "stage20_ispy2_pcr.csv")


if __name__ == "__main__":
    main()
