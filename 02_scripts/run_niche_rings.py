"""只重算"读数 x 距离环"的平均评分剖面 —— 用于画梯度曲线。

与 run_niche_spatial.py 的读数定义、基因过滤、标准化完全一致, 但不做 MSR
模拟(那一步很贵), 因此可以在几秒内跑完 30 个样本。

顺带修正: 最外环标签是 out_gt400(不是 r5_gt400), 之前写错导致该环缺失。

输出
----
01_data/processed/niche_ring_profiles.csv   样本 x 环 x 读数
03_results/tables/niche_ring_meta.csv       跨样本各环均值(直接画图)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
from run_msr_spatial import MIN_SPOT_FRAC             # noqa: E402
import niche_lib as NL                               # noqa: E402
from run_niche_spatial import build_blocks, RING_ORDER  # noqa: E402

DIR_H5 = CFG.DIR_PROC / "spot_level"
PROC = CFG.DIR_PROC
TAB = CFG.DIR_TAB


def sample_rings(path: Path, blocks):
    ad = sc.read_h5ad(path)
    obs = ad.obs
    if "ring_id" not in obs.columns:
        return None
    ring = obs["ring_id"].astype(str).to_numpy()
    tc = np.asarray(obs["total_counts"].to_numpy(float))
    X = ad.X
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)
    genes = np.asarray(ad.var_names)
    n = X.shape[0]
    del ad

    sel = (X > 0).mean(axis=0) >= MIN_SPOT_FRAC
    Xs = X[:, sel].astype(np.float64)
    genes = genes[sel]
    lib = tc.copy()
    lib[lib <= 0] = 1.0
    Xl = np.log1p(Xs / lib[:, None] * 1e4)
    del Xs
    Xl = (Xl - Xl.mean(axis=0)) / (Xl.std(axis=0) + 1e-12)

    S, meta, _ = NL.build_score_matrix(
        Xl, genes, blocks,
        min_frac_by_kind={"celltype": 0.5, "lr": 1.0, "pathway": 0.5})
    del Xl
    names = [m["name"] for m in meta]

    rows = []
    for rk in RING_ORDER:
        m = ring == rk
        if m.sum() < 5:
            continue
        row = dict(sample=path.stem, ring=rk, n_spot=int(m.sum()))
        row.update({nm: float(S[m, i].mean()) for i, nm in enumerate(names)})
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    blocks = build_blocks(None)
    print(f"[rings] celltype {len(blocks['celltype'])}, lr {len(blocks['lr'])}, "
          f"pathway {len(blocks['pathway'])}")
    files = sorted(DIR_H5.glob("*.h5ad"))
    out, t0 = [], time.time()
    for i, p in enumerate(files, 1):
        try:
            df = sample_rings(p, blocks)
        except Exception as e:
            print(f"  [skip] {p.stem}: {type(e).__name__}: {e}", flush=True)
            continue
        if df is None or not len(df):
            continue
        out.append(df)
        print(f"  [{i}/{len(files)}] {p.stem} {time.time()-t0:.0f}s", flush=True)
    R = pd.concat(out, ignore_index=True)
    R.to_csv(PROC / "niche_ring_profiles.csv", index=False)

    val = R.drop(columns=["sample", "n_spot"])
    M = val.groupby("ring").mean()
    M = M.reindex([r for r in RING_ORDER if r in M.index])
    M.insert(0, "n_spot", R.groupby("ring")["n_spot"].sum().reindex(M.index))
    M.to_csv(TAB / "niche_ring_meta.csv")
    print(f"\n[rings] {len(R)} 行, 环 {list(M.index)}")
    print(M.iloc[:, :8].round(3).to_string())
    return R, M


if __name__ == "__main__":
    main()
