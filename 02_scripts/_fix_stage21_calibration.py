# -*- coding: utf-8 -*-
"""修正 Stage 21 轴梯度 meta 的随机基线 (named-agg bug 复发).

与 N13.22 同款问题: R.groupby("name").agg(size=("rho","size")) 把每随机集的
样本行数当作基因集大小, 导致 bucket 全部落到伪桶 29、P 全在 floor。
本脚本从 per_sample 表的 name 前缀 RAND|s{gsize}| 提取真实基因集大小,
重新分桶计算 baseline / excess / z_vs_random / p_bucket / p_emp_rho。

输出: stage21_axis_gradient_meta.corrected.csv (覆盖入稿用)
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TAB = Path(r"D:/嗜神经浸润/03_results/tables/spatial_tlscaf")

# 实测读数基因集大小 (与 build_blocks 一致)
MEMBERS = {
    "MOD|CGRP_axis": 4, "MOD|CGRP_core": 2, "MOD|NGF_TrkA": 4,
    "MOD|Sensory_nf": 4,
    "GENE|CALCA": 1, "GENE|RAMP1": 1, "GENE|CALCRL": 1, "GENE|CALCR": 1,
    "GENE|NGF": 1, "GENE|NTRK1": 1, "GENE|NGFR": 1, "GENE|TRPV1": 1,
    "GENE|SCN9A": 1, "GENE|SCN10A": 1, "GENE|S100B": 1,
}


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[o[i]] * n / (i + 1))
        q[o[i]] = prev
    return np.minimum(q, 1.0)


def main():
    P = pd.read_csv(TAB / "stage21_axis_gradient_per_sample.csv.gz")
    M0 = pd.read_csv(TAB / "stage21_axis_gradient_meta.csv")
    M0 = M0.set_index("name")

    # ---- 随机集: 每个集合跨样本 Fisher-z 平均 ----
    R = P[P.kind == "random"].copy()
    R["gsize"] = R["name"].str.extract(r"RAND\|s(\d+)\|")[0].astype(int)
    R["z"] = np.arctanh(np.clip(R["rho"].astype(float), -0.999, 0.999))
    per_set = R.groupby(["gsize", "name"]).agg(
        z_mean=("z", "mean"), n=("z", "size")).reset_index()
    per_set["mean_rho"] = np.tanh(per_set["z_mean"])

    zr_all = []
    buckets = {}
    for s, g in per_set.groupby("gsize"):
        br = np.abs(g["mean_rho"].to_numpy(float))
        mu, sd = br.mean(), br.std(ddof=1)
        buckets[int(s)] = g
        zr_all.extend(((np.abs(g["mean_rho"]) - mu) / (sd + 1e-12)).tolist())
    zr_all = np.sort(np.asarray(zr_all))
    sizes = np.array(sorted(buckets))
    print(f"真实 bucket: {sizes.tolist()}, 总随机集 {len(per_set)}")

    def nearest(s):
        return int(sizes[np.argmin(np.abs(np.log(sizes) - np.log(max(1.0, s))))])

    rows = []
    for name, r0 in M0.iterrows():
        s = MEMBERS.get(name, int(r0["size"]))
        b = nearest(s)
        g = buckets[b]
        br = np.abs(g["mean_rho"].to_numpy(float))
        ar = abs(float(r0["rho_mean"]))
        zr = (ar - br.mean()) / (br.std(ddof=1) + 1e-12)
        # 方向匹配经验 P: 同方向随机集比例
        dirmatch = (np.sign(float(r0["rho_mean"])) == np.sign(g["mean_rho"])) & \
                   (np.abs(g["mean_rho"]) >= ar)
        rows.append(dict(
            name=name, size=s, n_samples=int(r0["n_samples"]),
            rho_mean=float(r0["rho_mean"]), rho_lo=float(r0["rho_lo"]),
            rho_hi=float(r0["rho_hi"]), p_rho=float(r0["p_rho"]),
            stouffer_Z=float(r0["stouffer_Z"]), p_stouffer=float(r0["p_stouffer"]),
            bucket=b, n_random=int(len(g)),
            baseline_rho=float(br.mean()), excess_rho=ar - float(br.mean()),
            z_vs_random=float(zr),
            p_emp_rho=float((1 + int((zr_all >= zr).sum())) / (1 + len(zr_all))),
            p_bucket=float((1 + int((br >= ar).sum())) / (1 + len(br))),
            p_dir=float((1 + int(dirmatch.sum())) / (1 + len(br))),
        ))
    M = pd.DataFrame(rows)
    M["fdr_rho"] = bh(M["p_rho"].to_numpy())
    M["fdr_emp_rho"] = bh(M["p_emp_rho"].to_numpy())
    M["fdr_bucket"] = bh(M["p_bucket"].to_numpy())
    M = M.sort_values("rho_mean")
    out = TAB / "stage21_axis_gradient_meta.corrected.csv"
    M.to_csv(out, index=False)
    print(M[["name", "size", "n_samples", "rho_mean", "rho_lo", "rho_hi",
             "baseline_rho", "p_bucket", "p_emp_rho", "p_dir"]].to_string(index=False))
    print(f"\nsaved -> {out}")


if __name__ == "__main__":
    main()
