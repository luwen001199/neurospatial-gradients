"""第三道门槛: 大小匹配随机基因集的经验校准。

为什么需要
----------
MSR 零分布回答的是"这个空间场是否比空间自相关随机的场更极端"。
但它**不回答**另一个问题: 一个同样大小、但毫无生物学意义的随机基因集,
是不是也能达到同样的效应量?

本队列的答案(见结果)是: 能, 而且很容易 —— 510 个随机基因集里 86.7% 同时
通过 "p<0.05 且 |rho|>=0.03"。原因是距离场 Moran's I≈0.97, 近乎确定性梯度,
且读数评分是成员 z 的均值 —— 基因越多, 均值场越平滑, 与平滑距离场的相关就
越高。这个"尺寸效应"与生物学无关。

因此只有**超过同尺寸随机集基线**的读数, 才能称为神经周围特异。

方法
----
对每个真实读数, 取其尺寸最近的随机集尺寸桶(30 个随机集), 用
    p_emp = (1 + #{random: |stat| >= |obs|}) / (1 + 30)
作为经验 p 值, 再跨读数做 BH-FDR。同时对 |mean_rho| 与 |Stouffer Z| 各做一遍。

输出
----
03_results/tables/niche_calibrated.csv  每个读数的基线、经验 p、FDR、是否通过
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
from nsgrad import config as CFG  # noqa: E402

TAB = CFG.DIR_TAB
FDR = 0.05


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
    R = pd.read_csv(TAB / "niche_random_null.csv")
    G = pd.read_csv(TAB / "niche_meta_summary.clean.csv", index_col=0)
    G["size"] = pd.to_numeric(G["size"], errors="coerce")

    buckets = {}
    for s, g in R.groupby("size"):
        buckets[int(s)] = dict(
            rho=np.sort(np.abs(g["mean_rho"].to_numpy(float))),
            Z=np.sort(np.abs(g["stouffer_Z"].to_numpy(float))))
    sizes = np.array(sorted(buckets))

    def nearest(s):
        if not np.isfinite(s):
            return None
        return int(sizes[np.argmin(np.abs(np.log(sizes) - np.log(max(1.0, s))))])

    # ---- 按尺寸标准化的合并零分布 ----
    # 每桶只有 30 个随机集, 直接比较会把经验 p 的下限卡在 1/31,
    # 无法分辨真正的强信号。改为: 桶内 z-score 后把 510 个随机集合并成
    # 一个零分布, 分辨率提高到 1/511。
    zr_all, zz_all = [], []
    for s, g in R.iterrows():
        b = int(g["size"])
        br = np.abs(R.loc[R["size"] == b, "mean_rho"].to_numpy(float))
        bz = np.abs(R.loc[R["size"] == b, "stouffer_Z"].to_numpy(float))
        zr_all.append((abs(g["mean_rho"]) - br.mean()) / (br.std(ddof=1) + 1e-12))
        zz_all.append((abs(g["stouffer_Z"]) - bz.mean()) / (bz.std(ddof=1) + 1e-12))
    zr_all = np.sort(np.asarray(zr_all))
    zz_all = np.sort(np.asarray(zz_all))

    rows = []
    for name, r in G.iterrows():
        b = nearest(float(r["size"])) if pd.notna(r["size"]) else None
        if b is None:
            continue
        br, bz = buckets[b]["rho"], buckets[b]["Z"]
        ar, az = abs(float(r["mean_rho"])), abs(float(r["stouffer_Z"]))
        zr = (ar - br.mean()) / (br.std(ddof=1) + 1e-12)
        zz = (az - bz.mean()) / (bz.std(ddof=1) + 1e-12)
        p_rho = float((1 + int((zr_all >= zr).sum())) / (1 + len(zr_all)))
        p_Z = float((1 + int((zz_all >= zz).sum())) / (1 + len(zz_all)))
        p_bucket = float((1 + int((br >= ar).sum())) / (1 + len(br)))
        rows.append(dict(
            name=name, kind=r["kind"], size=int(r["size"]), bucket=b,
            n_samples=int(r["n_samples"]),
            mean_rho=float(r["mean_rho"]),
            baseline_rho=float(np.mean(br)),
            excess_rho=ar - float(np.mean(br)),
            z_vs_random=float(zr),
            stouffer_Z=float(r["stouffer_Z"]),
            baseline_Z=float(np.mean(bz)),
            z_vs_random_Z=float(zz),
            p_emp_rho=p_rho, p_emp_Z=p_Z, p_bucket=p_bucket,
            circular=bool(r["circular"]),
            measurable=bool(r["measurable"]),
        ))
    C = pd.DataFrame(rows).set_index("name")
    C["fdr_emp_rho"] = bh(C["p_emp_rho"].to_numpy())
    C["fdr_emp_Z"] = bh(C["p_emp_Z"].to_numpy())
    C["pass_calibrated"] = C["measurable"] & (~C["circular"]) & \
        (C["fdr_emp_rho"] < FDR) & (C["fdr_emp_Z"] < FDR)
    C = C.sort_values("p_emp_rho")
    C.to_csv(TAB / "niche_calibrated.csv")

    print(f"读数 {len(C)} (celltype {int((C.kind=='celltype').sum())}, "
          f"lr {int((C.kind=='lr').sum())}, pathway {int((C.kind=='pathway').sum())})")
    print(f"\n通过三道门槛(可测 + 非循环 + 超随机基线): "
          f"{int(C.pass_calibrated.sum())} / {len(C)}")
    for kind in ["celltype", "lr", "pathway"]:
        s = C[C.kind == kind]
        print(f"  {kind:9s}: {int(s.pass_calibrated.sum())} / {len(s)}")
    print("\n=== 通过者(按经验 p 排序) ===")
    cols = ["kind", "size", "bucket", "n_samples", "mean_rho", "baseline_rho",
            "excess_rho", "z_vs_random", "stouffer_Z", "baseline_Z",
            "p_emp_rho", "fdr_emp_rho", "p_emp_Z", "fdr_emp_Z"]
    print(C[C.pass_calibrated][cols].round(4).to_string())
    print("\n=== 各类的最强 10 个(不论是否通过, 供正文挑选) ===")
    for kind in ["celltype", "lr", "pathway"]:
        s = C[(C.kind == kind) & C.measurable & (~C.circular)].sort_values("p_emp_rho")
        print(f"\n-- {kind} --")
        print(s.head(10)[cols].round(4).to_string())
    return C


if __name__ == "__main__":
    main()
