"""微环境读数的两道真实性门槛。

门槛 A —— 循环性 (circularity)
    神经 spot 是用 ALL_NERVE_GENES (UCHL1/TUBB3/PLP1/... ) 定义的,
    因此任何**成员里含神经标记基因**的读数, 都会因为定义本身而与
    "到最近神经的距离"负相关。这不是发现, 是同义反复。
    本脚本量化每个读数与神经定义基因的重叠, 并给出剔除后的结果。

门槛 B —— 可测性 (measurability)
    读数只在部分样本里可算(成员基因在该样本检出率不足)。
    n_samples 太少的读数(如 CALCA-RAMP1 仅 2 例)即便 Z 值大也不可信,
    按 MIN_SAMPLES 剔除。

输出
----
03_results/tables/niche_circularity.csv          每个读数的重叠统计
03_results/tables/niche_meta_summary.clean.csv    叠加两道门槛后的总表
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                      # noqa: E402
import niche_lib as NL                               # noqa: E402

TAB = CFG.DIR_TAB
MIN_SAMPLES = 10          # 读数至少要在 10/30 例里可算
FDR = 0.05
MIN_ABS_RHO = 0.03

# 神经定义基因 = 用于判定 nerve spot 的标记 + 细胞类型面板里的神经/施万细胞标记
NERVE_DEF = set(CFG.ALL_NERVE_GENES) | set(
    NL.CELLTYPE_MARKERS.get("Neuronal_Schwann", []))


def main():
    # ---- 读数成员 ----
    _, ct = NL.build_celltype_sets(None)
    lr = {f"LR|{lig}-{rec}|{axis}": [lig, rec] for lig, rec, axis, _ in NL.LR_PAIRS}
    pw = NL.load_msigdb(collections=("HALLMARK", "REACTOME", "KEGG"))
    blocks = {"celltype": ct, "lr": lr, "pathway": pw}

    rows = []
    for kind, d in blocks.items():
        for nm, gl in d.items():
            gs = set(gl)
            hit = sorted(gs & NERVE_DEF)
            rows.append(dict(name=nm, kind=kind, size=len(gs),
                             n_nerve_def=len(hit),
                             frac_nerve_def=len(hit) / max(1, len(gs)),
                             nerve_genes=";".join(hit)))
    C = pd.DataFrame(rows).set_index("name")

    G = pd.read_csv(TAB / "niche_meta_summary.csv", index_col=0)
    G = G.join(C[["size", "n_nerve_def", "frac_nerve_def", "nerve_genes"]],
               how="left")

    G["circular"] = G["n_nerve_def"].fillna(0) > 0
    G["measurable"] = G["n_samples"] >= MIN_SAMPLES
    G["pass_strict"] = (~G["circular"]) & G["measurable"] & \
        (G["fdr_msr"] < FDR) & (G["mean_rho"].abs() >= MIN_ABS_RHO)

    C.to_csv(TAB / "niche_circularity.csv")
    G.to_csv(TAB / "niche_meta_summary.clean.csv")

    print(f"读数 {len(G)}: celltype {int((G.kind=='celltype').sum())}, "
          f"lr {int((G.kind=='lr').sum())}, pathway {int((G.kind=='pathway').sum())}")
    print(f"\n含神经定义基因(循环性可疑): {int(G.circular.sum())} / {len(G)}")
    print(f"可测样本数 < {MIN_SAMPLES}: {int((~G.measurable).sum())}")

    for kind in ["celltype", "lr", "pathway"]:
        s = G[G.kind == kind]
        base = ((s.fdr_msr < FDR) & (s.mean_rho.abs() >= MIN_ABS_RHO)).sum()
        print(f"\n=== {kind}: 原始通过 {base}/{len(s)} -> 严格通过 "
              f"{int(s.pass_strict.sum())}/{len(s)}")
        sig = s[s.pass_strict].sort_values("stouffer_Z")
        cols = ["n_samples", "size", "mean_rho", "consistency", "stouffer_Z",
                "fdr_msr", "n_nerve_def"]
        print(sig.head(15)[cols].round(4).to_string())
        print("  反向(正向富集, 即远离神经升高) top5:")
        print(sig.tail(5)[cols].round(4).to_string())
    return G


if __name__ == "__main__":
    main()
