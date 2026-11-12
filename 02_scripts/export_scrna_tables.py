# -*- coding: utf-8 -*-
"""方案 2-A 收尾: 把 scRNA 锚定/解卷积结果导出为补充表 CSV（全量数据跑完后执行）。"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "01_data/processed"
SUPP = ROOT / "03_results/tables/supplementary"
SUPP.mkdir(parents=True, exist_ok=True)

# ---- S18: 模块基因细胞来源 ----
mo = pd.read_csv(PROC / "scrna_module_origin.csv")
cols = ["gene", "module", "status", "top_frac_type", "top_frac",
        "frac_myCAF", "frac_iCAF", "frac_CancerEpithelial", "frac_Tcells",
        "frac_Myeloid", "frac_Endothelial", "frac_PVL", "mean_cp10k_top",
        "enrichment_top", "total_frac_expr"]
mo[[c for c in cols if c in mo.columns]].to_csv(SUPP / "TableS18_module_cell_origin.csv", index=False)

# ---- S19: 解卷积汇总 ----
frames = []
gm = pd.read_csv(PROC / "deconv_gradient_meta.csv"); gm["analysis"] = "abundance~distance (partial Spearman)"; gm["panel"] = "primary"
frames.append(gm.assign(metric="rho_mean")[["analysis", "panel", "metric", "celltype", "k", "rho_mean", "p", "fdr_bh"]])
gm2 = pd.read_csv(PROC / "deconv_gradient_meta_nomodule.csv")
frames.append(gm2.assign(analysis="abundance~distance (partial Spearman)", panel="module-excluded",
                         metric="rho_mean")[["analysis", "panel", "metric", "celltype", "k", "rho_mean", "p", "fdr_bh"]])
nv = pd.read_csv(PROC / "deconv_nerve_vs_non_meta.csv"); nv["analysis"] = "nerve vs non-nerve spots"; nv["panel"] = "primary"
frames.append(nv.rename(columns={"mean_log2fc": "rho_mean", "p_stouffer": "p"})
              .assign(metric="log2FC mean")[["analysis", "panel", "metric", "celltype", "k", "rho_mean", "p", "fdr_bh"]])
cm = pd.read_csv(PROC / "deconv_module_coaggr_meta.csv")
cm1 = cm.assign(analysis="module score ~ abundance", panel="primary", metric="rho vs neural-intrinsic")
frames.append(cm1[["analysis", "panel", "metric", "celltype", "k", "rho_neural", "p_neural", "fdr_neural"]]
              .rename(columns={"rho_neural": "rho_mean", "p_neural": "p", "fdr_neural": "fdr_bh"}))
cm2 = cm.assign(analysis="module score ~ abundance", panel="primary", metric="rho vs perineural-stroma")
frames.append(cm2[["analysis", "panel", "metric", "celltype", "k", "rho_stroma", "p_stroma", "fdr_stroma"]]
              .rename(columns={"rho_stroma": "rho_mean", "p_stroma": "p", "fdr_stroma": "fdr_bh"}))
out = pd.concat(frames, ignore_index=True)
out.to_csv(SUPP / "TableS19_deconv_summary.csv", index=False)
print("wrote TableS18/S19:", len(mo), "genes;", len(out), "rows")
