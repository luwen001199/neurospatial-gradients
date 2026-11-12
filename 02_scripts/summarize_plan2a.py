"""汇总方案 2-A 全量结果的关键数字 (供 Results 3.6 / Note 11 回填)."""
from pathlib import Path

import pandas as pd

ROOT = Path("D:/嗜神经浸润")
PROC = ROOT / "01_data/processed"


def p(name):
    f = PROC / name
    return pd.read_csv(f) if f.exists() else None


def star(fdr):
    return "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else ""


print("=" * 72)
print("[1] 细胞来源归属 (GSE176078 全量参考)")
mod = p("scrna_module_origin.csv")
if mod is not None:
    for m in mod.module.unique():
        sub = mod[mod.module == m]
        print(f"\n  {m} (n={len(sub)} genes):")
        print("    top_frac_type 分布:", sub.top_frac_type.value_counts().head(8).to_dict())
        hi = sub.nlargest(5, "top_frac")[["gene", "top_frac_type", "top_frac"]]
        print("    top5:", ", ".join(f"{r.gene}→{r.top_frac_type}({r.top_frac:.2f})" for r in hi.itertuples()))
    # 代表基因
    for g in ["COL1A1", "SPARC", "POSTN", "LUM", "DCN", "ACTA2", "EPCAM", "TUBB3", "NEFL"]:
        row = mod[mod.gene == g]
        if len(row):
            r = row.iloc[0]
            print(f"    {g:8s} top={r.top_frac_type} frac={r.top_frac:.2f} "
                  f"(myCAF={r.frac_myCAF:.2f} CancerEpi={r.frac_CancerEpithelial if pd.notna(r.frac_CancerEpithelial) else float('nan'):.2f})")
sm = p("scrna_samples_meta.csv")
if sm is not None:
    print("\n  参考集: n_samples =", len(sm), "| total cells =", sm.n_cells.sum() if "n_cells" in sm else "?")

print("\n" + "=" * 72)
print("[2] 解卷积 QC (30 例 Visium)")
for sfx, tag in [("", "primary"), ("_nomodule", "module-excluded")]:
    qc = p(f"deconv_qc{sfx}.csv")
    if qc is not None:
        print(f"  {tag:16s} n={len(qc)}  median_rel_err: median={qc.median_rel_err.median():.3f} "
              f"IQR={qc.median_rel_err.quantile(.25):.3f}-{qc.median_rel_err.quantile(.75):.3f}  "
              f"panel_genes={int(qc.panel_genes_used.iloc[0])}")

print("\n" + "=" * 72)
print("[3] 共聚集: 细胞类型梯度 (偏 Spearman ~ 距神经距离, Fisher-z meta)")
for sfx, tag in [("", "primary"), ("_nomodule", "module-excluded")]:
    m = p(f"deconv_gradient_meta{sfx}.csv")
    if m is not None:
        print(f"\n  {tag}:")
        for r in m.itertuples():
            print(f"    {r.celltype:16s} rho={r.rho_mean:+.3f}  k={r.k}  p={r.p:.1e}  fdr={r.fdr_bh:.1e} {star(r.fdr_bh)}")

print("\n" + "=" * 72)
print("[4] 神经 spot vs 非神经 (每样本 mean log2FC, Stouffer meta)")
for sfx, tag in [("", "primary"), ("_nomodule", "module-excluded")]:
    nv = p(f"deconv_nerve_vs_non_meta{sfx}.csv")
    if nv is not None:
        print(f"\n  {tag}:")
        for r in nv.itertuples():
            print(f"    {r.celltype:16s} log2fc={r.mean_log2fc:+.3f}  k={r.k}  p={r.p_stouffer:.1e}  fdr={r.fdr_bh:.1e} {star(r.fdr_bh)}")

print("\n" + "=" * 72)
print("[5] 模块得分 ~ 细胞丰度 (每样本偏 Spearman, meta)")
for sfx, tag in [("", "primary"), ("_nomodule", "module-excluded")]:
    mo = p(f"deconv_module_coaggr_meta{sfx}.csv")
    if mo is not None:
        print(f"\n  {tag}:")
        for r in mo.itertuples():
            print(f"    {r.celltype:16s} neural={r.rho_neural:+.3f}(p={r.p_neural:.1e},fdr={r.fdr_neural:.1e}{star(r.fdr_neural)})  "
                  f"stroma={r.rho_stroma:+.3f}(p={r.p_stroma:.1e},fdr={r.fdr_stroma:.1e}{star(r.fdr_stroma)})  k={r.k}")

print("\n" + "=" * 72)
print("[6] GSE75688 交叉验证 (如已运行)")
cv = p("scrna_crossval_summary.csv")
if cv is not None:
    for m in cv.module.unique():
        sub = cv[cv.module == m]
        print(f"\n  {m} (n={len(sub)}): {sub.top_type.value_counts().head(6).to_dict()}")
else:
    print("  (crossval 未运行或未完成)")
print("\nDONE")
