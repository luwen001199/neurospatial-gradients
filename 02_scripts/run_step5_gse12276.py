# -*- coding: utf-8 -*-
"""Step5 —— GSE12276 脑转移竞争风险验证 (Fine-Gray)。

流程
----
1. 解析 GEO series matrix (gzip 文本): 样本元数据 + 表达矩阵。
2. `inspect` 模式: 打印所有 characteristics 字段的取值分布, 人工确认再映射。
3. `run` 模式: 映射脑转移终点 -> Cox + Fine-Gray(cause1=脑转移, cause2=其他部位转移)
   + NPP / 神经轴打分。

用法
----
python 02_scripts/run_step5_gse12276.py inspect
python 02_scripts/run_step5_gse12276.py run
"""
from __future__ import annotations
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG          # noqa: E402
from nsgrad import bulk as nb            # noqa: E402
from nsgrad import surv as nsurv         # noqa: E402

SMX = CFG.DIR_RAW / "gse12276" / "GSE12276_series_matrix.txt.gz"
TAB = CFG.DIR_TAB
TAB.mkdir(parents=True, exist_ok=True)


def parse_series_matrix(path):
    meta, expr_rows, in_expr = {}, [], False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_expr = True
                header = None
                continue
            if line.startswith("!series_matrix_table_end"):
                in_expr = False
                continue
            if in_expr:
                if header is None:
                    header = line.split("\t")
                    header[0] = header[0].lstrip('"').rstrip('"')
                else:
                    parts = line.split("\t")
                    expr_rows.append([p.strip().strip('"') for p in parts])
                continue
            if line.startswith("!"):
                key = line.split("\t")[0]
                vals = [v.strip().strip('"') for v in line.split("\t")[1:]]
                # 每一行单独存放: 同名 tag(如 Sample_characteristics_ch1)可能出现多次,
                # 各代表一个独立的临床字段, 合并会错位
                meta.setdefault(key, []).append(vals)
    expr = pd.DataFrame(expr_rows).set_index(0)
    expr.index.name = "ID_REF"
    acc = meta.get("!Sample_geo_accession")
    if acc:
        expr.columns = acc[0][: expr.shape[1]]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return meta, expr


def inspect():
    meta, expr = parse_series_matrix(SMX)
    print(f"expression matrix: {expr.shape[0]} probes x {expr.shape[1]} samples")
    for k in ["!Series_title", "!Series_summary", "!Series_platform_id",
              "!Series_overall_design"]:
        if k in meta:
            print(f"\n{k}:")
            for v in meta[k]:
                print("   ", v[:300])
    ch = [k for k in meta if k.startswith("!Sample_characteristics") or
          k in ("!Sample_title", "!Sample_source_name")]
    for k in ch:
        vals = pd.Series(meta[k]).astype(str)
        print(f"\n{k} (n_uniq={vals.nunique()}):")
        print(vals.value_counts().head(12).to_string())


ANNOT = CFG.DIR_RAW / "gse12276" / "GPL570.annot.gz"


def load_probe_map():
    """GPL570 annot.gz -> DataFrame(probe, symbol)。

    GEO 注释文件格式: `!platform_table_begin` 之后是制表符表, 列 ID / Gene symbol;
    一个 probe 的 symbol 可能是 `A///B///C`, 取第一段作为主符号。
    """
    rows = []
    in_tab = False
    with gzip.open(ANNOT, "rt", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_tab = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_tab:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                i_id = header.index("ID")
                i_sym = header.index("Gene symbol")
                continue
            sym = parts[i_sym].split("///")[0].strip()
            if sym:
                rows.append((parts[i_id].strip(), sym))
    m = pd.DataFrame(rows, columns=["probe", "symbol"]).drop_duplicates("probe")
    return m


def clinical_table(meta):
    """从 characteristics 字段解析 (site, time)。

    本 series 的两个字段都在 !Sample_characteristics_ch1 下:
      field A: 'site of relapse (brain or other): other/local/brain/brain+other'
      field B: 'survival time (months): <数值>'   —— 实为**无病至复发时间**
    """
    fields = meta.get("!Sample_characteristics_ch1", [])
    site = pd.Series(fields[0], dtype=str).str.strip()
    time_f = None
    for f in fields[1:]:
        s = pd.Series(f, dtype=str).str.strip()
        if s.str.contains("survival time", case=False).any():
            time_f = s
    assert time_f is not None, "找不到 survival time 字段"
    site_v = site.str.split(":", n=1).str[1].str.strip().str.lower()
    time_v = pd.to_numeric(time_f.str.split(":", n=1).str[1].str.strip(),
                           errors="coerce")
    d = pd.DataFrame({"site": site_v, "time": time_v})
    d.index = meta["!Sample_geo_accession"][0][: len(d)]
    return d


def zscore(s):
    return (s - s.mean()) / s.std(ddof=0)


def run():
    meta, expr = parse_series_matrix(SMX)
    print(f"[Step5] probes {expr.shape}")
    annot = load_probe_map()
    expr = expr[expr.index.isin(annot["probe"])]
    expr.index = annot.set_index("probe").loc[expr.index, "symbol"].values
    expr = expr.groupby(level=0).mean()          # 多 probe 同基因取平均
    print(f"[Step5] gene-level matrix {expr.shape}")

    clin = clinical_table(meta)
    print(clin["site"].value_counts().to_string())

    # ---- NPP 与神经轴打分 (与 Step4b 完全一致的签名定义) ----
    sig = json.load(open(CFG.DIR_PROC / "npp_signature.json", encoding="utf-8"))
    from run_step4b_nerve_axes import AXES     # 复用同一套轴定义
    # gene-wise z (跨样本), 与 Step4b 打分方式一致
    Zg = (expr.sub(expr.mean(axis=1), axis=0)
              .div(expr.std(axis=1).replace(0, np.nan), axis=0))

    def cov(genes):
        return [g for g in genes if g in Zg.index]

    up, dn = cov(sig["up"]), cov(sig["down"])
    print(f"[Step5] NPP coverage: up {len(up)}/{len(sig['up'])}, "
          f"down {len(dn)}/{len(sig['down'])}")
    S = pd.DataFrame(index=Zg.columns)
    S["NPPscore"] = Zg.loc[up].mean(axis=0) - (Zg.loc[dn].mean(axis=0)
                                               if dn else 0)
    for k, genes in AXES.items():
        g = cov(genes)
        S[k] = Zg.loc[g].mean(axis=0) if g else np.nan
        print(f"    axis {k:18s} coverage {len(g)}/{len(genes)}")

    # ---- 终点定义: 脑转移为 cause 1, 其他远处转移为 cause 2, 局部复发删失 ----
    d = clin.join(S, how="inner")
    d["status"] = np.where(d["site"].isin(["brain", "brain+other"]), 1,
                           np.where(d["site"] == "other", 2, 0))
    d["tmo"] = d["time"]
    print(f"\n[Step5] analytic set n={len(d)}; "
          f"brain events={int((d.status==1).sum())}, "
          f"other-site competing={int((d.status==2).sum())}, "
          f"censored(local)={int((d.status==0).sum())}")

    rows = []
    for var in ["NPPscore"] + list(AXES.keys()):
        if var not in d.columns:
            continue
        zz = zscore(d[var]).values
        # (a) 全因远处复发 Cox (cause1+cause2 合并)
        e_any = (d["status"] > 0).astype(float).values
        r1 = nsurv.cox(d["tmo"], e_any, zz, names=[var])[0]
        # (b) Fine-Gray: cause1=brain, cause2=other distant
        r2 = nsurv.finegray(d["tmo"], d["status"], zz, names=[var])
        r2 = r2[0] if r2 else nsurv.CoxResult(var, np.nan, np.nan, np.nan,
                                              np.nan, len(d), int((d.status == 1).sum()))
        rows.append(dict(cohort="GSE12276", variable=var,
                         cause="any_distant_met", **{k: getattr(r1, k)
                                                     for k in
                                                     ["hr", "lo", "hi", "p", "n", "events"]}))
        rows.append(dict(cohort="GSE12276", variable=var,
                         cause="brain_first_finegray",
                         **{k: getattr(r2, k)
                            for k in ["hr", "lo", "hi", "p", "n", "events"]}))
        print(f"    {var:18s} any-met HR={r1.hr:.3f} ({r1.lo:.2f}-{r1.hi:.2f}) "
              f"p={r1.p:.3g} | brain FG sHR={r2.hr:.3f} "
              f"({r2.lo:.2f}-{r2.hi:.2f}) p={r2.p:.3g}")

    out = pd.DataFrame(rows)

    # ---- 敏感性: 严格"首部位=仅脑"定义 (brain+other 归入竞争事件) ----
    d["status_strict"] = np.where(d["site"] == "brain", 1,
                                  np.where(d["site"].isin(["other", "brain+other"]),
                                           2, 0))
    for var in ["NPPscore", "CD8_cytotox"]:
        zz = zscore(d[var]).values
        r = nsurv.finegray(d["tmo"], d["status_strict"], zz, names=[var])
        r = r[0] if r else nsurv.CoxResult(var, np.nan, np.nan, np.nan, np.nan,
                                           len(d), int((d.status_strict == 1).sum()))
        out = pd.concat([out, pd.DataFrame([dict(
            cohort="GSE12276", variable=var, cause="brain_only_strict_finegray",
            **{k: getattr(r, k) for k in ["hr", "lo", "hi", "p", "n", "events"]})])],
            ignore_index=True)
        print(f"    [strict] {var:16s} FG sHR={r.hr:.3f} ({r.lo:.2f}-{r.hi:.2f}) "
              f"p={r.p:.3g}")

    # ---- NPP 三分位的 Aalen-Johansen CIF (供作图) ----
    d["NPP_tertile"] = pd.qcut(d["NPPscore"], 3, labels=["T1_low", "T2", "T3_high"])
    cif_rows = []
    for lev, gg in d.groupby("NPP_tertile", observed=True):
        t, c1, c2 = nsurv.cif(gg["tmo"], gg["status"])
        for tt, cc in zip(t, c1):
            cif_rows.append(dict(tertile=str(lev), time=tt, cif_brain=cc))
    pd.DataFrame(cif_rows).to_csv(TAB / "step5_gse12276_cif_tertiles.csv",
                                  index=False)
    # log-rank 型检验 (脑事件 vs 竞争, 用 T1 vs T3 的 FG 模型三分类)
    ind = (d["NPP_tertile"] == "T3_high").astype(float).values
    r3 = nsurv.finegray(d["tmo"], d["status"], ind, names=["T3_vs_T1"])
    r3 = r3[0]
    out = pd.concat([out, pd.DataFrame([dict(
        cohort="GSE12276", variable="NPP_tertile(T3 vs T1)",
        cause="brain_first_finegray",
        **{k: getattr(r3, k) for k in ["hr", "lo", "hi", "p", "n", "events"]})])],
        ignore_index=True)
    print(f"    [tertile] T3 vs T1 FG sHR={r3.hr:.3f} ({r3.lo:.2f}-{r3.hi:.2f}) "
          f"p={r3.p:.3g}")

    out.to_csv(TAB / "step5_gse12276_brainmet.csv", index=False)
    d.to_csv(CFG.DIR_PROC / "gse12276_analysis_table.csv")
    print(f"\n  -> {TAB/'step5_gse12276_brainmet.csv'}")
    return out, d


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "inspect"
    if mode == "inspect":
        inspect()
    else:
        run()
