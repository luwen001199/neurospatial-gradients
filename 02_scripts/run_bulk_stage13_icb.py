# -*- coding: utf-8 -*-
"""Stage 13 —— GSE91061 ICB 探索性分析（跨适应症，不计入主结论）。

队列：Riaz et al., Cell 2017 (PMID 29033130)，BMS038 试验
  nivolumab ± ipilimumab 治疗的晚期黑色素瘤，65 患者 109 RNA-seq 样本
  （51 pre-treatment / 58 on-treatment）。GEO 无治疗臂、无 PFS/OS ——
  终点只能做 RECIST 响应类别（PD / SD / PRCR）。

边界声明（与稿件 Note 13 一致）：
  * 跨适应症（黑色素瘤 vs 乳腺）、样本量小（pre-treatment 可评估 n≈49）
  * 假设生成性质，不计入主结论，不参与多重校正主分析
  * FPKM 平台与 bulk RNA-seq 队列不同源，模块评分在 Pre 子集内 z 化

输出（03_results/tables/bulk_extrapolation/）：
  stage13_icb_logistic.csv      模块 → ORR / benefit 两个终点的 logistic + Fisher
  stage13_random_baseline.csv   510 组 size-matched 随机集基线（ORR 终点）
  stage13_paired_delta.csv      pre→on 配对变化 × 响应
"""
from __future__ import annotations

import gzip
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from bulk_lib import ROOT, RAW, TAB, load_modules, score_mean, score_diff, fdr_bh  # noqa: E402

RAW91061 = RAW / "gse91061"
SEED = 20251109
N_RANDOM = 510

# ---------------------------------------------------------------- symbol→Entrez
def fetch_symbol_entrez(symbols: list[str]) -> dict[str, str]:
    """mygene.info 逐符号查询 symbol→entrez_id（GET 稳定，POST 在本网络 404）。
    结果缓存到本地 csv 保证可复现。返回 {symbol: entrez_id}。"""
    cache = RAW91061 / "symbol2entrez.csv"
    m: dict[str, str] = {}
    if cache.exists():
        c = pd.read_csv(cache, dtype=str)
        m = dict(zip(c["symbol"], c["entrez"]))
    todo = [s for s in symbols if s not in m]
    if todo:
        def one(s: str):
            import urllib.request, ssl, time, urllib.parse
            q = urllib.parse.quote(s)
            url = (f"https://mygene.info/v3/query?q={q}"
                   f"&fields=entrezgene&species=human&size=1")
            for _ in range(4):
                try:
                    ctx = ssl.create_default_context()
                    ctx.check_hostname = False
                    ctx.verify_mode = ssl.CERT_NONE
                    with urllib.request.urlopen(url, timeout=30, context=ctx) as r:
                        d = json.loads(r.read().decode())
                        hits = d.get("hits", [])
                        if hits:
                            eg = hits[0].get("entrezgene")
                            if eg:
                                return s, str(eg)
                        return s, ""
                except Exception:
                    time.sleep(2)
            return s, ""
        with ThreadPoolExecutor(8) as ex:
            for s, e in ex.map(one, todo):
                m[s] = e
        pd.DataFrame({"symbol": list(m), "entrez": [m[k] for k in m]}
                     ).to_csv(cache, index=False)
    return m


# ---------------------------------------------------------------- 数据装载
def load_gse91061() -> tuple[pd.DataFrame, pd.DataFrame]:
    """FPKM 矩阵 (sample × symbol) + 临床表 (index = sample 列名)。"""
    fp = RAW91061 / "GSE91061_fpkm.csv.gz"
    df = pd.read_csv(fp, index_col=0)
    df.index = df.index.astype(str).str.strip()
    # 行 = Entrez ID（GEO 声明 NCBI Entrez Gene IDs）
    mods = load_modules()
    krt = ["KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17", "KRT15"]
    need = sorted({g for v in mods.values() for g in v} | set(krt))
    idx_is_symbol = df.index.str.match(r"^[A-Za-z][A-Za-z0-9@._-]*$").mean() > 0.9
    if not idx_is_symbol:
        sym2ent = fetch_symbol_entrez(need)
        mapped = {}
        for sym, ent in sym2ent.items():
            if ent and ent in df.index:
                mapped[ent] = sym
        df = df.rename(index=mapped)
        print(f"[map] symbol→Entrez 命中 {len(mapped)}/{len(need)} 个模块基因")
    # 样本列名 = 完整 title（Pt#_Pre/On_AD编号），与 clinical.tsv 的 title 一致
    clin = pd.read_csv(RAW91061 / "clinical.tsv", sep="\t")
    title2row = {}
    for _, r in clin.iterrows():
        title2row[r["title"].strip()] = (r["visit"], r["response"], r["title"])
    cols = {}
    for c in df.columns:
        c0 = str(c).strip().strip('"')
        if c0 in title2row:
            cols[c0] = title2row[c0]
    if len(cols) < 50:
        print("[warn] title 列名匹配率低:", list(df.columns[:5]))
    df = df[list(cols)]
    clin_df = pd.DataFrame(
        {"visit": [cols[c][0] for c in df.columns],
         "response": [cols[c][1] for c in df.columns],
         "title": [cols[c][2] for c in df.columns]},
        index=df.columns)
    df = np.log2(df.astype(float) + 1.0).T   # sample × gene
    print(f"[load] {df.shape[0]} 样本 × {df.shape[1]} 基因; "
          f"Pre={int((clin_df.visit=='Pre').sum())}")
    return df, clin_df


# ---------------------------------------------------------------- logistic
def logistic_or(y: np.ndarray, x: np.ndarray) -> tuple[float, float, float, float]:
    """返回 beta, se, OR, p。样本量小, 用带 ridge epsilon 的 IRLS 防完全分离。"""
    X = np.column_stack([np.ones(len(x)), x.astype(float)])
    y = y.astype(float)
    b = np.zeros(2)
    for _ in range(60):
        eta = np.clip(X @ b, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-6, None) + 1e-9
        z = X @ b + (y - p) / W
        XtW = X.T * W
        A = XtW @ X + np.eye(2) * 1e-8
        b_new = np.linalg.solve(A, XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-10:
            b = b_new
            break
        b = b_new
    eta = np.clip(X @ b, -30, 30)
    p = 1.0 / (1.0 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-6, None) + 1e-9
    A = (X.T * W) @ X + np.eye(2) * 1e-8
    cov = np.linalg.inv(A)
    se = np.sqrt(np.diag(cov))
    from scipy.stats import norm
    zstat = b[1] / se[1]
    pval = 2 * norm.sf(abs(zstat))
    return b[1], se[1], float(np.exp(b[1])), float(pval)


def auc_score(y: np.ndarray, x: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, x))


# ---------------------------------------------------------------- 主流程
def main() -> None:
    df, clin = load_gse91061()
    mods = load_modules()
    pre = clin.index[clin["visit"] == "Pre"]
    expr_pre = df.loc[pre]
    resp = clin.loc[pre, "response"]

    # ---- 终点
    y_orr = (resp == "PRCR").astype(int)          # CR/PR vs SD+PD
    ev_orr = resp.isin(["PRCR", "SD", "PD"])
    y_ben = resp.replace({"PRCR": 1, "SD": 1, "PD": 0})  # 疾病控制 vs PD
    print(f"[endpoint] ORR: {y_orr[ev_orr].sum()}/{ev_orr.sum()} 响应; "
          f"benefit: {y_ben[ev_orr & y_ben.notna()].sum()}")

    # ---- 模块评分（Pre 子集内 z 化）
    scores = {}
    scores["NPP58"] = score_diff(expr_pre, mods["NPP_up"], mods["NPP_down"])
    for k in ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM", "LRRC15_myCAF",
              "CD8_cytotox", "TGFb_core", "TLS", "IFN_gamma", "M2_macro",
              "Treg", "SPP1_macro"]:
        scores[k] = score_mean(expr_pre, mods[k])
    S = pd.DataFrame(scores)

    # ---- TNBC-like（basal 角蛋白）可行性核查
    krt = ["KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17", "KRT15"]
    have = [g for g in krt if g in expr_pre.columns]
    krt_score = score_mean(expr_pre, have, min_genes=2)
    print(f"[tnbc-like] basal-keratin 基因检出 {len(have)}/{len(krt)}; "
          f"评分 IQR 宽度={np.nanpercentile(krt_score,75)-np.nanpercentile(krt_score,25):.2f}")

    # ---- 主分析: logistic per SD
    rows = []
    for k in S.columns:
        x = S[k][ev_orr]
        yy = y_orr[ev_orr]
        beta, se, orr, p = logistic_or(yy.to_numpy(), x.to_numpy())
        auc = auc_score(yy.to_numpy(), x.to_numpy())
        med = x > x.median()
        tbl = [[int(((yy == 1) & med).sum()), int(((yy == 1) & ~med).sum())],
               [int(((yy == 0) & med).sum()), int(((yy == 0) & ~med).sum())]]
        fish_p = stats.fisher_exact(tbl)[1]
        # benefit 终点
        evb = ev_orr & y_ben.notna() & (resp != "UNK")
        yb = y_ben[evb].astype(int)
        xb = S[k][evb]
        _, _, orr_b, p_b = logistic_or(yb.to_numpy(), xb.to_numpy())
        auc_b = auc_score(yb.to_numpy(), xb.to_numpy())
        rows.append(dict(module=k, n=int(ev_orr.sum()), events=int(yy.sum()),
                         OR_orr=round(orr, 3), OR_lo=round(np.exp(beta - 1.96 * se), 3),
                         OR_hi=round(np.exp(beta + 1.96 * se), 3),
                         AUC_orr=round(auc, 3), P_orr=round(p, 4),
                         P_fisher=round(fish_p, 4),
                         OR_benefit=round(orr_b, 3), AUC_benefit=round(auc_b, 3),
                         P_benefit=round(p_b, 4)))
    out = pd.DataFrame(rows)
    out["q_orr"] = fdr_bh(out["P_orr"].to_numpy())
    out.to_csv(TAB / "stage13_icb_logistic.csv", index=False)
    print(out.to_string(index=False))

    # ---- 随机集基线（ORR 终点, size-matched = stroma 模块大小）
    rng = np.random.default_rng(SEED)
    yy = y_orr[ev_orr].to_numpy()
    size = len(mods["Perineural_stroma"])
    det = [g for g in expr_pre.columns
           if (expr_pre[g] > 0).sum() >= 0.5 * len(expr_pre)]
    print(f"[null] 检出基因池 {len(det)}, 随机集大小 {size} × {N_RANDOM}")
    rec = []
    for i in range(N_RANDOM):
        gs = rng.choice(det, size=size, replace=False)
        x = score_mean(expr_pre, list(gs)).to_numpy()
        x = x[ev_orr]
        if np.std(x) == 0:
            continue
        _, _, orr, p = logistic_or(yy, x)
        rec.append((orr, auc_score(yy, x), p))
    rb = pd.DataFrame(rec, columns=["OR", "AUC", "P"])
    rb.to_csv(TAB / f"stage13_random_baseline.csv", index=False)

    for k in ["Perineural_stroma", "Neural_intrinsic", "NPP58", "CAF_ECM",
              "LRRC15_myCAF", "CD8_cytotox"]:
        obs_or = float(out.loc[out.module == k, "OR_orr"].iloc[0])
        obs_auc = float(out.loc[out.module == k, "AUC_orr"].iloc[0])
        # 双侧经验 P: 随机 OR 距 1 的偏差 >= 观测距 1 的偏差
        dev = abs(np.log(rb["OR"].to_numpy())) >= abs(np.log(obs_or))
        emp = min(1.0, 2 * min(dev.mean(), 1 - dev.mean() + 1 / len(rb)))
        auc_emp = min(1.0, 2 * min((rb["AUC"] >= obs_auc).mean(),
                                   (rb["AUC"] <= obs_auc).mean()))
        print(f"[null] {k}: 观测 OR={obs_or:.2f} AUC={obs_auc:.2f} | "
              f"随机 OR 中位={rb['OR'].median():.2f} | "
              f"经验P(OR)={emp:.3f} P(AUC)={auc_emp:.3f}")

    # ---- 配对 pre→on 变化 × 响应
    pair_rows = []
    by_pt: dict[str, list[str]] = {}
    for s in clin.index:
        pt = clin.loc[s, "title"].split("_")[0]
        by_pt.setdefault(pt, []).append(s)
    for pt, ss in by_pt.items():
        v = {clin.loc[s, "visit"]: s for s in ss}
        if "Pre" not in v or "On" not in v:
            continue
        r = clin.loc[v["Pre"], "response"]
        if r not in ("PRCR", "SD", "PD"):
            continue
        row = {"patient": pt, "response": r}
        for k in ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM",
                  "CD8_cytotox", "TGFb_core"]:
            # 全样本集合内评分保持单调（z 化用 Pre 的均值方差会破坏可比性,
            # 这里直接用 log-FPKM 的 rank-invariant mean-z 在 Pre+On 联合空间算分）
            pass
        pair_rows.append({"patient": pt, "response": r,
                          "gsm_pre": v["Pre"], "gsm_on": v["On"]})
    # 联合空间评分
    sel = ["Perineural_stroma", "Neural_intrinsic", "CAF_ECM", "CD8_cytotox",
           "TGFb_core", "NPP58"]
    Sj = pd.DataFrame({k: (score_diff(df, mods["NPP_up"], mods["NPP_down"])
                           if k == "NPP58" else score_mean(df, mods[k]))
                       for k in sel})
    pr = []
    for p in pair_rows:
        row = {"patient": p["patient"], "response": p["response"]}
        for k in sel:
            row[f"d_{k}"] = float(Sj.loc[p["gsm_on"], k] - Sj.loc[p["gsm_pre"], k])
        pr.append(row)
    dd = pd.DataFrame(pr)
    dd.to_csv(TAB / "stage13_paired_delta.csv", index=False)
    print(f"[paired] {len(dd)} 例配对 (Pre+On)")
    for k in sel:
        r = dd[dd.response == "PRCR"][f"d_{k}"]
        nr = dd[dd.response == "PD"][f"d_{k}"]
        u = stats.mannwhitneyu(r, nr, alternative="two-sided")
        print(f"  d_{k}: PRCR 中位 {np.median(r):+.2f} vs PD {np.median(nr):+.2f} "
              f"(P={u.pvalue:.3f})")

    print("\n[done] 探索性声明: 跨适应症(黑色素瘤)、n 小、GEO 无臂别/生存数据; "
          "结果仅假设生成, 不计入稿件主结论。")


if __name__ == "__main__":
    main()
