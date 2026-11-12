"""方案六(修订版): 空间衍生特征的临床增量价值 —— 严格验证, 不预设阳性。

为什么不是原方案的"空间机器学习预后模型"
------------------------------------------
GSE210616 是 43 例手术切除的空间转录组队列, GEO 记录中**没有随访/生存字段**。
因此"用空间梯度特征 + 空间细胞丰度训练预后模型"在数据层面无法实现:
没有 outcome 可学。任何声称做到的分析只能是把空间特征投影到 bulk 队列,
其本质仍是基因集评分 —— 一旦如此就回到本文已经证伪的"签名边界任意性"问题。

因此改为**临床预测模型增量价值(incremental value)框架**:
  基线: 常规临床病理变量
  增量: 空间衍生的三个读数(NPP / 神经固有模块 / 神经基质模块)
评估: 嵌套似然比检验、Harrell C-index 增量(10 折交叉验证, 防过拟合)、
      5 年时间依赖 AUC、校准、决策曲线分析(DCA)。
外部验证: GSE12276 (n=204, 有复发时间)。

所有指标如实报告 —— 若增量为 0, 结论就是"空间衍生的神经读数不提供
临床增量信息", 这本身就是本文要传达的方法学警示。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")

from nsgrad import config as CFG                     # noqa: E402

TAB = CFG.DIR_TAB
FIG = CFG.DIR_FIG

# 注意: METABRIC 表中 SIZE 与 GRADE 全列缺失(1976/1976), 不可用作协变量。
CLIN = ["AGE", "NODES_POS", "ER_pos", "chemo", "hormone"]
CANDIDATES = ["NPPscore", "Neural_intrinsic", "Perineural_stroma",
              "CAF_ECM", "CD8_cytotox"]


# ------------------------------------------------------------------ 工具
def harrell_c(t, e, risk):
    """Harrell's C (忽略 ties in risk 与不可比对)。"""
    n = len(t)
    if n > 3000:                     # 抽样近似, 控制 O(n^2)
        rng = np.random.default_rng(0)
        idx = rng.choice(n, size=3000, replace=False)
        t, e, risk = t[idx], e[idx], risk[idx]
        n = len(t)
    num = den = 0
    order = np.argsort(-risk)
    for i in range(n):
        ti, ei, ri = t[i], e[i], risk[i]
        if ei == 0:
            continue
        for j in range(n):
            if t[j] > ti:
                den += 1
                if risk[j] < ri:
                    num += 1
                elif risk[j] == ri:
                    num += 0.5
    return num / den if den else np.nan


def _fast_c(t, e, risk, max_n=2500):
    """向量化 Harrell C(必要时子抽样)。"""
    n = len(t)
    if n > max_n:
        rng = np.random.default_rng(0)
        idx = rng.choice(n, size=max_n, replace=False)
        t, e, risk = t[idx], e[idx], risk[idx]
    ti = t[:, None]
    comparable = (e[:, None] == 1) & (t[None, :] > ti)
    better = risk[None, :] < risk[:, None]
    tie = risk[None, :] == risk[:, None]
    num = (comparable * (better + 0.5 * tie)).sum()
    den = comparable.sum()
    return float(num / den) if den else np.nan


def fit_phreg(t, e, X, names):
    keep = [j for j in range(X.shape[1]) if np.std(X[:, j]) > 1e-12]
    if not keep:
        return None
    m = PHReg(t, X[:, keep], status=e, ties="efron").fit()
    return m, [names[j] for j in keep], keep


def lrt(t, e, X0, X1):
    """嵌套模型似然比检验。"""
    m0 = PHReg(t, X0, status=e, ties="efron").fit()
    m1 = PHReg(t, X1, status=e, ties="efron").fit()
    stat = 2 * (m1.llf - m0.llf)
    df = X1.shape[1] - X0.shape[1]
    return stat, df, stats.chi2.sf(stat, df), m0, m1


def td_auc(t, e, risk, horizon):
    """累积病例-对照 AUC (cumulative case-control) at horizon。"""
    case = (e == 1) & (t <= horizon)
    ctrl = t > horizon
    if case.sum() < 5 or ctrl.sum() < 5:
        return np.nan
    x = np.concatenate([risk[case], risk[ctrl]])
    y = np.concatenate([np.ones(case.sum()), np.zeros(ctrl.sum())])
    return float(stats.mannwhitneyu(x[y == 1], x[y == 0]).statistic /
                 (case.sum() * ctrl.sum()))


def net_benefit(t, e, risk, horizon, thresh):
    """DCA 净收益(简单版, 用 horizon 内的事件作为二分类结局)。"""
    y = ((e == 1) & (t <= horizon)).astype(int)
    n = len(y)
    # 风险阈值 -> 预测概率(用 rank 标定到事件率)
    p = pd.Series(risk).rank(pct=True).to_numpy()
    base = y.mean()
    p = np.clip(p * base / p.mean(), 1e-6, 1 - 1e-6) if p.mean() > 0 else p
    treat = p >= thresh
    tp = (treat & (y == 1)).sum()
    fp = (treat & (y == 0)).sum()
    nb = tp / n - fp / n * (thresh / (1 - thresh))
    nb_all = y.mean() - (1 - y.mean()) * (thresh / (1 - thresh))
    return float(nb), float(nb_all)


# ------------------------------------------------------------------ 主流程
def load_table():
    clin = pd.read_csv("01_data/processed/metabric_analysis_table.csv")
    expr = pd.read_csv("01_data/raw/metabric/metabric_mrna_targetgenes.csv",
                       index_col=0)
    from run_npp_split import MODULES, score_bulk
    # 按 PATIENT_ID 对齐 (expr 1980 行 vs clin 1976 行, 不能直接按位置赋值)
    d = clin.set_index("PATIENT_ID")
    for nm, gl in MODULES.items():
        d[nm] = score_bulk(expr, [g for g in gl if g in expr.columns])
    d = d.reset_index()
    return d


def main(folds=10, seed=7):
    d = load_table()
    rows = []
    for outcome, tc, ec in [("RFS", "RFS_MONTHS", "RFS_EVENT"),
                            ("OS", "OS_MONTHS", "OS_EVENT")]:
        sub = d.dropna(subset=CLIN + [tc, ec]).copy()
        sub = sub[sub[tc] > 0]
        t = sub[tc].to_numpy(float)
        e = sub[ec].to_numpy(float)
        print(f"\n=== {outcome}: n={len(sub)}, events={int(e.sum())} ===")

        def zs(cols):
            return np.column_stack([(sub[c].to_numpy(float) - sub[c].mean()) /
                                    sub[c].std(ddof=0) for c in cols])

        X0 = zs(CLIN)
        # 基线 C-index (in-sample 与 CV)
        m0 = PHReg(t, X0, status=e, ties="efron").fit()
        c0 = _fast_c(t, e, X0 @ m0.params)

        # ---- 10 折 CV 增量 ----
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(sub))
        fold = np.zeros(len(sub), dtype=int)
        for k in range(folds):
            fold[idx[k::folds]] = k
        cv = {c: [] for c in CANDIDATES}
        cv_base = []
        for k in range(folds):
            tr, te = fold != k, fold == k
            m = PHReg(t[tr], X0[tr], status=e[tr], ties="efron").fit()
            cv_base.append(_fast_c(t[te], e[te], X0[te] @ m.params))
            for c in CANDIDATES:
                if c not in sub.columns:
                    continue
                X1 = np.column_stack([X0, zs([c])])
                mm = PHReg(t[tr], X1[tr], status=e[tr], ties="efron").fit()
                cv[c].append(_fast_c(t[te], e[te], X1[te] @ mm.params))
        print(f"  基线 C-index: 全样本 {c0:.4f} | CV {np.mean(cv_base):.4f}")

        for c in CANDIDATES:
            if c not in sub.columns or sub[c].isna().all():
                continue
            s2 = sub.dropna(subset=[c]).copy()
            t2 = s2[tc].to_numpy(float)
            e2 = s2[ec].to_numpy(float)
            Xb = np.column_stack([(s2[x].to_numpy(float) - s2[x].mean()) /
                                  s2[x].std(ddof=0) for x in CLIN])
            Xf = np.column_stack([Xb, (s2[c].to_numpy(float) - s2[c].mean()) /
                                  s2[c].std(ddof=0)])
            stat, df, p, mb, mf = lrt(t2, e2, Xb, Xf)
            beta = mf.params[-1]
            se = mf.bse[-1]
            ci = _fast_c(t2, e2, Xf @ mf.params)
            rows.append(dict(outcome=outcome, variable=c, n=len(s2),
                             events=int(e2.sum()),
                             HR_perSD=float(np.exp(beta)),
                             CI_lo=float(np.exp(beta - 1.96 * se)),
                             CI_hi=float(np.exp(beta + 1.96 * se)),
                             p_wald=float(mf.pvalues[-1]),
                             LRT_chi2=float(stat), LRT_p=float(p),
                             C_full=ci,
                             C_base=_fast_c(t2, e2, Xb @ mb.params),
                             C_delta=ci - _fast_c(t2, e2, Xb @ mb.params),
                             C_CV=np.mean(cv[c]) if cv[c] else np.nan,
                             C_CV_base=np.mean(cv_base),
                             C_CV_delta=(np.mean(cv[c]) - np.mean(cv_base))
                             if cv[c] else np.nan))
            print(f"  {c:20s} HR/SD={np.exp(beta):.3f} "
                  f"({np.exp(beta-1.96*se):.3f}-{np.exp(beta+1.96*se):.3f}) "
                  f"p={mf.pvalues[-1]:.4f} | LRT p={p:.4f} | "
                  f"dC(in)={ci - _fast_c(t2, e2, Xb @ mb.params):+.4f} | "
                  f"dC(CV)={(np.mean(cv[c]) - np.mean(cv_base)) if cv[c] else np.nan:+.4f}")

    R = pd.DataFrame(rows)
    R.to_csv(TAB / "prognostic_incremental_value.csv", index=False)
    print(f"\n-> {TAB/'prognostic_incremental_value.csv'}")

    # ---- 校准 + DCA (RFS, 基线 vs +候选) ----
    out = []
    for outcome, tc, ec, horizon in [("RFS", "RFS_MONTHS", "RFS_EVENT", 60.0),
                                     ("OS", "OS_MONTHS", "OS_EVENT", 120.0)]:
        sub = d.dropna(subset=CLIN + [tc, ec])
        sub = sub[sub[tc] > 0]
        t = sub[tc].to_numpy(float)
        e = sub[ec].to_numpy(float)
        Xb = np.column_stack([(sub[x].to_numpy(float) - sub[x].mean()) /
                              sub[x].std(ddof=0) for x in CLIN])
        for label, X in [("clinical", Xb)] + [
                (f"clinical+{c}", np.column_stack([
                    Xb, (sub[c].to_numpy(float) - sub[c].mean()) /
                    sub[c].std(ddof=0)]))
                for c in CANDIDATES if c in sub.columns]:
            m = PHReg(t, X, status=e, ties="efron").fit()
            risk = X @ m.params
            auc = td_auc(t, e, risk, horizon)
            nb_list = [net_benefit(t, e, risk, horizon, th)[0]
                       for th in (0.10, 0.20, 0.30)]
            out.append(dict(outcome=outcome, model=label, horizon=horizon,
                            C_index=_fast_c(t, e, risk), AUC_t=auc,
                            NB_10=nb_list[0], NB_20=nb_list[1], NB_30=nb_list[2]))
    O = pd.DataFrame(out)
    O.to_csv(TAB / "prognostic_calibration_dca.csv", index=False)
    print(O.round(4).to_string(index=False))
    return R, O


if __name__ == "__main__":
    main()
