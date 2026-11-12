# -*- coding: utf-8 -*-
"""生存分析工具: Cox (statsmodels PHReg) / Kaplan-Meier / log-rank /
Fine-Gray 竞争风险亚分布风险模型 / CIF。

不依赖 R 和 lifeline(环境里装不上), 全部用 numpy + statsmodels 实现。
"""
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg


# --------------------------------------------------------------------------- #
# Cox
# --------------------------------------------------------------------------- #
@dataclass
class CoxResult:
    name: str
    hr: float
    lo: float
    hi: float
    p: float
    n: int
    events: int
    coef: float = np.nan
    se: float = np.nan


def cox(time, event, X, names=None, strata=None, ties="efron"):
    """Cox 比例风险。X: (n, k) 或 (n,)。

    X 只有 1 列时即单变量模型; 多列时是**联合**多变量模型
    (早期版本对每列单独拟合, 导致"多变量"结果实际等于单变量 —— 已修正)。
    返回 list[CoxResult], 每个变量一项。
    """
    t = np.asarray(time, dtype=float)
    e = np.asarray(event, dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    m = np.isfinite(t) & np.isfinite(e) & np.isfinite(X).all(axis=1)
    t, e, X = t[m], e[m], X[m]
    if names is None:
        names = [f"x{i}" for i in range(X.shape[1])]
    nv = len(names)
    res = [CoxResult(nm, np.nan, np.nan, np.nan, np.nan, int(len(t)), int(e.sum()))
           for nm in names]
    try:
        fit = PHReg(t, X, status=e, ties=ties,
                    strata=strata[m] if strata is not None else None).fit()
        b = np.asarray(fit.params, dtype=float).ravel()
        se = np.asarray(fit.bse, dtype=float).ravel()
        for j in range(min(nv, len(b))):
            z = b[j] / se[j] if se[j] > 0 else np.nan
            p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
            res[j] = CoxResult(names[j], float(np.exp(b[j])),
                               float(np.exp(b[j] - 1.96 * se[j])),
                               float(np.exp(b[j] + 1.96 * se[j])), p,
                               int(len(t)), int(e.sum()), float(b[j]), float(se[j]))
    except Exception:
        pass
    return res


def cox_table(rows, title=""):
    df = pd.DataFrame([{"variable": r.name, "HR": round(r.hr, 3),
                        "CI95_low": round(r.lo, 3), "CI95_high": round(r.hi, 3),
                        "p": f"{r.p:.3g}" if np.isfinite(r.p) else "NA",
                        "n": r.n, "events": r.events} for r in rows])
    if title:
        print(f"\n  {title}")
    print(df.to_string(index=False))
    return df


# --------------------------------------------------------------------------- #
# Kaplan-Meier / log-rank
# --------------------------------------------------------------------------- #
def km(time, event):
    """返回 (t, S, n_risk, n_event) 阶跃曲线; t[0]=0, S[0]=1。"""
    t = np.asarray(time, dtype=float)
    e = np.asarray(event, dtype=float)
    m = np.isfinite(t) & np.isfinite(e)
    t, e = t[m], e[m]
    ut = np.unique(t[e > 0])
    ts = [0.0]; ss = [1.0]; nr = [len(t)]; ne = [0]
    S = 1.0
    for tt in ut:
        at_risk = (t >= tt).sum()
        d = ((t == tt) & (e > 0)).sum()
        if at_risk == 0:
            break
        S *= (1 - d / at_risk)
        ts.append(tt); ss.append(S); nr.append(int(at_risk)); ne.append(int(d))
    return np.array(ts), np.array(ss), np.array(nr), np.array(ne)


def logrank(t1, e1, t2, e2):
    t = np.concatenate([t1, t2]); e = np.concatenate([e1, e2])
    g = np.concatenate([np.zeros(len(t1)), np.ones(len(t2))])
    m = np.isfinite(t) & np.isfinite(e)
    t, e, g = t[m], e[m], g[m]
    ut = np.unique(t[e > 0])
    O1 = E1 = V = 0.0
    for tt in ut:
        at = (t >= tt)
        n = at.sum(); d = ((t == tt) & (e > 0)).sum()
        if n < 2:
            continue
        n1 = (at & (g == 0)).sum()
        d1 = ((t == tt) & (e > 0) & (g == 0)).sum()
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1) if n > 1 else 0.0
    if V <= 0:
        return np.nan, np.nan
    chi = (O1 - E1) ** 2 / V
    return float(chi), float(stats.chi2.sf(chi, 1))


# --------------------------------------------------------------------------- #
# 竞争风险: CIF + Fine-Gray
# --------------------------------------------------------------------------- #
def cif(time, status, cause=1):
    """Aalen-Johansen 累积发生率函数。status: 0 删失 / 1 目标事件 / 2 竞争事件。"""
    t = np.asarray(time, dtype=float); s = np.asarray(status, dtype=float)
    m = np.isfinite(t) & np.isfinite(s)
    t, s = t[m], s[m]
    ut = np.unique(t[s > 0])
    S = 1.0
    ts = [0.0]; c1 = [0.0]; c2 = [0.0]
    for tt in ut:
        at_risk = (t >= tt).sum()
        d1 = ((t == tt) & (s == 1)).sum()
        d2 = ((t == tt) & (s == 2)).sum()
        if at_risk == 0:
            break
        c1.append(c1[-1] + S * d1 / at_risk)
        c2.append(c2[-1] + S * d2 / at_risk)
        S *= (1 - (d1 + d2) / at_risk)
        ts.append(tt)
    return np.array(ts), np.array(c1), np.array(c2)


def _km_censoring(time, status):
    """删失分布生存函数 G(t) = P(C > t); 事件指示 = I(status==0)。"""
    t = np.asarray(time, dtype=float); s = np.asarray(status, dtype=float)
    c = (s == 0).astype(float)
    ut = np.unique(t[c > 0])
    order = np.argsort(t)
    G = 1.0
    grid_t = [0.0]; grid_G = [1.0]
    for tt in ut:
        at_risk = (t >= tt).sum()
        d = ((t == tt) & (c > 0)).sum()
        if at_risk == 0:
            break
        G *= (1 - d / at_risk)
        grid_t.append(tt); grid_G.append(G)
    gt = np.array(grid_t); gg = np.array(grid_G)

    def Gfun(x):
        # 右连续阶梯: G(t) = P(C > t) = 取最后一个 grid_t <= x 对应的 G
        idx = np.searchsorted(gt, x, side="right") - 1
        idx = np.clip(idx, 0, len(gg) - 1)
        return gg[idx]
    return Gfun


def finegray(time, status, X, names=None, max_iter=100, tol=1e-8):
    """Fine-Gray 亚分布风险模型 (cause 1 vs cause 2)。

    实现要点
    --------
    风险集 R(t):  {j : T_j >= t}  ∪  {j : T_j < t 且 status_j == 2}
        —— 发生竞争事件的个体永久留在风险集中(这正是"亚分布"的含义)。
    权重  w_j(t) = G(t) / G(T_j ∧ t)   (删失的 IPCW 校正)
    Sk(t) = Σ_{j∈R(t)} w_j e^{β'Z_j} Z_j^{⊗k}
    U(β)  = Σ_{i: status_i=1} [ Z_i - S1(T_i)/S0(T_i) ]
    I(β)  = Σ_{i: status_i=1} [ S2/S0 - (S1/S0)(S1/S0)' ]
    """
    t = np.asarray(time, dtype=float); s = np.asarray(status, dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    m = np.isfinite(t) & np.isfinite(s) & np.isfinite(X).all(axis=1)
    t, s, X = t[m], s[m], X[m]
    n, k = X.shape
    if names is None:
        names = [f"x{i}" for i in range(k)]
    G = _km_censoring(t, s)
    Gt = G(t)
    ev1 = np.where(s == 1)[0]
    if len(ev1) < 5:
        return None

    beta = np.zeros(k)
    for _ in range(max_iter):
        U = np.zeros(k); I = np.zeros((k, k))
        eta = X @ beta
        for i in ev1:
            ti = t[i]
            # 风险集
            r = (t >= ti) | ((t < ti) & (s == 2))
            if r.sum() < 2:
                continue
            tj = t[r]; zj = X[r]; gj = Gt[r]
            w = G(ti) / np.where(tj < ti, gj, G(ti))
            w = np.clip(w, 0, 1e6)
            e = np.exp(zj @ beta)
            wwe = w * e
            S0 = wwe.sum()
            S1 = (wwe[:, None] * zj).sum(axis=0)
            S2 = (wwe[:, None, None] * zj[:, :, None] * zj[:, None, :]).sum(axis=0)
            zbar = S1 / S0
            U += X[i] - zbar
            I += S2 / S0 - np.outer(zbar, zbar)
        try:
            step = np.linalg.solve(I + 1e-10 * np.eye(k), U)
        except np.linalg.LinAlgError:
            break
        step = np.clip(step, -5, 5)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    try:
        cov = np.linalg.inv(I + 1e-10 * np.eye(k))
        se = np.sqrt(np.diag(cov))
    except Exception:
        se = np.full(k, np.nan)
    out = []
    for j, nm in enumerate(names):
        z = beta[j] / se[j] if np.isfinite(se[j]) and se[j] > 0 else np.nan
        p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
        out.append(CoxResult(nm, float(np.exp(beta[j])),
                             float(np.exp(beta[j] - 1.96 * se[j])) if np.isfinite(se[j]) else np.nan,
                             float(np.exp(beta[j] + 1.96 * se[j])) if np.isfinite(se[j]) else np.nan,
                             p, n, int((s == 1).sum()), float(beta[j]), float(se[j])))
    return out
