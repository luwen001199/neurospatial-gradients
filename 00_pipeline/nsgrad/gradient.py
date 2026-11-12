# -*- coding: utf-8 -*-
"""环形缓冲区构建 + 指数衰减梯度拟合 —— 本研究的灵魂。

模型:  y(d) = a * exp(-b * d) + c

  a  近端初始信号强度 (距离神经源点趋近 0 时的读数)
  b  衰减常数 -- 核心指标, 越大表示信号随距离衰减越快
  c  远端本底 (没有神经影响时的基线)
  1/b  拐点距离 / 有效作用半径

为什么用这个模型而不是线性回归:
  指数衰减是扩散型信号(神经肽通过组织间隙扩散)的自然形式,
  线性模型无法给出"作用半径"这个有物理意义的参数, 而这正是本文概念创新的落点。

稳健性:
  * 拟合失败时返回 NaN 而不是硬报错, 由上层决定降级为距离分箱(近/中/远),
    对应方案"阴性风险控制矩阵"里的 b≈0 情形。
  * 同时输出 R^2 与 RMSE, 用于过滤拟合质量过差的样本。
"""
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree

from .config import RING_EDGES_UM, RING_LABELS


def build_distance_fields(ad, nerve_mask, max_dist_um=None):
    """对每个 spot 计算到最近神经源点的欧氏距离(um)与所属环带。

    参数
    ----
    ad : AnnData, 需要提供 obs['x_um'], obs['y_um']
    nerve_mask : shape (n_spot,) bool, 神经源点
    返回
    ----
    DataFrame: distance_to_nearest_nerve, nearest_nerve_idx, ring_id
    """
    xy = np.c_[ad.obs["x_um"].values, ad.obs["y_um"].values]
    n_nerve = int(nerve_mask.sum())
    out = pd.DataFrame(index=ad.obs_names)
    if n_nerve == 0:
        out["distance_to_nearest_nerve"] = np.nan
        out["nearest_nerve_idx"] = -1
        out["ring_id"] = "no_nerve"
        return out, (0, -1)

    src = xy[nerve_mask]
    tree = cKDTree(src)
    dist, idx = tree.query(xy, k=1)
    out["distance_to_nearest_nerve"] = dist
    out["nearest_nerve_idx"] = idx
    out["ring_id"] = pd.cut(dist, bins=RING_EDGES_UM + [np.inf],
                            labels=RING_LABELS, include_lowest=True).astype(str)
    # 神经 spot 自身距离为 0
    out.loc[nerve_mask, "distance_to_nearest_nerve"] = 0.0
    out.loc[nerve_mask, "ring_id"] = "r1_0_50"
    return out, (n_nerve, tree)


def ring_profile(ad, scores: Dict[str, np.ndarray], dist_field,
                 exclude_nerve_spots=True, nerve_mask=None):
    """每个环带内的各指标均值/中位数。返回 tidy DataFrame。"""
    df = pd.DataFrame(scores, index=ad.obs_names)
    df = df.join(dist_field[["distance_to_nearest_nerve", "ring_id"]])
    if exclude_nerve_spots and nerve_mask is not None:
        df = df[~np.asarray(nerve_mask)]
    rows = []
    for ring, sub in df.groupby("ring_id", observed=False):
        if len(sub) == 0:
            continue
        rec = {"ring": ring, "n_spots": len(sub),
               "mean_dist_um": float(sub["distance_to_nearest_nerve"].mean())}
        for c in scores:
            rec[f"{c}_mean"] = float(sub[c].mean())
            rec[f"{c}_median"] = float(sub[c].median())
        rows.append(rec)
    # 保证环带顺序
    order = {lbl: i for i, lbl in enumerate(RING_LABELS)}
    prof = pd.DataFrame(rows)
    if len(prof):
        prof["_o"] = prof["ring"].map(order).fillna(99)
        prof = prof.sort_values("_o").drop(columns="_o")
    return prof


# b 的物理约束: 有效作用半径限定在 25~2000 um 之间。
# 低于 25um 意味着信号只在一个 spot 内、超出分辨率;
# 高于 2000um 意味着在整个组织切片尺度上几乎不衰减, 与"局灶神经信号"的生物学前提矛盾,
# 且此时 exp(-b*d) 趋近常数, a 与 c 不可辨识(表现为 a 撞到天文数字), 属退化解。
B_BOUNDS = (1.0 / 2000.0, 1.0 / 25.0)


def _exp_model(d, a, b, c):
    return a * np.exp(-b * d) + c


def _fit_core(x, y, w=None, y_vals=None, n_restarts=6):
    """带多起点重启的指数拟合。x,y 为聚合后的曲线点(而非 raw spot)。"""
    yv = np.asarray(y, float) if y_vals is None else np.asarray(y_vals, float)
    yrange = float(np.ptp(yv)) if len(yv) else 1.0
    yrange = max(yrange, 1e-3)
    ylo, yhi = float(np.min(yv)) - yrange, float(np.max(yv)) + yrange
    ab = max(5 * yrange, 1e-3)
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    w = np.maximum(w, 1e-9)

    best, best_sse = None, np.inf
    # 多起点: 假设作用半径分别为 50/80/120/200/400/800 um
    starts = [1.0 / r for r in (50, 80, 120, 200, 400, 800)][:max(n_restarts, 1)]
    for b0 in starts:
        for a0 in (yrange * 0.5, -yrange * 0.5, yrange * 0.05):
            try:
                popt, _ = curve_fit(_exp_model, x, y, sigma=1.0 / np.sqrt(w),
                                    absolute_sigma=False,
                                    p0=[a0, b0, float(np.median(y))],
                                    bounds=([-ab, B_BOUNDS[0], ylo],
                                            [ab, B_BOUNDS[1], yhi]),
                                    maxfev=40000)
                sse = float((((y - _exp_model(x, *popt)) ** 2) * w).sum())
                if sse < best_sse and np.isfinite(sse):
                    best_sse, best = sse, popt
            except Exception:
                continue
    if best is None:
        return None
    a, b, c = (float(v) for v in best)
    pred = _exp_model(x, *best)
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {"a": a, "b": b, "c": c,
            "r2": float(1 - ss_res / ss_tot) if ss_tot > 1e-12 else np.nan,
            "rmse": float(np.sqrt(ss_res / max(len(x), 1))),
            "decay_points_um": float(1.0 / b) if b > 1e-9 else np.nan}


def aggregate_distance_curve(dist, y, n_bins=12, min_per_bin=5):
    """把 raw spot 聚合成距离-表达曲线。

    这是**主**拟合路径: 直接拿上千个 spot 的散点去拟合 3 参数非线性模型,
    噪声会完全压过信号(实测 R2<0.05 且 a 退化为天文数字)。
    按距离分位数分箱后取加权均值, 既压噪又保留曲线形状, 是空间梯度分析的通行做法。
    """
    dist = np.asarray(dist, float)
    y = np.asarray(y, float)
    m = np.isfinite(dist) & np.isfinite(y)
    dist, y = dist[m], y[m]
    if len(dist) < 30:
        return None
    qs = np.quantile(dist, np.linspace(0, 1, n_bins + 1))
    qs = np.unique(qs)
    if len(qs) < 4:
        return None
    qs[0] -= 1e-6
    idx = np.digitize(dist, qs) - 1
    xs, ys, ws = [], [], []
    for b in range(len(qs) - 1):
        sel = idx == b
        if sel.sum() >= min_per_bin:
            xs.append(np.median(dist[sel]))
            ys.append(np.mean(y[sel]))
            ws.append(int(sel.sum()))
    if len(xs) < 4:
        # 分位数分箱失败时退化为等宽分箱
        ed = np.linspace(dist.min(), dist.max(), n_bins + 1)
        idx = np.digitize(dist, ed) - 1
        xs, ys, ws = [], [], []
        for b in range(n_bins):
            sel = idx == b
            if sel.sum() >= min_per_bin:
                xs.append(np.median(dist[sel]))
                ys.append(np.mean(y[sel]))
                ws.append(int(sel.sum()))
    if len(xs) < 4:
        return None
    return np.asarray(xs), np.asarray(ys), np.asarray(ws, float)


def fit_decay_aggregated(dist, y, n_bins=12, min_per_bin=5):
    """主拟合入口: 距离聚合曲线 -> 多起点加权非线性拟合。"""
    agg = aggregate_distance_curve(dist, y, n_bins=n_bins, min_per_bin=min_per_bin)
    if agg is None:
        return {"a": np.nan, "b": np.nan, "c": np.nan, "r2": np.nan,
                "rmse": np.nan, "n_points": int(np.isfinite(dist).sum()),
                "ok": False, "note": "insufficient for aggregation"}
    x, yy, w = agg
    r = _fit_core(x, yy, w, y_vals=y)
    if r is None:
        return {"a": np.nan, "b": np.nan, "c": np.nan, "r2": np.nan,
                "rmse": np.nan, "n_points": len(y), "ok": False,
                "note": "no convergent restart"}
    r.update({"rmse": r["rmse"], "ok": True, "note": "aggregated",
              "n_points": int(np.isfinite(dist).sum()) ,
              "n_bins": len(x)})
    return r


def fit_decay_binned(dist, y, n_bins=10, min_per_bin=5):
    """降级方案: 用环带/分箱均值而非 raw spot 做拟合, 抗噪更强。

    当 spot 级拟合不收敛时先走这里, 对应风险控制矩阵里的 "梯度拟合失败" 路径。
    """
    dist = np.asarray(dist, float)
    y = np.asarray(y, float)
    m = np.isfinite(dist) & np.isfinite(y) & (dist > 0)
    dist, y = dist[m], y[m]
    if len(dist) < 15:
        return {"a": np.nan, "b": np.nan, "c": np.nan, "r2": np.nan,
                "rmse": np.nan, "n_points": len(dist), "ok": False,
                "note": "insufficient points"}
    qs = np.quantile(dist, np.linspace(0, 1, n_bins + 1))
    qs[0] -= 1e-6
    idx = np.digitize(dist, qs) - 1
    xs, ys, ws = [], [], []
    for b in range(n_bins):
        sel = idx == b
        if sel.sum() >= min_per_bin:
            xs.append(dist[sel].mean())
            ys.append(y[sel].mean())
            ws.append(sel.sum())
    if len(xs) < 4:
        return {"a": np.nan, "b": np.nan, "c": np.nan, "r2": np.nan,
                "rmse": np.nan, "n_points": len(dist), "ok": False,
                "note": "too few bins"}
    xs, ys, ws = np.array(xs), np.array(ys), np.array(ws, float)
    try:
        popt, _ = curve_fit(_exp_model, xs, ys, w=np.sqrt(ws),
                            p0=[ys.max() - ys.min(), 0.01, ys.min()],
                            bounds=([-np.inf, 1e-6, -np.inf], [np.inf, 5.0, np.inf]),
                            maxfev=20000)
        pred = _exp_model(xs, *popt)
        ss_res = float(((ys - pred) ** 2).sum())
        ss_tot = float(((ys - ys.mean()) ** 2).sum())
        return {"a": float(popt[0]), "b": float(popt[1]), "c": float(popt[2]),
                "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
                "rmse": float(np.sqrt(ss_res / len(xs))),
                "n_points": len(dist), "ok": True, "note": "binned",
                "decay_points_um": float(1.0 / popt[1]) if popt[1] > 1e-9 else np.nan}
    except Exception as e:
        return {"a": np.nan, "b": np.nan, "c": np.nan, "r2": np.nan,
                "rmse": np.nan, "n_points": len(dist), "ok": False,
                "note": f"{type(e).__name__}"}


def fit_sample_gradients(ad, scores_dict, nerve_mask, metrics=None, verbose=True):
    """对单个样本的全部指标做梯度拟合, 返回 tidy DataFrame(每个指标一行)。"""
    dist_field, meta = build_distance_fields(ad, nerve_mask)
    dist = dist_field["distance_to_nearest_nerve"].values
    rows = []
    metrics = metrics or list(scores_dict.keys())
    for mname in metrics:
        y = np.asarray(scores_dict[mname], dtype=float)
        # 主路径: 距离聚合曲线 + 多起点加权拟合
        r = fit_decay_aggregated(dist, y)
        r = dict(r)
        if not r["ok"]:
            # 二次退路: 放宽最小每箱点数再试一次
            r2b = fit_decay_aggregated(dist, y, n_bins=6, min_per_bin=3)
            if r2b["ok"]:
                r = dict(r2b)
        r.setdefault("fallback_used", False)
        r["metric"] = mname
        rows.append(r)
        if verbose:
            print(f"      {mname:<16} b={r['b']:.5f}  1/b={r.get('decay_points_um', float('nan')):>8.1f}um  "
                  f"R2={r['r2']:.3f}  ok={r['ok']}{'  [binned fallback]' if r['fallback_used'] else ''}")
    df = pd.DataFrame(rows)
    return df, dist_field


def discrete_ring_features(ad, scores_dict, nerve_mask, metrics=None, near_max_um=150, far_min_um=400):
    """离散降级特征: 近神经 / 远神经分箱比值, 当连续梯度不成立时使用。"""
    dist_field, _ = build_distance_fields(ad, nerve_mask)
    d = dist_field["distance_to_nearest_nerve"].values
    near = d <= near_max_um
    far = (d >= far_min_um)
    out = {}
    for k, v in scores_dict.items():
        y = np.asarray(v, dtype=float)
        if near.sum() >= 5 and far.sum() >= 5:
            out[f"{k}_near"] = float(np.mean(y[near]))
            out[f"{k}_far"] = float(np.mean(y[far]))
            out[f"{k}_ratio"] = float(np.mean(y[near]) / (np.mean(y[far]) + 1e-9))
        else:
            out[f"{k}_near"] = np.nan
            out[f"{k}_far"] = np.nan
            out[f"{k}_ratio"] = np.nan
    return out, int(near.sum()), int(far.sum())
