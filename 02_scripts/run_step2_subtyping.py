# -*- coding: utf-8 -*-
"""Step2 空间神经-免疫-基质架构分型.

设计原则(读方法学时可直接引用):
1. 特征矩阵构建在**样本级**, 每个 spot 级指标通过中位数/比值聚合,
   避免样本间 spot 数不同导致大样本主导聚类。
2. 梯度特征一律用 **近/远比值** (near<=150um / far>=400um) 而不是拟合的 b,
   因为 25 个样本里近半数的 b 撞在参数边界上(见 diagnose_step1),
   不可辨识的 b 会把聚类拉向噪声。比值特征对全部样本都有定义, 更稳健。
3. 共识聚类 (Monti 2003) 自己实现, K=2..6 遍历, 用 PAC + 轮廓系数定 K;
   1000 次子抽样(80%), 层次聚类 ward。
4. 分型不是终点 —— 会进入 Step4 与 METABRIC/GSE12276 的临床终点做独立验证。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
from nsgrad import config as CFG

PROC = CFG.DIR_PROC
NEAR_MAX_UM = 150.0
FAR_MIN_UM = 400.0
N_SUBSAMPLE = 0.8
N_BOOT = 1000
SEED = 42


def build_feature_matrix(ss: pd.DataFrame):
    """从 sample_summary 构建分型特征.

    关键修订(相对第一版):
    1. 用 **near - far 差值** 而非 near/far 比值。原因: 签名得分是 z-score 的均值,
       可为负值, 作分母后比值的符号与量级都不再可解释, 且会跨 4 个数量级,
       单个极端样本就能把聚类拉成 21 vs 1 的失衡。差值有明确物理意义:
       "神经邻域比远端平均高/低多少个 z 单位"。
    2. 神经密度与 b 做 log1p / log 变换, 压偏态, 避免离群样本主导距离度量。
    """
    feats = pd.DataFrame(index=ss["sample"].values)
    feats["nerve_density_log"] = np.log1p(ss["pct_nerve_spots"].values / 100.0)
    for key, label in [("Nervesignal", "nerve_sig_diff"),
                       ("CollagenIndex", "collagen_diff"),
                       ("myCAFscore", "myCAF_diff"),
                       ("CD8density", "CD8_diff"),
                       ("TLSscore", "TLS_diff"),
                       ("IFNsignature", "IFN_diff")]:
        near, far = ss.get(f"{key}_near"), ss.get(f"{key}_far")
        if near is not None and far is not None:
            feats[label] = near.values - far.values
    # b 仅在可辨识时使用; log 压缩量级
    bcols = [c for c in ss.columns if c.startswith("b_")]
    bs = ss[bcols].astype(float).copy()
    ok = (bs.values > 1 / 2000 * 1.001) & (bs.values < 1 / 25 * 0.999)
    bs = bs.where(pd.DataFrame(ok, index=bs.index, columns=bs.columns))
    feats["b_median_log"] = np.log(bs.median(axis=1).values
                                    + 1 / 2000.0)       # log(0) 保护
    feats["b_nsignal_log"] = np.log(ss["b_Nervesignal"].astype(float).values
                                     + 1 / 2000.0)
    feats.index.name = "sample"
    return feats


def znorm(X):
    mu = X.mean(0)
    sd = X.std(0).replace(0, 1.0)
    return (X - mu) / sd


def consensus_matrix(X, k, n_boot=N_BOOT, frac=N_SUBSAMPLE, seed=SEED):
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    M = np.zeros((n, n))
    C = np.zeros((n, n))
    L = linkage(X, method="ward")        # 全数据参考层次结构
    for _ in range(n_boot):
        idx = np.sort(rng.choice(n, int(n * frac), replace=False))
        sub = X[idx]
        if len(np.unique(sub, axis=0)) < k:
            continue
        labels = fcluster(linkage(sub, method="ward"), t=k, criterion="maxclust")
        M[np.ix_(idx, idx)] += 1
        for a in range(len(idx)):
            for b_ in range(a + 1, len(idx)):       # 严格 a<b_, 避免对角被计 2 次
                if labels[a] == labels[b_]:
                    C[idx[a], idx[b_]] += 1
                    C[idx[b_], idx[a]] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(M > 0, C / np.maximum(M, 1), 0)
    return P


def pac_score(P, lo=0.1, hi=0.9):
    """Proportion of Ambiguous Clustering: 共识矩阵值落在 (lo,hi) 的比例, 越小越稳定。"""
    triu = P[np.triu_indices_from(P, k=1)]
    return float(((triu > lo) & (triu < hi)).mean())


def pick_k(X, k_range=(2, 3, 4), seed=SEED, n_boot=N_BOOT):
    """K 选择: 用 **共识矩阵自身的 silhouette**(Monti 2003 推荐) 而非原始数据的 silhouette.

    原始 silhouette 在小样本上会随 K 单调上升(本次实测 K=6 最高但每群只剩 3-4 个),
    是典型偏差; 共识矩阵的 silhouette 衡量"同一 K 重复抽样下划分的可复现程度",
    与样本量约束一起使用更稳健。K 上限设为 4: 22 个样本分 5 群以上没有生物学意义。
    """
    res = {}
    for k in k_range:
        P = consensus_matrix(X, k, n_boot=n_boot, seed=seed)
        D = 1 - P
        np.fill_diagonal(D, 0)
        D = (D + D.T) / 2
        D = np.clip(D, 0, None)
        try:
            labels = fcluster(linkage(squareform(D, checks=False), method="average"),
                              t=k, criterion="maxclust")
        except Exception:
            continue
        # 共识矩阵的 silhouette: sklearn 的 precomputed 需要方形距离矩阵, 不是压缩向量
        try:
            sil_consensus = silhouette_score(D, labels, metric="precomputed")
        except Exception:
            sil_consensus = np.nan
        try:
            sil_data = silhouette_score(X, labels)
        except Exception:
            sil_data = np.nan
        res[k] = {"PAC": pac_score(P),
                  "silhouette_consensus": sil_consensus,
                  "silhouette_data": sil_data,
                  "labels": labels, "consensus": P}
    return res


def assign_subtype_names(features, labels):
    """按生物学含义给聚类命名, 避免 'cluster 0/1/2' 的黑箱感。

    K>=3 时:
      Neural-Barriered : 高神经密度 + CD8 近/远比低(免疫被排斥在神经邻域外)
      Neural-Inflamed  : 高神经密度 + CD8 近/远比高(免疫可浸润神经邻域)
      Neural-Cold      : 低神经密度
    K=2 时退化为二分:
      Neural-Cold      : 低神经密度
      Neural-Active    : 高神经密度(其内部是否再分 barriered/inflamed 由 K 选择决定)
    """
    F = features.copy()
    F["cluster"] = labels
    prof = F.groupby("cluster").mean()
    cold = int(prof["nerve_density_log"].idxmin())
    rest = [c for c in prof.index if c != cold]
    if len(rest) == 1:
        return {cold: "Neural-Cold", rest[0]: "Neural-Active"}, prof
    barriered = int(prof.loc[rest, "CD8_diff"].idxmin())
    inflamed = [c for c in rest if c != barriered][0]
    naming = {cold: "Neural-Cold", barriered: "Neural-Barriered",
              inflamed: "Neural-Inflamed"}
    return {c: naming[c] for c in prof.index}, prof


def main():
    ss = pd.read_csv(PROC / "sample_summary.csv")
    F = build_feature_matrix(ss)
    F = F.dropna(axis=1, how="all")
    F = F.fillna(F.median(numeric_only=True))
    Fz = znorm(F)
    print(f"[Step2] features: {Fz.shape[0]} samples x {Fz.shape[1]} features")
    print(Fz.describe().T[["mean", "std", "min", "max"]].round(2).to_string())

    print("\n[Step2] consensus clustering over K=2..4 ...")
    res = pick_k(Fz.values)
    K = pd.DataFrame([{"k": k, "PAC": v["PAC"],
                       "sil_consensus": round(v["silhouette_consensus"], 3),
                       "sil_data": round(v["silhouette_data"], 3)}
                      for k, v in res.items()]).sort_values("k")
    print(K.to_string(index=False))
    # 选 K: 优先共识矩阵 silhouette 最大, 平手取 PAC 较小
    Kc = K.dropna(subset=["sil_consensus"])
    best_k = int(Kc.sort_values(["sil_consensus", "PAC"],
                                ascending=[False, True]).iloc[0]["k"])
    print(f"      best K = {best_k}")

    labels = res[best_k]["labels"]
    naming, prof = assign_subtype_names(F, labels)
    out = pd.DataFrame({"sample": F.index, "cluster": labels,
                        "subtype": [naming[c] for c in labels]})
    out = out.merge(F.reset_index().rename(columns={"index": "sample"}), on="sample")
    out.to_csv(PROC / "sample_subtypes.csv", index=False)

    # 保存一致性矩阵和统计
    np.save(PROC / "consensus_matrix.npy", res[best_k]["consensus"])
    K.to_csv(PROC / "consensus_k_selection.csv", index=False)
    prof.to_csv(PROC / "subtype_profiles.csv")

    print("\n[Step2] subtype profile (feature means by subtype):")
    print(prof.round(3).to_string())
    print("\n[Step2] sample counts per subtype:")
    print(out["subtype"].value_counts().to_string())
    print(f"\n[Step2] saved -> {PROC / 'sample_subtypes.csv'}")


if __name__ == "__main__":
    main()
