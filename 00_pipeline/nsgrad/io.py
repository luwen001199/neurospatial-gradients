# -*- coding: utf-8 -*-
"""Visium 数据读取 —— 从 SpaceRanger 原始输出构建带空间坐标的 AnnData.

针对 GSE210616 的实际文件形态做了适配:
  * h5 为标准 CellRanger/Feature-Barcode 矩阵 (CSC 稀疏)
  * tissue_positions_list 落盘时被去掉了 .gz 后缀但内容仍是 gzip,
    因此统一用 magic bytes 判定压缩, 不依赖扩展名
  * 通过 scalefactors_json.spot_diameter_fullres 把像素坐标换算为微米,
    保证后续环形缓冲区距离是真实物理距离而非像素数
"""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
import h5py

from .config import SPOT_DIAMETER_UM, QC


def _open_text(p):
    """按 magic bytes 自动识别 gzip, 返回 (opener, mode)."""
    with open(p, "rb") as fb:
        gz = fb.read(2) == b"\x1f\x8b"
    return (gzip.open, "rt") if gz else (open, "rt")


def load_scalefactors(p):
    op, _ = _open_text(p)
    try:
        with op(p, "rt", errors="replace") as f:
            return json.load(f)
    except Exception:
        return {}


def load_positions(p):
    """返回 DataFrame: barcode, in_tissue, array_row, array_col, x(um), y(um)"""
    op, _ = _open_text(p)
    with op(p, "rt", errors="replace") as f:
        rows = [l.strip().split(",") for l in f if l.strip()]
    if not rows:
        raise ValueError(f"empty positions: {p}")
    header = rows[0]
    if header[0] == "barcode":                # SpaceRanger v2+ 带表头
        rows = rows[1:]
    nf = len(rows[0])
    if nf >= 6:
        cols = ["barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres"]
    elif nf == 4:
        cols = ["barcode", "in_tissue", "array_row", "array_col"]
    else:
        raise ValueError(f"unexpected positions columns={nf}: {p}")
    df = pd.DataFrame(rows, columns=cols[:nf]).set_index("barcode")
    return df


def read_visium_sample(sample_dir, sample_id=None, apply_qc=True, verbose=False):
    """把一个 GSM 目录读成 AnnData, obs=spot, var=gene, obsm['spatial']=原始像素坐标,
    额外在 obs 里写入以微米为单位的 x_um / y_um。"""
    sample_dir = Path(sample_dir)
    h5f = sorted(sample_dir.glob("*_filtered_feature_bc_matrix.h5"))
    posf = sorted(sample_dir.glob("*tissue_positions*"))
    scf = sorted(sample_dir.glob("*scalefactors*"))
    if not h5f or not posf:
        raise FileNotFoundError(f"{sample_dir.name}: missing h5 or positions")

    h5_path, pos_path = h5f[0], posf[0]
    if sample_id is None:
        sample_id = sample_dir.name

    # ---------- matrix ----------
    with h5py.File(h5_path, "r") as h:
        grp = h["matrix"] if "matrix" in h else h
        barcodes = [x.decode() if isinstance(x, bytes) else str(x)
                    for x in grp["barcodes"][:]]
        ft = grp["features"]
        gname = [x.decode() if isinstance(x, bytes) else str(x)
                 for x in ft["name"][:]]
        gid = [x.decode() if isinstance(x, bytes) else str(x)
               for x in ft["id"][:]]
        shape = tuple(int(v) for v in grp["shape"][:])          # (genes, spots)
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        # CellRanger 存的是 genes x spots (CSC), AnnData 需要 spots x genes
        m = sp.csc_matrix((data, indices, indptr), shape=shape)
        X = m.T.tocsr()

    var = pd.DataFrame({"gene_id": gid}, index=gname)
    if var.index.duplicated().any():          # 同 symbol 多探针 -> 保留首个
        var = var[~var.index.duplicated(keep="first")]
        keep = [gname.index(g) for g in var.index]
        X = X[:, keep]
    ad = AnnData(X=X.tocsr(), obs=pd.DataFrame(index=barcodes), var=var)
    ad.var_names = ad.var_names.astype(str)

    # ---------- coordinates ----------
    pos = load_positions(pos_path)
    pos = pos.reindex(ad.obs_names)
    sf = load_scalefactors(scf[0]) if scf else {}

    # 像素 -> 微米: Visium spot 直径恒为 55um
    dia_px = sf.get("spot_diameter_fullres", None)
    if dia_px:
        um_per_px = SPOT_DIAMETER_UM / float(dia_px)
    else:
        um_per_px = SPOT_DIAMETER_UM / 45.0     # fallback 经验值
    if "pxl_col_in_fullres" in pos.columns:
        full_x = pos["pxl_col_in_fullres"].astype(float) * um_per_px
        full_y = pos["pxl_row_in_fullres"].astype(float) * um_per_px
    else:
        full_x = pos["array_col"].astype(float) * SPOT_DIAMETER_UM * 2
        full_y = pos["array_row"].astype(float) * SPOT_DIAMETER_UM * 2

    ad.obs["in_tissue"] = pos["in_tissue"].astype(str).eq("1").values
    ad.obs["array_row"] = pos["array_row"].astype(float).values
    ad.obs["array_col"] = pos["array_col"].astype(float).values
    ad.obs["x_um"] = full_x.values
    ad.obs["y_um"] = full_y.values
    ad.obsm["spatial"] = np.c_[
        pos["pxl_col_in_fullres"].astype(float).values
        if "pxl_col_in_fullres" in pos.columns else full_x.values,
        pos["pxl_row_in_fullres"].astype(float).values
        if "pxl_row_in_fullres" in pos.columns else full_y.values,
    ]
    ad.uns["spatial_um_per_px"] = um_per_px
    ad.uns["sample_id"] = sample_id
    ad.uns["source_dir"] = str(sample_dir)

    # ---------- basic metrics ----------
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    ad.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    ad.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
    mt_mask = ad.var["mt"].values
    if mt_mask.any():
        ad.obs["pct_mt"] = np.asarray(X[:, mt_mask].sum(axis=1)).ravel() \
                            / np.maximum(ad.obs["total_counts"].values, 1) * 100
    else:
        ad.obs["pct_mt"] = 0.0

    n0 = ad.n_obs
    if apply_qc:
        # pandas>=3 的列可能是只读视图, 必须 .copy() 后才能原地 &= 运算
        ok = np.asarray(ad.obs["in_tissue"].values).copy().astype(bool)
        ok &= np.asarray(ad.obs["total_counts"].values) >= QC["min_umi_per_spot"]
        ok &= np.asarray(ad.obs["n_genes"].values) >= QC["min_genes_per_spot"]
        ok &= np.asarray(ad.obs["pct_mt"].values) <= QC["mito_max_pct"]
        ad = ad[ok].copy()
    if verbose:
        print(f"    {sample_id}: spots {n0} -> {ad.n_obs} after QC, "
              f"genes={ad.n_vars}, um/px={um_per_px:.4f}")
    return ad


def list_ready_samples(manifest_csv):
    """从完整性审计表里取出 status=='ok' 的样本目录."""
    import pandas as pd
    df = pd.read_csv(manifest_csv)
    ok = df[df["status"] == "ok"]["sample"].tolist()
    return sorted(ok)
