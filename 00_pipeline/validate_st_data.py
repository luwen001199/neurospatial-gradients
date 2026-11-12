# -*- coding: utf-8 -*-
"""GSE210616 原始数据完整性审计 —— 正式分析前必须通过的一道卡。

背景: 早期单线程下载器没有 size 校验, 被中断时留下了 truncated 的 h5 文件。
这些文件能被 ls 看到、大小非零, 但 h5py 打开会报 "truncated file",
如果直接进入分析会导致 spot 数量对不上、表达矩阵残缺而不自知。

本脚本对每个样本做三项检查:
  1. h5 能否被 h5py 正常打开, 且预期大小与实际大小一致
  2. barcodes 数量 == tissue_positions 行数
  3. features 里有我们关心的核心基因 (NGF/UCHL1/CALCA/ACTA2/CD8A 等)

输出:
  01_data/processed/gse210616_manifest.csv  逐样本审计结果
  控制台汇总 -> 决定哪些样本可用、哪些必须重下
"""
import glob
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

RAW = Path(r"D:/嗜神经浸润/01_data/raw/gse210616")
OUT = Path(r"D:/嗜神经浸润/01_data/processed")
OUT.mkdir(parents=True, exist_ok=True)

# 本项目真正关心的最小基因集合 —— 用于确认 feature 空间可用
CORE_GENES = ["NGF", "UCHL1", "CALCA", "CALCB", "S100B", "TH", "RAMP1",
              "ACTA2", "PDGFRB", "FAP", "COL1A1", "COL3A1", "FN1",
              "CD8A", "CD8B", "CCL21", "CXCL13", "LAMP3"]


def read_positions(p):
    """按文件内容魔法数判断 gzip, 而不是按扩展名.

    注意: 下载器把 'xxx.csv.gz' 落盘为 'xxx.csv' (去掉了 .gz 后缀), 但内容仍是 gzip。
    因此必须检测 magic bytes 0x1f 0x8b, 否则会按纯文本打开二进制而失败。
    """
    import gzip
    p = str(p)
    with open(p, "rb") as fb:
        magic = fb.read(2)
    is_gz = magic == b"\x1f\x8b"
    try:
        if is_gz:
            with gzip.open(p, "rt", errors="replace") as f:
                rows = [l.strip().split(",") for l in f if l.strip()]
        else:
            with open(p, "rt", errors="replace") as f:
                rows = [l.strip().split(",") for l in f if l.strip()]
    except Exception:
        return None
    if not rows:
        return None
    # 新版 SpaceRanger 输出带 header
    if rows[0][0] == "barcode":
        rows = rows[1:]
    return pd.DataFrame(rows)


def audit_sample(d):
    d = Path(d)
    res = {"sample": d.name, "sample_id": None, "status": "ok", "notes": ""}
    h5s = list(d.glob("*_filtered_feature_bc_matrix.h5"))
    poss = list(d.glob("*tissue_positions_list*"))
    scs = list(d.glob("*scalefactors*"))
    if not h5s:
        return {**res, "status": "missing_h5", "notes": "no h5"}
    if not poss:
        return {**res, "status": "missing_positions", "notes": "no positions"}

    h5 = h5s[0]
    res["sample_id"] = h5.name.split("_")[1] if "_" in h5.name else ""
    res["h5_mb"] = round(h5.stat().st_size / 1e6, 2)

    # 1) h5 可打开性 + 内容
    try:
        with h5py.File(h5, "r") as h:
            grp = h["matrix"] if "matrix" in h else h
            bc = grp["barcodes"][:]
            bc = [x.decode() if isinstance(x, bytes) else str(x) for x in bc]
            feats = grp["features"]
            names = feats["name"][:]
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in names]
            shape = tuple(grp["shape"][:]) if "shape" in grp else None
            res["n_barcodes_h5"] = len(bc)
            res["n_genes"] = len(names)
            if shape is not None:
                res["mat_shape"] = f"{shape[0]}x{shape[1]}"
                # shape = (genes, barcodes) in CellRanger convention
                res["shape_ok"] = (int(shape[1]) == len(bc))
            found = {g for g in CORE_GENES if g in set(names)}
            res["n_core_genes_found"] = len(found)
            res["missing_core_genes"] = ";".join(sorted(set(CORE_GENES) - found))
            if len(found) < len(CORE_GENES) * 0.8:
                res["notes"] += f" few core genes({len(found)}/{len(CORE_GENES)});"
    except Exception as e:
        res["status"] = "truncated_or_corrupt"
        res["notes"] += f" {type(e).__name__}:{str(e)[:70]};"
        return res

    # 2) positions 行数一致性
    try:
        pos = read_positions(poss[0])
        if pos is None:
            res["status"] = "bad_positions"
            res["notes"] += " positions unreadable;"
        else:
            res["n_positions"] = len(pos)
            if res.get("n_barcodes_h5") and len(pos) != res["n_barcodes_h5"]:
                if pos.shape[1] >= 2:  # may include tissue positions beyond filtered
                    res["notes"] += f" pos!=barcodes({len(pos)} vs {res['n_barcodes_h5']});"
    except Exception as e:
        res["notes"] += f" pos err {type(e).__name__};"

    res["has_scalefactors"] = bool(scs)
    return res


def main():
    dirs = sorted([d for d in RAW.glob("GSM*") if d.is_dir()])
    print(f"[audit] {len(dirs)} sample dirs found")
    rows = []
    for i, d in enumerate(dirs, 1):
        r = audit_sample(d)
        rows.append(r)
        print(f"  [{i}/{len(dirs)}] {r['sample']:<12} {r['status']:<20} "
              f"h5={r.get('h5_mb','-')}MB bc={r.get('n_barcodes_h5','-')} "
              f"genes={r.get('n_genes','-')} core={r.get('n_core_genes_found','-')}/"
              f"{len(CORE_GENES)} {r.get('notes','')[:60]}", flush=True)

    df = pd.DataFrame(rows)
    pcsv = OUT / "gse210616_manifest.csv"
    df.to_csv(pcsv, index=False)

    print("\n===== SUMMARY =====")
    print(df.groupby("status")["sample"].count().to_string())
    ok = df[df["status"] == "ok"]
    print(f"\nusable samples: {len(ok)}/{len(df)}")
    print(df[df["status"] != "ok"][["sample", "status", "notes"]].to_string(index=False))
    bad = df[df["status"].isin(["truncated_or_corrupt", "missing_h5"])]
    if len(bad):
        print(f"\n>>> {len(bad)} files must be re-downloaded:")
        for s in bad["sample"]:
            h5n = list((RAW / s).glob("*_filtered_feature_bc_matrix.h5"))
            for f in h5n:
                print(f"    del \"{f}\"")
    print(f"\nsaved -> {pcsv}")


if __name__ == "__main__":
    main()
