# -*- coding: utf-8 -*-
"""把 GSE213688_RAW.tar 解包并整理成 nsgrad 期望的 per-GSM 目录结构。

扁平文件名形如:
    GSM6592048_M1_filtered_feature_bc_matrix.h5
    GSM6592048_M1_tissue_positions_list.csv.gz
    GSM6592048_M1_scalefactors_json.json.gz
    GSM6592048_M1_pathologist_annotations.csv.gz
整理后:
    01_data/raw/gse213688/GSM6592048/<同上, 保留原文件名>
"""
import re
import shutil
import tarfile
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
RAW = ROOT / "01_data" / "raw" / "gse213688"
TAR = RAW / "GSE213688_RAW.tar"
FLAT = RAW / "_flat"

KEEP = ("_filtered_feature_bc_matrix.h5", "_tissue_positions", "_scalefactors",
        "_pathologist_annotations")


def main():
    if not TAR.exists():
        raise SystemExit(f"missing {TAR}")
    FLAT.mkdir(exist_ok=True)
    print(f"[prep] extracting {TAR.stat().st_size/1e6:.0f} MB ...")
    with tarfile.open(TAR) as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        print(f"[prep] {len(members)} files in archive")
        tf.extractall(FLAT)
    n_ok = 0
    for f in sorted(FLAT.iterdir()):
        if not f.is_file():
            continue
        m = re.match(r"(GSM\d+)_", f.name)
        if not m or not any(k in f.name for k in KEEP):
            continue
        d = RAW / m.group(1)
        d.mkdir(exist_ok=True)
        shutil.move(str(f), str(d / f.name))
        n_ok += 1
    print(f"[prep] moved {n_ok} files into {len(list(RAW.glob('GSM*')))} sample dirs")
    for d in sorted(RAW.glob("GSM*")):
        have = sorted(p.name.split("_", 2)[-1] for p in d.iterdir())
        print(f"  {d.name}: {len(have)} files | h5={'Y' if any('h5' in h for h in have) else 'N'} "
              f"pos={'Y' if any('positions' in h for h in have) else 'N'}")


if __name__ == "__main__":
    main()
