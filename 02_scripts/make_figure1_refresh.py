# -*- coding: utf-8 -*-
"""只重跑 Fig1 中不依赖原始 h5 的面板: A (流程示意) + G (拟合质量) + H (丰度/效应量)。

用于修正:
  * Fig1A 旧的 K=2..6 / 三型标签 -> K=2, Neural-Active / Neural-Cold
  * Fig1H 旧数据 (median=15.0%) -> 当前 sample_summary.csv / nerve_validation.csv
B–F 为空间图与拟合曲线, 由原始 h5 直接生成, 不受本次修订影响, 不重跑。
"""
import sys
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
sys.path.insert(0, str(ROOT / "00_pipeline"))
sys.path.insert(0, str(ROOT / "02_scripts"))

import make_figure1 as F1  # noqa: E402
import pandas as pd  # noqa: E402

PROC = F1.PROC

if __name__ == "__main__":
    gp = pd.read_csv(PROC / "gradient_params.csv")
    nrv = pd.read_csv(PROC / "nerve_validation.csv")
    prev = pd.read_csv(PROC / "sample_summary.csv")

    F1.fig1a();  print("Fig1A done (K=2 labels)")
    F1.fig1g(gp); print("Fig1G done")
    F1.fig1h(prev, nrv); print("Fig1H done")
    print("median pct_nerve_spots =", prev["pct_nerve_spots"].median(),
          "| range", prev["pct_nerve_spots"].min(), "-", prev["pct_nerve_spots"].max())
