import sys, numpy as np, pandas as pd
sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
import nsgrad.surv as nsurv
d = pd.read_csv(r"D:/嗜神经浸润/01_data/processed/metabric_analysis_table.csv")
print("cols sample:", [c for c in d.columns][:60])
print()
for c in ["NPPscore","CAF_ECM","NPP_resid_stroma","AGE","nodes_pos","NODES_POS","ER_pos","chemo","hormone"]:
    if c in d.columns:
        print(f"  {c:20s} dtype={d[c].dtype}  nonnull={d[c].notna().sum()}  mean={pd.to_numeric(d[c],errors='coerce').mean():.4f}")
