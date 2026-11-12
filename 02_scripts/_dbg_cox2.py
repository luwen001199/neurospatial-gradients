import sys, numpy as np, pandas as pd
sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
import nsgrad.surv as nsurv
d = pd.read_csv(r"D:/嗜神经浸润/01_data/processed/metabric_analysis_table.csv")
def z(s):
    s = pd.to_numeric(s, errors="coerce")
    return (s - s.mean()) / s.std(ddof=0)
for outcome, tcol, ecol in [("RFS","RFS_MONTHS","RFS_EVENT"),("OS","OS_MONTHS","OS_EVENT")]:
    sub = d.dropna(subset=[tcol, ecol]); sub = sub[sub[tcol] > 0]
    for var in ["NPPscore","NPP_resid_stroma"]:
        covs = [var,"AGE","nodes_pos","ER_pos","chemo","hormone","CAF_ECM"]
        s2 = sub.dropna(subset=covs)
        Zc = pd.DataFrame({c: z(s2[c]) for c in covs})
        keep = [c for c in covs if np.isfinite(Zc[c]).all() and Zc[c].std(ddof=0) > 1e-9]
        r_uni = nsurv.cox(s2[tcol], s2[ecol], Zc[var].values, names=[var])[0]
        r_mul = [r for r in nsurv.cox(s2[tcol], s2[ecol], Zc[keep].values, names=keep) if r.name==var][0]
        print(f"{outcome} {var:18s} n={len(s2)} ev={int(s2[ecol].sum())} | uni HR={r_uni.hr:.3f} ({r_uni.lo:.3f}-{r_uni.hi:.3f}) p={r_uni.p:.3g} | adj HR={r_mul.hr:.3f} ({r_mul.lo:.3f}-{r_mul.hi:.3f}) p={r_mul.p:.3g}")
        print(f"      keep={keep}")
