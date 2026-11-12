# -*- coding: utf-8 -*-
"""Supplementary Figure S6 — functional readouts of the perineural niche.

A  per-sample gradient distributions for the measurable cell-type panels
B  size-matched random baseline: |mean rho| versus gene-set size (510 random sets)
C  ligand-receptor measurability: |mean rho| versus number of assessable samples
D  pathway volcano: excess over random baseline versus consistency z
E  empirical CDF of |mean rho| by readout family versus random sets
F  ring profiles for the measurable cell-type panels
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "01_data" / "processed"
TAB = ROOT / "03_results" / "tables"
SUP = TAB / "supplementary"
FIG = ROOT / "03_results" / "figures"

C_MAIN = "#2E5F8A"
C_ALT = "#B4453C"
C_GREY = "#8A8F98"
C_LIGHT = "#C6CBD1"

plt.rcParams.update({
    "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.labelsize": 7, "axes.titlesize": 7.5, "legend.fontsize": 6,
    "font.family": "DejaVu Sans",
})


def load():
    per = pd.read_csv(PROC / "niche_per_sample.csv.gz")
    rnd = pd.read_csv(PROC / "niche_per_sample.random.csv.gz")
    cal = pd.read_csv(SUP / "TableS4_readout_measurability_and_circularity.csv")
    ring = pd.read_csv(SUP / "TableS13_readout_ring_profiles.csv")
    return per, rnd, cal, ring


def panel_a(per, cal, ax):
    keep = cal[(cal.kind == "celltype") & cal.measurable]["name"].tolist()
    g = per[per.name.isin(keep)]
    order = g.groupby("name")["rho"].mean().sort_values().index.tolist()
    data = [g.loc[g.name == n, "rho"].values for n in order]
    bp = ax.boxplot(data, vert=False, widths=0.6, patch_artist=True,
                    medianprops=dict(color="black", lw=0.8),
                    flierprops=dict(marker="o", ms=1.6, mfc=C_GREY,
                                    mec="none", alpha=0.6),
                    whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    for b in bp["boxes"]:
        b.set(facecolor="#DCE6F1", edgecolor=C_MAIN, lw=0.5)
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.set_yticklabels([n.replace("_", " ") for n in order], fontsize=5.6)
    ax.set_xlabel("partial ρ (per sample); negative = enriched near nerve")
    ax.set_title(f"A  Per-sample cell-type gradients\n({len(keep)} measurable panels)",
                 loc="left", fontsize=7)


def panel_b(rnd, cal, ax):
    rr = rnd.groupby("name").agg(size=("name", "size"), rho=("rho", "mean")).reset_index()
    rr["size"] = [int(n.split("|")[1][1:]) for n in rr["name"]]
    rr["arho"] = rr["rho"].abs()
    ax.scatter(rr["size"], rr["arho"], s=7, color=C_LIGHT, alpha=0.7,
               edgecolors="none", label="random gene sets (n=510)")
    sz = np.sort(rr["size"].unique())
    mu = [rr.loc[rr["size"] == s, "arho"].mean() for s in sz]
    sd = [rr.loc[rr["size"] == s, "arho"].std(ddof=1) for s in sz]
    mu, sd = np.array(mu), np.array(sd)
    ax.plot(sz, mu, "-", color="black", lw=1.0, label="baseline mean")
    ax.fill_between(sz, mu - sd, mu + sd, color="black", alpha=0.12,
                    label="±1 SD")
    for kind, col, mk in [("celltype", C_MAIN, "o"), ("lr", C_ALT, "s"),
                          ("pathway", "#4E8B57", "^")]:
        s = cal[cal.kind == kind]
        ax.scatter(s["size"], s["mean_rho"].abs(), s=9, color=col, marker=mk,
                   alpha=0.65, edgecolors="none", label=f"{kind} readouts")
    ax.set_xscale("log")
    ax.set_xlabel("gene-set size (log scale)")
    ax.set_ylabel("|mean ρ|")
    ax.set_title("B  Calibration baseline: random sets vs real readouts",
                 loc="left")
    ax.legend(frameon=False, fontsize=5.2, loc="upper left")


def panel_c(cal, ax):
    lr = cal[cal.kind == "lr"].copy()
    ok = lr["measurable"].to_numpy(bool)
    ax.scatter(lr.loc[ok, "n_samples"], lr.loc[ok, "mean_rho"].abs(), s=14,
               color=C_MAIN, edgecolors="white", lw=0.4,
               label=f"measurable (n={int(ok.sum())})")
    ax.scatter(lr.loc[~ok, "n_samples"], lr.loc[~ok, "mean_rho"].abs(), s=14,
               color=C_LIGHT, edgecolors="white", lw=0.4,
               label=f"below threshold (n={int((~ok).sum())})")
    ax.axvline(10, color=C_ALT, lw=0.8, ls="--")
    ax.text(10.6, 0.015, "measurability\nthreshold", fontsize=5.0,
            color=C_ALT, va="bottom")
    hi = lr[~ok].nlargest(4, "mean_rho")
    for _, r in hi.iterrows():
        nm = r["name"].split("|")[1] if "|" in r["name"] else r["name"]
        ax.annotate(nm, (r["n_samples"], abs(r["mean_rho"])), fontsize=4.6,
                    xytext=(3, -4), textcoords="offset points", color="#555555")
    ax.set_xlabel("samples in which the pair was assessable (of 30)")
    ax.set_ylabel("|mean ρ|")
    ax.set_title("C  Ligand–receptor measurability audit: large effects "
                 "concentrate at small n", loc="left")
    ax.legend(frameon=False, fontsize=5.2, loc="upper right")


def panel_d(cal, ax):
    pw = cal[cal.kind == "pathway"]
    pas = pw["pass_calibrated"].to_numpy(bool)
    ax.scatter(pw.loc[~pas, "z_vs_random"], pw.loc[~pas, "z_vs_random_Z"],
               s=5, color=C_LIGHT, alpha=0.55, edgecolors="none",
               label=f"not passing (n={int((~pas).sum())})")
    ax.scatter(pw.loc[pas, "z_vs_random"], pw.loc[pas, "z_vs_random_Z"],
               s=22, color=C_MAIN, edgecolors="white", lw=0.4,
               label=f"passing all thresholds (n={int(pas.sum())})")
    ax.axvline(0, color="black", lw=0.4)
    ax.axhline(0, color="black", lw=0.4)
    for i, (_, r) in enumerate(pw[pas].nlargest(6, "z_vs_random").iterrows()):
        nm = (r["name"].replace("REACTOME_", "").replace("HALLMARK_", "")
              .replace("KEGG_", "").replace("_", " "))
        dy = (i % 3) * 0.55 - 0.2
        ax.annotate(nm[:30], (r["z_vs_random"], r["z_vs_random_Z"] + dy),
                    fontsize=4.3, xytext=(4, 0), textcoords="offset points")
    ax.set_xlabel("excess over size-matched random baseline (effect size)")
    ax.set_ylabel("excess over random baseline (consistency)")
    ax.set_title("D  Pathways: effect size vs between-sample consistency",
                 loc="left")
    ax.legend(frameon=False, fontsize=5.2, loc="upper left")


def panel_e(per, rnd, cal, ax):
    rnd_m = rnd.groupby("name")["rho"].mean().abs()
    x = np.linspace(0, 0.20, 400)
    ax.plot(x, [np.mean(rnd_m <= v) for v in x], color="black", lw=1.2,
            label="random gene sets (n=510)")
    for kind, col in [("celltype", C_MAIN), ("lr", C_ALT),
                      ("pathway", "#4E8B57")]:
        s = cal[cal.kind == kind]["mean_rho"].abs()
        ax.plot(np.sort(s), np.arange(1, len(s) + 1) / len(s), lw=1.0,
                color=col, label=f"{kind} (n={len(s)})")
    ax.set_xlabel("|mean partial ρ|")
    ax.set_ylabel("empirical cumulative fraction")
    ax.set_title("E  Empirical CDFs: real readouts vs random baseline",
                 loc="left")
    ax.legend(frameon=False, fontsize=5.2, loc="lower right")


def panel_f(ring, cal, ax):
    keep = cal[(cal.kind == "celltype") & cal.measurable]["name"].tolist()
    d = ring.set_index("ring")[keep]
    mapping = {"r1_0_50": "0–50", "r2_50_100": "50–100",
               "r3_100_200": "100–200", "r4_200_400": "200–400",
               "r5_gt400": ">400"}
    order = [o for o in mapping if o in d.index] or list(d.index)
    labels = [mapping.get(o, o) for o in order]
    d = d.loc[order]
    for col in d.columns:
        ax.plot(range(len(d)), d[col].values, "-o", ms=3, lw=1.0,
                label=col.replace("_", " "))
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels(labels, fontsize=5.6)
    ax.set_xlabel("distance to nearest nerve source")
    ax.set_ylabel("mean score (z, cross-sample mean)")
    ax.set_title("F  Ring profiles: measurable cell-type panels", loc="left")
    ax.legend(frameon=False, fontsize=4.8, ncol=2, loc="center left",
              bbox_to_anchor=(1.0, 0.5))


def main():
    per, rnd, cal, ring = load()
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.2))
    panel_a(per, cal, axes[0, 0])
    panel_b(rnd, cal, axes[0, 1])
    panel_c(cal, axes[0, 2])
    panel_d(cal, axes[1, 0])
    panel_e(per, rnd, cal, axes[1, 1])
    panel_f(ring, cal, axes[1, 2])
    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"SuppFigS6_functional_readouts.{ext}",
                    dpi=400 if ext == "png" else None,
                    bbox_inches="tight")
    print("  -> SuppFigS6_functional_readouts.png/.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
