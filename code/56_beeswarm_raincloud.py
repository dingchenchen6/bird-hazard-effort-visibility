#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
美观蜂群图 + 雨林图 / Publication beeswarm + raincloud of the climate x effort
interaction HR, across 4 effort specs x 3 datasets (v2 concurrent, v2 lagged t-1,
v3 relaxed). 自举 N(beta, se) -> exp -> HR draws; 直观展示交互在所有规格/数据集
上稳定 >1(含滞后内生性检验)。
输入: results/tables/{table_province_v2_coefs,table_effort_lag_refit,
       table_province_v3_all_specs_coefs}.csv
输出: figures/main/Figure_R1_interaction_raincloud.{png,pdf,svg}
      figures/main/Figure_R2_interaction_beeswarm.{png,pdf,svg}
运行: python3 code/56_beeswarm_raincloud.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

ROOT = Path(".").resolve()
T = ROOT / "results" / "tables"
OUTD = ROOT / "figures" / "main"
OUTD.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.7, "axes.edgecolor": "#333333",
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 150,
})

SPECS = ["spec_A", "spec_B", "spec_C", "spec_D"]
SPEC_LBL = {"spec_A": "A\nrecords", "spec_B": "B\nvisits",
            "spec_C": "C\nPCA", "spec_D": "D\nbirding-days"}
DSETS = ["v2", "lag", "v3"]
DS_LBL = {"v2": "v2 conservative (concurrent)",
          "lag": "v2 lagged effort (t–1)",
          "v3": "v3 relaxed (concurrent)"}
COL = {"v2": "#0072B2", "lag": "#E69F00", "v3": "#009E73"}  # Okabe-Ito

def interaction_betas():
    """return dict[(spec,dset)] = (beta, se)."""
    out = {}
    v2 = pd.read_csv(T/"table_province_v2_coefs.csv")
    s = v2[(v2.model == "M4") & (v2.term == "climate_z:effort_z")]
    for _, r in s.iterrows():
        out[(r.spec_id, "v2")] = (r.beta, r.se)
    lag = pd.read_csv(T/"table_effort_lag_refit.csv")
    for _, r in lag.iterrows():
        out[(r.spec_id, "lag")] = (r.beta, r.se)
    v3 = pd.read_csv(T/"table_province_v3_all_specs_coefs.csv")
    s = v3[(v3.model == "M4") & (v3.term == "climate_z:effort_z")]
    for _, r in s.iterrows():
        out[(r.spec_id, "v3")] = (r.beta, r.se)
    return out

B = interaction_betas()
NDRAW = 4000
draws = {k: np.exp(rng.normal(b, se, NDRAW)) for k, (b, se) in B.items()}

DODGE = {"v2": -0.27, "lag": 0.0, "v3": 0.27}
XBASE = {sp: i for i, sp in enumerate(SPECS)}

def base_axes(ax, title):
    ax.axhline(1.0, ls="--", lw=0.8, color="#888888", zorder=1)
    ax.set_xticks(list(XBASE.values()))
    ax.set_xticklabels([SPEC_LBL[s] for s in SPECS])
    ax.set_ylabel("Climate × effort interaction hazard ratio")
    ax.set_yscale("log")
    ax.set_yticks([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    ax.set_yticklabels(["0.9", "1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6"])
    ax.set_ylim(0.85, 1.75)
    ax.set_xlim(-0.6, 3.6)
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=8)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    leg = [Patch(facecolor=COL[d], edgecolor="none", label=DS_LBL[d]) for d in DSETS]
    ax.legend(handles=leg, frameon=False, fontsize=7.5, loc="upper right",
              ncol=1, handlelength=1.2)

# ---------------- Figure R1 — RAINCLOUD ----------------
fig, ax = plt.subplots(figsize=(8.6, 5.4))
base_axes(ax, "Climate × effort interaction — raincloud across effort specs and risk sets")
for sp in SPECS:
    for d in DSETS:
        x0 = XBASE[sp] + DODGE[d]
        v = draws[(sp, d)]
        col = COL[d]
        # half-violin (KDE) to the RIGHT of x0
        kde = gaussian_kde(v)
        ys = np.linspace(v.min(), v.max(), 200)
        dens = kde(ys); dens = dens / dens.max() * 0.11
        ax.fill_betweenx(ys, x0, x0 + dens, color=col, alpha=0.45, lw=0)
        ax.plot(x0 + dens, ys, color=col, lw=0.8, alpha=0.9)
        # rain (jitter) to the LEFT
        jit = rng.uniform(-0.10, -0.01, size=400)
        idx = rng.choice(len(v), 400, replace=False)
        ax.scatter(x0 + jit, v[idx], s=2.2, color=col, alpha=0.25, lw=0, zorder=2)
        # median + 95% interval
        med = np.median(v); lo, hi = np.percentile(v, [2.5, 97.5])
        ax.plot([x0, x0], [lo, hi], color="#222222", lw=1.2, zorder=4)
        ax.scatter([x0], [med], s=22, color="white", edgecolor=col,
                   linewidth=1.6, zorder=5)
ax.text(0.005, -0.14, "Dashed line = no effect (HR 1). Point = median; "
        "black line = 95% interval; cloud = bootstrap density (4000 draws ~ N(β, SE)).",
        transform=ax.transAxes, fontsize=6.8, color="#555555")
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUTD/f"Figure_R1_interaction_raincloud.{ext}",
                dpi=600 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("[56] wrote Figure_R1_interaction_raincloud.{png,pdf,svg}")

# ---------------- Figure R2 — BEESWARM ----------------
def beeswarm_x(values, x0, width=0.11, nbins=34):
    """simple beeswarm: bin by y, spread within bin symmetrically."""
    order = np.argsort(values)
    yv = values[order]
    xs = np.zeros_like(yv)
    bins = np.linspace(yv.min(), yv.max() + 1e-9, nbins + 1)
    which = np.digitize(yv, bins)
    for b in np.unique(which):
        m = which == b
        n = m.sum()
        offs = (np.arange(n) - (n - 1) / 2.0)
        if n > 1:
            offs = offs / (np.abs(offs).max()) * width
        xs[m] = x0 + offs
    out = np.zeros_like(xs); out[order] = xs
    return out

fig, ax = plt.subplots(figsize=(8.6, 5.4))
base_axes(ax, "Climate × effort interaction — beeswarm across effort specs and risk sets")
for sp in SPECS:
    for d in DSETS:
        x0 = XBASE[sp] + DODGE[d]
        v = draws[(sp, d)]
        idx = rng.choice(len(v), 600, replace=False)
        vv = v[idx]
        xs = beeswarm_x(vv, x0, width=0.11)
        ax.scatter(xs, vv, s=3.0, color=COL[d], alpha=0.35, lw=0, zorder=2)
        med = np.median(v); lo, hi = np.percentile(v, [2.5, 97.5])
        ax.plot([x0 - 0.13, x0 + 0.13], [med, med], color="#222222", lw=1.6, zorder=5)
        ax.plot([x0, x0], [lo, hi], color="#222222", lw=1.0, zorder=4)
ax.text(0.005, -0.14, "Each swarm = 600 bootstrap draws of the interaction HR; "
        "black bar = median; vertical line = 95% interval. All swarms sit above HR 1.",
        transform=ax.transAxes, fontsize=6.8, color="#555555")
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUTD/f"Figure_R2_interaction_beeswarm.{ext}",
                dpi=600 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("[56] wrote Figure_R2_interaction_beeswarm.{png,pdf,svg}")
print("[56] DONE")
