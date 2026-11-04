"""
figure1c.py

Raincloud plot — Weight PRE vs POST (Allen et al. 2019).

Three layers per condition:
  - Half-violin (KDE): PRE opens leftward, POST opens rightward
  - Jittered raw points with patient connecting lines
  - Thin boxplot overlay (IQR + whiskers + median)

Physical size matches figure1a/1b so all three panels align when assembled.
BH correction uses the same all-variable denominator as figure1a/1b.

Output:
  results/figures/figure1c.png
"""

import os
import yaml
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = os.path.join(BASE_DIR, "config/config.yaml")
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT      = os.path.join(BASE_DIR, "results/figures/figure1c.png")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

COL_PRE  = "#4A0080"
COL_POST = "#007A87"
COL_LINE = "#aaaaaa"

EXCLUDE = {"hip", "MetS_score", "waist", "kg fat mass"}

# ---------------------------------------------------------------------------
with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

df = pd.read_csv(DATA, na_values="")
pre  = df[df["timepoint"] == "PRE"].set_index("pz")
post = df[df["timepoint"] == "POST"].set_index("pz")
shared = pre.index.intersection(post.index)
pre  = pre.loc[shared]
post = post.loc[shared]

# ---------------------------------------------------------------------------
# BH correction over all continuous vars — same denominator as figure1a/1b
# ---------------------------------------------------------------------------
var_list = []
pval_raw = {}

for section, variables in cfg["clinical_descriptive"].items():
    for var, meta in variables.items():
        if isinstance(meta, dict):
            vtype   = meta["type"]
            compare = bool(meta.get("compare", True))
        else:
            vtype, compare = meta, True
        if vtype != "continuous" or not compare:
            continue
        if var in EXCLUDE or var not in pre.columns:
            continue
        both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
        if len(both) < 5:
            continue
        try:
            d = both["post"] - both["pre"]
            _, p = stats.wilcoxon(both["pre"], both["post"], alternative="two-sided") \
                   if not (d == 0).all() else (None, np.nan)
        except Exception:
            p = np.nan
        var_list.append((section, var))
        pval_raw[var] = p

vars_ordered = [v for _, v in var_list]
pvals = np.array([pval_raw[v] for v in vars_ordered])
valid = ~np.isnan(pvals)
padj  = np.full(len(pvals), np.nan)
if valid.sum() > 1:
    _, padj_v, _, _ = multipletests(pvals[valid], method="fdr_bh")
    padj[valid] = padj_v
padj_dict = {v: padj[i] for i, v in enumerate(vars_ordered)}


def sig_stars(p):
    if np.isnan(p): return ""
    if p < 0.001:   return "***"
    if p < 0.01:    return "**"
    if p < 0.05:    return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
weight_both = pd.DataFrame({"pre": pre["weight"], "post": post["weight"]}).dropna()
pre_w  = weight_both["pre"].values
post_w = weight_both["post"].values
n_w    = len(weight_both)

# ---------------------------------------------------------------------------
# Plot — same physical size as figure1a
# ---------------------------------------------------------------------------
N_A_VARS = 11
fig_h    = N_A_VARS * 0.105 + 0.28

fig, ax = plt.subplots(figsize=(1.3, fig_h))
ax.set_facecolor("white")
fig.patch.set_facecolor("white")

X_PRE, X_POST = 0.0, 1.0
VIO_W = 0.30    # max half-violin width in x-axis units
bw    = 0.35

# KDE y-grid covering both distributions
y_lo    = min(pre_w.min(), post_w.min())
y_hi    = max(pre_w.max(), post_w.max())
y_span  = y_hi - y_lo
y_grid  = np.linspace(y_lo - y_span * 0.08, y_hi + y_span * 0.08, 400)

# Half-violins: PRE opens LEFT, POST opens RIGHT
for vals, xc, col, side in [
    (pre_w,  X_PRE,  COL_PRE,  "left"),
    (post_w, X_POST, COL_POST, "right"),
]:
    kde_obj = gaussian_kde(vals, bw_method=bw)
    dens    = kde_obj(y_grid)
    data_sd = vals.std()
    # clip tails beyond data ± 0.5 SD (same clipping as ridge plots)
    mask = (y_grid < vals.min() - 0.5*data_sd) | (y_grid > vals.max() + 0.5*data_sd)
    dens[mask] = 0
    dens = dens / dens.max() * VIO_W if dens.max() > 0 else dens

    if side == "left":
        ax.fill_betweenx(y_grid, xc, xc - dens, alpha=0.55, color=col, linewidth=0)
        ax.plot(xc - dens, y_grid, color=col, linewidth=0.3, alpha=0.6)
    else:
        ax.fill_betweenx(y_grid, xc, xc + dens, alpha=0.55, color=col, linewidth=0)
        ax.plot(xc + dens, y_grid, color=col, linewidth=0.3, alpha=0.6)

# Patient connecting lines (drawn before points so points sit on top)
rng = np.random.default_rng(42)
jit = rng.uniform(-0.05, 0.05, n_w)
for j in range(n_w):
    ax.plot([X_PRE + jit[j], X_POST + jit[j]],
            [pre_w[j], post_w[j]],
            color=COL_LINE, linewidth=0.3, alpha=0.35, zorder=2)

# Jittered raw points
ax.scatter(X_PRE  + jit, pre_w,  color=COL_PRE,  s=3, zorder=3, edgecolors="none", alpha=0.85)
ax.scatter(X_POST + jit, post_w, color=COL_POST, s=3, zorder=3, edgecolors="none", alpha=0.85)

# Thin boxplot overlay
box_w = 0.06
for vals, xc, col in [(pre_w, X_PRE, COL_PRE), (post_w, X_POST, COL_POST)]:
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    iqr = q3 - q1
    lo  = max(vals.min(), q1 - 1.5 * iqr)
    hi  = min(vals.max(), q3 + 1.5 * iqr)
    rect = mpatches.FancyBboxPatch(
        (xc - box_w/2, q1), box_w, iqr,
        boxstyle="square,pad=0",
        linewidth=0.5, edgecolor=col, facecolor="white", zorder=4, alpha=0.9
    )
    ax.add_patch(rect)
    ax.plot([xc - box_w/2, xc + box_w/2], [med, med],
            color=col, linewidth=0.8, zorder=5)
    ax.plot([xc, xc], [lo, q1], color=col, linewidth=0.4, zorder=4)
    ax.plot([xc, xc], [q3, hi], color=col, linewidth=0.4, zorder=4)

# Significance bracket
p_w     = padj_dict.get("weight", np.nan)
stars_w = sig_stars(p_w)
y_sig   = y_hi + y_span * 0.09
ax.plot([X_PRE, X_POST], [y_sig, y_sig], color="#777", linewidth=0.4)
ax.text(0.5, y_sig + y_span * 0.015, stars_w,
        ha="center", va="bottom",
        fontsize=4.5 if stars_w != "ns" else 3.5,
        color="#9B2226" if stars_w not in ("ns", "") else "#aaa",
        fontweight="bold" if stars_w != "ns" else "normal")

ax.set_xticks([X_PRE, X_POST])
ax.set_xticklabels(["PRE", "POST"], fontsize=3.5)
ax.set_xlim(-VIO_W - 0.05, 1.0 + VIO_W + 0.05)
ax.set_ylabel("Weight (kg)", fontsize=3.5)
ax.set_ylim(y_lo - y_span * 0.05, y_sig + y_span * 0.15)
ax.tick_params(axis="y", labelsize=3.5, width=0.4, color="#888", pad=1)
ax.tick_params(axis="x", width=0.4, color="#888")
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_linewidth(0.4)
ax.spines["left"].set_color("#888")
ax.spines["bottom"].set_color("#888")

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
