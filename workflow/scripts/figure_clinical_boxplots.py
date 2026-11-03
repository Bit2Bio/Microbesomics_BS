"""
figure_clinical_boxplots.py

Paired boxplots (PRE vs POST) for key clinical variables.
Mirrors Olsson et al. (Nature Metabolism 2026) Figure 1b.

Each panel shows:
  - Thin grey lines connecting each patient PRE→POST
  - Jittered dots (PRE left, POST right)
  - Notched-box overlay (median, IQR, 1.5×IQR whiskers)
  - BH-corrected Wilcoxon p-value

Variables selected as primary metabolic endpoints from bariatric surgery:
  Weight, BMI, Waist, HOMA-IR, HbA1c, FBG, Triglycerides, HDL

Output:
  results/figures/figure_clinical_boxplots.png
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

BASE_DIR = "/home/lorenzo/Microbesomics"
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT      = os.path.join(BASE_DIR, "results/figures/figure_clinical_boxplots.png")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

COL_PRE  = "#2A5C45"   # deep forest green
COL_POST = "#8B2635"   # deep burgundy
COL_LINE = "#999999"   # patient trajectory lines

# ---------------------------------------------------------------------------
# Variables to show and their display labels / units
# ---------------------------------------------------------------------------
PANELS = [
    ("weight",        "Weight",            "kg"),
    ("bmi",           "BMI",               "kg/m²"),
    ("waist",         "Waist",             "cm"),
    ("homa",          "HOMA-IR",           ""),
    ("hba1c",         "HbA1c",             "mmol/mol"),
    ("glucose 0'",    "Fasting glucose",   "mg/dL"),
    ("triglycerides", "Triglycerides",     "mg/dL"),
    ("hdl",           "HDL",               "mg/dL"),
]

# ---------------------------------------------------------------------------
df = pd.read_csv(DATA, na_values="")
pre  = df[df["timepoint"] == "PRE"].set_index("pz")
post = df[df["timepoint"] == "POST"].set_index("pz")
shared = pre.index.intersection(post.index)
pre  = pre.loc[shared]
post = post.loc[shared]

# ---------------------------------------------------------------------------
# Wilcoxon signed-rank + BH correction across the panel variables
# ---------------------------------------------------------------------------
pvals_raw = {}
for var, label, unit in PANELS:
    both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
    if len(both) < 5:
        pvals_raw[var] = np.nan
        continue
    try:
        d = both["post"] - both["pre"]
        if (d == 0).all():
            pvals_raw[var] = np.nan
        else:
            _, p = stats.wilcoxon(both["pre"], both["post"], alternative="two-sided")
            pvals_raw[var] = p
    except Exception:
        pvals_raw[var] = np.nan

vars_ordered = [v for v, _, _ in PANELS]
praw  = np.array([pvals_raw[v] for v in vars_ordered])
valid = ~np.isnan(praw)
padj  = np.full(len(praw), np.nan)
if valid.sum() > 1:
    _, padj_v, _, _ = multipletests(praw[valid], method="fdr_bh")
    padj[valid] = padj_v
padj_dict = {v: padj[i] for i, v in enumerate(vars_ordered)}


def sig_label(p: float) -> str:
    if np.isnan(p): return ""
    if p < 0.001:   return "***"
    if p < 0.01:    return "**"
    if p < 0.05:    return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
N_COLS = 4
N_ROWS = int(np.ceil(len(PANELS) / N_COLS))
rng = np.random.default_rng(42)   # reproducible jitter

fig, axes = plt.subplots(N_ROWS, N_COLS,
                          figsize=(N_COLS * 2.2, N_ROWS * 2.8),
                          constrained_layout=True)
axes = np.atleast_2d(axes).flatten()
fig.patch.set_facecolor("white")

X_PRE  = 0.0
X_POST = 1.0
JIT    = 0.07   # horizontal jitter amplitude

for idx, (var, label, unit) in enumerate(PANELS):
    ax = axes[idx]
    ax.set_facecolor("white")

    both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
    if both.empty:
        ax.axis("off")
        continue

    pre_vals  = both["pre"].values
    post_vals = both["post"].values
    n         = len(both)
    jit       = rng.uniform(-JIT, JIT, n)

    # patient trajectory lines — drawn first so dots sit on top
    for j in range(n):
        ax.plot([X_PRE + jit[j], X_POST + jit[j]],
                [pre_vals[j], post_vals[j]],
                color=COL_LINE, linewidth=0.55, alpha=0.45, zorder=1)

    # jittered dots
    ax.scatter(X_PRE  + jit, pre_vals,  color=COL_PRE,  s=18, zorder=3,
               edgecolors="white", linewidths=0.4, alpha=0.9)
    ax.scatter(X_POST + jit, post_vals, color=COL_POST, s=18, zorder=3,
               edgecolors="white", linewidths=0.4, alpha=0.9)

    # box (IQR + whiskers) — plotted on top without outlier markers because we
    # already show all individual points
    for vals, xc, col in [(pre_vals, X_PRE, COL_PRE), (post_vals, X_POST, COL_POST)]:
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        lo  = max(vals.min(), q1 - 1.5 * iqr)
        hi  = min(vals.max(), q3 + 1.5 * iqr)
        box_w = 0.18
        rect = mpatches.FancyBboxPatch(
            (xc - box_w / 2, q1), box_w, iqr,
            boxstyle="square,pad=0",
            linewidth=1.2, edgecolor=col, facecolor="none", zorder=4
        )
        ax.add_patch(rect)
        # median line
        ax.plot([xc - box_w / 2, xc + box_w / 2], [med, med],
                color=col, linewidth=2.0, zorder=5)
        # whiskers
        ax.plot([xc, xc], [lo, q1], color=col, linewidth=1.0, zorder=4)
        ax.plot([xc, xc], [q3, hi], color=col, linewidth=1.0, zorder=4)

    # significance annotation
    p_adj = padj_dict[var]
    stars = sig_label(p_adj)
    y_max = max(pre_vals.max(), post_vals.max())
    y_rng = y_max - min(pre_vals.min(), post_vals.min())
    y_sig = y_max + y_rng * 0.08
    star_col = "#9B2226" if stars not in ("ns", "") else "#aaa"
    ax.plot([X_PRE, X_POST], [y_sig, y_sig], color="#777", linewidth=0.8)
    ax.text((X_PRE + X_POST) / 2, y_sig + y_rng * 0.01, stars,
            ha="center", va="bottom", fontsize=9 if stars != "ns" else 7,
            color=star_col, fontweight="bold" if stars != "ns" else "normal")

    # axes styling
    ax.set_xticks([X_PRE, X_POST])
    ax.set_xticklabels(["PRE", "POST"], fontsize=8)
    ax.set_xlim(-0.4, 1.4)
    title = f"{label}" + (f"\n({unit})" if unit else "")
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=4)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.8)

# hide unused panels
for idx in range(len(PANELS), len(axes)):
    axes[idx].axis("off")

# shared legend
handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_PRE,
           markersize=7, label="PRE"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_POST,
           markersize=7, label="POST"),
    Line2D([0], [0], color=COL_LINE, linewidth=1, alpha=0.7,
           label="Individual patient"),
]
fig.legend(handles=handles, fontsize=7.5, loc="upper right",
           bbox_to_anchor=(1.0, 1.0), framealpha=0.9, edgecolor="#ccc",
           ncol=1)

fig.suptitle("Key clinical changes — PRE vs POST bariatric surgery (n=33)",
             fontsize=9, fontweight="bold")

plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
