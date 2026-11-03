"""
figure1_combined.py

Combined Figure 1 — mirrors Olsson et al. (Nature Metabolism 2026) Figure 1.

Portrait layout:
  Left column  — a (top): ridge plot, anthropometrics / vitals / lipids / liver
                 b (bottom): ridge plot, glycaemia + OGTT
  Right column — c (top 1/3): paired boxplot for weight PRE vs POST
                 (bottom 2/3): empty

Output:
  results/figures/figure1_combined.png
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
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = os.path.join(BASE_DIR, "config/config.yaml")
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT      = os.path.join(BASE_DIR, "results/figures/figure1_combined.png")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

COL_PRE  = "#2A5C45"
COL_POST = "#8B2635"
COL_LINE = "#aaaaaa"

# Block B = glucose homeostasis + lipids (12 vars); Block A = body/vitals/liver (11 vars)
GLYCEMIA_VARS = {
    "glucose 0'", "hba1c", "insulin 0'", "homa",
    "matsuda", "quicki", "auc glucose ogtt", "auc insulin ogtt",
    "totale cholesterol", "hdl", "ldl", "triglycerides",
}
EXCLUDE = {"hip", "MetS_score", "waist", "kg fat mass"}

LABELS = {
    "weight": "Weight", "bmi": "BMI",
    "sbp": "SBP", "dbp": "DBP", "pulse": "Pulse",
    "% fat mass": "Fat mass (%)",
    "glucose 0'": "FBG", "hba1c": "HbA1c", "insulin 0'": "Fasting insulin",
    "homa": "HOMA-IR", "matsuda": "Matsuda", "quicki": "QUICKI",
    "auc glucose ogtt": "AUC glucose", "auc insulin ogtt": "AUC insulin",
    "totale cholesterol": "Tot. cholesterol", "hdl": "HDL", "ldl": "LDL",
    "triglycerides": "Triglycerides",
    "got": "AST", "gpt": "ALT", "ggt": "GGT", "fib4": "FIB-4",
    "steatosis (saf)": "Steatosis (SAF)",
}

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
# Wilcoxon + BH
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


other_list    = [(s, v) for s, v in var_list if v not in GLYCEMIA_VARS]
glycemia_list = [(s, v) for s, v in var_list if v in GLYCEMIA_VARS]
n_other    = len(other_list)
n_glycemia = len(glycemia_list)

# ---------------------------------------------------------------------------
# Figure — portrait, two stacked ridge blocks on left, boxplot upper-right
# ---------------------------------------------------------------------------
kde_h  = 0.40
x_lim  = 2.0
x_grid = np.linspace(-x_lim - 0.3, x_lim + 0.3, 500)
bw     = 0.35

fig_h = n_other * 0.13 + 0.55   # short
fig_w = 7.5                       # wide — full page

fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

gs = gridspec.GridSpec(
    1, 3,
    width_ratios=[2.5, 2.0, 0.8],
    wspace=0.12,
    left=0.08, right=0.98, top=0.93, bottom=0.14,
)

ax_a = fig.add_subplot(gs[0])   # ridge — other
ax_b = fig.add_subplot(gs[1])   # ridge — glycemia
ax_c = fig.add_subplot(gs[2])   # weight boxplot

for ax in (ax_a, ax_b, ax_c):
    ax.set_facecolor("white")


# ---------------------------------------------------------------------------
# Helper: draw ridge block
# ---------------------------------------------------------------------------
def draw_ridge(ax, vlist, panel_label, n_ref=None):
    n       = len(vlist)
    # stretch spacing so shorter block fills same y-span as the reference block
    n_ref   = n_ref or n
    spacing = (n_ref - 1) / (n - 1) if n > 1 else 1.0
    y_ticks        = []
    y_labels_list  = []
    section_breaks = []
    prev_section   = None

    for i, (section, var) in enumerate(vlist):
        y = (n - 1 - i) * spacing

        if section != prev_section:
            if prev_section is not None:
                section_breaks.append(y + spacing * 0.5)
            prev_section = section

        both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
        mu, sd = both["pre"].mean(), both["pre"].std()
        if sd == 0:
            continue

        z_pre  = (both["pre"]  - mu) / sd
        z_post = (both["post"] - mu) / sd

        for z_vals, color, alpha in [(z_pre, COL_PRE, 0.55), (z_post, COL_POST, 0.55)]:
            vals = z_vals[np.isfinite(z_vals)]
            if len(vals) < 3 or vals.std() == 0:
                continue
            kde_obj = gaussian_kde(vals, bw_method=bw)
            dens    = kde_obj(x_grid)
            data_sd = vals.std()
            mask = (x_grid < vals.min() - 0.5*data_sd) | (x_grid > vals.max() + 0.5*data_sd)
            dens[mask] = 0
            dens = dens / dens.max() * kde_h if dens.max() > 0 else dens
            ax.fill_between(x_grid, y, y + dens, alpha=alpha, color=color, linewidth=0)
            ax.plot(x_grid, y + dens, color=color, linewidth=0.7, alpha=0.9)

        for z_vals, color in [(z_pre, COL_PRE), (z_post, COL_POST)]:
            med = np.median(z_vals[np.isfinite(z_vals)])
            ax.vlines(med, y, y + kde_h, color=color, linewidth=1.2,
                      linestyle="--", alpha=0.9)

        stars = sig_stars(padj_dict[var])
        star_color = "#9B2226" if stars not in ("ns", "") else "#aaa"
        ax.text(x_lim + 0.06, y + kde_h * 0.45, stars,
                va="center", ha="left", fontsize=4.5,
                color=star_color, fontweight="bold" if stars != "ns" else "normal")

        y_ticks.append(y + kde_h * 0.45)
        y_labels_list.append(LABELS.get(var, var))

    for yb in section_breaks:
        ax.axhline(yb, color="#ddd", linewidth=0.5, linestyle="-")

    ax.axvline(0, color="#555", linewidth=0.7, linestyle="-", alpha=0.3, zorder=0)
    ax.set_xlim(-x_lim - 0.1, x_lim + 0.3)
    ax.set_ylim(-0.4, (n - 1) * spacing + kde_h + 0.3)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels_list, fontsize=5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=5)
    ax.text(-0.16, 1.02, panel_label, transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="bottom")


draw_ridge(ax_a, other_list,    "a", n_ref=n_other)
draw_ridge(ax_b, glycemia_list, "b", n_ref=n_other)

ax_a.set_xlabel("PRE SD units", fontsize=6.5)
ax_b.set_xlabel("PRE SD units", fontsize=6.5)

# legend on ax_a
patch_pre  = mpatches.Patch(color=COL_PRE,  alpha=0.7, label="PRE")
patch_post = mpatches.Patch(color=COL_POST, alpha=0.7, label="POST")
ax_a.legend(handles=[patch_pre, patch_post], fontsize=5, loc="lower right",
            framealpha=0.9, edgecolor="#ddd", handlelength=1.0)
ax_b.text(x_lim + 0.06, -0.35,
          "*<0.05  **<0.01  ***<0.001",
          fontsize=4, color="#888", va="bottom", ha="left")

# ---------------------------------------------------------------------------
# Panel c — weight boxplot
# ---------------------------------------------------------------------------
weight_both = pd.DataFrame({"pre": pre["weight"], "post": post["weight"]}).dropna()
pre_w  = weight_both["pre"].values
post_w = weight_both["post"].values
n_w    = len(weight_both)

rng = np.random.default_rng(42)
jit = rng.uniform(-0.07, 0.07, n_w)
X_PRE, X_POST = 0.0, 1.0

for j in range(n_w):
    ax_c.plot([X_PRE + jit[j], X_POST + jit[j]],
              [pre_w[j], post_w[j]],
              color=COL_LINE, linewidth=0.5, alpha=0.4, zorder=1)

ax_c.scatter(X_PRE  + jit, pre_w,  color=COL_PRE,  s=12, zorder=3,
             edgecolors="white", linewidths=0.3, alpha=0.9)
ax_c.scatter(X_POST + jit, post_w, color=COL_POST, s=12, zorder=3,
             edgecolors="white", linewidths=0.3, alpha=0.9)

box_w = 0.15
for vals, xc, col in [(pre_w, X_PRE, COL_PRE), (post_w, X_POST, COL_POST)]:
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    iqr = q3 - q1
    lo  = max(vals.min(), q1 - 1.5 * iqr)
    hi  = min(vals.max(), q3 + 1.5 * iqr)
    rect = mpatches.FancyBboxPatch(
        (xc - box_w/2, q1), box_w, iqr,
        boxstyle="square,pad=0",
        linewidth=1.1, edgecolor=col, facecolor="none", zorder=4
    )
    ax_c.add_patch(rect)
    ax_c.plot([xc - box_w/2, xc + box_w/2], [med, med],
              color=col, linewidth=1.8, zorder=5)
    ax_c.plot([xc, xc], [lo, q1], color=col, linewidth=0.8, zorder=4)
    ax_c.plot([xc, xc], [q3, hi], color=col, linewidth=0.8, zorder=4)

p_w     = padj_dict.get("weight", np.nan)
stars_w = sig_stars(p_w)
y_max_w = max(pre_w.max(), post_w.max())
y_min_w = min(pre_w.min(), post_w.min())
y_rng_w = y_max_w - y_min_w
y_sig_w = y_max_w + y_rng_w * 0.07
ax_c.plot([X_PRE, X_POST], [y_sig_w, y_sig_w], color="#777", linewidth=0.7)
ax_c.text(0.5, y_sig_w + y_rng_w * 0.01, stars_w,
          ha="center", va="bottom",
          fontsize=8 if stars_w != "ns" else 6,
          color="#9B2226" if stars_w not in ("ns", "") else "#aaa",
          fontweight="bold" if stars_w != "ns" else "normal")

ax_c.set_xticks([X_PRE, X_POST])
ax_c.set_xticklabels(["PRE", "POST"], fontsize=5)
ax_c.set_xlim(-0.4, 1.4)
ax_c.set_ylabel("Weight (kg)", fontsize=5)
ax_c.set_ylim(y_min_w - y_rng_w * 0.05, y_sig_w + y_rng_w * 0.18)
ax_c.tick_params(axis="y", labelsize=5)
ax_c.spines[["top", "right"]].set_visible(False)
ax_c.spines[["left", "bottom"]].set_linewidth(0.6)
ax_c.text(-0.25, 1.02, "c", transform=ax_c.transAxes,
          fontsize=8, fontweight="bold", va="bottom")
patch_pre  = mpatches.Patch(color=COL_PRE,  alpha=0.7, label="PRE")
patch_post = mpatches.Patch(color=COL_POST, alpha=0.7, label="POST")
ax_c.legend(handles=[patch_pre, patch_post], fontsize=4.5, loc="lower left",
            framealpha=0.9, edgecolor="#ddd", handlelength=0.8)

# ---------------------------------------------------------------------------
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
