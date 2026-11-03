"""
figure_lfc_clinical.py

Ridge plot: PRE vs POST distributions for continuous clinical variables.
Mirrors Olsson et al. (Nature Metabolism 2026) Figure 1a.

For each variable:
  - x-axis: value standardised to PRE mean/SD (so PRE ≈ N(0,1), POST shifted)
  - PRE distribution in blue, POST in coral
  - Wilcoxon signed-rank BH-adjusted p-value as significance stars

Output:
  results/figures/figure_lfc_clinical.png
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

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = os.path.join(BASE_DIR, "config/config.yaml")
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT      = os.path.join(BASE_DIR, "results/figures/figure_lfc_clinical.png")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

COL_PRE  = "#2A5C45"   # deep forest green
COL_POST = "#8B2635"   # deep burgundy

# ---------------------------------------------------------------------------
# Load config and data
# ---------------------------------------------------------------------------
with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

sections = cfg["clinical_descriptive"]

df = pd.read_csv(DATA, na_values="")
if "height" in df.columns:
    df["height"] = pd.to_numeric(df["height"], errors="coerce")

pre  = df[df["timepoint"] == "PRE"].set_index("pz")
post = df[df["timepoint"] == "POST"].set_index("pz")
shared = pre.index.intersection(post.index)
pre  = pre.loc[shared]
post = post.loc[shared]

# Variables excluded from this figure (redundant or non-independent)
EXCLUDE = {"hip", "MetS_score"}

LABELS = {
    "weight": "Weight", "bmi": "BMI", "waist": "Waist",
    "sbp": "SBP", "dbp": "DBP", "pulse": "Pulse",
    "% fat mass": "Fat mass (%)", "kg fat mass": "Fat mass (kg)",
    "glucose 0'": "FBG", "hba1c": "HbA1c", "insulin 0'": "Fasting insulin",
    "homa": "HOMA-IR", "matsuda": "Matsuda", "quicki": "QUICKI",
    "auc glucose ogtt": "AUC gluc (OGTT)", "auc insulin ogtt": "AUC ins (OGTT)",
    "totale cholesterol": "Tot. cholesterol", "hdl": "HDL", "ldl": "LDL",
    "triglycerides": "Triglycerides",
    "got": "AST", "gpt": "ALT", "ggt": "GGT", "fib4": "FIB-4",
    "steatosis (saf)": "Steatosis (SAF)", "MetS_score": "MetS score",
}

# ---------------------------------------------------------------------------
# Collect variables and compute Wilcoxon p-values
# ---------------------------------------------------------------------------
var_list = []
pval_raw = {}

for section, variables in sections.items():
    for var, meta in variables.items():
        if isinstance(meta, dict):
            vtype   = meta["type"]
            compare = bool(meta.get("compare", True))
        else:
            vtype, compare = meta, True

        if vtype != "continuous" or not compare:
            continue
        if var in EXCLUDE:
            continue
        if var not in pre.columns:
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

# BH correction
vars_ordered = [v for _, v in var_list]
pvals  = np.array([pval_raw[v] for v in vars_ordered])
valid  = ~np.isnan(pvals)
padj   = np.full(len(pvals), np.nan)
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
# Plot
# ---------------------------------------------------------------------------
n_vars  = len(var_list)
spacing = 1.0
x_lim   = 4.5     # axis in PRE SD units; wider to reduce tail clipping
x_grid  = np.linspace(-x_lim, x_lim, 500)
kde_h   = 0.65    # max KDE height per row

fig, ax = plt.subplots(figsize=(5.0, n_vars * 0.33 + 1.2))
ax.set_facecolor("white")
fig.patch.set_facecolor("white")

section_breaks = []   # y positions where section changes
prev_section   = None

y_ticks  = []
y_labels = []

for i, (section, var) in enumerate(var_list):
    y = (n_vars - 1 - i) * spacing

    if section != prev_section:
        if prev_section is not None:
            section_breaks.append(y + spacing * 0.5)
        prev_section = section

    both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
    mu  = both["pre"].mean()
    sd  = both["pre"].std()
    if sd == 0:
        continue

    z_pre  = (both["pre"]  - mu) / sd
    z_post = (both["post"] - mu) / sd

    for z_vals, color, alpha in [(z_pre, COL_PRE, 0.55), (z_post, COL_POST, 0.55)]:
        vals = z_vals[np.isfinite(z_vals)]
        if len(vals) < 3 or vals.std() == 0:
            continue
        # narrow bandwidth to avoid inflated tails with small n
        bw = 0.35
        kde  = gaussian_kde(vals, bw_method=bw)
        dens = kde(x_grid)
        # zero out density outside [min-0.5sd, max+0.5sd] of actual data
        data_sd = vals.std()
        mask = (x_grid < vals.min() - 0.5 * data_sd) | (x_grid > vals.max() + 0.5 * data_sd)
        dens[mask] = 0
        dens = dens / dens.max() * kde_h if dens.max() > 0 else dens
        ax.fill_between(x_grid, y, y + dens, alpha=alpha, color=color, linewidth=0)
        ax.plot(x_grid, y + dens, color=color, linewidth=0.7, alpha=0.9)

    # median markers
    for z_vals, color in [(z_pre, COL_PRE), (z_post, COL_POST)]:
        med = np.median(z_vals[np.isfinite(z_vals)])
        ax.vlines(med, y, y + kde_h, color=color, linewidth=1.4,
                  linestyle="--", alpha=0.95)

    # significance
    stars = sig_stars(padj_dict[var])
    star_color = "#c0392b" if stars not in ("ns", "") else "#999"
    ax.text(x_lim + 0.1, y + kde_h * 0.45, stars,
            va="center", ha="left", fontsize=7,
            color=star_color, fontweight="bold" if stars != "ns" else "normal")

    y_ticks.append(y + kde_h * 0.45)
    y_labels.append(LABELS.get(var, var))

# section separator lines
for yb in section_breaks:
    ax.axhline(yb, color="#ccc", linewidth=0.6, linestyle="-")

# x=0 reference
ax.axvline(0, color="#333", linewidth=0.8, linestyle="-", alpha=0.4, zorder=0)

# axes formatting
ax.set_xlim(-x_lim, x_lim + 0.7)
ax.set_ylim(-0.5, n_vars * spacing + 0.3)
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontsize=7)
ax.set_xlabel("Standardised value (PRE SD units)", fontsize=7.5)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.tick_params(axis="y", length=0, pad=3)
ax.tick_params(axis="x", labelsize=7)

# legend
from matplotlib.patches import Patch
handles = [Patch(color=COL_PRE, alpha=0.7, label="PRE"),
           Patch(color=COL_POST, alpha=0.7, label="POST")]
ax.legend(handles=handles, fontsize=7, loc="lower right",
          framealpha=0.9, edgecolor="#ccc")

ax.text(x_lim + 0.1, -0.3,
        "* p-adj<0.05\n** p-adj<0.01\n*** p-adj<0.001",
        fontsize=6, color="#777", va="bottom", ha="left")

ax.set_title("Clinical changes — PRE vs POST bariatric surgery (n=33)",
             fontsize=8, pad=6, loc="left")

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}")
