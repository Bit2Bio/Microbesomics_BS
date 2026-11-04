"""
figure1a.py

Ridge plot — PRE vs POST distributions for:
  Anthropometrics, Vitals, Body composition, Liver enzymes

Standalone panel to be combined later with figure1b.py and figure1c.py.

Output:
  results/figures/figure1a.png
"""

import os
import yaml
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR  = "/home/lorenzo/Microbesomics"
CONFIG    = os.path.join(BASE_DIR, "config/config.yaml")
DATA      = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
STATS_CSV = os.path.join(BASE_DIR, "results/tables/descriptive_baseline_stats.csv")
OUT       = os.path.join(BASE_DIR, "results/figures/figure1a.png")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

COL_PRE  = "#4A0080"
COL_POST = "#007A87"

# Variables NOT in this panel (glycemia/lipids go to figure1b)
BLOCK_B = {
    "glucose 0'", "hba1c", "insulin 0'", "homa",
    "matsuda", "quicki", "auc glucose ogtt", "auc insulin ogtt",
    "totale cholesterol", "hdl", "ldl", "triglycerides",
}
EXCLUDE = {"hip", "MetS_score", "waist", "kg fat mass"}

LABELS = {
    "weight": "Weight", "bmi": "BMI",
    "sbp": "SBP", "dbp": "DBP", "pulse": "Pulse",
    "% fat mass": "Fat mass (%)",
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
# Load pre-computed p-adj from descriptive_baseline.csv (produced by 06a)
# ---------------------------------------------------------------------------
stats_df  = pd.read_csv(STATS_CSV)
padj_dict = dict(zip(stats_df["Variable"], stats_df["p_adj_num"]))

# Collect block A variables in config order
var_list = []
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
        var_list.append((section, var))

# keep only block A vars, preserving config order
panel_vars = [(s, v) for s, v in var_list if v not in BLOCK_B]
n_vars = len(panel_vars)


def sig_stars(p):
    if np.isnan(p): return ""
    if p < 0.001:   return "***"
    if p < 0.01:    return "**"
    if p < 0.05:    return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
spacing = 1.0
x_lim   = 3.5
x_grid  = np.linspace(-x_lim - 0.2, x_lim + 0.2, 500)
kde_h   = 0.50
bw      = 0.35

fig, ax = plt.subplots(figsize=(1.3, n_vars * 0.105 + 0.28))
ax.set_facecolor("white")
fig.patch.set_facecolor("white")

y_ticks        = []
y_labels_list  = []
section_breaks = []
prev_section   = None

for i, (section, var) in enumerate(panel_vars):
    y = (n_vars - 1 - i) * spacing

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
        ax.plot(x_grid, y + dens, color=color, linewidth=0.3, alpha=0.6)

    for z_vals, color in [(z_pre, COL_PRE), (z_post, COL_POST)]:
        med = np.median(z_vals[np.isfinite(z_vals)])
        ax.vlines(med, y, y + kde_h, color=color, linewidth=0.5,
                  linestyle="--", alpha=0.7)

    stars = sig_stars(padj_dict[var])
    star_color = "#9B2226" if stars not in ("ns", "") else "#aaa"
    ax.text(x_lim + 0.06, y + kde_h * 0.45, stars,
            va="center", ha="left", fontsize=3.5,
            color=star_color, fontweight="bold" if stars != "ns" else "normal")

    y_ticks.append(y + kde_h * 0.45)
    y_labels_list.append(LABELS.get(var, var))

    # section break lines removed — too cluttered at small size

ax.axvline(0, color="#555", linewidth=0.4, linestyle="-", alpha=0.25, zorder=0)

ax.set_xlim(-x_lim - 0.1, x_lim + 0.4)
ax.set_ylim(-0.4, (n_vars - 1) * spacing + kde_h + 0.3)
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels_list, fontsize=3.5)
ax.set_xlabel("z-score", fontsize=3.5)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.spines["bottom"].set_linewidth(0.4)
ax.spines["bottom"].set_color("#888")
ax.tick_params(axis="y", length=0, pad=1)
ax.tick_params(axis="x", labelsize=3.5, width=0.4, color="#888")

patch_pre  = mpatches.Patch(color=COL_PRE,  alpha=0.7, label="PRE")
patch_post = mpatches.Patch(color=COL_POST, alpha=0.7, label="POST")
ax.legend(handles=[patch_pre, patch_post], fontsize=2.8, loc="upper right",
          framealpha=0.9, edgecolor="#ddd", handlelength=0.6, borderpad=0.25,
          handletextpad=0.3, labelspacing=0.2)

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close()
print(f"Written: {OUT}  ({n_vars} variables)")
