"""
figure1a_palette_test.py

Generates 7 colour-palette variants of figure1a for visual comparison.
Run manually — not in Snakefile. Output: results/figures/palette_test/
"""

import os, yaml
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
OUT_DIR  = os.path.join(BASE_DIR, "results/figures/palette_test")
os.makedirs(OUT_DIR, exist_ok=True)

PALETTES = [
    # batch 4 — genuinely different hue families
    ("28_black_coral",       "#1A1A1A", "#E05A40"),   # near-black + vivid coral
    ("29_ultramarine_amber", "#1B2CC1", "#E8920A"),   # ultramarine + amber
    ("30_teal_crimson",      "#006D77", "#C1121F"),   # medical teal + crimson
    ("31_violet_lime",       "#5A189A", "#80B918"),   # deep violet + lime
    ("32_navy_scarlet",      "#14213D", "#E63946"),   # dark navy + pure scarlet
    ("33_purple_teal",       "#4A0080", "#007A87"),   # purple + bright teal
    ("34_black_blue",        "#1A1A2E", "#4CC9F0"),   # near-black + electric sky
    ("35_indigo_orange",     "#3D405B", "#F2A65A"),   # muted indigo + warm orange
    ("36_charcoal_magenta",  "#2D2D2D", "#C2185B"),   # charcoal + magenta
    ("37_slate_cyan",        "#334155", "#06B6D4"),   # dark slate + cyan
]

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

with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

df   = pd.read_csv(DATA, na_values="")
pre  = df[df["timepoint"] == "PRE"].set_index("pz")
post = df[df["timepoint"] == "POST"].set_index("pz")
shared = pre.index.intersection(post.index)
pre  = pre.loc[shared]
post = post.loc[shared]

# collect vars + BH once
var_list, pval_raw = [], {}
for section, variables in cfg["clinical_descriptive"].items():
    for var, meta in variables.items():
        if isinstance(meta, dict):
            vtype, compare = meta["type"], bool(meta.get("compare", True))
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

panel_vars = [(s, v) for s, v in var_list if v not in BLOCK_B]
n_vars     = len(panel_vars)


def sig_stars(p):
    if np.isnan(p): return ""
    if p < 0.001:   return "***"
    if p < 0.01:    return "**"
    if p < 0.05:    return "*"
    return "ns"


def make_figure(col_pre, col_post, out_path):
    spacing = 1.0
    x_lim   = 3.5
    x_grid  = np.linspace(-x_lim - 0.2, x_lim + 0.2, 500)
    kde_h   = 0.50
    bw      = 0.35

    fig, ax = plt.subplots(figsize=(1.3, n_vars * 0.105 + 0.28))
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    y_ticks, y_labels_list = [], []

    for i, (section, var) in enumerate(panel_vars):
        y = (n_vars - 1 - i) * spacing

        both = pd.DataFrame({"pre": pre[var], "post": post[var]}).dropna()
        mu, sd = both["pre"].mean(), both["pre"].std()
        if sd == 0:
            continue

        z_pre  = (both["pre"]  - mu) / sd
        z_post = (both["post"] - mu) / sd

        for z_vals, color, alpha in [(z_pre, col_pre, 0.55), (z_post, col_post, 0.55)]:
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

        for z_vals, color in [(z_pre, col_pre), (z_post, col_post)]:
            med = np.median(z_vals[np.isfinite(z_vals)])
            ax.vlines(med, y, y + kde_h, color=color, linewidth=0.5,
                      linestyle="--", alpha=0.7)

        stars = sig_stars(padj_dict[var])
        star_color = "#9B2226" if stars not in ("ns", "") else "#aaa"
        ax.text(x_lim + 0.08, y + kde_h * 0.45, stars,
                va="center", ha="left", fontsize=3.5,
                color=star_color, fontweight="bold" if stars != "ns" else "normal")

        y_ticks.append(y + kde_h * 0.45)
        y_labels_list.append(LABELS.get(var, var))

    ax.axvline(0, color="#555", linewidth=0.4, linestyle="-", alpha=0.25, zorder=0)
    ax.set_xlim(-x_lim - 0.1, x_lim + 0.4)
    ax.set_ylim(-0.4, (n_vars - 1) * spacing + kde_h + 0.3)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels_list, fontsize=3.5)
    ax.set_xlabel("PRE SD units", fontsize=3.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["bottom"].set_color("#888")
    ax.tick_params(axis="y", length=0, pad=1)
    ax.tick_params(axis="x", labelsize=3.5, width=0.4, color="#888")

    patch_pre  = mpatches.Patch(color=col_pre,  alpha=0.7, label="PRE")
    patch_post = mpatches.Patch(color=col_post, alpha=0.7, label="POST")
    ax.legend(handles=[patch_pre, patch_post], fontsize=3.5, loc="lower right",
              framealpha=0.9, edgecolor="#ddd", handlelength=0.8, borderpad=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


for name, col_pre, col_post in PALETTES:
    out = os.path.join(OUT_DIR, f"{name}.png")
    make_figure(col_pre, col_post, out)
    print(f"Written: {out}")
