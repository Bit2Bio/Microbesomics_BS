#!/usr/bin/env python3
"""
clinical_continuous_diagnostics.py

QC grid: raw vs log-transformed distribution for every continuous clinical
variable at PRE timepoint. Use this to decide which variables benefit from
log-normalization before KDE ridge plotting.

Two panels per variable:
  LEFT  — raw histogram (light blue)
  RIGHT — log(x) histogram (light green)

Red dashed lines mark individual patients with |z| > OUTLIER_Z in that scale.
Skewed variables (e.g. HOMA, insulin, TG, GGT) will look more symmetric on the
log panel; if both look roughly symmetric, raw scale is fine.

Standalone QC script — not wired into the Snakefile. Run manually:
    python3 workflow/scripts/clinical_continuous_diagnostics.py

Output:
    results/figures/qc/continuous_distributions_pre.png
"""

import sys
import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = f"{BASE_DIR}/config/config.yaml"
DATA     = f"{BASE_DIR}/data/processed/bariatric_clinical_66.csv"
OUT      = Path(f"{BASE_DIR}/results/figures/qc/continuous_distributions_pre.png")
OUTLIER_Z = 3.0

OUT.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

df_all = pd.read_csv(DATA, na_values="")
pre = df_all[df_all["timepoint"] == "PRE"].set_index("pz")

# Collect continuous variables with compare: yes — these are the ones that go
# into figure_lfc_clinical.py, so they're the ones worth inspecting.
EXCLUDE = {"hip", "MetS_score"}

cont_vars = []
for section, variables in cfg["clinical_descriptive"].items():
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
            print(f"  WARNING: '{var}' not in PRE data — skipping", file=sys.stderr)
            continue
        cont_vars.append((section, var))

n_vars   = len(cont_vars)
n_cols   = 3          # variables per row (each occupies 2 panels: raw + log)
n_rows   = int(np.ceil(n_vars / n_cols))

fig, axes = plt.subplots(n_rows, n_cols * 2,
                          figsize=(n_cols * 2 * 3.0, n_rows * 2.4))
axes = np.atleast_2d(axes)


def plot_panel(ax, values: pd.Series, title: str, color: str) -> None:
    """Histogram with outlier markers."""
    valid = values.dropna()
    if len(valid) < 3:
        ax.text(0.5, 0.5, "n<3", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="#999")
        ax.set_title(title, fontsize=7)
        return

    ax.hist(valid, bins=15, color=color, edgecolor="white", linewidth=0.5, alpha=0.85)

    # mark outliers (|z| > threshold)
    z = (valid - valid.mean()) / valid.std()
    for idx, zi in z.items():
        if abs(zi) > OUTLIER_Z:
            v = valid.loc[idx]
            ax.axvline(v, color="#e74c3c", linewidth=1.2, linestyle="--", alpha=0.9)
            ax.text(v, ax.get_ylim()[1] * 0.92, str(idx),
                    rotation=90, fontsize=5.5, color="#e74c3c",
                    ha="right", va="top")

    # annotate skewness — helps quantify how lopsided the distribution is
    sk = float(valid.skew())
    ax.text(0.97, 0.95, f"skew={sk:+.2f}", transform=ax.transAxes,
            fontsize=6, ha="right", va="top", color="#555")

    ax.set_title(title, fontsize=7.5)
    ax.tick_params(labelsize=6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)


for i, (section, var) in enumerate(cont_vars):
    row   = i // n_cols
    block = i % n_cols
    ax_raw = axes[row, block * 2]
    ax_log = axes[row, block * 2 + 1]

    s = pre[var]

    # raw panel
    plot_panel(ax_raw, s, f"{var}\nraw", color="#5b9bd5")

    # log panel — only for strictly positive values; shift by tiny epsilon if needed
    s_pos = s.copy()
    min_pos = s_pos[s_pos > 0].min() if (s_pos > 0).any() else np.nan
    if np.isnan(min_pos):
        ax_log.text(0.5, 0.5, "all ≤ 0", ha="center", va="center",
                    transform=ax_log.transAxes, fontsize=8, color="#999")
        ax_log.set_title(f"{var}\nlog (N/A)", fontsize=7.5)
    else:
        # floor near-zero values at half the minimum positive observed value
        s_pos = s_pos.clip(lower=min_pos * 0.5)
        plot_panel(ax_log, np.log(s_pos), f"{var}\nlog", color="#6dbf8c")

# hide unused panels
for i in range(n_vars, n_rows * n_cols):
    row   = i // n_cols
    block = i % n_cols
    axes[row, block * 2].axis("off")
    axes[row, block * 2 + 1].axis("off")

fig.suptitle(
    f"Continuous clinical variables at PRE — raw vs log  "
    f"(red dashed = |z|>{OUTLIER_Z:.0f}; n=33)",
    fontsize=10, fontweight="bold"
)
plt.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
