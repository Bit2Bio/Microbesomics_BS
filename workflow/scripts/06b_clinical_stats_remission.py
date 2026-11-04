"""
descriptive_remission_stratified.py

Clinical characteristics at PRE stratified by remission status:
remitters (n=22) vs non_remitters (n=6). Only MetS+ subjects included (n=28).
MetS− (mets_neg, n=5) excluded — remission is undefined for them.

Tests (unpaired):
  Continuous → Mann-Whitney U (scipy.stats.mannwhitneyu, two-sided)
  Binary     → Fisher's exact test (scipy.stats.fisher_exact)
  Categorical → no p-value (cells too sparse with n=6 in non_remitter group)

BH correction applied separately per test family.

Note: non_remitters n=6 — tests are underpowered. p-values reported for
completeness but interpret with caution.

Output:
  results/tables/descriptive_remission_stratified.csv
  results/tables/descriptive_remission_stratified.html
"""

import sys
import os
import yaml
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = os.path.join(BASE_DIR, "config/config.yaml")
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT_CSV  = os.path.join(BASE_DIR, "results/tables/descriptive_remission_stratified.csv")
OUT_HTML = os.path.join(BASE_DIR, "results/tables/descriptive_remission_stratified.html")

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

# ---------------------------------------------------------------------------
# Load config and data
# ---------------------------------------------------------------------------
with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

sections = cfg["clinical_descriptive"]

df_all = pd.read_csv(DATA, na_values="")
pre = df_all[df_all["timepoint"] == "PRE"].copy()

if "height" in pre.columns:
    pre["height"] = pd.to_numeric(pre["height"], errors="coerce")

rem     = pre[pre["remission"] == "remitter"]
non_rem = pre[pre["remission"] == "non_remitter"]
n_rem   = len(rem)
n_non   = len(non_rem)
print(f"Remitters n={n_rem}, Non-remitters n={n_non}")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def coverage_str(col: pd.Series, n_total: int) -> str:
    n_obs     = col.notna().sum()
    n_missing = n_total - n_obs
    return f"{n_obs}/{n_total}" if n_missing == 0 else f"{n_obs}/{n_total} ({n_missing} miss.)"


def stat_continuous(col: pd.Series) -> str:
    q1, med, q3 = col.quantile([0.25, 0.50, 0.75])
    return f"{med:.1f} [{q1:.1f}–{q3:.1f}]"


def stat_binary(col: pd.Series) -> str:
    n_obs = col.notna().sum()
    n_pos = (col == 1).sum()
    pct   = 100 * n_pos / n_obs if n_obs > 0 else np.nan
    return f"{n_pos} ({pct:.0f}%)"


def pval_str(p: float) -> str:
    if np.isnan(p):
        return "—"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def mannwhitney(col_a: pd.Series, col_b: pd.Series) -> float:
    a = col_a.dropna()
    b = col_b.dropna()
    if len(a) < 3 or len(b) < 2:
        return np.nan
    try:
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return p
    except Exception:
        return np.nan


def fisher_exact(col_a: pd.Series, col_b: pd.Series) -> float:
    n1   = col_a.notna().sum()
    n2   = col_b.notna().sum()
    pos1 = (col_a == 1).sum()
    pos2 = (col_b == 1).sum()
    table = [[pos1, n1 - pos1],
             [pos2, n2 - pos2]]
    try:
        _, p = stats.fisher_exact(table, alternative="two-sided")
        return p
    except Exception:
        return np.nan


# ---------------------------------------------------------------------------
# Compute rows
# ---------------------------------------------------------------------------
rows = []

for section, variables in sections.items():
    for var, meta in variables.items():
        if isinstance(meta, dict):
            vtype = meta["type"]
        else:
            vtype = meta

        # skip remission itself and MetS-derived variables (used for grouping)
        if var in ("remission", "MetS", "MetS_score",
                   "HyperglicemiaMetS", "HypertensionMetS",
                   "HDL_MetS", "TrygliceridesMetS"):
            continue

        if var not in pre.columns:
            print(f"  WARNING: '{var}' not in data — skipping", file=sys.stderr)
            continue

        col_rem = rem[var]
        col_non = non_rem[var]

        if vtype == "continuous":
            rows.append({
                "Section":       section,
                "Variable":      var,
                "Remitters":     stat_continuous(col_rem),
                "n Remitters":   coverage_str(col_rem, n_rem),
                "Non-remitters": stat_continuous(col_non),
                "n Non-rem":     coverage_str(col_non, n_non),
                "_p_raw":        mannwhitney(col_rem, col_non),
                "_test":         "mannwhitney",
            })

        elif vtype == "binary":
            rows.append({
                "Section":       section,
                "Variable":      var,
                "Remitters":     stat_binary(col_rem),
                "n Remitters":   coverage_str(col_rem, n_rem),
                "Non-remitters": stat_binary(col_non),
                "n Non-rem":     coverage_str(col_non, n_non),
                "_p_raw":        fisher_exact(col_rem, col_non),
                "_test":         "fisher",
            })

        elif vtype == "categorical":
            levels    = sorted(pre[var].dropna().unique())
            n_obs_rem = col_rem.notna().sum()
            n_obs_non = col_non.notna().sum()
            for i, lvl in enumerate(levels):
                n_r   = (col_rem == lvl).sum()
                n_n   = (col_non == lvl).sum()
                pct_r = 100 * n_r / n_obs_rem if n_obs_rem > 0 else np.nan
                pct_n = 100 * n_n / n_obs_non if n_obs_non > 0 else np.nan
                rows.append({
                    "Section":       section,
                    "Variable":      f"{var} = {lvl}",
                    "Remitters":     f"{n_r} ({pct_r:.0f}%)",
                    "n Remitters":   coverage_str(col_rem, n_rem) if i == 0 else "",
                    "Non-remitters": f"{n_n} ({pct_n:.0f}%)",
                    "n Non-rem":     coverage_str(col_non, n_non) if i == 0 else "",
                    "_p_raw":        np.nan,
                    "_test":         None,
                })

result = pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# BH correction per family
# ---------------------------------------------------------------------------
for test_family in ("mannwhitney", "fisher"):
    mask  = result["_test"] == test_family
    pvals = result.loc[mask, "_p_raw"].values
    valid = ~np.isnan(pvals)
    if valid.sum() > 1:
        _, p_adj, _, _ = multipletests(pvals[valid], method="fdr_bh")
        adj_full = np.full(len(pvals), np.nan)
        adj_full[valid] = p_adj
        result.loc[mask, "_p_adj"] = adj_full
    elif valid.sum() == 1:
        result.loc[mask & result["_p_raw"].notna(), "_p_adj"] = pvals[valid]

if "_p_adj" not in result.columns:
    result["_p_adj"] = np.nan

result["p-value"] = result["_p_raw"].apply(pval_str)
result["p-adj"]   = result["_p_adj"].apply(pval_str)
result = result.drop(columns=["_p_raw", "_p_adj", "_test"])

# ---------------------------------------------------------------------------
# Save CSV
# ---------------------------------------------------------------------------
result.to_csv(OUT_CSV, index=False)
print(f"Written: {OUT_CSV}")

# ---------------------------------------------------------------------------
# Save HTML
# ---------------------------------------------------------------------------
def _build_html(df_in: pd.DataFrame) -> str:

    def _cov(val):
        val = str(val) if pd.notna(val) and str(val).strip() not in ("", "nan") else ""
        if not val:
            return ""
        style = ("color:#e67e22;font-weight:600;" if "miss." in val
                 else "color:#888;font-size:0.82em;")
        return f'<span style="{style}">{val}</span>'

    def _p(val):
        val = str(val) if pd.notna(val) and str(val) != "nan" else "—"
        if val == "<0.001":
            return f'<strong style="color:#c0392b;">{val}</strong>'
        try:
            if float(val) < 0.05:
                return f'<strong style="color:#c0392b;">{val}</strong>'
        except ValueError:
            pass
        return val

    lines = [
        '<table class="desc-table">',
        "<thead><tr>"
        f"<th>Variable</th>"
        f"<th>Remitters (n={n_rem})</th><th>n</th>"
        f"<th>Non-remitters (n={n_non})</th><th>n</th>"
        "<th>p-value</th><th>p-adj (BH)</th>"
        "</tr></thead>",
        "<tbody>",
    ]
    prev_section = None
    for _, row in df_in.iterrows():
        if row["Section"] != prev_section:
            lines.append(
                f'<tr style="background:var(--accent-bg,#eff6ff);font-weight:600;">'
                f'<td colspan="7">{row["Section"]}</td></tr>'
            )
            prev_section = row["Section"]
        lines.append(
            f'<tr><td>{row["Variable"]}</td>'
            f'<td>{row["Remitters"]}</td><td>{_cov(row["n Remitters"])}</td>'
            f'<td>{row["Non-remitters"]}</td><td>{_cov(row["n Non-rem"])}</td>'
            f'<td>{_p(row["p-value"])}</td><td>{_p(row["p-adj"])}</td></tr>'
        )
    lines.append("</tbody></table>")
    return "\n".join(lines)


html_style = """<style>
.desc-table { border-collapse: collapse; width: 100%; font-size: 0.88em; }
.desc-table th, .desc-table td { padding: 5px 9px; border: 1px solid #ddd; }
.desc-table thead { background: #2c3e50; color: #fff; }
.desc-table tr:nth-child(even) { background: #fafbfc; }
</style>"""

with open(OUT_HTML, "w") as f:
    f.write(f"<!DOCTYPE html><html><body>{html_style}{_build_html(result)}</body></html>")
print(f"Written: {OUT_HTML}")

# ---------------------------------------------------------------------------
# Print to stdout
# ---------------------------------------------------------------------------
print(f"\nRemitters (n={n_rem}) vs Non-remitters (n={n_non}) at PRE\n")
prev = None
for _, row in result.iterrows():
    if row["Section"] != prev:
        print(f"\n  [{row['Section']}]")
        prev = row["Section"]
    print(f"    {row['Variable']:<38} Rem: {row['Remitters']:<28} Non-rem: {row['Non-remitters']:<22} p={row['p-value']:<10} p-adj={row['p-adj']}")
