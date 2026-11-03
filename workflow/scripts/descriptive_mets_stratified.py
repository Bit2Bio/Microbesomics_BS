"""
descriptive_mets_stratified.py

Clinical characteristics at PRE stratified by MetS status (MetS+ vs MetS−).
MetS+ n=28, MetS− n=5. Unpaired comparison — patients are different individuals.

Tests:
  Continuous → Mann-Whitney U (scipy.stats.mannwhitneyu, two-sided)
  Binary     → Fisher's exact test (scipy.stats.fisher_exact)
  Categorical → no p-value (cells too sparse with n=5 in MetS− group)

BH correction applied separately for each test family (Mann-Whitney / Fisher).

Note: with MetS− n=5, tests are severely underpowered. p-values are reported
for completeness but should be interpreted with caution.

Output:
  results/tables/descriptive_mets_stratified.csv
  results/tables/descriptive_mets_stratified.html
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
OUT_CSV  = os.path.join(BASE_DIR, "results/tables/descriptive_mets_stratified.csv")
OUT_HTML = os.path.join(BASE_DIR, "results/tables/descriptive_mets_stratified.html")

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

mets_pos = pre[pre["MetS"] == 1]
mets_neg = pre[pre["MetS"] == 0]
n_pos    = len(mets_pos)
n_neg    = len(mets_neg)
assert n_pos + n_neg == 33, f"Expected 33 PRE rows, got {n_pos + n_neg}"
print(f"MetS+ n={n_pos}, MetS− n={n_neg}")


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


def mannwhitney(col_pos: pd.Series, col_neg: pd.Series) -> float:
    a = col_pos.dropna()
    b = col_neg.dropna()
    if len(a) < 3 or len(b) < 2:
        return np.nan
    try:
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return p
    except Exception:
        return np.nan


def fisher_exact(col_pos: pd.Series, col_neg: pd.Series) -> float:
    """Fisher's exact test for binary variable."""
    n1   = col_pos.notna().sum()
    n2   = col_neg.notna().sum()
    pos1 = (col_pos == 1).sum()
    pos2 = (col_neg == 1).sum()
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

        # MetS itself and remission are the stratification variables — skip
        if var in ("MetS", "MetS_score", "remission",
                   "HyperglicemiaMetS", "HypertensionMetS",
                   "HDL_MetS", "TrygliceridesMetS"):
            continue

        if var not in pre.columns:
            print(f"  WARNING: '{var}' not in data — skipping", file=sys.stderr)
            continue

        col_pos = mets_pos[var]
        col_neg = mets_neg[var]

        if vtype == "continuous":
            rows.append({
                "Section":  section,
                "Variable": var,
                "MetS+":    stat_continuous(col_pos),
                "n MetS+":  coverage_str(col_pos, n_pos),
                "MetS−":    stat_continuous(col_neg),
                "n MetS−":  coverage_str(col_neg, n_neg),
                "_p_raw":   mannwhitney(col_pos, col_neg),
                "_test":    "mannwhitney",
            })

        elif vtype == "binary":
            rows.append({
                "Section":  section,
                "Variable": var,
                "MetS+":    stat_binary(col_pos),
                "n MetS+":  coverage_str(col_pos, n_pos),
                "MetS−":    stat_binary(col_neg),
                "n MetS−":  coverage_str(col_neg, n_neg),
                "_p_raw":   fisher_exact(col_pos, col_neg),
                "_test":    "fisher",
            })

        elif vtype == "categorical":
            levels   = sorted(pre[var].dropna().unique())
            n_obs_pos = col_pos.notna().sum()
            n_obs_neg = col_neg.notna().sum()
            for i, lvl in enumerate(levels):
                n_p   = (col_pos == lvl).sum()
                n_n   = (col_neg == lvl).sum()
                pct_p = 100 * n_p / n_obs_pos if n_obs_pos > 0 else np.nan
                pct_n = 100 * n_n / n_obs_neg if n_obs_neg > 0 else np.nan
                rows.append({
                    "Section":  section,
                    "Variable": f"{var} = {lvl}",
                    "MetS+":    f"{n_p} ({pct_p:.0f}%)",
                    "n MetS+":  coverage_str(col_pos, n_pos) if i == 0 else "",
                    "MetS−":    f"{n_n} ({pct_n:.0f}%)",
                    "n MetS−":  coverage_str(col_neg, n_neg) if i == 0 else "",
                    "_p_raw":   np.nan,
                    "_test":    None,
                })

result = pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# BH correction — separately for each test family
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
        f"<th>MetS+ (n={n_pos})</th><th>n</th>"
        f"<th>MetS− (n={n_neg})</th><th>n</th>"
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
            f'<td>{row["MetS+"]}</td><td>{_cov(row["n MetS+"])}</td>'
            f'<td>{row["MetS−"]}</td><td>{_cov(row["n MetS−"])}</td>'
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
print(f"\nMetS+ (n={n_pos}) vs MetS− (n={n_neg}) at PRE\n")
prev = None
for _, row in result.iterrows():
    if row["Section"] != prev:
        print(f"\n  [{row['Section']}]")
        prev = row["Section"]
    print(f"    {row['Variable']:<38} MetS+: {row['MetS+']:<28} MetS−: {row['MetS−']:<22} p={row['p-value']:<10} p-adj={row['p-adj']}")
