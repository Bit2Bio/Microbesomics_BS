"""
descriptive_clinical.py

PRE vs POST descriptive statistics table for the bariatric surgery cohort (n=33 pairs).
Variable selection and grouping driven by config/config.yaml (clinical_descriptive).

Columns:
  Variable | PRE | n PRE | POST | n POST | p-value | p-adj (BH)

Tests (paired, n=33):
  Continuous → Wilcoxon signed-rank (scipy.stats.wilcoxon, two-sided)
  Binary     → McNemar (statsmodels, with continuity correction)
  Categorical → no p-value

BH correction applied separately for each test family:
  - all Wilcoxon p-values corrected together
  - all McNemar p-values corrected together

Output:
  results/tables/descriptive_baseline.csv
  results/tables/descriptive_baseline.html
"""

import sys
import os
import yaml
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests

BASE_DIR = "/home/lorenzo/Microbesomics"
CONFIG   = os.path.join(BASE_DIR, "config/config.yaml")
DATA     = os.path.join(BASE_DIR, "data/processed/bariatric_clinical_66.csv")
OUT_CSV  = os.path.join(BASE_DIR, "results/tables/descriptive_baseline.csv")
OUT_HTML = os.path.join(BASE_DIR, "results/tables/descriptive_baseline.html")

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

# ---------------------------------------------------------------------------
# Load config and data
# ---------------------------------------------------------------------------
with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

sections = cfg["clinical_descriptive"]

df_all = pd.read_csv(DATA, na_values="")

# height stored as character in source — coerce to numeric
for tp in ("PRE", "POST"):
    mask = df_all["timepoint"] == tp
    if "height" in df_all.columns:
        df_all.loc[mask, "height"] = pd.to_numeric(df_all.loc[mask, "height"], errors="coerce")

pre  = df_all[df_all["timepoint"] == "PRE"].set_index("pz")
post = df_all[df_all["timepoint"] == "POST"].set_index("pz")

# Align on the same patients (same 33 pz in both)
shared_pz = pre.index.intersection(post.index)
pre  = pre.loc[shared_pz]
post = post.loc[shared_pz]
n_pairs = len(shared_pz)
assert n_pairs == 33, f"Expected 33 pairs, got {n_pairs}"


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


def wilcoxon_paired(col_pre: pd.Series, col_post: pd.Series) -> float:
    """Wilcoxon signed-rank on complete pairs only."""
    both = pd.DataFrame({"pre": col_pre, "post": col_post}).dropna()
    if len(both) < 5:
        return np.nan
    d = both["post"] - both["pre"]
    if (d == 0).all():
        return np.nan
    try:
        _, p = stats.wilcoxon(both["pre"], both["post"], alternative="two-sided")
        return p
    except Exception:
        return np.nan


def mcnemar_paired(col_pre: pd.Series, col_post: pd.Series) -> float:
    """McNemar test on complete pairs."""
    both = pd.DataFrame({"pre": col_pre, "post": col_post}).dropna()
    if len(both) < 5:
        return np.nan
    # 2×2 contingency: rows=PRE, cols=POST
    ct = pd.crosstab(both["pre"], both["post"])
    # Ensure 2×2
    if ct.shape != (2, 2):
        return np.nan
    try:
        result = mcnemar(ct.values, exact=False, correction=True)
        return result.pvalue
    except Exception:
        return np.nan


# ---------------------------------------------------------------------------
# Compute rows (raw p-values stored as float; BH correction applied after)
# ---------------------------------------------------------------------------
rows = []  # each row carries _p_raw (float or NaN) and _test ("wilcoxon"/"mcnemar"/None)

for section, variables in sections.items():
    for var, meta in variables.items():
        if isinstance(meta, dict):
            vtype   = meta["type"]
            compare = bool(meta.get("compare", True))
        else:
            vtype   = meta
            compare = True

        if var not in pre.columns:
            print(f"  WARNING: '{var}' not in data — skipping", file=sys.stderr)
            continue

        col_pre  = pre[var]
        col_post = post[var]

        if vtype == "continuous":
            stat_p = stat_continuous(col_pre)
            cov_p  = coverage_str(col_pre, n_pairs)
            if compare:
                stat_po = stat_continuous(col_post)
                cov_po  = coverage_str(col_post, n_pairs)
                p_raw   = wilcoxon_paired(col_pre, col_post)
                test    = "wilcoxon"
            else:
                stat_po = cov_po = "—"
                p_raw   = np.nan
                test    = None
            rows.append({"Section": section, "Variable": var,
                         "PRE": stat_p, "n PRE": cov_p,
                         "POST": stat_po, "n POST": cov_po,
                         "_p_raw": p_raw, "_test": test})

        elif vtype == "binary":
            stat_p = stat_binary(col_pre)
            cov_p  = coverage_str(col_pre, n_pairs)
            if compare:
                stat_po = stat_binary(col_post)
                cov_po  = coverage_str(col_post, n_pairs)
                p_raw   = mcnemar_paired(col_pre, col_post)
                test    = "mcnemar"
            else:
                stat_po = cov_po = "—"
                p_raw   = np.nan
                test    = None
            rows.append({"Section": section, "Variable": var,
                         "PRE": stat_p, "n PRE": cov_p,
                         "POST": stat_po, "n POST": cov_po,
                         "_p_raw": p_raw, "_test": test})

        elif vtype == "categorical":
            levels   = sorted(pd.concat([col_pre, col_post]).dropna().unique())
            n_obs_p  = col_pre.notna().sum()
            n_obs_po = col_post.notna().sum()
            for i, lvl in enumerate(levels):
                n_p   = (col_pre == lvl).sum()
                pct_p = 100 * n_p / n_obs_p if n_obs_p > 0 else np.nan
                label = f"{var} = {lvl}"
                if compare:
                    n_po      = (col_post == lvl).sum()
                    pct_po    = 100 * n_po / n_obs_po if n_obs_po > 0 else np.nan
                    post_stat = f"{n_po} ({pct_po:.0f}%)"
                    post_cov  = coverage_str(col_post, n_pairs) if i == 0 else ""
                else:
                    post_stat = post_cov = "—"
                rows.append({"Section": section, "Variable": label,
                             "PRE": f"{n_p} ({pct_p:.0f}%)",
                             "n PRE": coverage_str(col_pre, n_pairs) if i == 0 else "",
                             "POST": post_stat, "n POST": post_cov,
                             "_p_raw": np.nan, "_test": None})

result = pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# BH correction — separately for Wilcoxon and McNemar families
# ---------------------------------------------------------------------------
for test_family in ("wilcoxon", "mcnemar"):
    mask = result["_test"] == test_family
    p_vals = result.loc[mask, "_p_raw"].values
    valid  = ~np.isnan(p_vals)
    if valid.sum() > 1:
        _, p_adj, _, _ = multipletests(p_vals[valid], method="fdr_bh")
        adj_full = np.full(len(p_vals), np.nan)
        adj_full[valid] = p_adj
        result.loc[mask, "_p_adj"] = adj_full
    elif valid.sum() == 1:
        result.loc[mask & result["_p_raw"].notna(), "_p_adj"] = p_vals[valid]

if "_p_adj" not in result.columns:
    result["_p_adj"] = np.nan

# ---------------------------------------------------------------------------
# Format p-values as strings for display
# ---------------------------------------------------------------------------
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
    lines = ['<table class="desc-table">',
             "<thead><tr>"
             "<th>Variable</th>"
             "<th>PRE</th><th>n PRE</th>"
             "<th>POST</th><th>n POST</th>"
             "<th>p-value</th><th>p-adj (BH)</th>"
             "</tr></thead>",
             "<tbody>"]
    prev_section = None
    for _, row in df_in.iterrows():
        if row["Section"] != prev_section:
            lines.append(
                f'<tr style="background:var(--accent-bg,#eff6ff);font-weight:600;">'
                f'<td colspan="7">{row["Section"]}</td></tr>'
            )
            prev_section = row["Section"]

        def _cov(val):
            val = str(val) if pd.notna(val) else ""
            if not val:
                return ""
            style = "color:#e67e22;font-weight:600;" if "miss." in val else "color:#888;font-size:0.82em;"
            return f'<span style="{style}">{val}</span>'

        def _p(val):
            val = str(val) if pd.notna(val) else "—"
            if val not in ("—", "") and val != "<0.001":
                try:
                    fval = float(val)
                    if fval < 0.05:
                        return f'<strong style="color:#c0392b;">{val}</strong>'
                except ValueError:
                    pass
            if val == "<0.001":
                return f'<strong style="color:#c0392b;">{val}</strong>'
            return val

        lines.append(
            f'<tr>'
            f'<td>{row["Variable"]}</td>'
            f'<td>{row["PRE"]}</td><td>{_cov(row["n PRE"])}</td>'
            f'<td>{row["POST"]}</td><td>{_cov(row["n POST"])}</td>'
            f'<td>{_p(row["p-value"])}</td>'
            f'<td>{_p(row["p-adj"])}</td>'
            f'</tr>'
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
print(f"\nPRE vs POST (n={n_pairs} pairs)\n")
prev = None
for _, row in result.iterrows():
    if row["Section"] != prev:
        print(f"\n  [{row['Section']}]")
        prev = row["Section"]
    print(f"    {row['Variable']:<38} PRE: {row['PRE']:<28} POST: {row['POST']:<28} p={row['p-value']:<10} p-adj={row['p-adj']}")
