#!/usr/bin/env python3
"""
make_html_report.py — Microbesomics_BS HTML report

Translates the content of the retired generate_pptx.py into a self-contained
HTML report following the Atheromics_TVAS style: base64-encoded images,
sortable tables, dark-mode support, no external network calls.

Report sections (mirrors the 7-slide PPTX):
  1. Study overview  — cohort description + ASV filtering pipeline
  2. MetS criteria   — NCEP ATP III table + implementation notes
  3. Volcano plots   — MetS baseline | Remission | PRE vs POST
  4. Heatmap         — log2FC across comparisons (intersection)
  5. Barplots Phylum — relative abundance, 3 comparisons
  6. Barplots Family
  7. Barplots Genus

Arguments:
    --out   output HTML path (passed by Snakemake)
"""

import argparse
import base64
import io
import os
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True)
args = parser.parse_args()

BASE_DIR = "/home/lorenzo/Microbesomics"


def img_tag(path: str, caption: str = "", style: str = "width:100%;height:auto;border-radius:8px;") -> str:
    """Return an <img> tag with base64-encoded PNG. Returns empty string if file missing."""
    p = Path(path)
    if not p.exists():
        return f'<div class="missing">Image not found: {path}</div>'
    data = base64.b64encode(p.read_bytes()).decode("ascii")
    src  = f"data:image/png;base64,{data}"
    cap  = f'<div class="meta">{caption}</div>' if caption else ""
    return f'<img src="{src}" style="{style}">{cap}'


def three_panel(images: list, labels: list) -> str:
    """Lay out three images side-by-side with labels underneath."""
    cells = ""
    for path, label in zip(images, labels):
        cells += f"""
        <div style="flex:1;min-width:0;text-align:center;">
          {img_tag(path, style="width:100%;height:auto;border-radius:6px;")}
          <div class="panel-label">{label}</div>
        </div>"""
    return f'<div style="display:flex;gap:16px;align-items:flex-start;">{cells}</div>'


def xlsx_download_button(df: pd.DataFrame, filename: str, label: str) -> str:
    """Encode a DataFrame as xlsx base64 and return an HTML download button."""
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="clinical_data")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    return (
        f'<a href="data:{mime};base64,{b64}" download="{filename}" '
        f'style="display:inline-block;padding:8px 18px;background:var(--accent);'
        f'color:#fff;border-radius:8px;text-decoration:none;font-size:0.88rem;'
        f'font-weight:600;margin-top:10px;">{label}</a>'
    )


def df_to_table(df: pd.DataFrame, table_id: str) -> str:
    thead = "".join(
        f'<th onclick="sortTable(\'{table_id}\',{i})">{c}</th>'
        for i, c in enumerate(df.columns)
    )
    rows = []
    for _, r in df.iterrows():
        cells = "".join(f"<td>{'' if pd.isna(v) else v}</td>" for v in r)
        rows.append(f"<tr>{cells}</tr>")
    tbody = "".join(rows)
    return f'<table id="{table_id}"><thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table>'


# --------------------------------------------------------------------------
# Clinical data
# --------------------------------------------------------------------------
clinical_df = pd.read_csv(f"{BASE_DIR}/data/processed/clinical_all.csv", na_values="")

# --------------------------------------------------------------------------
# Descriptive baseline table
# --------------------------------------------------------------------------
_desc_path = f"{BASE_DIR}/results/tables/descriptive_baseline.csv"

def _build_desc_html() -> str:
    if not Path(_desc_path).exists():
        return "<p class='missing'>descriptive_baseline.csv not found — run descriptive_clinical.py first.</p>"
    desc = pd.read_csv(_desc_path)

    def _cov(val):
        val = str(val) if pd.notna(val) and str(val).strip() not in ("", "nan") else ""
        if not val:
            return ""
        style = ("color:#e67e22;font-weight:600;" if "miss." in val
                 else "color:var(--muted);font-size:0.82em;")
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

    lines = ['<table id="tbl-desc"><thead><tr>'
             '<th>Variable</th>'
             '<th>PRE</th><th>n PRE</th>'
             '<th>POST</th><th>n POST</th>'
             '<th>p-value</th><th>p-adj (BH)</th>'
             '</tr></thead><tbody>']
    prev_section = None
    for _, row in desc.iterrows():
        if row["Section"] != prev_section:
            lines.append(
                f'<tr style="background:var(--accent-bg);font-weight:600;">'
                f'<td colspan="7">{row["Section"]}</td></tr>'
            )
            prev_section = row["Section"]
        lines.append(
            f'<tr><td>{row["Variable"]}</td>'
            f'<td>{row["PRE"]}</td><td>{_cov(row["n PRE"])}</td>'
            f'<td>{row["POST"]}</td><td>{_cov(row["n POST"])}</td>'
            f'<td>{_p(row["p-value"])}</td><td>{_p(row["p-adj"])}</td></tr>'
        )
    lines.append("</tbody></table>")
    return "\n".join(lines)

desc_table_html = _build_desc_html()

# --------------------------------------------------------------------------
# MetS-stratified descriptive table
# --------------------------------------------------------------------------
_mets_path = f"{BASE_DIR}/results/tables/descriptive_mets_stratified.csv"

def _build_mets_html() -> str:
    if not Path(_mets_path).exists():
        return "<p class='missing'>descriptive_mets_stratified.csv not found — run descriptive_mets_stratified.py first.</p>"
    desc = pd.read_csv(_mets_path)

    def _cov(val):
        val = str(val) if pd.notna(val) and str(val).strip() not in ("", "nan") else ""
        if not val:
            return ""
        style = ("color:#e67e22;font-weight:600;" if "miss." in val
                 else "color:var(--muted);font-size:0.82em;")
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

    lines = ['<table id="tbl-mets-strat"><thead><tr>'
             '<th>Variable</th>'
             '<th>MetS+ (n=28)</th><th>n</th>'
             '<th>MetS− (n=5)</th><th>n</th>'
             '<th>p-value</th><th>p-adj (BH)</th>'
             '</tr></thead><tbody>']
    prev_section = None
    for _, row in desc.iterrows():
        if row["Section"] != prev_section:
            lines.append(
                f'<tr style="background:var(--accent-bg);font-weight:600;">'
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

mets_strat_html = _build_mets_html()

# --------------------------------------------------------------------------
# Remission-stratified descriptive table
# --------------------------------------------------------------------------
_rem_path = f"{BASE_DIR}/results/tables/descriptive_remission_stratified.csv"

def _build_rem_html() -> str:
    if not Path(_rem_path).exists():
        return "<p class='missing'>descriptive_remission_stratified.csv not found — run descriptive_remission_stratified.py first.</p>"
    desc = pd.read_csv(_rem_path)

    def _cov(val):
        val = str(val) if pd.notna(val) and str(val).strip() not in ("", "nan") else ""
        if not val:
            return ""
        style = ("color:#e67e22;font-weight:600;" if "miss." in val
                 else "color:var(--muted);font-size:0.82em;")
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

    lines = ['<table id="tbl-rem-strat"><thead><tr>'
             '<th>Variable</th>'
             '<th>Remitters (n=22)</th><th>n</th>'
             '<th>Non-remitters (n=6)</th><th>n</th>'
             '<th>p-value</th><th>p-adj (BH)</th>'
             '</tr></thead><tbody>']
    prev_section = None
    for _, row in desc.iterrows():
        if row["Section"] != prev_section:
            lines.append(
                f'<tr style="background:var(--accent-bg);font-weight:600;">'
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

rem_strat_html = _build_rem_html()

# --------------------------------------------------------------------------
# MetS criteria table
# --------------------------------------------------------------------------
CRITERIA = [
    ("Abdominal Obesity",    "waist ≥ 102 cm (M) / ≥ 88 cm (F)",                    "AbdominalObesityMetS"),
    ("Hyperglycemia",        "fasting glucose ≥ 100 mg/dL  OR  diabetes == 1",       "HyperglicemiaMetS"),
    ("Hypertension",         "SBP ≥ 130  OR  DBP ≥ 85  OR  hypertension == 1",      "HypertensionMetS"),
    ("Low HDL",              "HDL &lt; 40 mg/dL (M) / &lt; 50 mg/dL (F)",           "HDL_MetS"),
    ("Hypertriglyceridemia", "TG ≥ 150 mg/dL  OR  TG_pharmacotherapy == 1",         "TrygliceridesMetS"),
]

NOTES = [
    "Positive: ≥ 3 of 5 criteria met. Negative: impossible to reach 3 even accounting for missing data. NA: missing data is decisive.",
    "sex coding: 0 = Male, 1 = Female (inverted from the intuitive convention).",
    "TG_pharmacotherapy = 1 for 3 samples on TG-specific drugs (fenofibrate / omega-3): S85 (C19, grp3 PRE), S200 (CO-71, grp5 PRE), S215 (DA-157, grp6 POST). These always count as positive for criterion 5 regardless of measured TG.",
    "dyslipidaemia flag NOT used for HDL/TG: in this cohort it primarily reflects high LDL; using it would add up to 2 false MetS points per patient.",
    "hypertension flag IS used: cross-check with drug columns confirmed 0 cases of antihypertensive treatment with hypertension == 0.",
]

mets_df = pd.DataFrame(CRITERIA, columns=["Criterion", "Threshold / Rule", "Variable"])

# --------------------------------------------------------------------------
# DESeq2 result counts (load from xlsx if present)
# --------------------------------------------------------------------------
def load_deseq2_summary(xlsx_path: str) -> str:
    p = Path(xlsx_path)
    if not p.exists():
        return "<p class='missing'>Results file not found.</p>"
    df = pd.read_excel(p)
    # Expect columns: padj, log2FoldChange (or log2FC)
    fc_col  = next((c for c in df.columns if "log2" in c.lower()), None)
    adj_col = next((c for c in df.columns if "padj" in c.lower()), None)
    if adj_col is None:
        return "<p class='meta'>padj column not found.</p>"
    sig = df[df[adj_col] < 0.05]
    total = len(df)
    n_sig = len(sig)
    up = dn = "n/a"
    if fc_col:
        up = len(sig[sig[fc_col] > 0])
        dn = len(sig[sig[fc_col] < 0])
    return f"""<p class="meta">
        Total tested: <strong>{total}</strong> ASVs &nbsp;|&nbsp;
        Significant (FDR &lt; 0.05): <strong>{n_sig}</strong> &nbsp;|&nbsp;
        Enriched (positive log2FC): <strong>{up}</strong> &nbsp;|&nbsp;
        Depleted (negative log2FC): <strong>{dn}</strong>
    </p>"""


# --------------------------------------------------------------------------
# Beta diversity PERMANOVA table
# --------------------------------------------------------------------------
_perm_path = f"{BASE_DIR}/results/tables/beta_diversity_permanova.csv"

def _build_permanova_html() -> str:
    if not Path(_perm_path).exists():
        return "<p class='missing'>beta_diversity_permanova.csv not found.</p>"
    df = pd.read_csv(_perm_path)
    def _p(val):
        try:
            v = float(val)
            s = "<0.001" if v < 0.001 else f"{v:.3f}"
            color = ";color:#c0392b;font-weight:600" if v < 0.05 else ""
            return f'<span style="font-variant-numeric:tabular-nums{color}">{s}</span>'
        except Exception:
            return str(val)
    cols = ["distance", "n_samples", "n_pairs", "F_statistic", "R2", "p_value"]
    heads = ["Distance", "N samples", "N pairs", "F", "R²", "p-value"]
    rows = []
    for _, r in df.iterrows():
        cells = "".join(
            f"<td>{_p(r[c]) if c in ('p_value',) else round(r[c], 3) if isinstance(r[c], float) else r[c]}</td>"
            for c in cols
        )
        rows.append(f"<tr>{cells}</tr>")
    thead = "".join(f"<th>{h}</th>" for h in heads)
    return f'<table><thead><tr>{thead}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


# --------------------------------------------------------------------------
# Build HTML
# --------------------------------------------------------------------------
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Microbesomics_BS — Analysis Report</title>
<style>
  :root {{
    --bg: #f7f8fa; --card: #ffffff; --text: #1a1d23; --muted: #6b7280;
    --border: #e5e7eb; --accent: #2563eb; --accent-bg: #eff6ff;
    --stripe: #fafbfc; --badge-up: #dcfce7; --badge-dn: #fee2e2;
    --badge-up-text: #166534; --badge-dn-text: #991b1b;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #0f1115; --card: #171a21; --text: #e6e8eb; --muted: #9aa3af;
      --border: #2a2e37; --accent: #5b9bff; --accent-bg: #14213a;
      --stripe: #1b1e26; --badge-up: #14532d; --badge-dn: #450a0a;
      --badge-up-text: #86efac; --badge-dn-text: #fca5a5;
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 32px 16px; background: var(--bg); color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  }}
  .wrap {{ max-width: 1200px; margin: 0 auto; }}
  h1 {{ font-size: 1.7rem; margin-bottom: 4px; }}
  .subtitle {{ color: var(--muted); margin-bottom: 32px; font-size: 0.95rem; }}
  .card {{
    background: var(--card); border: 1px solid var(--border); border-radius: 12px;
    padding: 20px 24px; margin-bottom: 28px; box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  }}
  h2 {{ font-size: 1.15rem; margin-top: 0; margin-bottom: 14px; color: var(--text); }}
  h3 {{ font-size: 1rem; margin: 20px 0 8px; color: var(--accent); }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); }}
  th {{
    cursor: pointer; user-select: none; color: var(--muted); font-weight: 600;
    position: sticky; top: 0; background: var(--card);
  }}
  th:hover {{ color: var(--accent); }}
  th::after {{ content: " \\21C5"; font-size: 0.7em; opacity: 0.4; }}
  tbody tr:nth-child(even) {{ background: var(--stripe); }}
  tbody tr:hover {{ background: var(--accent-bg); }}
  .table-scroll {{ overflow-x: auto; max-height: 400px; overflow-y: auto; }}
  .meta {{ color: var(--muted); font-size: 0.82rem; margin-top: 8px; line-height: 1.6; }}
  .panel-label {{ color: var(--muted); font-size: 0.82rem; margin-top: 6px; }}
  .missing {{ color: #ef4444; font-size: 0.85rem; padding: 8px; background: #fee2e2; border-radius: 6px; }}
  .flow {{
    display: flex; flex-direction: column; align-items: center; gap: 0; max-width: 480px;
  }}
  .flow-box {{
    width: 100%; padding: 8px 14px; border-radius: 8px; text-align: center;
    font-size: 0.85rem; border: 1px solid var(--border);
    background: var(--accent-bg); color: var(--text);
  }}
  .flow-box.dark {{ background: var(--accent); color: #fff; font-weight: 600; }}
  .flow-arrow {{ color: var(--accent); font-size: 1.2rem; line-height: 1; padding: 2px 0; }}
  .flow-detail {{ font-size: 0.78rem; color: var(--muted); }}
  .two-col {{ display: flex; gap: 32px; }}
  .two-col > div {{ flex: 1; min-width: 0; }}
  ul.notes {{ padding-left: 18px; margin: 0; }}
  ul.notes li {{ font-size: 0.84rem; color: var(--muted); margin-bottom: 6px; line-height: 1.5; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>Microbesomics_BS — Analysis Report</h1>
  <div class="subtitle">Bariatric surgery cohort · Gut microbiome · MetS · 16S rRNA amplicon sequencing</div>

  <!-- ======================================================= -->
  <!-- 1. Study Overview                                        -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>1. Study Overview</h2>
    <div class="two-col">

      <div>
        <h3>Cohort</h3>
        <ul class="notes" style="color:var(--text);">
          <li>n = 33 patients, complete PRE / POST bariatric surgery pairs (66 samples)</li>
          <li>MetS at PRE: <strong>28 MetS+</strong> (85%) · <strong>5 MetS−</strong> (15%)</li>
          <li>MetS at POST: <strong>6 MetS+</strong> (18%) · <strong>27 MetS−</strong> (82%)</li>
        </ul>

        <h3>Remission classification (PRE → POST)</h3>
        <ul class="notes" style="color:var(--text);">
          <li><strong>Remitters</strong>: MetS+ → MetS− &nbsp; n = 22 (67%)</li>
          <li><strong>Non-remitters</strong>: MetS+ → MetS+ &nbsp; n = 6 (18%)</li>
          <li><strong>MetS−</strong>: MetS− at PRE &nbsp; n = 5 (15%)</li>
        </ul>

        <h3>Omics platform</h3>
        <ul class="notes" style="color:var(--text);">
          <li>16S rRNA amplicon sequencing (gut microbiome)</li>
          <li>QIIME2 processing — ASVs, SILVA taxonomy, FastTree phylogeny</li>
        </ul>

        <h3>Comparisons (DESeq2)</h3>
        <ul class="notes" style="color:var(--text);">
          <li>PRE vs POST surgery — paired design (n = 33)</li>
          <li>Remitters vs Non-remitters at PRE (22 vs 6)</li>
          <li>MetS+ vs MetS− at PRE (28 vs 5)</li>
          <li>Remitters vs MetS− at PRE (22 vs 5)</li>
        </ul>
      </div>

      <div>
        <h3>ASV filtering pipeline</h3>
        <div class="flow">
          <div class="flow-box dark">Raw phyloseq<br><span class="flow-detail" style="color:#d1d5db;">247 samples · 14 741 ASVs</span></div>
          <div class="flow-arrow">▼</div>
          <div class="flow-box">Subset: 33 complete pairs<br><span class="flow-detail">66 samples · 14 741 ASVs</span></div>
          <div class="flow-arrow">▼</div>
          <div class="flow-box">Prevalence filter<br><span class="flow-detail">≥ 4 reads in ≥ 10% of samples (≥ 7 / 66)</span></div>
          <div class="flow-arrow">▼</div>
          <div class="flow-box">Variance filter<br><span class="flow-detail">bottom 10% by IQR removed</span></div>
          <div class="flow-arrow">▼</div>
          <div class="flow-box dark">Final dataset<br><span class="flow-detail" style="color:#d1d5db;">66 samples · 707 ASVs</span></div>
        </div>
      </div>

    </div>
  </div>

  <!-- ======================================================= -->
  <!-- QC — Rarefaction curves                                  -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>QC — Rarefaction curves</h2>
    <p class="meta" style="margin-bottom:14px;">
      Metodi: rarefazione manuale con <code>rrarefy()</code> (vegan) su ps_filt_5vs6 (66 campioni, ASV filtrati).
      Viola = PRE · Teal = POST. Linee tratteggiate: minimo (8,069 reads — outlier SS-04),
      10° percentile (~123,923), 20° percentile (~144,977). La quasi totalità delle curve
      raggiunge il plateau prima di 100,000 reads; soglia di rarefazione scelta = 75,000 reads.
      SS-04 (PRE, 8,069 reads) escluso dalle analisi del microbioma.
    </p>
    {img_tag(f"{BASE_DIR}/results/figures/qc/rarefaction_curves.png",
             style="width:80%;display:block;margin:0 auto;")}
  </div>

  <!-- ======================================================= -->
  <!-- Clinical data + descriptive tables                       -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>Dati clinici — caratteristiche baseline</h2>
    <p class="meta" style="margin-bottom:14px;">
      Metodi: variabili continue mediana [IQR], Wilcoxon signed-rank paired con correzione BH.
      Binarie: n (%), McNemar. Criteri MetS: NCEP ATP III (revisione AHA/NHLBI 2005).
    </p>
    <div class="table-scroll">{desc_table_html}</div>
    {xlsx_download_button(clinical_df, "clinical_all.xlsx", "⬇ Scarica clinical_all.xlsx")}
  </div>

  <div class="card">
    <h2>MetS+ vs MetS− al baseline (n=28 vs n=5)</h2>
    <p class="meta" style="margin-bottom:14px;">
      Metodi: Mann-Whitney U per continui, Fisher exact per binari, correzione BH per famiglia.
      <strong>Nota:</strong> MetS− n=5 — sottopotenziato, p-value puramente esplorativi.
    </p>
    <div class="table-scroll">{mets_strat_html}</div>
  </div>

  <div class="card">
    <h2>Remitters vs Non-remitters al baseline (n=22 vs n=6)</h2>
    <p class="meta" style="margin-bottom:14px;">
      Solo MetS+ (n=28). Metodi: Mann-Whitney U, Fisher exact, correzione BH per famiglia.
      <strong>Nota:</strong> Non-remitters n=6 — risultati esplorativi.
    </p>
    <div class="table-scroll">{rem_strat_html}</div>
  </div>

  <!-- ======================================================= -->
  <!-- Figure 1 — Clinical characterization                     -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>Figure 1 — Caratterizzazione clinica PRE → POST</h2>
    <p class="meta" style="margin-bottom:14px;">
      Metodi: Wilcoxon signed-rank paired, correzione BH su tutte le variabili continue congiuntamente.
      <strong>a</strong> Ridge plot — antropometria, segni vitali, composizione corporea, fegato (11 variabili, unità SD del PRE).
      <strong>b</strong> Ridge plot — glicemia, OGTT, lipidi (12 variabili).
      <strong>c</strong> Raincloud — peso (kg) PRE vs POST.
      <strong>d</strong> Alluvial — flusso MetS PRE→POST: 22/28 MetS+ raggiungono la remissione.
      Risultati: la chirurgia bariatrica induce miglioramenti significativi (p-adj &lt;0.05) in quasi tutte le variabili metaboliche e antropometriche.
    </p>
    {img_tag(f"{BASE_DIR}/results/figures/figure1_combined.png")}
  </div>

  <!-- ======================================================= -->
  <!-- Alpha diversity                                          -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>Alpha diversity PRE → POST</h2>
    <p class="meta" style="margin-bottom:14px;">
      Metodi: phyloseq non filtrato (tutti gli ASVs), SS-04 escluso, rarefazione a 75,000 reads
      (n=32 coppie paired). Wilcoxon signed-rank paired. Metriche: Observed ASVs e Shannon index.
      Risultati: nessuna differenza significativa tra PRE e POST
      (Observed: p=0.78 · Shannon: p=0.25). La richness e la diversità alpha
      non cambiano in modo rilevante dopo chirurgia bariatrica in questa coorte.
    </p>
    {img_tag(f"{BASE_DIR}/results/figures/figure_alpha_diversity.png",
             style="width:40%;display:block;margin:0 auto;")}
  </div>

  <!-- ======================================================= -->
  <!-- 2. MetS Classification                                   -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>2. MetS Classification — NCEP ATP III (AHA/NHLBI 2005 revision)</h2>
    <p class="meta" style="margin-bottom:14px;">
      Positive = ≥ 3 of 5 criteria met. Waist thresholds: ATP III 2001 revised 2005.
      Glucose threshold 100 mg/dL (lowered from 110 in 2005 revision).
    </p>
    <div class="table-scroll">
      {df_to_table(mets_df, "tbl-mets")}
    </div>
    <h3 style="margin-top:20px;">Implementation notes</h3>
    <ul class="notes">
      {''.join(f"<li>{n}</li>" for n in NOTES)}
    </ul>
  </div>

  <!-- ======================================================= -->
  <!-- 3. Volcano plots                                         -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>3. Differential abundance — DESeq2 volcano plots</h2>

    <h3>3.1 MetS+ vs MetS− at baseline (PRE)</h3>
    {load_deseq2_summary(f"{BASE_DIR}/results/deseq2/deseq2_mets_baseline_results.xlsx")}
    {img_tag(f"{BASE_DIR}/results/deseq2/volcano_mets_baseline.png", "Positive log2FC = enriched in MetS+")}

    <h3>3.2 Remitters vs Non-remitters at baseline (PRE)</h3>
    {load_deseq2_summary(f"{BASE_DIR}/results/deseq2/deseq2_remission_results.xlsx")}
    {img_tag(f"{BASE_DIR}/results/deseq2/volcano_remission.png", "Positive log2FC = enriched in non-remitters")}

    <h3>3.3 PRE vs POST bariatric surgery (paired)</h3>
    {load_deseq2_summary(f"{BASE_DIR}/results/deseq2/deseq2_5vs6_results.xlsx")}
    {img_tag(f"{BASE_DIR}/results/deseq2/volcano_5vs6.png", "Positive log2FC = enriched POST-surgery (group 6)")}

    <h3>3.4 MetS+ remitters vs MetS− at baseline (PRE)</h3>
    {load_deseq2_summary(f"{BASE_DIR}/results/deseq2/deseq2_remitters_vs_metsneg_results.xlsx")}
    {img_tag(f"{BASE_DIR}/results/deseq2/volcano_remitters_vs_metsneg.png", "Positive log2FC = enriched in remitters vs MetS−")}
  </div>

  <!-- ======================================================= -->
  <!-- 4. Heatmap intersection                                  -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>4. ASV log2FC across comparisons</h2>
    <p class="meta" style="margin-bottom:14px;">
      ASVs significant (FDR &lt; 0.05) in at least one of: remission or MetS baseline comparisons.
      Columns = ASVs (clustered by log2FC profile). Rows = 3 comparisons. Fill = log2FC; * = FDR &lt; 0.05.
    </p>
    {img_tag(f"{BASE_DIR}/results/deseq2/heatmap_intersection.png")}
  </div>

  <!-- ======================================================= -->
  <!-- 5–7. Taxonomy barplots                                   -->
  <!-- ======================================================= -->
  <div class="card">
    <h2>5. Taxonomic composition — Phylum level (top 15)</h2>
    {three_panel(
        [f"{BASE_DIR}/results/metagenomics/barplot_Phylum_mets.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Phylum_remission.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Phylum_5vs6.png"],
        ["MetS+ vs MetS−", "Remitters vs Non-remitters", "PRE vs POST"]
    )}
  </div>

  <div class="card">
    <h2>6. Taxonomic composition — Family level (top 15)</h2>
    {three_panel(
        [f"{BASE_DIR}/results/metagenomics/barplot_Family_mets.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Family_remission.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Family_5vs6.png"],
        ["MetS+ vs MetS−", "Remitters vs Non-remitters", "PRE vs POST"]
    )}
  </div>

  <div class="card">
    <h2>7. Taxonomic composition — Genus level (top 15)</h2>
    {three_panel(
        [f"{BASE_DIR}/results/metagenomics/barplot_Genus_mets.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Genus_remission.png",
         f"{BASE_DIR}/results/metagenomics/barplot_Genus_5vs6.png"],
        ["MetS+ vs MetS−", "Remitters vs Non-remitters", "PRE vs POST"]
    )}
  </div>

</div>

<script>
function sortTable(tableId, colIdx) {{
  const table = document.getElementById(tableId);
  const tbody = table.tBodies[0];
  const rows  = Array.from(tbody.rows);
  const asc   = table.dataset.sortCol == colIdx && table.dataset.sortDir !== 'asc';
  rows.sort((a, b) => {{
    let x = a.cells[colIdx].innerText, y = b.cells[colIdx].innerText;
    const fx = parseFloat(x), fy = parseFloat(y);
    if (!isNaN(fx) && !isNaN(fy)) {{ x = fx; y = fy; }}
    if (x < y) return asc ? -1 : 1;
    if (x > y) return asc ? 1  : -1;
    return 0;
  }});
  rows.forEach(r => tbody.appendChild(r));
  table.dataset.sortCol = colIdx;
  table.dataset.sortDir = asc ? 'asc' : 'desc';
}}
</script>
</body>
</html>
"""

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(html, encoding="utf-8")
print(f"Saved: {out}")
