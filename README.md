# Microbesomics_MetS

Multi-omics analysis of gut microbiome signatures associated with Metabolic Syndrome (MetS) and its individual components, including pre/post intervention comparisons (dietary intervention and bariatric surgery) and external cohort validation.

## Study Design

- **Main cohort**: n=177 cross-sectional samples, MetS classification based on NCEP ATP III criteria (≥3 of 5)
- **Diet intervention**: paired pre/post comparison (Groups 1 vs 2)
- **Bariatric surgery**: paired pre/post comparison (Groups 5 vs 6), n=33 complete pairs with metagenomics at both timepoints
- **External validation**: Floromidia cohort (independent replication)

### MetS Criteria (NCEP ATP III — AHA/NHLBI 2005 revision)

MetS is defined as the presence of **≥3 of the following 5 criteria**:

| Criterion | Threshold | Notes |
|-----------|-----------|-------|
| Abdominal Obesity | Waist ≥ 102 cm (M) / ≥ 88 cm (F) | ATP III 2001 used strict >; 2005 revision uses ≥ |
| Hyperglycemia | Fasting glucose ≥ 100 mg/dL **or** `diabetes == 1` | 2005 revision lowered threshold from 110 to 100 |
| Hypertension | SBP ≥ 130 mmHg **or** DBP ≥ 85 mmHg **or** `hypertension == 1` | Binary flag confirmed to fully capture treated patients |
| Low HDL | HDL < 40 mg/dL (M) / < 50 mg/dL (F) | Lab values only — see drug note below |
| Hypertriglyceridemia | Triglycerides ≥ 150 mg/dL | Lab values only — see drug note below |

> **Note on waist thresholds**: IDF criteria (94/80 cm) were considered but are not used. ATP III thresholds (102/88 cm) are applied throughout, as the cohort is entirely obese and the IDF lower bounds would classify virtually all subjects as positive.

### Drug Column Handling

The dataset merges two cohorts with different data entry formats for drug columns:

- **Mingrone cohort** (`M*` codes): drug columns encoded as `"0"` (not on drug) or drug name string
- **Bariatric cohort** (`XX-NNN` codes): drug columns are `NA` (absent) or drug name string

**Why drug columns are NOT used for HDL and triglycerides criteria:**  
Inspection of all lipid-lowering drug names in the dataset shows they are predominantly statins and ezetimibe, which target LDL and are not relevant for the HDL or triglycerides MetS criteria. Three samples are on TG-specific therapy at specific timepoints and are flagged with `TG_pharmacotherapy = 1`; all others are 0.

| sample_id | pz | group | drug |
|-----------|-----|-------|------|
| S85 | C19 | 3 (PRE) | Fenofibrate |
| S200 | CO-71 | 5 (PRE) | Fulcosupra (fenofibrate) |
| S215 | DA-157 | 6 (POST) | Eskimo (omega-3 ethyl esters) |

These samples count as positive for the hypertriglyceridemia criterion regardless of their measured TG value, as the drug may artificially suppress TG below the 150 mg/dL threshold.

**Why `dyslipidaemia == 1` is NOT used for HDL and triglycerides criteria:**  
Manual inspection confirmed that `dyslipidaemia == 1` in this cohort largely reflects high LDL. Several patients with `dyslipidaemia == 1` had normal HDL and triglycerides. Using this flag would add up to 2 false MetS points per patient.

**Why `hypertension == 1` IS used for the hypertension criterion:**  
Cross-checking with HYP drug columns confirmed 0 cases where a patient was on antihypertensive therapy but had `hypertension == 0`. The binary flag is specific and complete.

## Computational Environment & Reproducibility

Package management is handled via **`renv`**. To restore the exact environment used for this analysis:

```r
install.packages("renv")
renv::restore()   # installs all packages at the versions recorded in renv.lock
```

All scripts use absolute paths anchored to the project root (`/home/lorenzo/Microbesomics/`). If running on a different machine, update `BASE_DIR` at the top of each script.

**Core R packages**: `phyloseq`, `DESeq2`, `vegan`, `ggplot2`, `ComplexUpset`, `microbiomeMarker`, `openxlsx`, `readxl`, `dplyr`

**Python** (called via `reticulate`): `pandas`, `numpy`, `statsmodels`, `pingouin` — environment managed separately via `envs/`.

## Repository Structure

```
Microbesomics/
├── scripts/                          # Active analysis scripts (run in order)
├── data/
│   ├── raw/                          # Raw input files (not tracked by git)
│   │   ├── Microobesomics_clinical data_dec2022.xlsx
│   │   └── qza/                      # QIIME2 artifacts
│   └── processed/                    # Clean output datasets
│       └── bariatric_clinical.csv
├── results/                          # Figures and tables (not tracked by git)
└── _old/                             # Archived scripts from previous analysis
    └── R-scripts/
```

> Raw data and results are excluded from version control (patient data). Processed outputs in `data/processed/` are included as they contain only derived/aggregated variables.

## Scripts

### `scripts/01_import_qiime2.R`

**Purpose**: Imports QIIME2 artifacts (QZA files) and builds a bare phyloseq object with 247 samples and 14741 ASVs. No sample metadata is attached at this stage. Feature filtering is NOT applied here — it is handled group-specifically by `R/filter_phyloseq.R` in each analysis script.

**Input**: `data/raw/qza/` — `asv_table.qza`, `taxonomy.qza`, `fasttree_tree_rooted.qza`, `rep_seq.qza`

**Output**:
- `data/processed/ps_raw.rds` — unfiltered phyloseq object (247 samples, 14741 ASVs)

**How to run**:
```r
Rscript scripts/01_import_qiime2.R
```

**Note on sample naming**: QZA internal IDs (`ID2988-N-N-boxK`) are renamed to `S1`..`S247` by extracting N via regex. These correspond to CODE IGA in all clinical metadata files.

---

### `scripts/02_preprocess_clinical.R`

**Purpose**: Reads all four clinical sheets from the IGA Excel file, normalises column names, and stacks them into a single dataset with one row per sample. Columns present in some sheets but not others are filled with NA.

**Input**: `data/raw/Copia di Microobesomics_clinical data IGA 260324.xlsx`
- Sheet `group1vs2` — dietary intervention, pre/post
- Sheet `group3vs4` — Mingrone cohort, pre/post
- Sheet `group5vs6` — bariatric surgery, pre/post

> `clinical data merged - BASAL` is excluded: it contains only the baseline timepoint for the same patients already present in the group sheets, and would introduce duplicates.

**Output**:
- `data/processed/clinical_all.csv` — 247 rows × 169 columns, one row per QZA sample
- `data/processed/clinical_all.xlsx` — same, single sheet

**How to run**:
```r
Rscript scripts/02_preprocess_clinical.R
```

**Column name harmonisation**: The four sheets use inconsistent naming for the same variables. All names are lowercased; known synonyms are unified:

| Raw name (some sheets) | Unified name |
|------------------------|--------------|
| `DIABETE/IFG` | `diabetes` |
| `OSAS` | `saos` |
| `Weight pre` | `weight` |
| `BMI pre` | `bmi` |

Variables unique to one sheet (e.g. BIA body composition and PREDIMED diet score in `group1vs2`) are kept as-is with NA in all other sheets.

**Sample ID**: `sample_id` column is added as `"S" + CODE IGA`, linking each clinical row to the phyloseq object (`S1`..`S247`).

---


### `R/filter_phyloseq.R`

Reusable filtering function sourced by all group-specific scripts. Applies MicrobiomeAnalyst-equivalent filters to a phyloseq object already subsetted to the samples of interest:

1. **Low count filter**: ASV must have ≥ `min_count` reads (default 4) in ≥ `min_prev_frac` of samples (default 10%)
2. **Low variance filter**: removes bottom `var_pct`% by IQR on raw counts; set `var_pct = 0` (default) to skip — recommended for DESeq2 workflows where independent filtering is applied internally

Filtering within each group subset (rather than globally) avoids excluding ASVs that are relevant for a specific comparison but rare in the rest of the cohort.

---

### `scripts/03_filter_5vs6.R`

**Purpose**: Subsets the raw phyloseq to bariatric surgery samples (groups 5 and 6), attaches clinical metadata, and applies feature filtering via `R/filter_phyloseq.R`.

**Input**:
- `data/processed/ps_raw.rds`
- `data/processed/clinical_all.csv`

**Output**: `data/processed/ps_filt_5vs6.rds` — 66 samples (33 complete pairs), 707 ASVs

**How to run**:
```r
Rscript scripts/03_filter_5vs6.R
```

**Filters applied**: ≥4 reads in ≥10% of 66 samples (= 7 samples); bottom 10% by IQR removed (variance filter).

---

### `scripts/04_deseq2_5vs6.R`

**Purpose**: Differential abundance analysis (DESeq2), POST (group 6) vs PRE (group 5) bariatric surgery. Paired design with patient as blocking factor.

**Input**: `data/processed/ps_filt_5vs6.rds` — 66 samples, 707 ASVs

**Output**:
- `results/deseq2/deseq2_5vs6_results.xlsx` — full results table with taxonomy
- `results/deseq2/volcano_5vs6.png` — volcano plot (top 25 ASVs labelled)

**How to run**:
```r
Rscript scripts/04_deseq2_5vs6.R
```

**Design**: `~ pz + groupIGA` — pz (patient ID) blocks between-patient variance; contrast is groupIGA 6 vs 5 (positive log2FC = enriched POST-surgery).

**Results**: 255 / 707 ASVs significant at FDR < 0.05 (162 enriched POST, 93 enriched PRE).

---

### `scripts/05_deseq2_remission.R`

**Purpose**: Differential abundance analysis (DESeq2), MetS remitters vs non-remitters at baseline (PRE, group 5).

**Input**: `data/processed/ps_filt_5vs6.rds` — 66 samples, 707 ASVs (subset to 28 PRE MetS+ samples at runtime)

**Output**:
- `results/deseq2/deseq2_remission_results.xlsx`
- `results/deseq2/volcano_remission.png`

**How to run**:
```r
Rscript scripts/05_deseq2_remission.R
```

**Design**: `~ remission` — remitter (MetS+ PRE → MetS− POST) vs non_remitter (MetS+ PRE → MetS+ POST); positive log2FC = enriched in remitters.

**Results**: 59 / 707 ASVs significant at FDR < 0.05 (18 enriched in remitters, 41 in non-remitters).

---

### `scripts/06_deseq2_mets_baseline.R`

**Purpose**: Differential abundance analysis (DESeq2), MetS+ vs MetS− at baseline (PRE, group 5) among the 33 complete pairs.

**Input**: `data/processed/ps_filt_5vs6.rds` — 66 samples, 707 ASVs (subset to 33 PRE samples at runtime)

**Output**:
- `results/deseq2/deseq2_mets_baseline_results.xlsx`
- `results/deseq2/volcano_mets_baseline.png`

**How to run**:
```r
Rscript scripts/06_deseq2_mets_baseline.R
```

**Design**: `~ MetS_group` — MetS_pos vs MetS_neg (reference); positive log2FC = enriched in MetS+.

**Results**: 60 / 707 ASVs significant at FDR < 0.05 (27 enriched in MetS+, 33 enriched in MetS−).
