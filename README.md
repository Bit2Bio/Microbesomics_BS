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
Inspection of all lipid-lowering drug names in the dataset shows they are predominantly statins and ezetimibe, which target LDL and are not relevant for the HDL or triglycerides MetS criteria. Three patients across the cohort are on TG-specific therapy (fenofibrate or omega-3 ethyl esters) at specific timepoints; these are captured via the `TG_pharmacotherapy` flag and counted as positive for the TG criterion regardless of their measured lab value (see script 02).

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

**Purpose**: Imports QIIME2 artifacts (QZA files) and builds a bare phyloseq object with 247 samples and 14741 ASVs. No sample metadata is attached at this stage.

**Input**: `data/raw/qza/` — `asv_table.qza`, `taxonomy.qza`, `fasttree_tree_rooted.qza`, `rep_seq.qza`

**Output**:
- `data/processed/ps_raw.rds` — phyloseq object (247 samples, 14741 ASVs, 7 taxonomic ranks)
- `data/processed/otu_table.csv` — ASV count matrix
- `data/processed/taxonomy.csv` — per-ASV taxonomic assignments
- `data/processed/ps_raw.xlsx` — same tables as separate sheets

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

### `scripts/03_preprocess_bariatric.R`

**Purpose**: Extracts the 40 pre/post bariatric surgery patient pairs from the raw clinical Excel file, computes MetS criteria per NCEP ATP III, and writes a clean CSV.

**Input**: `data/raw/Microobesomics_clinical data_dec2022.xlsx` — sheet `clinical data merged`

**Output**: `data/processed/bariatric_clinical.csv` — 80 rows (40 PRE + 40 POST), 35 columns

**How to run**:
```r
Rscript scripts/03_preprocess_bariatric.R
```

**Key decisions documented in script comments**:
- Pre/post pairing verified explicitly by patient ID extracted from the `"N FU"` pattern, not by row order
- Drug columns excluded from HDL and triglycerides criteria (predominantly statins targeting LDL)
- `dyslipidaemia == 1` excluded from lipid criteria (generic diagnosis, confirmed to reflect mainly high LDL in this cohort)
- MetS classification uses three-way logic: certain YES / certain NO / NA — avoids propagating NA when outcome is already determined by observed criteria

**Results**:

| Timepoint | MetS− | MetS+ |
|-----------|-------|-------|
| PRE | 8 (20%) | 32 (80%) |
| POST | 34 (85%) | 6 (15%) |

---

### `scripts/04_filter_bariatric_paired.R`

**Purpose**: Filters the bariatric clinical dataset to retain only patients with metagenomics data at **both** timepoints (PRE and POST). Of the 40 pairs, 7 are missing at least one QZA sample; this script removes them to produce a clean paired set for downstream omics analyses.

**Input**:
- `data/processed/bariatric_clinical.csv` — script 03 output (80 rows)
- `data/processed/clinical_all.csv` — script 02 output (QZA-linked samples)

**Output**: `data/processed/bariatric_clinical_paired.csv` — 64 rows (33 complete pairs)

**How to run**:
```r
Rscript scripts/04_filter_bariatric_paired.R
```

**Results (33 complete pairs, IGA file)**:

| Timepoint | MetS− | MetS+ |
|-----------|-------|-------|
| PRE | 5 | 28 |
| POST | 27 | 6 |

Of the 28 PRE MetS+ patients: **22 achieve remission** (MetS− at POST), **6 do not**. The key comparison for understanding remission drivers is **22 remitters vs 6 non-remitters at baseline**. 5 patients were MetS− at PRE and remained MetS− at POST; none worsened.
