# Microbesomics_MetS

Multi-omics analysis of gut microbiome signatures associated with Metabolic Syndrome (MetS) and its individual components, including pre/post intervention comparisons (dietary intervention and bariatric surgery) and external cohort validation.

## Study Design

- **Main cohort**: n=177 cross-sectional samples, MetS classification based on NCEP ATP III criteria (≥3 of 5)
- **Diet intervention**: paired pre/post comparison (Groups 1 vs 2)
- **Bariatric surgery**: paired pre/post comparison (Groups 5 vs 6), n=40 pairs
- **External validation**: Floromidia cohort (independent replication)

### Bariatric Surgery — MetS Prevalence

| Timepoint | MetS− | MetS+ |
|-----------|-------|-------|
| PRE (before surgery) | 8 (20%) | 32 (80%) |
| POST (after surgery) | 34 (85%) | 6 (15%) |

85% of subjects achieved MetS remission after bariatric surgery. The 6 subjects remaining MetS+ post-surgery represent partial responders; the group is too small (n=6) for independent multi-omics analysis.

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
Inspection of all lipid-lowering drug names in the dataset shows they are predominantly statins and ezetimibe, which target LDL and are not relevant for the HDL or triglycerides MetS criteria. Only 2 bariatric patients use TG/HDL-specific drugs (1 fenofibrate, 1 omega-3 ethyl esters); both are covered by their measured lab values, which are complete for all 40 pairs.

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
├── R-scripts/          # Analysis scripts (run in numerical order)
├── data/
│   ├── clinics/        # Clinical metadata (.rds, .xlsx)
│   ├── metagenomics/   # Phyloseq objects (.rds) and QIIME2 inputs
│   ├── metabolomics/   # Metabolomic data (.rds)
│   └── ValidazioneFloromidia/  # External cohort data
└── results/            # Output figures and tables (organized by analysis)
```

## Pipeline

Scripts must be run in the order below. Steps 0.x generate `.rds` objects used by all downstream scripts.

### 0. Data Import *(run once to generate .rds objects)*

| Script | Output |
|---|---|
| `0_Import.R` | `ps_filt.rds`, `sampleDt.rds`, `mtbDt.rds` |
| `0_ImportPrePostDiet.R` | `ps_filt1vs2.rds`, `sampleDt1vs2.rds`, `mtbDt1vs2.rds` |
| `0_ImportPrePostBariatricSurgery.R` | `ps_filt5vs6.rds`, `sampleDt5vs6.rds`, `mtbDt5vs6.rds` |

### 1. Clinical Analysis

| Script | Description |
|---|---|
| `1_MetS_RiskFactors.Rmd` | Define MetS and component variables from clinical criteria |
| `2_0_AnalisiDescrittivaClinica.Rmd` | Descriptive statistics of clinical variables |

### 2. Metagenomic Analysis (Main Cohort)

| Script | Description |
|---|---|
| `2_1_Barplots.Rmd` | Taxonomic composition (Phylum → Family level) |
| `2_2_DEseq2.Rmd` | Differential abundance by MetS component (DESeq2) |
| `3_InterSections_upsetplot_Heatmap.Rmd` | Shared ASVs across MetS components (UpSet plot) |
| `5_BetaDiversity.Rmd` | Community structure — Bray-Curtis PCoA + ADONIS2 |
| `6_AlphaDiversity.Rmd` | Within-sample diversity (Observed, Fisher, ACE) |

### 3. Metabolomics & Cross-Omics

| Script | Description |
|---|---|
| `4_0_Metabolomics_CrossOmics.Rmd` | Metabolomic PCA, ADONIS2 by MetS component |
| `4_1_Metabolomics_CrossOmics1_vs_2.Rmd` | Metabolomic changes — diet intervention |
| `4_2_Metabolomics_CrossOmics5_vs_6.Rmd` | Metabolomic changes — bariatric surgery |

### 4. Intervention Analyses

| Script | Description |
|---|---|
| `7_1_DEseq2_1vs2_Diet.Rmd` | Microbiome changes — diet intervention (paired DESeq2) |
| `7_2_DEseq2_5vs6_Bariatric_Surgery.Rmd` | Microbiome changes — bariatric surgery (paired DESeq2) |
| `7_3_prevotella_9.Rmd` | Targeted analysis of *Prevotella copri* (Prevotella_9) |

### 5. Functional Prediction

| Script | Description |
|---|---|
| `8_Picrust.Rmd` | PICRUSt2 pathway inference — main cohort |
| `8_Picrust_BS.Rmd` | PICRUSt2 pathway inference — bariatric surgery |

### 6. External Validation

| Script | Description |
|---|---|
| `ValidazioneFloromidia.Rmd` | Replication in Floromidia cohort |
| `ValidazioneFloromidia_Updated.Rmd` | Updated validation analysis |
| `ValidazioneFloromidiaCrossOmics.Rmd` | Cross-omics validation (metagenomics + metabolomics + transcriptomics) |

## Key Outputs

| Analysis | Output file |
|---|---|
| Differential abundance (MetS) | `results/metagenomics/DEseq2/AnalisiComparativaMetagenomica.xlsx` |
| Intersection biomarkers | `results/metagenomics/intersections/IntersezioniUpsetPlotMicrobesomics.xlsx` |
| Functional pathways | `results/Picrust2_MetS_vs_BS.xlsx` |
| Clinical description | `results/clinics/Descrittiva Clinica.xlsx` |
| MetS classification | `results/MetS_RiskFactors/MetS_RisksFactors.xlsx` |
