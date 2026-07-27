# =============================================================================
# 02_preprocess_clinical.R
#
# Reads all four clinical sheets from the IGA Excel file, normalises column
# names, stacks them into a single dataset, and computes MetS criteria
# (NCEP ATP III, AHA/NHLBI 2005 revision) for all 247 samples.
#
# Input:  data/raw/Copia di Microobesomics_clinical data IGA 260324.xlsx
#   Sheets used:
#     "group1vs2"  — dietary intervention, pre/post
#     "group3vs4"  — Mingrone cohort, pre/post
#     "group5vs6"  — bariatric surgery, pre/post
#
#   "clinical data merged - BASAL" is intentionally excluded: it contains the
#   same patients as the group sheets (baseline timepoint only) and would
#   introduce duplicates.
#
# Output: data/processed/
#   clinical_all.csv   — stacked dataset with MetS columns, one row per sample
#   clinical_all.xlsx  — same, single sheet
#
# Column name strategy:
#   All names are lowercased.  Known synonyms across sheets are unified:
#     DIABETE/IFG  → diabetes
#     OSAS         → saos
#     Weight pre   → weight     (group1vs2 uses "pre" suffix for baseline weight)
#     BMI pre      → bmi
#   Everything else is kept as-is after lowercasing.  Variables present in
#   some sheets but not others are filled with NA (bind_rows behaviour).
#
# Sample ID:
#   CODE IGA is the integer key linking each clinical row to the phyloseq
#   object (sample names S1..S247).  A sample_id column ("S" + CODE IGA) is
#   added as the primary identifier.
#
# MetS criteria (NCEP ATP III, AHA/NHLBI 2005 revision — >= 3 of 5):
#   1. AbdominalObesityMetS : waist >= 102 cm (M) / >= 88 cm (F)
#   2. HyperglicemiaMetS    : glucose 0' >= 100 mg/dL OR diabetes == 1
#   3. HypertensionMetS     : sbp >= 130 OR dbp >= 85 OR hypertension == 1
#   4. HDL_MetS             : HDL < 40 mg/dL (M) / < 50 mg/dL (F)
#   5. TrygliceridesMetS    : triglycerides >= 150 mg/dL OR TG_pharmacotherapy == 1
#
# TG_pharmacotherapy:
#   3 patients are on TG-specific drugs (fenofibrate or omega-3) at specific
#   timepoints; their measured TG may be artificially reduced by treatment.
#   These rows are flagged with TG_pharmacotherapy = 1 and count as positive
#   for criterion 5 regardless of their lab value.
#     S200 (CO-71,  group 5 PRE)  — Fulcosupra (fenofibrate)
#     S215 (DA-157, group 6 POST) — Eskimo (omega-3 ethyl esters)
#     S85  (C19,    group 3 PRE)  — Fenofibrate
#
# Classification uses three-way logic:
#   Certain YES (>= 3 confirmed positive) / Certain NO / NA (missing decisive)
# =============================================================================

BASE_DIR  <- "/home/lorenzo/Microbesomics"
IGA_FILE  <- file.path(BASE_DIR, "data/raw/Copia di Microobesomics_clinical data IGA 260324.xlsx")
PROC_DIR  <- file.path(BASE_DIR, "data/processed")

library(readxl)
library(dplyr)
library(openxlsx)

# --- Column name normalisation -----------------------------------------------

# Synonym map: raw lowercased name → unified name.
# Add entries here whenever a new synonym is discovered.
SYNONYMS <- c(
  "diabete/ifg"  = "diabetes",
  "osas"         = "saos",
  "weight pre"   = "weight",   # group1vs2 uses "pre" suffix for baseline
  "bmi pre"      = "bmi"
)

normalise_names <- function(df) {
  nms <- tolower(colnames(df))
  nms <- ifelse(nms %in% names(SYNONYMS), SYNONYMS[nms], nms)
  colnames(df) <- nms
  df
}

# --- Read and clean each sheet -----------------------------------------------

read_sheet <- function(sheet_name) {
  dt <- read_excel(IGA_FILE, sheet = sheet_name)
  dt <- normalise_names(dt)

  # Drop any footer/summary rows where CODE IGA is not a parseable integer
  # (e.g. repeated header rows or totals rows that appear at the bottom)
  dt <- dt[!is.na(suppressWarnings(as.integer(dt[["code iga"]]))), ]

  # Coerce CODE IGA to integer
  dt[["code iga"]] <- as.integer(dt[["code iga"]])

  # Add sample_id as primary key linking to phyloseq sample names S1..S247
  dt[["sample_id"]] <- paste0("S", dt[["code iga"]])

  # Move sample_id, code iga, group iga, pz to the front
  id_cols <- c("sample_id", "code iga", "group iga", "pz")
  other   <- setdiff(colnames(dt), id_cols)
  dt[, c(id_cols, other)]
}

sheets <- list(
  "group1vs2",
  "group3vs4",
  "group5vs6"
)

sheet_list <- lapply(sheets, read_sheet)
names(sheet_list) <- c("1vs2", "3vs4", "5vs6")

cat("Rows per sheet:\n")
for (nm in names(sheet_list))
  cat(" ", nm, ":", nrow(sheet_list[[nm]]), "\n")

# --- Stack all sheets --------------------------------------------------------

# bind_rows fills missing columns with NA automatically
clinical_all <- bind_rows(sheet_list)

cat("\nTotal rows in merged dataset:", nrow(clinical_all), "\n")
cat("Total columns:", ncol(clinical_all), "\n")
cat("GROUP IGA distribution:\n")
print(table(clinical_all[["group iga"]], useNA = "ifany"))

# --- TG pharmacotherapy flag -------------------------------------------------

# Patients on TG-specific drugs at specific timepoints only.
# Flagged per-row because drug use does not always span both PRE and POST.
TG_DRUG_SAMPLES <- c("S200", "S215", "S85")
clinical_all[["TG_pharmacotherapy"]] <- as.integer(
  clinical_all[["sample_id"]] %in% TG_DRUG_SAMPLES
)
cat("\nTG_pharmacotherapy == 1:\n")
print(clinical_all[clinical_all[["TG_pharmacotherapy"]] == 1,
                   c("sample_id", "pz", "group iga")])

# --- MetS criteria (NCEP ATP III, AHA/NHLBI 2005 revision) ------------------

# Coerce variables used in criteria to numeric
num_vars <- c("sex", "waist", "glucose 0'", "diabetes",
              "sbp", "dbp", "hypertension", "hdl", "triglycerides")
clinical_all[num_vars] <- lapply(clinical_all[num_vars], as.numeric)

# 1. Abdominal Obesity
clinical_all$AbdominalObesityMetS <- case_when(
  clinical_all$sex == 0 & clinical_all$waist >= 102 ~ 1L,
  clinical_all$sex == 0 & clinical_all$waist <  102 ~ 0L,
  clinical_all$sex == 1 & clinical_all$waist >= 88  ~ 1L,
  clinical_all$sex == 1 & clinical_all$waist <  88  ~ 0L,
  TRUE ~ NA_integer_
)

# 2. Hyperglycemia
clinical_all$HyperglicemiaMetS <- case_when(
  clinical_all[["glucose 0'"]] >= 100 | clinical_all$diabetes == 1 ~ 1L,
  clinical_all[["glucose 0'"]] <  100 & clinical_all$diabetes == 0 ~ 0L,
  is.na(clinical_all[["glucose 0'"]]) & clinical_all$diabetes == 1 ~ 1L,
  is.na(clinical_all[["glucose 0'"]]) & clinical_all$diabetes == 0 ~ 0L,
  TRUE ~ NA_integer_
)

# 3. Hypertension — hypertension == 1 captures all treated patients
clinical_all$HypertensionMetS <- case_when(
  clinical_all$sbp >= 130 | clinical_all$dbp >= 85 | clinical_all$hypertension == 1 ~ 1L,
  clinical_all$sbp <  130 & clinical_all$dbp <  85 & clinical_all$hypertension == 0 ~ 0L,
  TRUE ~ NA_integer_
)

# 4. Low HDL
clinical_all$HDL_MetS <- case_when(
  clinical_all$sex == 0 & clinical_all$hdl <  40 ~ 1L,
  clinical_all$sex == 0 & clinical_all$hdl >= 40 ~ 0L,
  clinical_all$sex == 1 & clinical_all$hdl <  50 ~ 1L,
  clinical_all$sex == 1 & clinical_all$hdl >= 50 ~ 0L,
  TRUE ~ NA_integer_
)

# 5. Hypertriglyceridemia — lab value OR on TG-specific pharmacotherapy
clinical_all$TrygliceridesMetS <- case_when(
  clinical_all$TG_pharmacotherapy == 1              ~ 1L,
  clinical_all$triglycerides >= 150                 ~ 1L,
  clinical_all$triglycerides <  150                 ~ 0L,
  TRUE ~ NA_integer_
)

# --- MetS composite score ----------------------------------------------------

mets_criteria <- c("AbdominalObesityMetS", "HyperglicemiaMetS",
                   "HypertensionMetS", "HDL_MetS", "TrygliceridesMetS")

n_pos <- rowSums(clinical_all[mets_criteria] == 1, na.rm = TRUE)
n_na  <- rowSums(is.na(clinical_all[mets_criteria]))

clinical_all$MetS_score <- n_pos
clinical_all$MetS <- case_when(
  n_pos >= 3        ~ 1L,
  n_pos + n_na < 3  ~ 0L,
  TRUE              ~ NA_integer_
)
clinical_all$MetS_incomplete <- as.integer(is.na(clinical_all$MetS))

cat("\nMetS distribution by GROUP IGA:\n")
print(table(GROUP = clinical_all[["group iga"]], MetS = clinical_all$MetS, useNA = "ifany"))

# --- Write outputs -----------------------------------------------------------

write.csv(clinical_all,
          file = file.path(PROC_DIR, "clinical_all.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(PROC_DIR, "clinical_all.csv"))

write.xlsx(clinical_all,
           file = file.path(PROC_DIR, "clinical_all.xlsx"),
           overwrite = TRUE)
message("Written: ", file.path(PROC_DIR, "clinical_all.xlsx"))
