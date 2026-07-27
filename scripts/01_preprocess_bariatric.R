# =============================================================================
# 01_preprocess_bariatric.R
#
# Extracts the 40 pre/post bariatric surgery patient pairs from the raw
# clinical Excel file, computes MetS criteria per NCEP ATP III, and writes
# a clean CSV to data/processed/.
#
# Input:  data/raw/Microobesomics_clinical data_dec2022.xlsx
# Output: data/processed/bariatric_clinical.csv
#
# MetS definition (NCEP ATP III, AHA/NHLBI 2005 revision): >= 3 of 5 criteria
#   1. Abdominal Obesity    : waist >= 102 cm (M) / >= 88 cm (F)
#   2. Hyperglycemia        : fasting glucose >= 100 mg/dL OR diabetes == 1
#   3. Hypertension         : SBP >= 130 OR DBP >= 85 OR hypertension == 1
#   4. Low HDL              : HDL < 40 mg/dL (M) / < 50 mg/dL (F)
#   5. Hypertriglyceridemia : triglycerides >= 150 mg/dL
#
# All thresholds use >= / <= per the 2005 AHA/NHLBI harmonised revision
# (original ATP III 2001 used strict > for waist; revision aligned to >=).
#
# Why no drug columns or dyslipidaemia flag for criteria 4 and 5:
#   - Drug columns in the bariatric cohort contain free-text drug names.
#     Inspection revealed that among the 6 patients on lipid-lowering therapy,
#     4 are on statins or ezetimibe (LDL-targeted, not relevant for MetS HDL/TG
#     criteria), 1 on fenofibrate and 1 on omega-3 ethyl esters. These 2 cases
#     with TG/HDL-specific therapy are covered by their measured lab values,
#     which are available for all 40 pairs — so no proxy is needed.
#   - dyslipidaemia==1 is a generic diagnosis that in this cohort primarily
#     reflects high LDL (hypercholesterolaemia). Manual inspection confirmed
#     that several patients with dyslipidaemia==1 had normal HDL and
#     triglycerides, making the flag too unspecific for criteria 4 and 5.
#
# Why hypertension==1 IS used for criterion 3:
#   - hypertension==1 is diagnosis-specific (unlike dyslipidaemia).
#   - Cross-checking with HYP drug columns showed 0 cases where a patient
#     was on antihypertensive drugs but had hypertension==0, confirming the
#     binary flag fully captures treated patients without adding false positives.
# =============================================================================

# --- Environment --------------------------------------------------------------

BASE_DIR <- "/home/lorenzo/Microbesomics"

# renv::restore() will reinstall all packages at the exact versions in renv.lock
# Run this once after cloning on a new machine:
# renv::restore()

library(readxl)
library(dplyr)

# --- Load raw data ------------------------------------------------------------

excel_path <- file.path(BASE_DIR, "data/raw/Microobesomics_clinical data_dec2022.xlsx")

# The merged sheet has TWO header rows:
#   Row 1: machine-readable column names (used as colnames)
#   Row 2: human-readable labels (skipped)
#   Row 3+: patient data
raw <- read_excel(excel_path, sheet = "clinical data merged", col_names = FALSE)

col_names <- unlist(raw[1, ])   # grab row 1 as column names
data_rows <- raw[3:nrow(raw), ] # skip both header rows
colnames(data_rows) <- col_names

# --- Identify bariatric patients ----------------------------------------------

# Bariatric patients have codes in the format "XX-NNN" (e.g. "CA-89", "DQ-64").
# Their follow-up rows have no code; the patient column contains "N FU" (e.g. "12 FU").
is_bar_pre  <- grepl("-", data_rows[["code"]], fixed = TRUE) & !is.na(data_rows[["code"]])
is_bar_post <- grepl("FU", data_rows[["patient"]])

bar_pre  <- data_rows[is_bar_pre,  ]
bar_post <- data_rows[is_bar_post, ]

stopifnot(nrow(bar_pre) == 40, nrow(bar_post) == 40)

# Extract numeric patient ID from both sets:
#   PRE  rows: patient column contains the bare number ("1", "2", ...)
#   POST rows: patient column contains "N FU" — strip " FU" to get the number
bar_pre[["patient_id"]]  <- as.integer(bar_pre[["patient"]])
bar_post[["patient_id"]] <- as.integer(gsub("\\s*FU.*", "", bar_post[["patient"]]))

# Verify one-to-one correspondence — stop with an informative error if any ID
# is present in one set but not the other
missing_in_post <- setdiff(bar_pre[["patient_id"]], bar_post[["patient_id"]])
missing_in_pre  <- setdiff(bar_post[["patient_id"]], bar_pre[["patient_id"]])
if (length(missing_in_post) > 0)
  stop("PRE patient IDs with no matching POST: ", paste(missing_in_post, collapse = ", "))
if (length(missing_in_pre) > 0)
  stop("POST patient IDs with no matching PRE: ", paste(missing_in_pre, collapse = ", "))

# Sort both by patient_id before copying the code, so the assignment is
# guaranteed to be correct regardless of row order in the Excel sheet
bar_pre  <- bar_pre[order(bar_pre[["patient_id"]]),  ]
bar_post <- bar_post[order(bar_post[["patient_id"]]), ]

# Copy the patient code (XX-NNN) into POST rows, which have NA in that column
bar_post[["code"]] <- bar_pre[["code"]]

# Add timepoint labels
bar_pre[["timepoint"]]  <- "PRE"
bar_post[["timepoint"]] <- "POST"

# Stack into a single 80-row data frame, sorted by patient_id then timepoint
bar <- bind_rows(bar_pre, bar_post)
bar <- bar[order(bar[["patient_id"]], bar[["timepoint"]]), ]

# --- Select and coerce relevant columns --------------------------------------

numeric_cols <- c(
  "age", "bmi", "weight", "height", "waist", "hip",
  "sbp", "dbp",
  "glucose 0'", "HOMA", "hba1c",
  "HDL", "LDL", "triglycerides", "totale cholesterol",
  "uric acid", "hsCRP"
)

# dyslipidaemia is kept as a descriptive variable only — it is NOT used in
# MetS criteria (see header note above)
binary_cols <- c("sex", "diabetes", "hypertension", "dyslipidaemia", "Nafld", "saos")

keep_cols <- c("patient", "code", "timepoint", "VISIT DATE", binary_cols, numeric_cols)

bar <- bar[, keep_cols]

# Coerce numeric columns (may have been read as character due to Excel mixed types)
bar[numeric_cols] <- lapply(bar[numeric_cols], as.numeric)

# Coerce binary columns to integer (0/1)
bar[binary_cols] <- lapply(bar[binary_cols], function(x) as.integer(as.character(x)))

# --- Compute MetS criteria ---------------------------------------------------

# 1. Abdominal Obesity — sex-specific waist thresholds per AHA/NHLBI 2005 revision
bar$AbdominalObesityMetS <- case_when(
  bar$sex == 0 & bar$waist >= 102 ~ 1L,
  bar$sex == 0 & bar$waist <  102 ~ 0L,
  bar$sex == 1 & bar$waist >= 88  ~ 1L,
  bar$sex == 1 & bar$waist <  88  ~ 0L,
  TRUE ~ NA_integer_
)

# 2. Hyperglycemia
# Using the 2004 AHA/NHLBI revision threshold (>= 100 mg/dL vs original 110).
# diabetes==1 alone qualifies regardless of measured glucose.
bar$HyperglicemiaMetS <- case_when(
  bar[["glucose 0'"]] >= 100 | bar$diabetes == 1 ~ 1L,
  bar[["glucose 0'"]] <  100 & bar$diabetes == 0 ~ 0L,
  is.na(bar[["glucose 0'"]]) & bar$diabetes == 1 ~ 1L,
  is.na(bar[["glucose 0'"]]) & bar$diabetes == 0 ~ 0L,
  TRUE ~ NA_integer_
)

# 3. Hypertension
# hypertension==1 captures both diagnosed-untreated and treated patients.
# Cross-check confirmed that all patients on antihypertensive drugs also have
# hypertension==1, so the drug columns add no additional cases.
bar$HypertensionMetS <- case_when(
  bar$sbp >= 130 | bar$dbp >= 85 | bar$hypertension == 1 ~ 1L,
  bar$sbp <  130 & bar$dbp <  85 & bar$hypertension == 0 ~ 0L,
  TRUE ~ NA_integer_
)

# 4. Low HDL — sex-specific thresholds per NCEP ATP III, lab values only
bar$HDL_MetS <- case_when(
  bar$sex == 0 & bar$HDL <  40 ~ 1L,
  bar$sex == 0 & bar$HDL >= 40 ~ 0L,
  bar$sex == 1 & bar$HDL <  50 ~ 1L,
  bar$sex == 1 & bar$HDL >= 50 ~ 0L,
  TRUE ~ NA_integer_
)

# 5. Hypertriglyceridemia — lab values only
bar$TrygliceridesMetS <- case_when(
  bar$triglycerides >= 150 ~ 1L,
  bar$triglycerides <  150 ~ 0L,
  TRUE ~ NA_integer_
)

# --- Compute MetS composite score --------------------------------------------

mets_criteria <- c("AbdominalObesityMetS", "HyperglicemiaMetS",
                   "HypertensionMetS", "HDL_MetS", "TrygliceridesMetS")

# Count confirmed-positive criteria and missing criteria per patient
n_pos <- rowSums(bar[mets_criteria] == 1, na.rm = TRUE)
n_na  <- rowSums(is.na(bar[mets_criteria]))

# MetS_score: sum of confirmed-positive criteria (NA criteria excluded)
bar$MetS_score <- n_pos

# Determine MetS status with three-way logic:
#   - Certain YES : confirmed positives alone reach the threshold (>= 3)
#   - Certain NO  : even if all NAs were positive, threshold cannot be reached
#   - Uncertain   : NA (outcome depends on the missing values)
bar$MetS <- case_when(
  n_pos >= 3             ~ 1L,   # certain MetS positive
  n_pos + n_na < 3      ~ 0L,   # certain MetS negative
  TRUE                   ~ NA_integer_  # uncertain — missing data is decisive
)

# Flag rows where the MetS call is uncertain due to missing criteria
bar$MetS_incomplete <- as.integer(is.na(bar$MetS))

# --- Reorder columns for readability -----------------------------------------

id_cols   <- c("patient", "code", "timepoint", "VISIT DATE")
demo_cols <- c("age", "sex", "bmi", "weight", "height", "waist", "hip")
bp_cols   <- c("sbp", "dbp")
lab_cols  <- c("glucose 0'", "HOMA", "hba1c", "HDL", "LDL",
               "triglycerides", "totale cholesterol", "uric acid", "hsCRP")
comorbid  <- c("diabetes", "hypertension", "dyslipidaemia", "Nafld", "saos")
mets_cols <- c(mets_criteria, "MetS_score", "MetS", "MetS_incomplete")

bar <- bar[, c(id_cols, demo_cols, bp_cols, lab_cols, comorbid, mets_cols)]

# --- Write output ------------------------------------------------------------

out_path <- file.path(BASE_DIR, "data/processed/bariatric_clinical.csv")
write.csv(bar, file = out_path, row.names = FALSE, na = "")

message("Written: ", out_path)
message("Rows: ", nrow(bar), " | Cols: ", ncol(bar))

# Quick sanity check: MetS prevalence pre vs post
cat("\n--- MetS prevalence by timepoint ---\n")
print(table(Timepoint = bar$timepoint, MetS = bar$MetS, useNA = "ifany"))
