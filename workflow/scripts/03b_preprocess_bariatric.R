# =============================================================================
# 03b_preprocess_bariatric.R
#
# Extracts the 66 samples (33 complete PRE/POST pairs) from clinical_all.csv
# and produces a clean bariatric-specific clinical file with all available
# variables. This is the canonical clinical dataset for all downstream
# bariatric surgery analyses (Tabella 1, clinical correlations, remission
# classification, etc.).
#
# "Complete pair" means the patient has metagenomics data at BOTH group 5
# (PRE) and group 6 (POST) — same criterion used in 03_filter_5vs6.R to
# build ps_filt_5vs6.rds.
#
# Input:
#   data/processed/clinical_all.csv — 247 rows, all groups
#
# Output:
#   data/processed/bariatric_clinical_66.csv  — 66 rows (long format, PRE+POST)
#   data/processed/bariatric_clinical_66.xlsx — same, Excel
#
# Format: one row per sample (long format). Columns include:
#   - Identifiers : sample_id, pz, group iga, timepoint (PRE / POST)
#   - All numeric clinical variables available in group 5 / 6
#   - MetS criteria and composite score
#   - Remission label (assigned at patient level from PRE→POST MetS change)
#
# Remission classification (PRE samples only):
#   remitter     : MetS == 1 at PRE AND MetS == 0 at POST
#   non_remitter : MetS == 1 at PRE AND MetS == 1 at POST
#   mets_neg     : MetS == 0 at PRE (already metabolically healthy at baseline)
#   NA           : MetS status undetermined at PRE or POST
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"
PROC_DIR <- file.path(BASE_DIR, "data/processed")

library(dplyr)
library(openxlsx)

clin <- read.csv(file.path(PROC_DIR, "clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

# --- Subset to bariatric surgery groups (5 = PRE, 6 = POST) -----------------

bari <- clin[clin[["group iga"]] %in% c("5", "6"), ]

# Identify patients with metagenomics at BOTH timepoints
pz_pre  <- bari[bari[["group iga"]] == "5", "pz"]
pz_post <- bari[bari[["group iga"]] == "6", "pz"]
complete_pz <- intersect(pz_pre, pz_post)

cat("Patients with PRE sample:", length(pz_pre), "\n")
cat("Patients with POST sample:", length(pz_post), "\n")
cat("Complete pairs:", length(complete_pz), "\n\n")

bari <- bari[bari[["pz"]] %in% complete_pz, ]
cat("Rows retained:", nrow(bari), "\n")

# --- Add timepoint label and reorder columns ---------------------------------

bari[["timepoint"]] <- ifelse(bari[["group iga"]] == "5", "PRE", "POST")

# Drop columns that are entirely NA in this subset or non-informative
drop_cols <- c(
  "code iga", "date of birth", "visit date", "ultrasound date",
  "...33",                          # ghost column from Excel import
  "fatty liver index...50",         # duplicate
  "fatty liver index...51",         # duplicate
  "eosinoph",                       # duplicate of eosinophil
  "stiffness pre", "iqr pre", "cap  pre",  # PRE-only ultrasound columns
  "us steatosis pre", "cca prox dx pre", "cca prox sn pre",
  "bca dx  pre", "bca sn pre", "ica dx pre", "ica sn pre",
  "fib-4 index pre", "fatty liver index pre", "nafld fibrosis score pre",
  "mets-ir pre", "predimed score pre", "predimed categoria pre"
)
bari <- bari[, !colnames(bari) %in% drop_cols]

# Move identifier columns to the front
id_cols  <- c("sample_id", "pz", "group iga", "timepoint")
other    <- setdiff(colnames(bari), id_cols)
bari     <- bari[, c(id_cols, other)]

# Sort: PRE before POST within each patient, patients in consistent order
bari <- bari[order(match(bari[["pz"]], complete_pz), bari[["group iga"]]), ]

# --- Assign remission label --------------------------------------------------

# Remission is a patient-level label derived from MetS status at both timepoints.
# Applied to all rows of that patient (useful for filtering either timepoint).
mets_pre  <- setNames(
  bari[bari[["timepoint"]] == "PRE",  c("pz", "MetS")],
  c("pz", "MetS_pre")
)
mets_post <- setNames(
  bari[bari[["timepoint"]] == "POST", c("pz", "MetS")],
  c("pz", "MetS_post")
)
mets_pair <- merge(mets_pre, mets_post, by = "pz", all = TRUE)

mets_pair[["remission"]] <- case_when(
  mets_pair$MetS_pre == 1 & mets_pair$MetS_post == 0 ~ "remitter",
  mets_pair$MetS_pre == 1 & mets_pair$MetS_post == 1 ~ "non_remitter",
  mets_pair$MetS_pre == 0                             ~ "mets_neg",
  TRUE                                                ~ NA_character_
)

bari <- merge(bari, mets_pair[, c("pz", "remission")], by = "pz", all.x = TRUE)

# Re-sort after merge
bari <- bari[order(match(bari[["pz"]], complete_pz), bari[["group iga"]]), ]

cat("\nRemission distribution (at patient level):\n")
remission_summary <- unique(bari[, c("pz", "remission")])
print(table(remission_summary[["remission"]], useNA = "ifany"))

cat("\nMetS distribution by timepoint:\n")
print(table(timepoint = bari[["timepoint"]], MetS = bari[["MetS"]], useNA = "ifany"))

# --- Write outputs -----------------------------------------------------------

write.csv(bari,
          file = file.path(PROC_DIR, "bariatric_clinical_66.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(PROC_DIR, "bariatric_clinical_66.csv"))

write.xlsx(bari,
           file = file.path(PROC_DIR, "bariatric_clinical_66.xlsx"),
           overwrite = TRUE)
message("Written: ", file.path(PROC_DIR, "bariatric_clinical_66.xlsx"))
