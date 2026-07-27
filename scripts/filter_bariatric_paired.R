# =============================================================================
# 04_filter_bariatric_paired.R
#
# Filters the bariatric clinical dataset (script 03 output) to retain only
# patients with metagenomics data at BOTH timepoints (PRE and POST).
#
# Of the 40 bariatric pairs, 7 are missing at least one QZA sample. This
# script removes those incomplete pairs so that downstream paired analyses
# (DESeq2, etc.) work on a clean matched set.
#
# Input:
#   data/processed/bariatric_clinical.csv  — 80 rows (script 03 output)
#   data/processed/clinical_all.csv        — 247 QZA-linked samples (script 02)
#
# Output:
#   data/processed/bariatric_clinical_paired.csv — complete pairs only
#
# Logic:
#   A patient is kept if and only if their code (pz) appears in group 5 AND
#   in group 6 in clinical_all.csv (i.e. metagenomics was sequenced at both
#   PRE and POST).
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(dplyr)

bar  <- read.csv(file.path(BASE_DIR, "data/processed/bariatric_clinical.csv"),
                 check.names = FALSE, na.strings = "")
clin <- read.csv(file.path(BASE_DIR, "data/processed/clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

# Patients with QZA at PRE (group 5) and POST (group 6)
pz_pre  <- clin[clin[["group iga"]] == "5", "pz"]
pz_post <- clin[clin[["group iga"]] == "6", "pz"]
complete_pz <- intersect(pz_pre, pz_post)

cat("Pazienti con QZA al PRE:", length(pz_pre), "\n")
cat("Pazienti con QZA al POST:", length(pz_post), "\n")
cat("Coppie complete (entrambi i timepoint):", length(complete_pz), "\n\n")

bar_paired <- bar[bar[["code"]] %in% complete_pz, ]

cat("Righe prima del filtro:", nrow(bar), "\n")
cat("Righe dopo il filtro:", nrow(bar_paired), "\n\n")

cat("--- MetS prevalence (coppie complete) ---\n")
print(table(Timepoint = bar_paired$timepoint, MetS = bar_paired$MetS, useNA = "ifany"))

out_path <- file.path(BASE_DIR, "data/processed/bariatric_clinical_paired.csv")
write.csv(bar_paired, file = out_path, row.names = FALSE, na = "")
message("Written: ", out_path)
