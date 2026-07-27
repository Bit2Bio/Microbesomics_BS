# =============================================================================
# 02_preprocess_clinical.R
#
# Reads all four clinical sheets from the IGA Excel file, normalises column
# names, and stacks them into a single dataset.  Columns that exist in some
# sheets but not others are filled with NA.
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
#   clinical_all.csv   — stacked dataset, one row per sample
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

# --- Write outputs -----------------------------------------------------------

write.csv(clinical_all,
          file = file.path(PROC_DIR, "clinical_all.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(PROC_DIR, "clinical_all.csv"))

write.xlsx(clinical_all,
           file = file.path(PROC_DIR, "clinical_all.xlsx"),
           overwrite = TRUE)
message("Written: ", file.path(PROC_DIR, "clinical_all.xlsx"))
