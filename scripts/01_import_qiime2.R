# =============================================================================
# 02_import_qiime2.R
#
# Imports QIIME2 artifacts (QZA files) and builds a bare phyloseq object
# (OTU table + taxonomy + phylogenetic tree + reference sequences).
# Sample metadata is NOT attached here — that is done in downstream scripts.
#
# Input:  data/raw/qza/
#           asv_table.qza            — ASV count table
#           taxonomy.qza             — taxonomic assignments
#           fasttree_tree_rooted.qza — rooted phylogenetic tree
#           rep_seq.qza              — representative sequences
#
# Output: data/processed/
#           ps_raw.rds        — phyloseq object, 247 samples, no metadata
#           otu_table.csv     — ASV count matrix (ASVs × samples)
#           taxonomy.csv      — per-ASV taxonomic assignments
#           ps_raw.xlsx       — same tables as separate sheets
#
# Note on sample naming:
#   QZA sample IDs have the form "ID2988-N-N-boxK". N is the CODE IGA integer
#   used in all clinical metadata files, and the IDs come out of the BIOM file
#   already sorted 1..247, so a direct 1:247 rename is safe and keeps
#   downstream joins simple.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

# renv::restore()  # run once on a new machine to restore package versions

library(microbiomeMarker)
library(phyloseq)
library(openxlsx)

qza_dir  <- file.path(BASE_DIR, "data/raw/qza")
proc_dir <- file.path(BASE_DIR, "data/processed")

# Import all four QIIME2 artifacts into a single phyloseq object.
# import_qiime2() handles format conversion internally.
ps <- import_qiime2(
  otu_qza    = file.path(qza_dir, "asv_table.qza"),
  taxa_qza   = file.path(qza_dir, "taxonomy.qza"),
  sam_tab    = NULL,   # no sample metadata at this stage
  tree_qza   = file.path(qza_dir, "fasttree_tree_rooted.qza"),
  refseq_qza = file.path(qza_dir, "rep_seq.qza")
)

# QZA IDs have the form "ID2988-N-N-boxK". Extract N and prefix with "S"
# so sample names match CODE IGA in all clinical metadata files while
# avoiding bare integers as identifiers.
new_names <- paste0("S", sub("^ID[0-9]+-([0-9]+)-[0-9]+-.*$", "\\1", sample_names(ps)))
if (any(new_names == sample_names(ps)))
  stop("Regex did not match all sample IDs — check QZA naming format")
sample_names(ps) <- new_names

cat("Phyloseq object created:\n")
print(ps)

# --- Save RDS ----------------------------------------------------------------

saveRDS(ps, file = file.path(proc_dir, "ps_raw.rds"))
message("Written: ", file.path(proc_dir, "ps_raw.rds"))

# --- Extract flat tables -----------------------------------------------------

# OTU table: rows = ASVs, columns = samples (samples as columns is standard)
otu_df <- as.data.frame(otu_table(ps))
# Add ASV ID as an explicit column for readability in spreadsheet tools
otu_df <- cbind(ASV_ID = rownames(otu_df), otu_df)

# Taxonomy table
tax_df <- as.data.frame(tax_table(ps))
tax_df <- cbind(ASV_ID = rownames(tax_df), tax_df)

# --- Write CSVs --------------------------------------------------------------

write.csv(otu_df, file = file.path(proc_dir, "otu_table.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(proc_dir, "otu_table.csv"))

write.csv(tax_df, file = file.path(proc_dir, "taxonomy.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(proc_dir, "taxonomy.csv"))

# --- Write XLSX (one sheet per table) ----------------------------------------

write.xlsx(list(otu_table = otu_df, taxonomy = tax_df),
           file = file.path(proc_dir, "ps_raw.xlsx"), overwrite = TRUE)
message("Written: ", file.path(proc_dir, "ps_raw.xlsx"))
