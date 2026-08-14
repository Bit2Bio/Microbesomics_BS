# =============================================================================
# 01_import_qiime2.R
#
# Imports QIIME2 artifacts (QZA files) and builds a bare phyloseq object.
# Sample metadata is NOT attached here — that is done in downstream scripts.
# Feature filtering is NOT applied here — it is handled group-specifically
# by R/filter_phyloseq.R in each analysis script.
#
# Input:  data/raw/qza/
#           asv_table.qza            — ASV count table
#           taxonomy.qza             — taxonomic assignments
#           fasttree_tree_rooted.qza — rooted phylogenetic tree
#           rep_seq.qza              — representative sequences
#
# Output: data/processed/
#           ps_raw.rds — unfiltered phyloseq object (247 samples, 14741 ASVs)
#
# Note on sample naming:
#   QZA sample IDs have the form "ID2988-N-N-boxK". N is the CODE IGA integer
#   used in all clinical metadata files. Renamed to S1..S247.
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

saveRDS(ps, file = file.path(proc_dir, "ps_raw.rds"))
message("Written: ", file.path(proc_dir, "ps_raw.rds"))
