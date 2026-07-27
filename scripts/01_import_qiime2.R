# =============================================================================
# 01_import_qiime2.R
#
# Imports QIIME2 artifacts (QZA files), builds a bare phyloseq object, and
# applies MicrobiomeAnalyst-equivalent feature filtering.
# Sample metadata is NOT attached here — that is done in downstream scripts.
#
# Input:  data/raw/qza/
#           asv_table.qza            — ASV count table
#           taxonomy.qza             — taxonomic assignments
#           fasttree_tree_rooted.qza — rooted phylogenetic tree
#           rep_seq.qza              — representative sequences
#
# Output: data/processed/
#           ps_raw.rds        — unfiltered phyloseq object, 247 samples
#           ps_filt.rds       — filtered phyloseq object
#           otu_table.csv     — filtered ASV count matrix (ASVs × samples)
#           taxonomy.csv      — per-ASV taxonomic assignments (filtered)
#
# Note on sample naming:
#   QZA sample IDs have the form "ID2988-N-N-boxK". N is the CODE IGA integer
#   used in all clinical metadata files. Renamed to S1..S247.
#
# Filtering (replicates MicrobiomeAnalyst defaults):
#   1. Low count filter  : ASV must have >= 4 reads in >= 20% of samples
#                          (>= 50 out of 247 samples)
#   2. Low variance filter: ASVs in the bottom 10% by IQR are removed
#                           (applied to log2-transformed counts)
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

cat("Phyloseq object created (unfiltered):\n")
print(ps)

# --- Save raw RDS ------------------------------------------------------------

saveRDS(ps, file = file.path(proc_dir, "ps_raw.rds"))
message("Written: ", file.path(proc_dir, "ps_raw.rds"))

# --- Feature filtering -------------------------------------------------------

n_samples <- nsamples(ps)

# 1. Low count filter: >= 4 reads in >= 20% of samples
min_count     <- 4
min_prev_frac <- 0.20
min_samples   <- ceiling(n_samples * min_prev_frac)

prevalence <- apply(otu_table(ps), 1, function(x) sum(x >= min_count))
ps_filt <- prune_taxa(prevalence >= min_samples, ps)

cat("\nAfter low count filter (>= ", min_count, " reads in >= ",
    min_prev_frac * 100, "% of samples = ", min_samples, " samples):\n", sep = "")
cat("ASVs retained:", ntaxa(ps_filt), "of", ntaxa(ps), "\n")

# 2. Low variance filter: remove bottom 10% by IQR on raw filtered counts.
# MicrobiomeAnalystR applies this filter before any normalisation/transformation,
# so raw counts are used here to replicate that behaviour faithfully.
var_pct <- 0.10
iqr_vals <- apply(as.matrix(otu_table(ps_filt)), 1, IQR)
iqr_threshold <- quantile(iqr_vals, var_pct)
ps_filt <- prune_taxa(iqr_vals > iqr_threshold, ps_filt)

cat("After low variance filter (bottom ", var_pct * 100, "% by IQR removed):\n", sep = "")
cat("ASVs retained:", ntaxa(ps_filt), "\n")

cat("\nFiltered phyloseq object:\n")
print(ps_filt)

saveRDS(ps_filt, file = file.path(proc_dir, "ps_filt.rds"))
message("Written: ", file.path(proc_dir, "ps_filt.rds"))

# --- Export filtered tables --------------------------------------------------

otu_df <- as.data.frame(otu_table(ps_filt))
otu_df <- cbind(ASV_ID = rownames(otu_df), otu_df)

tax_df <- as.data.frame(tax_table(ps_filt))
tax_df <- cbind(ASV_ID = rownames(tax_df), tax_df)

write.csv(otu_df, file = file.path(proc_dir, "otu_table.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(proc_dir, "otu_table.csv"))

write.csv(tax_df, file = file.path(proc_dir, "taxonomy.csv"),
          row.names = FALSE, na = "")
message("Written: ", file.path(proc_dir, "taxonomy.csv"))
