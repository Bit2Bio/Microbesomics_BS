# =============================================================================
# 03_filter_5vs6.R
#
# Subsets the raw phyloseq object to bariatric surgery samples (groups 5 and 6),
# attaches clinical metadata, and applies feature filtering.
#
# Input:
#   data/processed/ps_raw.rds       — full unfiltered phyloseq (247 samples)
#   data/processed/clinical_all.csv — clinical metadata with sample_id S1..S247
#
# Output:
#   data/processed/ps_filt_5vs6.rds — filtered phyloseq, groups 5+6 only
#
# Filtering (via R/filter_phyloseq.R):
#   Low count filter: >= 4 reads in >= 10% of samples
#   For 66 samples (38 PRE + 28 POST after metagenomics availability):
#   threshold = ceiling(66 * 0.10) = 7 samples
#   Variance filter: skipped (var_pct = 0) — DESeq2 applies independent
#   filtering internally.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)

source(file.path(BASE_DIR, "R/filter_phyloseq.R"))

ps  <- readRDS(file.path(BASE_DIR, "data/processed/ps_raw.rds"))
clin <- read.csv(file.path(BASE_DIR, "data/processed/clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

# Attach clinical metadata to phyloseq
rownames(clin) <- clin[["sample_id"]]
sample_data(ps) <- sample_data(clin[sample_names(ps), ])

# Subset to bariatric surgery groups (5 = PRE, 6 = POST)
# subset_samples() uses NSE and doesn't handle column names with spaces;
# filter manually via sample_data index instead
keep <- sample_data(ps)[["group iga"]] %in% c("5", "6")
ps_5vs6 <- prune_samples(keep, ps)

cat("Samples in groups 5+6:", nsamples(ps_5vs6), "\n")
cat("ASVs before filtering:", ntaxa(ps_5vs6), "\n\n")

# Apply feature filtering on the subset
ps_filt_5vs6 <- filter_ps(ps_5vs6, min_count = 4, min_prev_frac = 0.10, var_pct = 0)

cat("\nFiltered phyloseq object:\n")
print(ps_filt_5vs6)

saveRDS(ps_filt_5vs6, file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
message("Written: ", file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
