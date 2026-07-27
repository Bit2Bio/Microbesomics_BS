# =============================================================================
# 03_filter_5vs6.R
#
# Subsets the raw phyloseq object to the 33 complete bariatric surgery pairs —
# patients with metagenomics at BOTH PRE (group 5) and POST (group 6) — and
# applies feature filtering.
#
# Input:
#   data/processed/ps_raw.rds       — full unfiltered phyloseq (247 samples)
#   data/processed/clinical_all.csv — clinical metadata with sample_id S1..S247
#
# Output:
#   data/processed/ps_filt_5vs6.rds — filtered phyloseq, 66 samples (33 pairs)
#
# Filtering (via R/filter_phyloseq.R):
#   Low count filter: >= 4 reads in >= 10% of 66 samples = 7 samples
#   Variance filter: skipped (var_pct = 0) — DESeq2 applies independent
#   filtering internally.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)

source(file.path(BASE_DIR, "R/filter_phyloseq.R"))

ps   <- readRDS(file.path(BASE_DIR, "data/processed/ps_raw.rds"))
clin <- read.csv(file.path(BASE_DIR, "data/processed/clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

# Attach clinical metadata to phyloseq
rownames(clin) <- clin[["sample_id"]]
sample_data(ps) <- sample_data(clin[sample_names(ps), ])

# Identify the 33 patients with metagenomics at BOTH timepoints:
# a patient is kept only if their pz code appears in group 5 AND in group 6
pz_pre  <- clin[clin[["group iga"]] == "5", "pz"]
pz_post <- clin[clin[["group iga"]] == "6", "pz"]
complete_pz <- intersect(pz_pre, pz_post)

cat("Patients with QZA at PRE:", length(pz_pre), "\n")
cat("Patients with QZA at POST:", length(pz_post), "\n")
cat("Complete pairs (both timepoints):", length(complete_pz), "\n\n")

# Subset phyloseq to the 66 samples belonging to complete pairs
keep <- sample_data(ps)[["group iga"]] %in% c("5", "6") &
        sample_data(ps)[["pz"]] %in% complete_pz
ps_5vs6 <- prune_samples(keep, ps)

cat("Samples retained:", nsamples(ps_5vs6), "\n")
cat("ASVs before filtering:", ntaxa(ps_5vs6), "\n\n")

# Apply feature filtering on the 66-sample subset
ps_filt_5vs6 <- filter_ps(ps_5vs6, min_count = 4, min_prev_frac = 0.10, var_pct = 0)

cat("\nFiltered phyloseq object:\n")
print(ps_filt_5vs6)

saveRDS(ps_filt_5vs6, file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
message("Written: ", file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
