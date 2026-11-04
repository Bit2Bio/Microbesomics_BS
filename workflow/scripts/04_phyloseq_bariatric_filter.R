# =============================================================================
# 03_filter_5vs6.R
#
# Subsets the raw phyloseq object to the 33 complete bariatric surgery pairs —
# patients with metagenomics at BOTH PRE (group 5) and POST (group 6) — and
# applies feature filtering.
#
# Input:
#   data/processed/ps_raw.rds                — full unfiltered phyloseq (247 samples)
#   data/processed/clinical_all.csv          — clinical metadata with sample_id S1..S247
#   data/processed/bariatric_clinical_66.csv — enriched 66-row clinical file with
#                                              remission label (produced by 03b)
#
# Output:
#   data/processed/ps_filt_5vs6.rds — filtered phyloseq, 66 samples (33 pairs)
#                                      sample_data includes remission column
#
# Filtering (via filter_phyloseq.R):
#   Low count filter : >= 4 reads in >= 10% of 66 samples = 7 samples
#   Variance filter  : remove bottom 10% by IQR (raw counts)
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)

source(file.path(BASE_DIR, "workflow/scripts/filter_phyloseq.R"))

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
# var_pct = 0.10 removes the bottom 10% of ASVs by IQR (applied on raw counts)
ps_filt_5vs6 <- filter_ps(ps_5vs6, min_count = 4, min_prev_frac = 0.10, var_pct = 0.10)

# Add remission label from bariatric_clinical_66.csv (computed by 03b).
# Avoids each downstream DESeq2 script having to recompute it from clinical_all.csv.
bari <- read.csv(file.path(BASE_DIR, "data/processed/bariatric_clinical_66.csv"),
                 check.names = FALSE, na.strings = "")
rem  <- unique(bari[, c("pz", "remission")])
sd   <- as(sample_data(ps_filt_5vs6), "data.frame")
sd   <- merge(sd, rem, by = "pz", all.x = TRUE)
rownames(sd) <- sd[["sample_id"]]
sd   <- sd[sample_names(ps_filt_5vs6), ]   # restore original sample order
sample_data(ps_filt_5vs6) <- sample_data(sd)

cat("\nFiltered phyloseq object:\n")
print(ps_filt_5vs6)
cat("\nRemission distribution in sample_data:\n")
print(table(sample_data(ps_filt_5vs6)[["remission"]], useNA = "ifany"))

saveRDS(ps_filt_5vs6, file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
message("Written: ", file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
