# =============================================================================
# R/filter_phyloseq.R
#
# Reusable feature filtering function for phyloseq objects.
# Replicates the MicrobiomeAnalyst filtering pipeline applied to a subset:
#   1. Low count filter  : ASV must have >= min_count reads in >= min_prev_frac
#                          of the samples in the object passed in.
#   2. Low variance filter (optional): ASVs in the bottom var_pct by IQR are
#                          removed. Applied on raw (untransformed) counts.
#                          Set var_pct = 0 to skip (default, recommended for
#                          DESeq2 workflows where independent filtering is
#                          applied internally).
#
# Note on IQR ties: with sparse microbiome data many ASVs have IQR = 0,
# which means quantile(iqr_vals, var_pct) may also be 0, causing > 10% of
# ASVs to be removed. The guard `if (var_pct > 0)` prevents unintended
# removal when var_pct = 0.
#
# Usage:
#   source("R/filter_phyloseq.R")
#   ps_filt <- filter_ps(ps_subset)                        # count filter only
#   ps_filt <- filter_ps(ps_subset, var_pct = 0.10)        # + variance filter
#
# Parameters:
#   ps           — phyloseq object (already subsetted to the samples of interest)
#   min_count    — minimum reads to consider a sample "present" (default 4)
#   min_prev_frac— fraction of samples that must be "present" (default 0.10)
#   var_pct      — fraction of lowest-IQR ASVs to remove; 0 = skip (default 0)
#   verbose      — print filtering summary (default TRUE)
# =============================================================================

filter_ps <- function(ps,
                      min_count     = 4,
                      min_prev_frac = 0.10,
                      var_pct       = 0,
                      verbose       = TRUE) {

  n_start   <- ntaxa(ps)
  n_samples <- nsamples(ps)
  min_samp  <- ceiling(n_samples * min_prev_frac)

  # 1. Low count filter
  prevalence <- apply(otu_table(ps), 1, function(x) sum(x >= min_count))
  ps <- prune_taxa(prevalence >= min_samp, ps)

  if (verbose)
    message("Low count filter (>= ", min_count, " reads in >= ",
            round(min_prev_frac * 100), "% of ", n_samples, " samples = ",
            min_samp, " samples): ", ntaxa(ps), " / ", n_start, " ASVs retained")

  # 2. Low variance filter (optional, skipped when var_pct = 0)
  if (var_pct > 0) {
    iqr_vals      <- apply(as.matrix(otu_table(ps)), 1, IQR)
    iqr_threshold <- quantile(iqr_vals, var_pct)
    ps <- prune_taxa(iqr_vals > iqr_threshold, ps)

    if (verbose)
      message("Low variance filter (bottom ", round(var_pct * 100),
              "% by IQR removed): ", ntaxa(ps), " ASVs retained")
  }

  ps
}
