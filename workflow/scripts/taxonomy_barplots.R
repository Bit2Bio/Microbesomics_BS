# =============================================================================
# taxonomy_barplots.R
#
# Generates a single stacked relative abundance barplot for a given taxonomic
# rank and comparison. Called by the Snakemake rule taxonomy_barplots once
# per (rank × comparison) combination — 9 runs total.
#
# This script replaces run_barplots.sh (retired), which embedded inline R
# code in a bash heredoc. The logic is identical; the interface is cleaner.
#
# Arguments (all required):
#   --rank         Taxonomic rank: Phylum | Family | Genus
#   --comparison   One of: 5vs6 | remission | mets
#   --top_n        Number of top taxa shown individually (rest → "Other")
#   --width        Output width in cm
#   --height       Output height in cm
#   --dpi          Output resolution
#   --ps           Path to ps_filt_5vs6.rds
#   --clinical     Path to clinical_all.csv
#   --out          Output PNG path
#
# Comparison logic:
#   5vs6       — all 66 samples, group.iga == "5" (PRE) vs "6" (POST)
#   remission  — PRE samples only, remitters vs non-remitters
#                (remitter = MetS 1→0 after surgery; excluded = MetS 0 at PRE)
#   mets       — PRE samples only, MetS+ vs MetS-
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

source(file.path(BASE_DIR, "workflow/scripts/plot_taxon_barplot.R"))

# --- Parse command-line arguments --------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) stop("Missing required argument: ", flag)
  args[idx + 1]
}

rank       <- get_arg("--rank")
comparison <- get_arg("--comparison")
top_n      <- as.integer(get_arg("--top_n"))
width_cm   <- as.numeric(get_arg("--width"))
height_cm  <- as.numeric(get_arg("--height"))
dpi        <- as.integer(get_arg("--dpi"))
ps_path    <- get_arg("--ps")
clin_path  <- get_arg("--clinical")
out_path   <- get_arg("--out")

stopifnot(rank       %in% c("Phylum", "Family", "Genus"))
stopifnot(comparison %in% c("5vs6", "remission", "mets"))

# --- Load data ---------------------------------------------------------------

ps   <- readRDS(ps_path)
clin <- read.csv(clin_path, check.names = FALSE, na.strings = "")

# --- Subset phyloseq and assign grouping variable ----------------------------

if (comparison == "5vs6") {

  # All 66 samples, grouped by timepoint
  p <- plot_taxon_barplot(
    ps, rank = rank, group_var = "group.iga", top_n = top_n,
    group_labels = c("5" = "PRE", "6" = "POST")
  )

} else if (comparison == "remission") {

  # PRE samples only; assign remission label per patient using MetS at both
  # timepoints. Patients with MetS == 0 at PRE are excluded (not informative
  # for remission: they had nothing to remit).
  pre  <- clin[clin[["group iga"]] == "5", c("pz", "MetS")]
  post <- clin[clin[["group iga"]] == "6", c("pz", "MetS")]
  complete_pz <- intersect(pre$pz, post$pz)

  pre  <- pre[pre$pz  %in% complete_pz, ]
  post <- post[post$pz %in% complete_pz, ]

  merged <- merge(pre, post, by = "pz", suffixes = c("_pre", "_post"))
  merged$remission <- ifelse(
    merged$MetS_pre == 1 & merged$MetS_post == 0, "remitter",
    ifelse(merged$MetS_pre == 1 & merged$MetS_post == 1, "non_remitter", NA)
  )
  keep_pz <- merged[!is.na(merged$remission), c("pz", "remission")]

  sd_df   <- data.frame(sample_data(ps), check.names = FALSE)
  keep    <- sd_df[["group.iga"]] == "5" & sd_df[["pz"]] %in% keep_pz$pz
  ps_sub  <- prune_samples(keep, ps)

  sd_sub  <- data.frame(sample_data(ps_sub), check.names = FALSE)
  sd_sub  <- merge(sd_sub, keep_pz, by = "pz", all.x = TRUE)
  rownames(sd_sub) <- sd_sub$sample_id
  sample_data(ps_sub) <- sample_data(sd_sub)

  p <- plot_taxon_barplot(
    ps_sub, rank = rank, group_var = "remission", top_n = top_n,
    group_labels = c("remitter" = "Remitters", "non_remitter" = "Non-remitters")
  )

} else if (comparison == "mets") {

  # PRE samples only; MetS+ vs MetS-
  pz_pre      <- clin[clin[["group iga"]] == "5", "pz"]
  pz_post     <- clin[clin[["group iga"]] == "6", "pz"]
  complete_pz <- intersect(pz_pre, pz_post)

  sd_df   <- data.frame(sample_data(ps), check.names = FALSE)
  keep    <- sd_df[["group.iga"]] == "5" & sd_df[["pz"]] %in% complete_pz
  ps_sub  <- prune_samples(keep, ps)

  sd_sub  <- data.frame(sample_data(ps_sub), check.names = FALSE)
  sd_sub$MetS_group <- factor(
    ifelse(sd_sub$MetS == 1, "MetS_pos", "MetS_neg"),
    levels = c("MetS_neg", "MetS_pos")
  )
  sample_data(ps_sub) <- sample_data(sd_sub)

  p <- plot_taxon_barplot(
    ps_sub, rank = rank, group_var = "MetS_group", top_n = top_n,
    group_labels = c("MetS_pos" = "MetS+", "MetS_neg" = "MetS−")
  )
}

# --- Save output -------------------------------------------------------------

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path, p, width = width_cm, height = height_cm,
       units = "cm", dpi = dpi)
message("Written: ", out_path)
