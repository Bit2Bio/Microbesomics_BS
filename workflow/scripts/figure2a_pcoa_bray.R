library(phyloseq)
library(vegan)
library(ggplot2)

PS_PATH   <- "/home/lorenzo/Microbesomics/data/processed/ps_rarefied_beta.rds"
OUT_BC    <- "/home/lorenzo/Microbesomics/results/figures/figure2a_pcoa_braycurtis.png"
OUT_BC_D  <- "/home/lorenzo/Microbesomics/results/figures/figure2b_bc_distances.png"
OUT_UUF   <- "/home/lorenzo/Microbesomics/results/figures/figure2c_pcoa_unifrac.png"
OUT_UUF_D <- "/home/lorenzo/Microbesomics/results/figures/figure2d_uuf_distances.png"
OUT_CSV   <- "/home/lorenzo/Microbesomics/results/tables/beta_diversity_permanova.csv"

for (d in c(dirname(OUT_BC), dirname(OUT_CSV)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

COL_PRE  <- "#4A0080"
COL_POST <- "#007A87"
FIG_W    <- 1.5
FIG_H    <- 1.5

ps   <- readRDS(PS_PATH)
meta <- as(sample_data(ps), "data.frame")
meta$timepoint <- factor(ifelse(meta$group.iga == 5, "PRE", "POST"),
                         levels = c("PRE", "POST"))

# meta_ord: sample rows aligned with OTU table, shared by both helpers
otu_mat  <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
meta_ord <- meta[rownames(otu_mat), ]

theme_beta <- function() {
  theme_classic(base_size = 5) +
    theme(
      axis.title       = element_text(size = 3.5),
      axis.text        = element_text(size = 3.5),
      axis.line        = element_line(linewidth = 0.3, colour = "#888"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "#888"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.2)
    )
}

# ── PCoA + paired PERMANOVA ───────────────────────────────────────────────────
run_pcoa <- function(d, dist_label, out_path) {
  pcoa <- ordinate(ps, method = "PCoA", distance = d)

  eig <- pcoa$values$Eigenvalues
  pct <- round(eig / sum(eig[eig > 0]) * 100, 1)
  ax_lab <- function(i) sprintf("PC%d (%.1f%%)", i, pct[i])

  coords <- as.data.frame(pcoa$vectors[, 1:2])
  colnames(coords) <- c("PC1", "PC2")
  coords$pz        <- meta[rownames(coords), "pz"]
  coords$timepoint <- meta[rownames(coords), "timepoint"]

  set.seed(42)
  perm <- adonis2(d ~ timepoint,
                  data         = meta_ord,
                  strata       = meta_ord$pz,
                  permutations = 999)

  r2 <- round(perm$R2[1], 3)
  pv <- perm$`Pr(>F)`[1]
  f  <- round(perm$F[1], 3)
  pv_label <- ifelse(pv < 0.001, "p<0.001",
              ifelse(pv < 0.01,  sprintf("p=%.3f", pv),
                                 sprintf("p=%.2f", pv)))
  cat(sprintf("%s — R²=%.3f, %s\n", dist_label, r2, pv_label))

  pairs_df <- merge(
    coords[coords$timepoint == "PRE",  c("pz", "PC1", "PC2")],
    coords[coords$timepoint == "POST", c("pz", "PC1", "PC2")],
    by = "pz", suffixes = c("_pre", "_post")
  )

  p <- ggplot(coords, aes(x = PC1, y = PC2, colour = timepoint, fill = timepoint)) +
    geom_segment(data = pairs_df,
                 aes(x = PC1_pre, y = PC2_pre, xend = PC1_post, yend = PC2_post),
                 colour = "#cccccc", linewidth = 0.25, alpha = 0.55,
                 inherit.aes = FALSE) +
    stat_ellipse(geom = "polygon", alpha = 0.10, linewidth = 0.3, level = 0.95) +
    stat_ellipse(linewidth = 0.4, level = 0.95) +
    geom_point(size = 0.8, alpha = 0.85, shape = 16) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("R²=%.3f\n%s", r2, pv_label),
             hjust = 1.05, vjust = 1.3, size = 1.5,
             colour = "#444", lineheight = 1.3) +
    scale_colour_manual(values = c("PRE" = COL_PRE, "POST" = COL_POST)) +
    scale_fill_manual(values   = c("PRE" = COL_PRE, "POST" = COL_POST)) +
    labs(x = ax_lab(1), y = ax_lab(2), colour = NULL, fill = NULL,
         title = dist_label) +
    theme_beta() +
    theme(
      legend.position   = c(0.15, 0.13),
      legend.text       = element_text(size = 3.5),
      legend.key.size   = unit(0.25, "cm"),
      legend.background = element_rect(fill = "white", colour = "#ddd", linewidth = 0.2),
      plot.title        = element_text(size = 4, face = "plain", colour = "#555",
                                       margin = margin(b = 2))
    )

  ggsave(out_path, p, width = FIG_W, height = FIG_H, units = "in", dpi = 300)
  cat(sprintf("Written: %s\n", out_path))

  list(dist_label = dist_label, r2 = r2, pv = pv, f = f,
       PC1_pct = pct[1], PC2_pct = pct[2])
}

# ── Within-patient vs between-patient pairwise distances ──────────────────────
# Within: d(PRE_i, POST_i) — same patient before/after surgery
# Between: d(sample_i, sample_j) for patients i ≠ j — baseline inter-individual spread
run_distances <- function(d, dist_label, out_path) {
  d_mat <- as.matrix(d)

  pre_samps  <- rownames(meta_ord)[meta_ord$timepoint == "PRE"]
  post_samps <- rownames(meta_ord)[meta_ord$timepoint == "POST"]

  pairs <- merge(
    data.frame(pz = meta_ord[pre_samps,  "pz"], pre  = pre_samps,  stringsAsFactors = FALSE),
    data.frame(pz = meta_ord[post_samps, "pz"], post = post_samps, stringsAsFactors = FALSE),
    by = "pz"
  )
  within_d <- mapply(function(pr, po) d_mat[pr, po], pairs$pre, pairs$post)

  idx   <- which(upper.tri(d_mat), arr.ind = TRUE)
  pz_i  <- meta_ord[rownames(d_mat)[idx[, 1]], "pz"]
  pz_j  <- meta_ord[rownames(d_mat)[idx[, 2]], "pz"]
  between_d <- d_mat[idx][pz_i != pz_j]

  wt <- wilcox.test(within_d, between_d)
  pv_label <- ifelse(wt$p.value < 0.001, "p<0.001",
              ifelse(wt$p.value < 0.01,  sprintf("p=%.3f", wt$p.value),
                                         sprintf("p=%.2f",  wt$p.value)))

  df <- data.frame(
    dist  = c(within_d, between_d),
    group = factor(c(rep("Within", length(within_d)),
                     rep("Between", length(between_d))),
                   levels = c("Within", "Between"))
  )

  y_max   <- max(df$dist)
  y_annot <- y_max * 1.08

  p <- ggplot(df, aes(x = group, y = dist)) +
    geom_boxplot(aes(fill = group), width = 0.45, outlier.size = 0.3,
                 linewidth = 0.3, alpha = 0.75) +
    annotate("segment",
             x = 1, xend = 2, y = y_annot, yend = y_annot,
             linewidth = 0.25, colour = "#555") +
    annotate("text",
             x = 1.5, y = y_annot * 1.02,
             label = pv_label, size = 1.5, colour = "#444", vjust = 0) +
    scale_fill_manual(values = c("Within"  = "#888888",
                                 "Between" = "#cccccc")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = NULL, y = dist_label) +
    theme_beta() +
    theme(legend.position = "none",
          axis.title.y    = element_text(size = 3.5))

  ggsave(out_path, p, width = FIG_W, height = FIG_H, units = "in", dpi = 300)
  cat(sprintf("Written: %s\n", out_path))
}

# ── Compute distances once, pass to both helpers ──────────────────────────────
d_bc  <- phyloseq::distance(ps, method = "bray")
d_uuf <- phyloseq::distance(ps, method = "unifrac")

res_bc  <- run_pcoa(d_bc,  "Bray-Curtis",        OUT_BC)
run_distances(d_bc,        "Bray-Curtis",         OUT_BC_D)

res_uuf <- run_pcoa(d_uuf, "Unweighted UniFrac",  OUT_UUF)
run_distances(d_uuf,       "Unweighted UniFrac",  OUT_UUF_D)

# ── PERMANOVA results table ───────────────────────────────────────────────────
results <- do.call(rbind, lapply(list(res_bc, res_uuf), function(r) {
  data.frame(
    comparison   = "PRE vs POST",
    distance     = r$dist_label,
    n_samples    = nsamples(ps),
    n_pairs      = nsamples(ps) / 2,
    F_statistic  = r$f,
    R2           = r$r2,
    p_value      = r$pv,
    permutations = 999,
    strata       = "patient (pz)",
    PC1_pct      = r$PC1_pct,
    PC2_pct      = r$PC2_pct
  )
}))
write.csv(results, OUT_CSV, row.names = FALSE)
cat(sprintf("Written: %s\n", OUT_CSV))
