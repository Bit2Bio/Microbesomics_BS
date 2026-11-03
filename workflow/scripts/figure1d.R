library(phyloseq)
library(ggplot2)
library(ggdist)
library(dplyr)

PS_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_rarefied_alpha.rds"
OUT     <- "/home/lorenzo/Microbesomics/results/figures/figure1d.png"

COL_PRE  <- "#4A0080"
COL_POST <- "#007A87"

ps   <- readRDS(PS_PATH)
meta <- as(sample_data(ps), "data.frame")

# ── 1. Calculate alpha diversity ─────────────────────────────────────────────
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha$sample_id <- rownames(alpha)
alpha$pz        <- meta[rownames(alpha), "pz"]
alpha$group     <- meta[rownames(alpha), "group.iga"]   # 5=PRE, 6=POST
alpha$timepoint <- ifelse(alpha$group == 5, "PRE", "POST")
alpha$timepoint <- factor(alpha$timepoint, levels = c("PRE", "POST"))

# ── 2. Paired Wilcoxon PRE vs POST + BH correction over both metrics ─────────
paired_test <- function(df, metric) {
  pre  <- df[df$timepoint == "PRE",  c("pz", metric)]
  post <- df[df$timepoint == "POST", c("pz", metric)]
  both <- merge(pre, post, by = "pz", suffixes = c("_pre", "_post"))
  wilcox.test(both[[paste0(metric, "_pre")]],
              both[[paste0(metric, "_post")]],
              paired = TRUE, alternative = "two.sided")$p.value
}

p_obs <- paired_test(alpha, "Observed")
p_sha <- paired_test(alpha, "Shannon")

# BH correction over the 2 tests
p_adj <- p.adjust(c(p_obs, p_sha), method = "BH")
names(p_adj) <- c("Observed", "Shannon")

cat(sprintf("Observed: p=%.4f, p-adj=%.4f\n", p_obs, p_adj["Observed"]))
cat(sprintf("Shannon:  p=%.4f, p-adj=%.4f\n", p_sha, p_adj["Shannon"]))

sig_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ── 3. Raincloud plot helper ──────────────────────────────────────────────────
# Same physical dimensions as figure1a/b/c: 1.3 × 1.435 inches
# Two metrics side by side → total width 2 × 1.3 = 2.6 inches

make_raincloud <- function(df, metric, y_label, stars) {

  d <- df[, c("pz", "timepoint", metric)]
  colnames(d)[3] <- "value"

  # summary stats for significance bracket
  y_max  <- max(d$value, na.rm = TRUE)
  y_min  <- min(d$value, na.rm = TRUE)
  y_span <- y_max - y_min
  y_sig  <- y_max + y_span * 0.08

  ggplot(d, aes(x = timepoint, y = value, fill = timepoint, colour = timepoint)) +

    # half-violin (rain cloud)
    ggdist::stat_halfeye(
      aes(fill = timepoint),
      adjust      = 0.6,
      width       = 0.4,
      .width      = 0,
      point_colour = NA,
      alpha       = 0.55,
      side        = "left"    # PRE left-side, POST left-side (mirrored by flip below)
    ) +

    # raw jittered points
    geom_point(
      position = position_jitter(width = 0.04, seed = 42),
      size = 0.6, alpha = 0.8, shape = 16
    ) +

    # thin boxplot overlay
    geom_boxplot(
      width = 0.07, outlier.shape = NA,
      colour = "white", fill = NA, linewidth = 0.4
    ) +
    geom_boxplot(
      width = 0.07, outlier.shape = NA,
      alpha = 0, linewidth = 0.4
    ) +

    # significance bracket
    annotate("segment",
             x = 1, xend = 2, y = y_sig, yend = y_sig,
             colour = "#777", linewidth = 0.3) +
    annotate("text",
             x = 1.5, y = y_sig + y_span * 0.02,
             label = stars,
             size  = ifelse(stars != "ns", 2.8, 2.2),
             colour = ifelse(stars %in% c("ns", ""), "#aaa", "#9B2226"),
             fontface = ifelse(stars != "ns", "bold", "plain"),
             vjust = 0) +

    scale_fill_manual(values   = c("PRE" = COL_PRE, "POST" = COL_POST)) +
    scale_colour_manual(values = c("PRE" = COL_PRE, "POST" = COL_POST)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 5) +
    theme(
      legend.position  = "none",
      axis.line        = element_line(linewidth = 0.3, colour = "#888"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "#888"),
      axis.text        = element_text(size = 3.5, colour = "#333"),
      axis.title.y     = element_text(size = 3.5, colour = "#333"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

p_observed <- make_raincloud(alpha, "Observed", "Observed ASVs",
                              sig_stars(p_adj["Observed"]))
p_shannon  <- make_raincloud(alpha, "Shannon",  "Shannon entropy",
                              sig_stars(p_adj["Shannon"]))

# ── 4. Combine and save ───────────────────────────────────────────────────────
library(patchwork)

combined <- p_observed + p_shannon +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(OUT, combined,
       width  = 2 * 1.3,
       height = 11 * 0.105 + 0.28,   # same height as figure1a
       units  = "in", dpi = 300)

cat(sprintf("Written: %s\n", OUT))
