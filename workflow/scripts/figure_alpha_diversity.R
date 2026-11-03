library(phyloseq)
library(ggplot2)
library(patchwork)

PS_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_rarefied_alpha.rds"
OUT_PNG <- "/home/lorenzo/Microbesomics/results/figures/figure_alpha_diversity.png"

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

COL_PRE  <- "#4A0080"
COL_POST <- "#007A87"

ps   <- readRDS(PS_PATH)
meta <- as(sample_data(ps), "data.frame")

alpha_df           <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha_df$pz        <- meta[rownames(alpha_df), "pz"]
alpha_df$timepoint <- factor(ifelse(meta[rownames(alpha_df), "group.iga"] == 5, "PRE", "POST"),
                              levels = c("PRE", "POST"))

paired_plot <- function(df, metric, y_label) {
  pairs <- merge(
    df[df$timepoint == "PRE",  c("pz", metric)],
    df[df$timepoint == "POST", c("pz", metric)],
    by = "pz", suffixes = c("_pre", "_post")
  )
  colnames(pairs)[2:3] <- c("pre", "post")

  long <- data.frame(
    pz        = rep(pairs$pz, 2),
    value     = c(pairs$pre, pairs$post),
    timepoint = factor(c(rep("PRE", nrow(pairs)), rep("POST", nrow(pairs))),
                       levels = c("PRE", "POST"))
  )

  wt <- wilcox.test(pairs$pre, pairs$post, paired = TRUE)
  pv_label <- ifelse(wt$p.value < 0.001, "p<0.001",
              ifelse(wt$p.value < 0.01,  sprintf("p=%.3f", wt$p.value),
                                         sprintf("p=%.2f",  wt$p.value)))
  cat(sprintf("%s: W=%g, p=%.4f\n", metric, wt$statistic, wt$p.value))

  y_rng   <- range(long$value)
  y_annot <- y_rng[2] + diff(y_rng) * 0.12

  ggplot(long, aes(x = timepoint, y = value, group = pz)) +
    geom_line(colour = "#cccccc", linewidth = 0.25, alpha = 0.8) +
    geom_point(aes(colour = timepoint), size = 0.8, alpha = 0.85, shape = 16) +
    annotate("segment", x = 1, xend = 2, y = y_annot, yend = y_annot,
             linewidth = 0.25, colour = "#555") +
    annotate("text", x = 1.5, y = y_annot + diff(y_rng) * 0.04,
             label = pv_label, size = 1.5, colour = "#444", vjust = 0) +
    scale_colour_manual(values = c("PRE" = COL_PRE, "POST" = COL_POST)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 5) +
    theme(
      legend.position  = "none",
      axis.title.y     = element_text(size = 3.5),
      axis.text        = element_text(size = 3.5),
      axis.line        = element_line(linewidth = 0.3, colour = "#888"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "#888"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.2)
    )
}

p_obs     <- paired_plot(alpha_df, "Observed", "Observed ASVs")
p_shannon <- paired_plot(alpha_df, "Shannon",  "Shannon index")

ggsave(OUT_PNG, p_obs + p_shannon,
       width = 2.4, height = 1.1, units = "in", dpi = 300)
cat(sprintf("Written: %s\n", OUT_PNG))
