library(phyloseq)
library(ggplot2)

PS_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_filt_5vs6.rds"
OUT_PNG <- "/home/lorenzo/Microbesomics/results/figures/figure1j_phylum_clr.png"

COL_PRE  <- "#4A0080"
COL_POST <- "#007A87"

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

ps   <- readRDS(PS_PATH)
ps_p <- tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_p, function(x) x / sum(x))
df   <- psmelt(ps_rel)
df$Phylum[is.na(df$Phylum) | df$Phylum == ""] <- "Unassigned"
df$timepoint <- ifelse(df$group.iga == 5, "PRE", "POST")

# Pseudocount = half the minimum non-zero relative abundance across all samples.
# Required to log-transform phyla that have zeros in some patients.
# Sensitivity analysis restricted to zero-free phyla (n=5) confirmed that the
# three significant results are not driven by the pseudocount (see Methods).
pseudo <- min(df$Abundance[df$Abundance > 0]) / 2

phyla <- unique(df$Phylum)
phyla <- phyla[phyla != "NA" & !is.na(phyla) & phyla != "Unassigned"]

res <- do.call(rbind, lapply(phyla, function(ph) {
  d <- aggregate(Abundance ~ pz + timepoint,
                 data = df[df$Phylum == ph, ], FUN = sum)
  paired <- merge(d[d$timepoint == "PRE",  c("pz", "Abundance")],
                  d[d$timepoint == "POST", c("pz", "Abundance")],
                  by = "pz", suffixes = c("_pre", "_post"))
  if (nrow(paired) < 5) return(NULL)

  # Geometric mean per sample across all phyla (with pseudocount)
  gm_pre  <- sapply(paired$pz, function(p)
    exp(mean(log(df$Abundance[df$timepoint == "PRE"  & df$pz == p] + pseudo))))
  gm_post <- sapply(paired$pz, function(p)
    exp(mean(log(df$Abundance[df$timepoint == "POST" & df$pz == p] + pseudo))))

  clr_pre  <- log(paired$Abundance_pre  + pseudo) - log(gm_pre)
  clr_post <- log(paired$Abundance_post + pseudo) - log(gm_post)

  tt <- t.test(clr_post, clr_pre, paired = TRUE)
  data.frame(Phylum   = ph,
             clr_diff = mean(clr_post - clr_pre),
             p_value  = tt$p.value)
}))

res$p_adj <- p.adjust(res$p_value, method = "BH")
res$sig   <- ifelse(res$p_adj < 0.001, "***",
             ifelse(res$p_adj < 0.01,  "**",
             ifelse(res$p_adj < 0.05,  "*", "")))
res$direction <- ifelse(res$clr_diff > 0, "up", "down")

res <- res[order(res$clr_diff), ]
res$Phylum <- factor(res$Phylum, levels = res$Phylum)

p <- ggplot(res, aes(x = clr_diff, y = Phylum, fill = direction)) +
  geom_col(width = 0.65, colour = NA) +
  geom_text(aes(label = sig,
                x = clr_diff + ifelse(clr_diff >= 0, 0.04, -0.04)),
            hjust = ifelse(res$clr_diff >= 0, 0, 1),
            size = 2.2, colour = "#333", vjust = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888") +
  scale_fill_manual(values = c("up" = COL_POST, "down" = COL_PRE)) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.18))) +
  labs(x = "CLR difference (POST − PRE)", y = NULL) +
  theme_classic(base_size = 5) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_text(size = 3.5, face = "italic"),
    axis.text.x        = element_text(size = 3.5),
    axis.title.x       = element_text(size = 3.5),
    axis.line          = element_line(linewidth = 0.3, colour = "#888"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "#888"),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_line(colour = "#f0f0f0", linewidth = 0.2)
  )

ggsave(OUT_PNG, p, width = 1.2, height = 1.0, units = "in", dpi = 300)
cat(sprintf("Written: %s\n", OUT_PNG))
