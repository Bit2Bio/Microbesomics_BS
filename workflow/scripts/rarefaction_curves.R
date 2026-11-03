library(phyloseq)
library(vegan)
library(ggplot2)

PS_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_filt_5vs6.rds"
OUT_DIR <- "/home/lorenzo/Microbesomics/results/figures/qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUT_DIR, "rarefaction_curves.png")

ps  <- readRDS(PS_PATH)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)   # samples as rows, taxa as columns

meta      <- as(sample_data(ps), "data.frame")
# group.iga: 5 = PRE, 6 = POST
grp       <- meta[rownames(otu), "group.iga"]
timepoint <- setNames(ifelse(grp == 5, "PRE", "POST"), rownames(otu))
lib_sizes <- rowSums(otu)

cat("Library sizes:\n")
print(summary(lib_sizes))
cat(sprintf("n samples: %d\n", nrow(otu)))

# depth grid: fine at low end, coarser at high end
max_depth <- max(lib_sizes)
depths <- unique(as.integer(c(
  seq(1000, 20000,    by = 2000),
  seq(20000, 100000,  by = 10000),
  seq(100000, max_depth, length.out = 30)
)))
depths <- sort(depths[depths >= 1])

set.seed(42)

# for each sample compute observed ASVs at each valid depth
results <- do.call(rbind, lapply(seq_len(nrow(otu)), function(i) {
  samp <- otu[i, ]
  n    <- sum(samp)
  vd   <- depths[depths <= n]
  if (length(vd) == 0) return(NULL)
  richness <- sapply(vd, function(d) {
    if (d == n) return(sum(samp > 0))
    sum(rrarefy(samp, d) > 0)   # single subsample per depth (fast)
  })
  data.frame(
    sample    = rownames(otu)[i],
    timepoint = timepoint[i],
    depth     = vd,
    richness  = richness,
    stringsAsFactors = FALSE
  )
}))

cat(sprintf("Curve data: %d rows, %d samples\n",
            nrow(results), length(unique(results$sample))))

# candidate thresholds
thresh <- data.frame(
  label = c("min", "10th pct", "20th pct"),
  depth = c(min(lib_sizes),
            as.integer(quantile(lib_sizes, 0.10)),
            as.integer(quantile(lib_sizes, 0.20))),
  col   = c("#e63946", "#f4a261", "#2a9d8f"),
  stringsAsFactors = FALSE
)

col_map <- c("PRE" = "#4A0080", "POST" = "#007A87")

p <- ggplot(results, aes(x = depth, y = richness,
                          group = sample, colour = timepoint)) +
  geom_line(alpha = 0.65, linewidth = 0.6) +
  geom_vline(data = thresh,
             aes(xintercept = depth, colour = label),
             linetype = "dashed", linewidth = 0.8,
             inherit.aes = FALSE) +
  scale_colour_manual(
    values = c(col_map,
               setNames(thresh$col, thresh$label)),
    breaks = c("PRE", "POST", thresh$label),
    labels = c("PRE", "POST",
               sprintf("%s (%s reads)",
                       thresh$label,
                       format(thresh$depth, big.mark = ",")))
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x      = "Sequencing depth (reads)",
    y      = "Observed ASVs",
    colour = NULL,
    title  = "Rarefaction curves — bariatric surgery cohort (n=33 PRE/POST)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 7),
    plot.title       = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.3)
  )

ggsave(OUT, p, width = 20, height = 13, units = "cm", dpi = 300)
cat(sprintf("Written: %s\n", OUT))
