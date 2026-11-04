library(ggplot2)
library(ggalluvial)

DATA <- "/home/lorenzo/Microbesomics/data/processed/bariatric_clinical_66.csv"
OUT  <- "/home/lorenzo/Microbesomics/results/figures/figure1d_mets_alluvial.png"

COL_POS <- "#C1121F"   # MetS+
COL_NEG <- "#2A9D8F"   # MetS-

clinical <- read.csv(DATA)
pre  <- clinical[clinical$timepoint == "PRE",  c("pz", "MetS")]
post <- clinical[clinical$timepoint == "POST", c("pz", "MetS")]
both <- merge(pre, post, by = "pz", suffixes = c("_pre", "_post"))
both <- both[!is.na(both$MetS_pre) & !is.na(both$MetS_post), ]

both$MetS_pre  <- ifelse(both$MetS_pre  == 1, "MetS+", "MetS−")
both$MetS_post <- ifelse(both$MetS_post == 1, "MetS+", "MetS−")
both$MetS_pre  <- factor(both$MetS_pre,  levels = c("MetS+", "MetS−"))
both$MetS_post <- factor(both$MetS_post, levels = c("MetS+", "MetS−"))

# counts per flow
flows <- as.data.frame(table(MetS_pre = both$MetS_pre, MetS_post = both$MetS_post))
flows <- flows[flows$Freq > 0, ]

# alluvial format: one row per patient
both$id <- seq_len(nrow(both))

p <- ggplot(both,
       aes(axis1 = MetS_pre, axis2 = MetS_post)) +

  geom_alluvium(aes(fill = MetS_post), alpha = 0.60, width = 0.18, curve_type = "sigmoid") +
  geom_stratum(width = 0.18, fill = "white", colour = "#888", linewidth = 0.25) +
  geom_text(stat = "stratum",
            aes(label = paste0(after_stat(stratum), "  n=", after_stat(n))),
            size = 1.6, fontface = "plain", colour = "#222", lineheight = 1.2) +

  scale_fill_manual(values = c("MetS+" = COL_POS, "MetS−" = COL_NEG)) +
  scale_x_discrete(limits = c("PRE", "POST"), expand = c(0.25, 0.25)) +

  theme_void(base_size = 4) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = 3.5, colour = "#333", margin = margin(t = 2)),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(OUT, p,
       width  = 1.3,
       height = 1.10,
       units  = "in", dpi = 300)
cat(sprintf("Written: %s\n", OUT))
