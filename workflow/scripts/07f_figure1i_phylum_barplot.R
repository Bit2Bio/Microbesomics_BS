library(phyloseq)
library(ggplot2)
library(dplyr)

PS_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_filt_5vs6.rds"
OUT_PNG <- "/home/lorenzo/Microbesomics/results/figures/figure1i_composition.png"
TOP_N   <- 10

COL_PRE  <- "#4A0080"
COL_POST <- "#007A87"

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

ps <- readRDS(PS_PATH)

ps_p   <- tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
ps_rel <- transform_sample_counts(ps_p, function(x) x / sum(x))
df     <- psmelt(ps_rel)
df$Phylum[is.na(df$Phylum) | df$Phylum == ""] <- "Unassigned"
df$timepoint <- factor(ifelse(df$group.iga == 5, "PRE", "POST"),
                       levels = c("PRE", "POST"))

top_phyla <- df %>%
  group_by(Phylum) %>%
  summarise(m = mean(Abundance), .groups = "drop") %>%
  arrange(desc(m)) %>% slice_head(n = TOP_N) %>% pull(Phylum)

df$Phylum <- ifelse(df$Phylum %in% top_phyla, df$Phylum, "Other")

df_samp <- df %>%
  group_by(Sample, pz, timepoint, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Order by Firmicutes abundance in PRE; same order applied to POST so paired
# samples occupy the same x-position in both facets. PRE shows a decreasing
# gradient; POST uses the same patient order to allow visual comparison.
firm_pre <- df_samp[df_samp$timepoint == "PRE" & df_samp$Phylum == top_phyla[1],
                    c("pz", "Abundance")]
firm_pre <- firm_pre[order(firm_pre$Abundance, decreasing = TRUE), ]
pz_order <- firm_pre$pz

df_samp$pz_rank <- match(as.character(df_samp$pz), pz_order)
df_samp <- df_samp[order(df_samp$timepoint, df_samp$pz_rank), ]
df_samp$Sample <- factor(df_samp$Sample,
  levels = unique(df_samp$Sample[order(df_samp$timepoint, df_samp$pz_rank)]))

phylum_order <- c(rev(top_phyla), "Other")
df_samp$Phylum <- factor(df_samp$Phylum, levels = phylum_order)

# Kelly-inspired palette: maximally distinct hues across the spectrum.
# The 4 CLR-significant phyla (Firmicutes, Actinobacteriota, Proteobacteria,
# and borderline Verrucomicrobiota) carry the most saturated colours.
pal <- c(
  Firmicutes        = "#E07B39",  # terracotta
  Bacteroidota      = "#2E7D6E",  # deep teal
  Actinobacteriota  = "#7B5EA7",  # violet
  Proteobacteria    = "#4A90C4",  # steel blue
  Verrucomicrobiota = "#F2C94C",  # golden yellow
  Synergistota      = "#BE0032",  # crimson
  Fusobacteriota    = "#8DB600",  # lime green
  Patescibacteria   = "#FF8C69",  # salmon
  Desulfobacterota  = "#A0522D",  # sienna
  Euryarchaeota     = "#E68FAC",  # rose pink
  Other             = "#CCCCCC"
)
pal <- pal[phylum_order]

p <- ggplot(df_samp, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 1, color = NA) +
  scale_fill_manual(values = pal,
                    guide = guide_legend(nrow = 2, byrow = TRUE,
                                         keywidth  = unit(0.2, "cm"),
                                         keyheight = unit(0.2, "cm"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.01)),
                     breaks = c(0, 0.5, 1)) +
  facet_wrap(~ timepoint, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "Relative\nabundance", fill = NULL) +
  theme_classic(base_size = 5) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line.x      = element_blank(),
    axis.line        = element_line(linewidth = 0.3, colour = "#888"),
    axis.ticks.y     = element_line(linewidth = 0.3, colour = "#888"),
    axis.title.y     = element_text(size = 3.5),
    axis.text.y      = element_text(size = 3.5),
    strip.background = element_blank(),
    strip.text       = element_text(size = 4, face = "bold", margin = margin(b = 2)),
    panel.spacing    = unit(0.08, "cm"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 3),
    legend.spacing.x = unit(0.05, "cm"),
    legend.margin    = margin(t = 1),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(OUT_PNG, p, width = 3.4, height = 1.2, units = "in", dpi = 300)
cat(sprintf("Written: %s\n", OUT_PNG))
