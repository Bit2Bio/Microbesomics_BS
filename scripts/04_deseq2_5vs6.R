# =============================================================================
# 05_deseq2_5vs6.R
#
# Differential abundance analysis (DESeq2) for bariatric surgery:
# POST (group 6) vs PRE (group 5), paired design blocking on patient.
#
# Input:
#   data/processed/ps_filt_5vs6.rds — 66 samples (33 complete pairs), 1039 ASVs
#
# Output:
#   results/deseq2/deseq2_5vs6_results.xlsx — full DESeq2 results + taxonomy
#   results/deseq2/volcano_5vs6.png          — volcano plot
#
# Design: ~ pz + groupIGA
#   pz (patient ID) is the blocking factor; groupIGA (5=PRE, 6=POST) is the
#   variable of interest. The paired design removes between-patient variance,
#   increasing power to detect pre/post changes.
#
# Contrast: groupIGA 6 vs 5 → positive log2FC = higher in POST (group 6).
#
# Pseudocount: +1 added to the raw count table before phyloseq_to_deseq2().
#   DESeq2 cannot handle exact zeros in the count matrix.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)

dir.create(file.path(BASE_DIR, "results/deseq2"), recursive = TRUE, showWarnings = FALSE)

# --- Load filtered phyloseq --------------------------------------------------

ps <- readRDS(file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
cat("Loaded:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# --- Prepare sample variables ------------------------------------------------

# groupIGA and pz must be factors for DESeq2 model matrix.
# phyloseq sample_data doesn't support $<- directly; modify via data.frame.
# phyloseq coerces column names: spaces → dots ("group iga" → "group.iga")
sd <- data.frame(sample_data(ps), check.names = FALSE)
sd$groupIGA <- factor(sd[["group.iga"]])
sd$pz       <- factor(sd[["pz"]])
sample_data(ps) <- sample_data(sd)

# Sanity check: every patient should appear exactly twice (PRE + POST)
pz_counts <- table(sample_data(ps)$pz)
if (any(pz_counts != 2))
  stop("Some patients do not have exactly 2 samples: ",
       paste(names(pz_counts[pz_counts != 2]), collapse = ", "))

cat("Patients:", nlevels(sample_data(ps)$pz), "\n")
cat("Groups:  ", levels(sample_data(ps)$groupIGA), "\n\n")

# --- DESeq2 ------------------------------------------------------------------

dds <- phyloseq_to_deseq2(ps, ~ pz + groupIGA)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Contrast: POST (6) vs PRE (5); positive log2FC = enriched POST-surgery
res <- results(dds,
               contrast      = c("groupIGA", "6", "5"),
               cooksCutoff   = FALSE,
               pAdjustMethod = "fdr")

cat("DESeq2 results summary:\n")
summary(res)

# --- Build results table with taxonomy ---------------------------------------

sigtab <- as.data.frame(res)
sigtab <- cbind(sigtab, as.data.frame(tax_table(ps)))

cat("\nSignificant ASVs (FDR < 0.05):", sum(sigtab$padj < 0.05, na.rm = TRUE), "\n")
cat("  Enriched POST (log2FC > 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Enriched PRE  (log2FC < 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange < 0, na.rm = TRUE), "\n\n")

out_xlsx <- file.path(BASE_DIR, "results/deseq2/deseq2_5vs6_results.xlsx")
write.xlsx(sigtab, out_xlsx, rowNames = TRUE, overwrite = TRUE)
message("Written: ", out_xlsx)

# --- Volcano plot ------------------------------------------------------------

alpha <- 0.05

df <- sigtab
df$ASV       <- rownames(df)
df$Phylum    <- ifelse(is.na(df$Phylum) | df$Phylum == "", "Unassigned", df$Phylum)
df$Genus     <- ifelse(is.na(df$Genus)  | df$Genus  == "", NA, df$Genus)
min_nonzero  <- min(df$padj[df$padj > 0], na.rm = TRUE)
df$padj_safe <- ifelse(is.na(df$padj) | df$padj <= 0, min_nonzero * 0.1, df$padj)

df$label <- ifelse(is.na(df$Genus),
                   paste0(df$ASV, "_", df$Phylum),
                   paste0(df$ASV, "_", df$Genus))

sig <- subset(df, padj_safe <= alpha & abs(log2FoldChange) > 1)
ns  <- subset(df, padj_safe >  alpha | abs(log2FoldChange) <= 1)

# Top 15 per side (by significance then |log2FC|) to label both directions
sig_up   <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]
sig_up   <- sig_up[order(sig_up$padj_safe,     -abs(sig_up$log2FoldChange)),   ]
sig_down <- sig_down[order(sig_down$padj_safe,  -abs(sig_down$log2FoldChange)), ]
sig_top  <- rbind(head(sig_up, 15), head(sig_down, 15))

all_phyla  <- sort(unique(df$Phylum))
phylum_pal <- setNames(scales::hue_pal()(length(all_phyla)), all_phyla)

p <- ggplot() +
  geom_point(data = ns,
             aes(x = log2FoldChange, y = -log10(padj_safe)),
             color = "grey80", alpha = 0.6, size = 2) +
  geom_point(data = sig,
             aes(x = log2FoldChange, y = -log10(padj_safe), color = Phylum),
             alpha = 0.9, size = 3) +
  geom_text_repel(data = sig_top,
                  aes(x = log2FoldChange, y = -log10(padj_safe),
                      label = label, color = Phylum),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  geom_vline(xintercept =  1, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  scale_color_manual(values = phylum_pal) +
  coord_cartesian(ylim = c(0, max(-log10(df$padj_safe), na.rm = TRUE) * 1.05)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x     = "log2 fold change (POST vs PRE bariatric surgery)",
       y     = "-log10 adjusted p-value",
       color = "Phylum")

out_png <- file.path(BASE_DIR, "results/deseq2/volcano_5vs6.png")
ggsave(out_png, p, width = 30, height = 25, units = "cm", dpi = 300)
message("Written: ", out_png)
