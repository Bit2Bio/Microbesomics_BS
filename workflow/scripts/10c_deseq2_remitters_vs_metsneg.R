# =============================================================================
# 08_deseq2_remitters_vs_metsneg.R
#
# Differential abundance analysis (DESeq2): MetS+ remitters vs MetS- at
# baseline (PRE, group 5) among the 33 complete bariatric surgery pairs.
#
# This comparison asks what distinguishes the gut microbiome of future
# remitters (MetS+ who will achieve remission after surgery) from subjects
# already metabolically healthy at baseline.
#
# Patient classification (PRE timepoint only):
#   remitter : MetS == 1 at PRE AND MetS == 0 at POST  (n = 22)
#   mets_neg : MetS == 0 at PRE                         (n =  5)
#   excluded : non-remitters (MetS == 1 at both PRE and POST, n = 6)
#
# Input:
#   data/processed/ps_filt_5vs6.rds  — 66 samples (33 pairs), 707 ASVs
#   data/processed/clinical_all.csv  — MetS status for all 247 samples
#
# Output:
#   results/deseq2/deseq2_remitters_vs_metsneg_results.xlsx
#   results/deseq2/volcano_remitters_vs_metsneg.png
#
# Design: ~ comp_group (two-group comparison at PRE timepoint)
#   Positive log2FC = enriched in remitters vs MetS-.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)

dir.create(file.path(BASE_DIR, "results/deseq2"), recursive = TRUE, showWarnings = FALSE)

# --- Load data ---------------------------------------------------------------

ps <- readRDS(file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
cat("Loaded phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# --- Classify patients using remission from sample_data ----------------------
# remission column added by 04_phyloseq_bariatric_filter.R:
#   remitter     = MetS+ at PRE, MetS- at POST
#   mets_neg     = MetS- at PRE
#   non_remitter = excluded from this comparison

sd <- data.frame(sample_data(ps), check.names = FALSE)

# Keep PRE samples of remitters and mets_neg only
keep   <- sd[["group.iga"]] == "5" & sd[["remission"]] %in% c("remitter", "mets_neg")
ps_pre <- prune_samples(keep, ps)

cat("\nSamples in comparison:", nsamples(ps_pre), "\n")

sd_pre <- data.frame(sample_data(ps_pre), check.names = FALSE)

# mets_neg as reference: positive log2FC = enriched in remitters
sd_pre$comp_group <- factor(
  ifelse(sd_pre$remission == "remitter", "remitter", "mets_neg"),
  levels = c("mets_neg", "remitter")
)
rownames(sd_pre) <- sd_pre$sample_id
sample_data(ps_pre) <- sample_data(sd_pre)

cat("Group sizes:\n")
print(table(sample_data(ps_pre)$comp_group))

# --- DESeq2 ------------------------------------------------------------------

dds <- phyloseq_to_deseq2(ps_pre, ~ comp_group)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Positive log2FC = enriched in remitters (reference = mets_neg)
res <- results(dds,
               contrast      = c("comp_group", "remitter", "mets_neg"),
               cooksCutoff   = FALSE,
               pAdjustMethod = "fdr")

cat("\nDESeq2 results summary:\n")
summary(res)

# --- Build results table with taxonomy ---------------------------------------

sigtab <- as.data.frame(res)
sigtab <- cbind(sigtab, as.data.frame(tax_table(ps_pre)))

cat("\nSignificant ASVs (FDR < 0.05):", sum(sigtab$padj < 0.05, na.rm = TRUE), "\n")
cat("  Enriched in remitters (log2FC > 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Enriched in MetS-     (log2FC < 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange < 0, na.rm = TRUE), "\n\n")

out_xlsx <- file.path(BASE_DIR, "results/deseq2/deseq2_remitters_vs_metsneg_results.xlsx")
write.xlsx(sigtab, out_xlsx, rowNames = TRUE, overwrite = TRUE)
message("Written: ", out_xlsx)

# --- Volcano plot ------------------------------------------------------------

alpha <- 0.05

df <- sigtab
df$ASV    <- rownames(df)
df$Phylum <- ifelse(is.na(df$Phylum) | df$Phylum == "", "Unassigned", df$Phylum)
df$Genus  <- ifelse(is.na(df$Genus)  | df$Genus  == "", NA, df$Genus)

df$pval_safe <- ifelse(is.na(df$pvalue) | df$pvalue <= 0, .Machine$double.eps, df$pvalue)

sig_mask           <- !is.na(df$padj) & df$padj < alpha
fdr_pval_threshold <- if (any(sig_mask)) max(df$pvalue[sig_mask], na.rm = TRUE) else alpha

df$label <- ifelse(is.na(df$Genus),
                   paste0(df$ASV, "_", df$Phylum),
                   paste0(df$ASV, "_", df$Genus))

sig <- subset(df, sig_mask)
ns  <- df[!sig_mask, ]

sig_up   <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]
sig_up   <- sig_up[order(sig_up$pval_safe,    -abs(sig_up$log2FoldChange)), ]
sig_down <- sig_down[order(sig_down$pval_safe, -abs(sig_down$log2FoldChange)), ]
sig_top  <- rbind(head(sig_up, 10), head(sig_down, 10))

all_phyla  <- sort(unique(df$Phylum))
phylum_pal <- setNames(scales::hue_pal()(length(all_phyla)), all_phyla)

p <- ggplot() +
  geom_point(data = ns,
             aes(x = log2FoldChange, y = -log10(pval_safe)),
             color = "grey80", alpha = 0.6, size = 2) +
  geom_point(data = sig,
             aes(x = log2FoldChange, y = -log10(pval_safe), color = Phylum),
             alpha = 0.9, size = 3) +
  geom_text_repel(data = sig_top,
                  aes(x = log2FoldChange, y = -log10(pval_safe),
                      label = label, color = Phylum),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  geom_vline(xintercept = 0,  linetype = "dashed") +
  geom_vline(xintercept =  1, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = -log10(fdr_pval_threshold), linetype = "dashed") +
  scale_color_manual(values = phylum_pal) +
  coord_cartesian(xlim = c(-15, 15),
                  ylim = c(0, max(-log10(df$pval_safe), na.rm = TRUE) * 1.05)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x     = "log2 fold change (remitters vs MetS-)",
       y     = "-log10 p-value (FDR threshold line at padj = 0.05)",
       color = "Phylum")

out_png <- file.path(BASE_DIR, "results/deseq2/volcano_remitters_vs_metsneg.png")
ggsave(out_png, p, width = 30, height = 25, units = "cm", dpi = 300)
message("Written: ", out_png)
