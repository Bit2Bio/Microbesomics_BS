# =============================================================================
# 06_deseq2_remission.R
#
# Differential abundance analysis (DESeq2): MetS remitters vs non-remitters
# at baseline (PRE, group 5) among the 33 complete bariatric surgery pairs.
#
# Remission is defined at the patient level using MetS status from
# clinical_all.csv:
#   remitter     : MetS == 1 at PRE (group 5) AND MetS == 0 at POST (group 6)
#   non-remitter : MetS == 1 at PRE (group 5) AND MetS == 1 at POST (group 6)
#   excluded     : MetS == 0 at PRE (group 5) — 5 patients already MetS-negative
#
# Input:
#   data/processed/ps_filt_5vs6.rds  — 66 samples (33 pairs), 707 ASVs
#                                       same filtering as script 03
#   data/processed/clinical_all.csv  — MetS status for all 247 samples
#
# Output:
#   results/deseq2/deseq2_remission_results.xlsx — DESeq2 results + taxonomy
#   results/deseq2/volcano_remission.png          — volcano plot
#
# Design: ~ remission (simple two-group comparison at the same timepoint)
#   Positive log2FC = enriched in non-remitters.
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

# --- Subset phyloseq to PRE samples of remitters/non-remitters ---------------
# remission column is already in sample_data (added by 04_phyloseq_bariatric_filter.R)

sd <- data.frame(sample_data(ps), check.names = FALSE)

cat("\nRemission status (complete pairs):\n")
print(table(sd[sd[["group.iga"]] == "5", "remission"], useNA = "ifany"))

# Keep only PRE samples of remitters and non-remitters (exclude mets_neg)
keep   <- sd[["group.iga"]] == "5" & sd[["remission"]] %in% c("remitter", "non_remitter")
ps_pre <- prune_samples(keep, ps)

cat("\nPRE samples after exclusion of MetS- at PRE:", nsamples(ps_pre), "\n")

sd_pre <- data.frame(sample_data(ps_pre), check.names = FALSE)
sd_pre$remission <- factor(sd_pre$remission, levels = c("non_remitter", "remitter"))
rownames(sd_pre) <- sd_pre$sample_id
sample_data(ps_pre) <- sample_data(sd_pre)

cat("Remission groups:\n")
print(table(sample_data(ps_pre)$remission))

# --- DESeq2 ------------------------------------------------------------------

dds <- phyloseq_to_deseq2(ps_pre, ~ remission)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Positive log2FC = enriched in non-remitters (reference = remitter)
res <- results(dds,
               contrast      = c("remission", "non_remitter", "remitter"),
               cooksCutoff   = FALSE,
               pAdjustMethod = "fdr")

cat("\nDESeq2 results summary:\n")
summary(res)

# --- Build results table with taxonomy ---------------------------------------

sigtab <- as.data.frame(res)
sigtab <- cbind(sigtab, as.data.frame(tax_table(ps_pre)))

cat("\nSignificant ASVs (FDR < 0.05):", sum(sigtab$padj < 0.05, na.rm = TRUE), "\n")
cat("  Enriched in non-remitters (log2FC > 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Enriched in remitters     (log2FC < 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange < 0, na.rm = TRUE), "\n\n")

out_xlsx <- file.path(BASE_DIR, "results/deseq2/deseq2_remission_results.xlsx")
write.xlsx(sigtab, out_xlsx, rowNames = TRUE, overwrite = TRUE)
message("Written: ", out_xlsx)

# --- Volcano plot ------------------------------------------------------------

alpha <- 0.05

df <- sigtab
df$ASV    <- rownames(df)
df$Phylum <- ifelse(is.na(df$Phylum) | df$Phylum == "", "Unassigned", df$Phylum)
df$Genus  <- ifelse(is.na(df$Genus)  | df$Genus  == "", NA, df$Genus)

# Use raw pvalue on y-axis: padj collapses many ASVs to the same floor value
# (when padj ≈ 0), losing vertical spread. pvalue preserves ranking.
# The FDR threshold line is drawn at the pvalue of the least significant
# hit that still passes FDR < 0.05.
df$pval_safe <- ifelse(is.na(df$pvalue) | df$pvalue <= 0, .Machine$double.eps, df$pvalue)

sig_mask    <- !is.na(df$padj) & df$padj < alpha
fdr_pval_threshold <- if (any(sig_mask)) max(df$pvalue[sig_mask], na.rm = TRUE) else alpha

df$label <- ifelse(is.na(df$Genus),
                   paste0(df$ASV, "_", df$Phylum),
                   paste0(df$ASV, "_", df$Genus))

sig <- subset(df, sig_mask)
ns  <- df[!sig_mask, ]

sig_lab  <- subset(sig, abs(log2FoldChange) > 1)
sig_up   <- sig_lab[sig_lab$log2FoldChange > 0, ]
sig_down <- sig_lab[sig_lab$log2FoldChange < 0, ]
sig_up   <- sig_up[order(sig_up$pval_safe,    -abs(sig_up$log2FoldChange)),   ]
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
  labs(x     = "log2 fold change (non-remitters vs remitters)",
       y     = "-log10 p-value (FDR threshold line at padj = 0.05)",
       color = "Phylum")

out_png <- file.path(BASE_DIR, "results/deseq2/volcano_remission.png")
ggsave(out_png, p, width = 30, height = 25, units = "cm", dpi = 300)
message("Written: ", out_png)
