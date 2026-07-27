# =============================================================================
# 07_deseq2_mets_baseline.R
#
# Differential abundance analysis (DESeq2): MetS+ vs MetS- at baseline
# (PRE, group 5) among the 33 complete bariatric surgery pairs.
#
#   MetS+ : 28 patients (MetS == 1 at PRE)
#   MetS- :  5 patients (MetS == 0 at PRE)
#
# Input:
#   data/processed/ps_filt_5vs6.rds  — 66 samples (33 pairs), 707 ASVs
#                                       same filtering as script 03
#   data/processed/clinical_all.csv  — MetS status for all 247 samples
#
# Output:
#   results/deseq2/deseq2_mets_baseline_results.xlsx — DESeq2 results + taxonomy
#   results/deseq2/volcano_mets_baseline.png          — volcano plot
#
# Design: ~ MetS (simple two-group comparison at the same timepoint)
#   Positive log2FC = enriched in MetS+ relative to MetS-.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)

dir.create(file.path(BASE_DIR, "results/deseq2"), recursive = TRUE, showWarnings = FALSE)

# --- Load data ---------------------------------------------------------------

ps   <- readRDS(file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
clin <- read.csv(file.path(BASE_DIR, "data/processed/clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

cat("Loaded phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# --- Identify complete pairs and subset to PRE --------------------------------

pz_pre  <- clin[clin[["group iga"]] == "5", "pz"]
pz_post <- clin[clin[["group iga"]] == "6", "pz"]
complete_pz <- intersect(pz_pre, pz_post)

sd <- data.frame(sample_data(ps), check.names = FALSE)
keep <- sd[["group.iga"]] == "5" & sd[["pz"]] %in% complete_pz
ps_pre <- prune_samples(keep, ps)

cat("PRE samples (33 complete pairs):", nsamples(ps_pre), "\n")

# --- Add MetS variable -------------------------------------------------------

sd_pre <- data.frame(sample_data(ps_pre), check.names = FALSE)

# MetS column is already in sample_data (computed in script 02)
# Recode to factor with MetS- as reference
sd_pre$MetS_group <- factor(
  ifelse(sd_pre$MetS == 1, "MetS_pos", "MetS_neg"),
  levels = c("MetS_neg", "MetS_pos")
)

sample_data(ps_pre) <- sample_data(sd_pre)

cat("\nMetS distribution at baseline:\n")
print(table(sd_pre$MetS_group, useNA = "ifany"))

# --- DESeq2 ------------------------------------------------------------------

ps_ds <- ps_pre
otu_table(ps_ds) <- otu_table(ps_ds) + 1

dds <- phyloseq_to_deseq2(ps_ds, ~ MetS_group)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Positive log2FC = enriched in MetS+ (reference = MetS_neg)
res <- results(dds,
               contrast      = c("MetS_group", "MetS_pos", "MetS_neg"),
               cooksCutoff   = FALSE,
               pAdjustMethod = "fdr")

cat("\nDESeq2 results summary:\n")
summary(res)

# --- Build results table with taxonomy ---------------------------------------

sigtab <- as.data.frame(res)
sigtab <- cbind(sigtab, as.data.frame(tax_table(ps_ds)))

cat("\nSignificant ASVs (FDR < 0.05):", sum(sigtab$padj < 0.05, na.rm = TRUE), "\n")
cat("  Enriched in MetS+  (log2FC > 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Enriched in MetS-  (log2FC < 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange < 0, na.rm = TRUE), "\n\n")

out_xlsx <- file.path(BASE_DIR, "results/deseq2/deseq2_mets_baseline_results.xlsx")
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

sig <- subset(df, padj_safe <= alpha)
ns  <- subset(df, padj_safe >  alpha)

sig_up   <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]
sig_up   <- sig_up[order(sig_up$padj_safe,    -abs(sig_up$log2FoldChange)),   ]
sig_down <- sig_down[order(sig_down$padj_safe, -abs(sig_down$log2FoldChange)), ]
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  scale_color_manual(values = phylum_pal) +
  coord_cartesian(ylim = c(0, max(-log10(df$padj_safe), na.rm = TRUE) * 1.05)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x     = "log2 fold change (MetS+ vs MetS-)",
       y     = "-log10 adjusted p-value",
       color = "Phylum")

out_png <- file.path(BASE_DIR, "results/deseq2/volcano_mets_baseline.png")
ggsave(out_png, p, width = 30, height = 25, units = "cm", dpi = 300)
message("Written: ", out_png)
