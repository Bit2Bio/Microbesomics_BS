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
#   Positive log2FC = enriched in remitters.
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"

library(phyloseq)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(openxlsx)

dir.create(file.path(BASE_DIR, "results/deseq2"), recursive = TRUE, showWarnings = FALSE)

# --- Load data ---------------------------------------------------------------

ps   <- readRDS(file.path(BASE_DIR, "data/processed/ps_filt_5vs6.rds"))
clin <- read.csv(file.path(BASE_DIR, "data/processed/clinical_all.csv"),
                 check.names = FALSE, na.strings = "")

cat("Loaded phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")

# --- Determine remission status per patient ----------------------------------

# MetS status at PRE (group 5) and POST (group 6) for the complete pairs
pre  <- clin[clin[["group iga"]] == "5", c("pz", "MetS")]
post <- clin[clin[["group iga"]] == "6", c("pz", "MetS")]

# Only patients with metagenomics at both timepoints (same filter as script 03)
complete_pz <- intersect(pre$pz, post$pz)
pre  <- pre[pre$pz  %in% complete_pz, ]
post <- post[post$pz %in% complete_pz, ]

merged <- merge(pre, post, by = "pz", suffixes = c("_pre", "_post"))

# Three-way classification
merged$remission <- ifelse(
  merged$MetS_pre == 1 & merged$MetS_post == 0, "remitter",
  ifelse(
    merged$MetS_pre == 1 & merged$MetS_post == 1, "non_remitter",
    NA  # MetS- at PRE: excluded
  )
)

cat("\nRemission status (complete pairs):\n")
print(table(merged$remission, useNA = "ifany"))

# Patients to include in the DESeq2 comparison (exclude MetS- at PRE)
keep_pz <- merged[!is.na(merged$remission), c("pz", "remission")]

# --- Subset phyloseq to PRE samples of remitters/non-remitters ---------------

sd <- data.frame(sample_data(ps), check.names = FALSE)

keep <- sd[["group.iga"]] == "5" & sd[["pz"]] %in% keep_pz$pz
ps_pre <- prune_samples(keep, ps)

cat("\nPRE samples after exclusion of MetS- at PRE:", nsamples(ps_pre), "\n")

# Add remission column to sample_data
sd_pre <- data.frame(sample_data(ps_pre), check.names = FALSE)
sd_pre <- merge(sd_pre, keep_pz, by = "pz", all.x = TRUE)
rownames(sd_pre) <- sd_pre$sample_id
sd_pre$remission <- factor(sd_pre$remission,
                           levels = c("non_remitter", "remitter"))
sample_data(ps_pre) <- sample_data(sd_pre)

cat("Remission groups:\n")
print(table(sample_data(ps_pre)$remission))

# --- DESeq2 ------------------------------------------------------------------

dds <- phyloseq_to_deseq2(ps_pre, ~ remission)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

# Positive log2FC = enriched in remitters (reference = non_remitter)
res <- results(dds,
               contrast      = c("remission", "remitter", "non_remitter"),
               cooksCutoff   = FALSE,
               pAdjustMethod = "fdr")

# Shrink log2FC estimates with apeglm to stabilise noisy fold changes.
# With only 6 non-remitters, raw log2FC for low-count ASVs can be very large
# but unreliable; apeglm pulls them toward zero proportionally to uncertainty,
# so the volcano correctly shows the most significant ASVs at the extremes.
res <- lfcShrink(dds,
                 coef = "remission_remitter_vs_non_remitter",
                 type = "apeglm",
                 res  = res)

cat("\nDESeq2 results summary (after lfcShrink):\n")
summary(res)

# --- Build results table with taxonomy ---------------------------------------

sigtab <- as.data.frame(res)
sigtab <- cbind(sigtab, as.data.frame(tax_table(ps_pre)))

cat("\nSignificant ASVs (FDR < 0.05):", sum(sigtab$padj < 0.05, na.rm = TRUE), "\n")
cat("  Enriched in remitters    (log2FC > 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Enriched in non-remitters (log2FC < 0):",
    sum(sigtab$padj < 0.05 & sigtab$log2FoldChange < 0, na.rm = TRUE), "\n\n")

out_xlsx <- file.path(BASE_DIR, "results/deseq2/deseq2_remission_results.xlsx")
write.xlsx(sigtab, out_xlsx, rowNames = TRUE, overwrite = TRUE)
message("Written: ", out_xlsx)

# --- Volcano plot ------------------------------------------------------------

alpha <- 0.05

df <- sigtab
df$ASV       <- rownames(df)
df$Phylum    <- ifelse(is.na(df$Phylum) | df$Phylum == "", "Unassigned", df$Phylum)
df$Genus     <- ifelse(is.na(df$Genus)  | df$Genus  == "", NA, df$Genus)
# Cap padj at the smallest observed non-zero value (preserves relative ordering
# and avoids -log10(0) = Inf collapsing all extreme points to the same y level)
min_nonzero <- min(df$padj[df$padj > 0], na.rm = TRUE)
df$padj_safe <- ifelse(is.na(df$padj) | df$padj <= 0, min_nonzero * 0.1, df$padj)

df$label <- ifelse(is.na(df$Genus),
                   paste0(df$ASV, "_", df$Phylum),
                   paste0(df$ASV, "_", df$Genus))

sig <- subset(df, padj_safe <= alpha)
ns  <- subset(df, padj_safe >  alpha)

# Select top 15 per side (by significance then |log2FC|) to ensure labels
# appear on both sides of the volcano even when many points share the same
# capped padj value
sig_up   <- sig[sig$log2FoldChange > 0, ]
sig_down <- sig[sig$log2FoldChange < 0, ]
sig_up   <- sig_up[order(sig_up$padj_safe,   -abs(sig_up$log2FoldChange)),   ]
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
  geom_vline(xintercept = 0,  linetype = "dashed") +
  geom_vline(xintercept =  1, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  scale_color_manual(values = phylum_pal) +
  coord_cartesian(ylim = c(0, max(-log10(df$padj_safe), na.rm = TRUE) * 1.05)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x     = "log2 fold change (remitters vs non-remitters)",
       y     = "-log10 adjusted p-value",
       color = "Phylum")

out_png <- file.path(BASE_DIR, "results/deseq2/volcano_remission.png")
ggsave(out_png, p, width = 30, height = 25, units = "cm", dpi = 300)
message("Written: ", out_png)
