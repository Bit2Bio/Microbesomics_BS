# =============================================================================
# 07_heatmap_intersection.R
#
# Heatmap of log2FC across the three DESeq2 analyses, restricted to ASVs
# significant (FDR < 0.05) in at least one of the two baseline comparisons:
#   - Remission  : remitters vs non-remitters at PRE (script 05)
#   - MetS       : MetS+ vs MetS- at PRE            (script 06)
#
# Layout:
#   Columns = selected ASVs (clustered by log2FC profile)
#   Rows    = 3 comparisons, fixed order:
#               1. Remission   (remitters vs non-remitters)
#               2. MetS        (MetS+ vs MetS-)
#               3. Surgery     (POST vs PRE)
#   Fill    = log2FC (diverging palette, blue = negative / red = positive)
#   Overlay = * where padj < 0.05, blank otherwise
#
# Input:
#   results/deseq2/deseq2_remission_results.xlsx
#   results/deseq2/deseq2_mets_baseline_results.xlsx
#   results/deseq2/deseq2_5vs6_results.xlsx
#
# Output:
#   results/deseq2/heatmap_intersection.png
# =============================================================================

BASE_DIR <- "/home/lorenzo/Microbesomics"
RES_DIR  <- file.path(BASE_DIR, "results/deseq2")

library(openxlsx)
library(ComplexHeatmap)
library(circlize)

# --- Load DESeq2 results -----------------------------------------------------

res_rem  <- read.xlsx(file.path(RES_DIR, "deseq2_remission_results.xlsx"),    rowNames = TRUE)
res_bas  <- read.xlsx(file.path(RES_DIR, "deseq2_mets_baseline_results.xlsx"), rowNames = TRUE)
res_surg <- read.xlsx(file.path(RES_DIR, "deseq2_5vs6_results.xlsx"),          rowNames = TRUE)

# --- Feature selection: significant in remission OR MetS baseline ------------

alpha <- 0.05

sig_rem <- rownames(res_rem)[!is.na(res_rem$padj)  & res_rem$padj  < alpha]
sig_bas <- rownames(res_bas)[!is.na(res_bas$padj)  & res_bas$padj  < alpha]
features <- union(sig_rem, sig_bas)

cat("Significant ASVs:\n")
cat("  Remission:    ", length(sig_rem), "\n")
cat("  MetS baseline:", length(sig_bas), "\n")
cat("  Union (cols): ", length(features), "\n\n")

# --- Build log2FC matrix (rows = comparisons, cols = ASVs) -------------------

row_names <- c(
  "Remission\n(non-remitters vs remitters)",
  "MetS baseline\n(MetS+ vs MetS-)",
  "Surgery\n(POST vs PRE)"
)

lfc_mat  <- matrix(NA_real_, nrow = 3, ncol = length(features),
                   dimnames = list(row_names, features))
padj_mat <- lfc_mat

for (asv in features) {
  if (asv %in% rownames(res_rem))  { lfc_mat[1, asv] <- res_rem[asv,  "log2FoldChange"]; padj_mat[1, asv] <- res_rem[asv,  "padj"] }
  if (asv %in% rownames(res_bas))  { lfc_mat[2, asv] <- res_bas[asv,  "log2FoldChange"]; padj_mat[2, asv] <- res_bas[asv,  "padj"] }
  if (asv %in% rownames(res_surg)) { lfc_mat[3, asv] <- res_surg[asv, "log2FoldChange"]; padj_mat[3, asv] <- res_surg[asv, "padj"] }
}

# --- Column labels: ASV_Genus (or ASV_Phylum if Genus is NA) ----------------

tax <- data.frame(
  Genus  = rep(NA_character_, length(features)),
  Phylum = rep(NA_character_, length(features)),
  row.names = features,
  stringsAsFactors = FALSE
)
for (asv in features) {
  for (res in list(res_rem, res_bas, res_surg)) {
    if (asv %in% rownames(res)) {
      if (is.na(tax[asv, "Genus"])  && !is.na(res[asv, "Genus"])  && res[asv, "Genus"]  != "") tax[asv, "Genus"]  <- res[asv, "Genus"]
      if (is.na(tax[asv, "Phylum"]) && !is.na(res[asv, "Phylum"]) && res[asv, "Phylum"] != "") tax[asv, "Phylum"] <- res[asv, "Phylum"]
    }
  }
}

col_labels <- ifelse(
  !is.na(tax$Genus) & tax$Genus != "",
  paste0(features, "_", tax$Genus),
  paste0(features, "_", tax$Phylum)
)

# --- Significance stars (padj < 0.05) ----------------------------------------

stars_mat <- ifelse(!is.na(padj_mat) & padj_mat < alpha, "*", "")

# --- Cap log2FC for visualization only (extreme values from unshrunken DESeq2
#     estimates distort the color scale; underlying data in lfc_mat unchanged) -

CAP <- 10
lfc_mat_plot <- pmin(pmax(lfc_mat, -CAP), CAP)

# --- Color scale (symmetric around 0) ----------------------------------------

lfc_max <- CAP

col_fun <- colorRamp2(
  c(-lfc_max, 0, lfc_max),
  c("#2166AC", "white", "#D6604D")
)

# --- Heatmap -----------------------------------------------------------------

cell_w_cm <- 0.28
cell_h_cm <- 1.2

ht <- Heatmap(
  lfc_mat_plot,
  name             = "log2FC",
  col              = col_fun,
  na_col           = "grey88",
  column_labels    = col_labels,
  column_names_gp  = gpar(fontsize = 7),
  column_names_rot = 45,
  row_names_gp     = gpar(fontsize = 9, fontface = "bold"),
  row_names_side   = "left",
  cluster_rows     = FALSE,
  cluster_columns  = TRUE,
  show_column_dend = FALSE,
  cell_fun         = function(j, i, x, y, width, height, fill) {
    if (stars_mat[i, j] != "")
      grid.text("*", x, y, gp = gpar(fontsize = 11, col = "black", fontface = "bold"))
  },
  heatmap_legend_param = list(
    title         = "log2FC",
    at            = c(-lfc_max, 0, lfc_max),
    legend_height = unit(3, "cm")
  ),
  width  = unit(cell_w_cm * length(features), "cm"),
  height = unit(cell_h_cm * 3,                "cm")
)

out_png <- file.path(RES_DIR, "heatmap_intersection.png")
png(out_png,
    width  = cell_w_cm * length(features) + 10,
    height = cell_h_cm * 3 + 5,
    units  = "cm", res = 300)
draw(ht, padding = unit(c(5, 5, 5, 40), "mm"))
dev.off()

message("Written: ", out_png)
cat("ASVs in heatmap:", length(features), "\n")
