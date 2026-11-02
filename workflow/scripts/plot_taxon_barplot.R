# =============================================================================
# R/plot_taxon_barplot.R
#
# plot_taxon_barplot(): stacked barplot of mean relative abundance per group
# at a given taxonomic rank.
#
# Arguments:
#   ps        — phyloseq object, already subsetted to samples of interest
#   rank      — taxonomic rank string: "Phylum", "Class", "Order",
#               "Family", "Genus"
#   group_var — name of the sample_data column to group by (e.g. "groupIGA",
#               "remission", "MetS_group")
#   top_n     — number of top taxa to show individually (default 15);
#               remaining taxa are collapsed into "Other"
#   group_labels — optional named character vector to rename group levels for
#               the plot (e.g. c("5"="PRE", "6"="POST"))
#   palette   — optional named or unnamed character vector of colours;
#               if NULL, uses a default palette
#
# Returns: a ggplot object (caller is responsible for saving)
# =============================================================================

plot_taxon_barplot <- function(ps,
                               rank,
                               group_var,
                               top_n       = 15,
                               group_labels = NULL,
                               palette      = NULL) {

  library(phyloseq)
  library(ggplot2)
  library(dplyr)

  # --- Aggregate at chosen rank -----------------------------------------------
  ps_glom <- tax_glom(ps, taxrank = rank, NArm = FALSE)

  # --- Relative abundance -----------------------------------------------------
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))

  # --- Melt to long format ----------------------------------------------------
  df <- psmelt(ps_rel)

  # Use the rank column as the taxon label; replace NA/empty with "Unassigned"
  df$Taxon <- as.character(df[[rank]])
  df$Taxon[is.na(df$Taxon) | df$Taxon == ""] <- "Unassigned"

  # Ensure group variable is a factor
  df$Group <- as.character(df[[group_var]])
  if (!is.null(group_labels)) {
    df$Group <- ifelse(df$Group %in% names(group_labels),
                       group_labels[df$Group],
                       df$Group)
  }

  # --- Select top N taxa by overall mean abundance ----------------------------
  taxon_means <- df %>%
    group_by(Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund))

  top_taxa <- head(taxon_means$Taxon, top_n)

  df$Taxon <- ifelse(df$Taxon %in% top_taxa, df$Taxon, "Other")

  # Reorder factor: top taxa by abundance, then Other at the bottom
  taxon_order <- c(rev(top_taxa), "Other")
  df$Taxon <- factor(df$Taxon, levels = taxon_order)

  # --- Mean relative abundance per group × taxon ------------------------------
  # Must sum within each sample first: psmelt gives one row per ASV, so "Other"
  # has many rows per sample. Averaging without summing first divides by
  # n_other_taxa × n_samples instead of just n_samples → bars don't reach 100%.
  df_sample <- df %>%
    group_by(Sample, Group, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

  df_mean <- df_sample %>%
    group_by(Group, Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

  # --- Colour palette ---------------------------------------------------------
  n_taxa <- length(taxon_order)
  if (is.null(palette)) {
    # Distinct palette: use a combination of Set1/Set2/Set3 or hue_pal
    pal_cols <- c(
      "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
      "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
      "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#B3B3B3"
    )
    # If more taxa than colours, recycle with hue_pal
    if (n_taxa > length(pal_cols)) {
      pal_cols <- c(pal_cols, scales::hue_pal()(n_taxa - length(pal_cols)))
    }
    palette <- setNames(pal_cols[seq_len(n_taxa)], taxon_order)
    # "Other" always grey
    palette["Other"] <- "#CCCCCC"
  }

  # --- Plot -------------------------------------------------------------------
  p <- ggplot(df_mean, aes(x = Group, y = mean_abund, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7,
             color = "black", linewidth = 0.25) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(x    = NULL,
         y    = "Mean relative abundance",
         fill = rank)

  return(p)
}
