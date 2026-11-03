library(phyloseq)
library(yaml)

CONFIG   <- "/home/lorenzo/Microbesomics/config/config.yaml"
PS_PATH  <- "/home/lorenzo/Microbesomics/data/processed/ps_filt_5vs6.rds"
OUT_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_rarefied_beta.rds"

cfg        <- yaml.load_file(CONFIG)
alpha_cfg  <- cfg$alpha_diversity   # reuse same exclusions and depth
excl_pairs <- alpha_cfg$exclude_pairs
rare_depth <- alpha_cfg$rarefaction_depth
seed       <- alpha_cfg$seed

cat(sprintf("Config: exclude_pairs=%s | rarefaction_depth=%d | seed=%d\n",
            paste(excl_pairs, collapse=", "), rare_depth, seed))

ps   <- readRDS(PS_PATH)
meta <- as(sample_data(ps), "data.frame")
cat(sprintf("ps_filt_5vs6: %d samples, %d ASVs\n", nsamples(ps), ntaxa(ps)))

# ── Exclude flagged pairs ────────────────────────────────────────────────────
excl_samples <- rownames(meta)[meta$pz %in% excl_pairs]
if (length(excl_samples) > 0) {
  cat(sprintf("Excluding: %s (pz: %s)\n",
              paste(excl_samples, collapse=", "),
              paste(excl_pairs,   collapse=", ")))
  ps <- prune_samples(!sample_names(ps) %in% excl_samples, ps)
}
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat(sprintf("After exclusion: %d samples, %d ASVs\n", nsamples(ps), ntaxa(ps)))

# ── Rarefy ───────────────────────────────────────────────────────────────────
lib_sizes <- sample_sums(ps)
cat(sprintf("Library sizes: min=%d  median=%d  max=%d\n",
            min(lib_sizes), as.integer(median(lib_sizes)), max(lib_sizes)))

n_below <- sum(lib_sizes < rare_depth)
if (n_below > 0)
  warning(sprintf("%d sample(s) below rarefaction depth — will be dropped.", n_below))

set.seed(seed)
ps_rare <- rarefy_even_depth(ps, sample.size = rare_depth,
                              rngseed = seed, replace = FALSE, verbose = TRUE)
cat(sprintf("After rarefaction: %d samples, %d ASVs\n",
            nsamples(ps_rare), ntaxa(ps_rare)))

saveRDS(ps_rare, OUT_PATH)
cat(sprintf("Written: %s\n", OUT_PATH))
