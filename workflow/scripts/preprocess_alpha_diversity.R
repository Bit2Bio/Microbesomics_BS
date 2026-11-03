library(phyloseq)
library(yaml)

CONFIG   <- "/home/lorenzo/Microbesomics/config/config.yaml"
PS_RAW   <- "/home/lorenzo/Microbesomics/data/processed/ps_raw.rds"
OUT_PATH <- "/home/lorenzo/Microbesomics/data/processed/ps_rarefied_alpha.rds"

cfg        <- yaml.load_file(CONFIG)
alpha_cfg  <- cfg$alpha_diversity
excl_pairs <- alpha_cfg$exclude_pairs     # patient IDs (pz column)
rare_depth <- alpha_cfg$rarefaction_depth
seed       <- alpha_cfg$seed

cat(sprintf("Config: exclude_pairs=%s | rarefaction_depth=%d | seed=%d\n",
            paste(excl_pairs, collapse=", "), rare_depth, seed))

ps <- readRDS(PS_RAW)
cat(sprintf("ps_raw: %d samples, %d ASVs\n", nsamples(ps), ntaxa(ps)))

# ── 1. Attach clinical metadata (ps_raw has no sample_data) ─────────────────
clinical <- read.csv("/home/lorenzo/Microbesomics/data/processed/clinical_all.csv",
                     na.strings = "")
rownames(clinical) <- clinical$sample_id
# keep only samples present in the phyloseq object
clinical <- clinical[intersect(rownames(clinical), sample_names(ps)), ]
sample_data(ps) <- sample_data(clinical)

# ── 2. Subset to the bariatric samples (group 5 = PRE, group 6 = POST) ──────
meta <- as(sample_data(ps), "data.frame")
keep_bariatric <- rownames(meta)[meta$group.iga %in% c(5, 6)]
ps <- prune_samples(keep_bariatric, ps)
cat(sprintf("After bariatric subset: %d samples\n", nsamples(ps)))

# ── 3. Exclude flagged pairs (both PRE and POST of each excluded patient) ────
meta <- as(sample_data(ps), "data.frame")
excl_samples <- rownames(meta)[meta$pz %in% excl_pairs]
if (length(excl_samples) > 0) {
  cat(sprintf("Excluding samples: %s (pair: %s)\n",
              paste(excl_samples, collapse=", "),
              paste(excl_pairs, collapse=", ")))
  ps <- prune_samples(!sample_names(ps) %in% excl_samples, ps)
}
cat(sprintf("After pair exclusion: %d samples\n", nsamples(ps)))

# ── 3. No ASV filtering — alpha diversity requires the full community ────────
# (filtering would artificially reduce observed richness and bias Shannon)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)   # only drop ASVs with zero counts across all retained samples
cat(sprintf("ASVs retained (sum > 0): %d\n", ntaxa(ps)))

# ── 4. Rarefy to even depth ──────────────────────────────────────────────────
lib_sizes <- sample_sums(ps)
cat(sprintf("Library sizes before rarefaction: min=%d, median=%d, max=%d\n",
            min(lib_sizes), as.integer(median(lib_sizes)), max(lib_sizes)))

n_below <- sum(lib_sizes < rare_depth)
if (n_below > 0) {
  warning(sprintf("%d sample(s) below rarefaction depth (%d) — they will be dropped by rarefy_even_depth.",
                  n_below, rare_depth))
  cat(sprintf("Samples below threshold: %s\n",
              paste(names(lib_sizes[lib_sizes < rare_depth]), collapse=", ")))
}

set.seed(seed)
ps_rare <- rarefy_even_depth(ps, sample.size = rare_depth,
                              rngseed = seed, replace = FALSE, verbose = TRUE)
cat(sprintf("After rarefaction: %d samples, %d ASVs\n",
            nsamples(ps_rare), ntaxa(ps_rare)))

saveRDS(ps_rare, OUT_PATH)
cat(sprintf("Written: %s\n", OUT_PATH))
