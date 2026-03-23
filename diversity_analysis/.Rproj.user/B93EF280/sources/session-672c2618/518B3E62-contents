#install.packages("BiocManager")
#BiocManager::install("biomformat")
#BiocManager::install("phyloseq")
#BiocManager::install("ALDEx2")
# ============================================================
# Assignment 3: Shotgun metagenomics diversity analysis
# Vegan (n = 3) vs Omnivore (n = 3)
# Based on BIOM output generated from Kraken2/Bracken pipeline
# ============================================================

# ----------------------------
# 1. Load libraries
# ----------------------------
library(biomformat)
library(phyloseq)
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(vegan)
library(ggrepel)
library(patchwork)

# Optional differential abundance package
# Install first if needed:
# install.packages("BiocManager")
# BiocManager::install("ALDEx2")
library(ALDEx2)

# ----------------------------
# 2. Set file paths
# ----------------------------
base_dir <- "~/Desktop/6110/assignment3"

biom_file <- file.path(base_dir, "results", "biom", "table.biom")
metadata_file <- file.path(base_dir, "data", "metadata.tsv")

fig_dir <- file.path(base_dir, "figs")
res_dir <- file.path(base_dir, "results", "R")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 3. Load BIOM and create phyloseq object
# ----------------------------
biom_data <- read_biom(biom_file)
ps <- import_biom(biom_data)

# Keep only Bracken species samples if both raw Kraken and Bracken samples exist
ps <- prune_samples(grepl("_bracken_species$", sample_names(ps)), ps)

# Clean sample names to match metadata
sample_names(ps) <- gsub("_bracken_species$", "", sample_names(ps))

# ----------------------------
# 4. Create / load metadata
# ----------------------------
# If metadata file does not exist, create it automatically
if (!file.exists(metadata_file)) {
  metadata_df <- data.frame(
    Sample_ID = c("SRR8146972", "SRR8146973", "SRR8146974",
                  "SRR8146975", "SRR8146976", "SRR8146977"),
    Diet = c("Omnivore", "Vegan", "Vegan",
             "Omnivore", "Omnivore", "Vegan")
  )
  
  write.table(
    metadata_df,
    metadata_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

metadata_df <- read_tsv(metadata_file, show_col_types = FALSE)
metadata_df <- as.data.frame(metadata_df)

rownames(metadata_df) <- metadata_df$Sample_ID
metadata_df$Sample_ID <- NULL

# Ensure sample names match
stopifnot(all(sample_names(ps) %in% rownames(metadata_df)))

# Reorder metadata to match phyloseq sample order
metadata_df <- metadata_df[sample_names(ps), , drop = FALSE]

# Attach metadata
sample_data(ps) <- metadata_df
sample_data(ps)$Diet <- as.factor(sample_data(ps)$Diet)

# ----------------------------
# 5. Clean taxonomy table
# ----------------------------
# ----------------------------
# 5. Clean taxonomy table
# ----------------------------
tax_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)

# Add missing rank columns if needed
expected_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

if (ncol(tax_df) < length(expected_names)) {
  for (i in (ncol(tax_df) + 1):length(expected_names)) {
    tax_df[[ expected_names[i] ]] <- NA
  }
}

colnames(tax_df)[1:7] <- expected_names

# Reorder to standard rank order
tax_df <- tax_df[, expected_names]

# Remove prefixes like d__, p__, c__, etc.
tax_df[] <- lapply(tax_df, function(x) gsub("^[a-z]__", "", x))

# Replace blanks / NA
tax_df[is.na(tax_df)] <- "Unknown"
tax_df[tax_df == ""] <- "Unknown"

# Put back into phyloseq
tax_table(ps) <- tax_table(as.matrix(tax_df))

# ----------------------------
# 6. Rarefaction
# ----------------------------
otu_df <- as.data.frame(t(otu_table(ps)))

png(file.path(fig_dir, "01_rarefaction_curve.png"), width = 3000, height = 2000, res = 300)
rarecurve(
  otu_df,
  step = 100,
  label = TRUE,
  xlab = "Sequencing Depth",
  ylab = "Number of Observed Taxa"
)
dev.off()

# ----------------------------
# 7. Relative abundance at Phylum level
# ----------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_phy <- tax_glom(ps_rel, taxrank = "Phylum")
df_phy <- psmelt(ps_phy)

# Top 10 phyla
top_phyla <- df_phy %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:10) %>%
  pull(Phylum)

df_phy$Phylum <- ifelse(df_phy$Phylum %in% top_phyla, df_phy$Phylum, "Other")

p_phylum <- ggplot(
  df_phy,
  aes(
    x = Sample,
    y = Abundance,
    fill = fct_reorder(Phylum, Abundance, .fun = sum, .desc = TRUE)
  )
) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    y = "Relative Abundance",
    x = "Sample",
    fill = "Phylum"
  ) +
  facet_wrap(~Diet, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fig_dir, "02_relative_abundance_phylum.png"),
  plot = p_phylum,
  width = 10,
  height = 6.5,
  dpi = 300
)

# ----------------------------
# 8. Alpha diversity
# ----------------------------
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)
alpha_div$Diet <- sample_data(ps)$Diet

write.csv(alpha_div, file.path(res_dir, "alpha_diversity.csv"), row.names = FALSE)

p_shannon <- ggplot(alpha_div, aes(x = Sample, y = Shannon, color = Diet, label = round(Shannon, 3))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  facet_wrap(~Diet, scales = "free_x") +
  labs(y = "Shannon Diversity", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_simpson <- ggplot(alpha_div, aes(x = Sample, y = Simpson, color = Diet, label = round(Simpson, 3))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  facet_wrap(~Diet, scales = "free_x") +
  labs(y = "Simpson Diversity", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fig_dir, "03_alpha_diversity_shannon.png"),
  plot = p_shannon,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "04_alpha_diversity_simpson.png"),
  plot = p_simpson,
  width = 10,
  height = 6,
  dpi = 300
)

# Optional group comparison
shannon_test <- wilcox.test(Shannon ~ Diet, data = alpha_div)
simpson_test <- wilcox.test(Simpson ~ Diet, data = alpha_div)

capture.output(shannon_test, file = file.path(res_dir, "shannon_wilcox.txt"))
capture.output(simpson_test, file = file.path(res_dir, "simpson_wilcox.txt"))

# ----------------------------
# 9. Beta diversity
# ----------------------------
ord_bray <- ordinate(ps, method = "PCoA", distance = "bray")
p_bray <- plot_ordination(ps, ord_bray, type = "samples", color = "Diet") +
  geom_point(size = 4) +
  labs(title = "Bray-Curtis PCoA") +
  theme_minimal()

ord_jaccard <- ordinate(ps, method = "PCoA", distance = "jaccard")
p_jaccard <- plot_ordination(ps, ord_jaccard, type = "samples", color = "Diet") +
  geom_point(size = 4) +
  labs(title = "Jaccard PCoA") +
  theme_minimal()

p_beta <- p_bray + p_jaccard + plot_layout(guides = "collect")

ggsave(
  file.path(fig_dir, "05_beta_diversity_pcoa.png"),
  plot = p_beta,
  width = 10,
  height = 6.5,
  dpi = 300
)

# PERMANOVA
metadata_beta <- as(sample_data(ps), "data.frame")

set.seed(123)
permanova_bray <- adonis2(
  phyloseq::distance(ps, method = "bray") ~ Diet,
  data = metadata_beta,
  permutations = 999
)

set.seed(123)
permanova_jaccard <- adonis2(
  phyloseq::distance(ps, method = "jaccard") ~ Diet,
  data = metadata_beta,
  permutations = 999
)

capture.output(permanova_bray, file = file.path(res_dir, "permanova_bray.txt"))
capture.output(permanova_jaccard, file = file.path(res_dir, "permanova_jaccard.txt"))

# ----------------------------
# 10. Differential abundance with ALDEx2
# ----------------------------
# Use raw counts
count_matrix <- as(otu_table(ps), "matrix")

if (!taxa_are_rows(ps)) {
  count_matrix <- t(count_matrix)
}

# Diet vector must align with columns (samples)
conds <- sample_data(ps)$Diet
conds <- as.character(conds)

set.seed(123)
aldex_out <- aldex.clr(count_matrix, conds, mc.samples = 128, denom = "all", verbose = FALSE)
aldex_tt <- aldex.ttest(aldex_out)
aldex_eff <- aldex.effect(aldex_out)

aldex_res <- cbind(
  data.frame(feature = rownames(aldex_tt)),
  aldex_tt,
  aldex_eff
)

# Attach taxonomy
tax_df <- as.data.frame(tax_table(ps))
tax_df$feature <- rownames(tax_df)

aldex_res <- left_join(aldex_res, tax_df, by = "feature")

write.csv(aldex_res, file.path(res_dir, "aldex2_results.csv"), row.names = FALSE)

# Volcano-style plot using effect vs adjusted p-value
aldex_res$significant <- aldex_res$wi.eBH < 0.05
aldex_res$label <- ifelse(
  aldex_res$significant & !is.na(aldex_res$Species),
  aldex_res$Species,
  ""
)

p_daa <- ggplot(aldex_res, aes(x = diff.btw, y = -log10(wi.eBH))) +
  geom_point(aes(color = significant), alpha = 0.8, size = 2.5) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    x = "Difference between groups",
    y = "-log10 adjusted p-value",
    color = "Significant"
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "06_differential_abundance_aldex2.png"),
  plot = p_daa,
  width = 9,
  height = 6,
  dpi = 300
)

# ----------------------------
# 11. Save phyloseq object
# ----------------------------
saveRDS(ps, file = file.path(res_dir, "phyloseq_object.rds"))

cat("Analysis complete.\n")