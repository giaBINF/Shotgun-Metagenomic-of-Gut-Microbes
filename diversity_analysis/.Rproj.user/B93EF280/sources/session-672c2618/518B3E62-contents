# Assignment 3: Shotgun metagenomics diversity analysis

# Load libraries
library(biomformat)
library(phyloseq)
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(vegan)
library(patchwork)
library(ALDEx2)

# File paths
base_dir  <- "~/Desktop/6110/assignment3"
biom_file <- file.path(base_dir, "results", "biom", "table.biom")
meta_file <- file.path(base_dir, "data", "metadata.tsv")
fig_dir   <- file.path(base_dir, "figs")
res_dir   <- file.path(base_dir, "results", "R")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
n_top_phyla    <- 6
diff_threshold <- 5

diet_colors <- c(
  "Omnivore" = "#D55E00",
  "Vegan"    = "#0072B2"
)

wong_palette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00"
)

# Load BIOM and metadata
biom_data <- read_biom(biom_file)
ps <- import_biom(biom_data)

# Keep only Bracken-derived samples
ps <- prune_samples(grepl("_bracken_species$", sample_names(ps)), ps)

# Clean sample names
sample_names(ps) <- gsub("_bracken_species$", "", sample_names(ps))

meta_df <- read_tsv(meta_file, show_col_types = FALSE) %>%
  as.data.frame()

rownames(meta_df) <- meta_df$Sample_ID
meta_df$Sample_ID <- NULL

sample_data(ps) <- meta_df
sample_data(ps)$Diet <- as.factor(sample_data(ps)$Diet)

# Clean taxonomy table
tax_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)

expected_names <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

# Add missing columns
if (ncol(tax_df) < length(expected_names)) {
  for (i in (ncol(tax_df) + 1):length(expected_names)) {
    tax_df[[expected_names[i]]] <- NA
  }
}

colnames(tax_df)[1:7] <- expected_names
tax_df <- tax_df[, expected_names]

tax_df[] <- lapply(tax_df, function(x) gsub("^[a-z]__", "", x))

tax_df[is.na(tax_df)] <- "Unknown"
tax_df[tax_df == ""] <- "Unknown"

tax_table(ps) <- tax_table(as.matrix(tax_df))

ps <- subset_taxa(ps, !Phylum %in% c("Chordata", "Unknown", ""))

base_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "grey92"),
    panel.grid.major = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

legend_theme <- theme(
  legend.position = "right",
  legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(fill = "white", color = NA)
)

# 1. Relative abundance at phylum level
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_phy <- tax_glom(ps_rel, taxrank = "Phylum")
df_phy <- psmelt(ps_phy)

top_phyla <- df_phy %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice(1:n_top_phyla) %>%
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
  scale_fill_manual(values = wong_palette) +
  labs(
    y = "Relative Abundance",
    x = "Sample",
    fill = "Phylum"
  ) +
  facet_wrap(~Diet, scales = "free_x") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  file.path(fig_dir, "01_relative_abundance_phylum.png"),
  plot = p_phylum,
  width = 14,
  height = 6.5,
  dpi = 300
)

# 2. Alpha diversity
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)
alpha_div$Diet   <- sample_data(ps)$Diet

write.csv(alpha_div, file.path(res_dir, "alpha_diversity.csv"), row.names = FALSE)

p_shannon <- ggplot(alpha_div, aes(x = Diet, y = Shannon, color = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.1, size = 3) +
  scale_color_manual(values = diet_colors, name = "Diet") +
  labs(title = "Shannon Diversity", x = "Diet", y = "Shannon Index") +
  base_theme +
  theme(legend.position = "none")

p_simpson <- ggplot(alpha_div, aes(x = Diet, y = Simpson, color = Diet)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.1, size = 3) +
  scale_color_manual(values = diet_colors, name = "Diet") +
  labs(title = "Simpson Diversity", x = "Diet", y = "Simpson Index") +
  base_theme +
  legend_theme

ggsave(
  file.path(fig_dir, "02_alpha_diversity.png"),
  plot = p_shannon + p_simpson + plot_layout(guides = "collect"),
  width = 10,
  height = 6,
  dpi = 300
)

capture.output(
  wilcox.test(Shannon ~ Diet, data = alpha_div),
  file = file.path(res_dir, "shannon_wilcox.txt")
)

capture.output(
  wilcox.test(Simpson ~ Diet, data = alpha_div),
  file = file.path(res_dir, "simpson_wilcox.txt")
)

# 3. Beta diversity
meta_plot_df <- as(sample_data(ps), "data.frame")
meta_plot_df$SampleID <- rownames(meta_plot_df)

make_pcoa_plot <- function(ord, title, meta_df) {
  coords <- as.data.frame(ord$vectors[, 1:2])
  colnames(coords) <- c("Axis1", "Axis2")
  coords$SampleID <- rownames(coords)
  coords <- left_join(coords, meta_df, by = "SampleID")
  
  pct <- round(ord$values$Relative_eig[1:2] * 100, 1)
  
  ggplot(coords, aes(x = Axis1, y = Axis2, color = Diet)) +
    geom_point(size = 4) +
    scale_color_manual(values = diet_colors, name = "Diet") +
    labs(
      title = title,
      x = paste0("PC1 [", pct[1], "%]"),
      y = paste0("PC2 [", pct[2], "%]")
    ) +
    base_theme
}

ord_bray <- ordinate(ps, method = "PCoA", distance = "bray")
ord_jaccard <- ordinate(ps, method = "PCoA", distance = "jaccard")

p_bray <- make_pcoa_plot(ord_bray, "Bray-Curtis PCoA", meta_plot_df) +
  theme(legend.position = "none")

p_jaccard <- make_pcoa_plot(ord_jaccard, "Jaccard PCoA", meta_plot_df) +
  legend_theme

ggsave(
  file.path(fig_dir, "03_beta_diversity_pcoa.png"),
  plot = p_bray + p_jaccard + plot_layout(guides = "collect"),
  width = 12,
  height = 6,
  dpi = 300
)

meta_beta <- as(sample_data(ps), "data.frame")

capture.output(
  adonis2(phyloseq::distance(ps, method = "bray") ~ Diet,
          data = meta_beta, permutations = 999),
  file = file.path(res_dir, "permanova_bray.txt")
)

capture.output(
  adonis2(phyloseq::distance(ps, method = "jaccard") ~ Diet,
          data = meta_beta, permutations = 999),
  file = file.path(res_dir, "permanova_jaccard.txt")
)

# 4. Differential abundance with ALDEx2
count_matrix <- as(otu_table(ps), "matrix")

if (!taxa_are_rows(ps)) {
  count_matrix <- t(count_matrix)
}

conds <- as.character(sample_data(ps)$Diet)

set.seed(123)

aldex_out <- aldex.clr(
  count_matrix,
  conds,
  mc.samples = 128,
  denom = "all",
  verbose = FALSE
)

aldex_tt  <- aldex.ttest(aldex_out)
aldex_eff <- aldex.effect(aldex_out)

aldex_res <- cbind(
  data.frame(feature = rownames(aldex_tt)),
  aldex_tt,
  aldex_eff
)

tax_df2 <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
tax_df2$feature <- rownames(tax_df2)

aldex_res <- left_join(aldex_res, tax_df2, by = "feature")

write.csv(aldex_res, file.path(res_dir, "aldex2_results.csv"), row.names = FALSE)

# 5. Differential abundance forest plot
forest_data <- aldex_res %>%
  filter(!is.na(Species), !is.na(diff.btw), abs(diff.btw) > diff_threshold) %>%
  group_by(Species) %>%
  slice_max(abs(diff.btw), n = 1) %>%
  ungroup() %>%
  arrange(desc(diff.btw)) %>%
  mutate(
    Scientific_Name = ifelse(
      grepl("Unknown|NA", Genus),
      Species,
      paste(Genus, Species)
    ),
    sig = ifelse(wi.eBH < 0.05, "< 0.05", "≥ 0.05"),
    Scientific_Name = factor(Scientific_Name, levels = unique(Scientific_Name))
  )

cat("Species in forest plot:", nrow(forest_data), "\n")

p_forest <- ggplot(forest_data, aes(x = diff.btw, y = Scientific_Name)) +
  geom_vline(xintercept = 0, color = "red", linewidth = 1) +
  geom_errorbarh(
    aes(xmin = diff.btw - diff.win, xmax = diff.btw + diff.win),
    height = 0.4,
    linewidth = 0.8
  ) +
  geom_point(aes(fill = sig), shape = 21, size = 4, color = "black") +
  scale_fill_manual(
    values = c("< 0.05" = "#D55E00", "≥ 0.05" = "grey60"),
    name = "Adjusted p-value"
  ) +
  labs(
    x = "Estimated CLR Abundance Difference (Vegan vs. Omnivore)",
    y = "Species"
  ) +
  base_theme +
  legend_theme +
  theme(
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(15, 20, 15, 15)
  )

ggsave(
  file.path(fig_dir, "04_differential_abundance_forest.png"),
  plot = p_forest,
  width = 10,
  height = max(6, nrow(forest_data) * 0.35),
  dpi = 300
)

# 6. Save phyloseq object
saveRDS(ps, file = file.path(res_dir, "phyloseq_object.rds"))

cat("Analysis complete.\n")