############### Gut microbial SNV ###################################

.libPaths(c("D:/2_Software/R/install/R-4.5.1/library_shuang",
            "D:/2_Software/R/install/R-4.5.1/library"))

# Set working directory
setwd("D:/3_Projects/2_Children_jaundice/3_R_analysis")


### Necessary Packages 
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(readr)
library(pheatmap)
library(tidytext)
library(forcats)
library(tibble)
library(vegan)
library(DirichletMultinomial)


### Load phyloseq data
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/Result4_R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Result4_R.RData")


############ Result: SNV feature-defined microbial genetic states track stepwise infant gut development #####

#### Metadata information 
meta_final <- read_excel("meta_final.xlsx")
meta_final #


#### Fig. 4a Selection of optimal cluster number using Laplace approximation ######

### Fecal samples ~ total SNV burden ~ DMM clustering

filtered_SNVs_coverage5_top5 <- read_tsv("filtered_SNVs_coverage5_top5.tsv", show_col_types = FALSE)
filtered_SNVs_coverage5_top5 #
dim(filtered_SNVs_coverage5_top5)


#### Create harmonized genome identifiers for bin-level mapping

setDT(filtered_SNVs_coverage5_top5)

filtered_SNVs_coverage5_top5[, Genome1 := sub("_[0-9]+$", "", scaffold)]
filtered_SNVs_coverage5_top5[, Genome2 := gsub("_k[0-9]+$", "", Genome1)]
filtered_SNVs_coverage5_top5[, Genome3 := sub("(_genomic).*", "_genomic", Genome2)]
filtered_SNVs_coverage5_top5[, Genome4 := sub("(_\\.\\d+.*)$", "", Genome3)]

filtered_SNVs_coverage5_top5[, Genome5 :=
                               fifelse(grepl("_genomic", Genome4),
                                       sub("(_genomic).*", "\\1", Genome4),
                                       sub("_ERS.*$", "", Genome4))
]

filtered_SNVs_coverage5_top5$genome <- filtered_SNVs_coverage5_top5$Genome1


#### Read genome-level inStrain information

all_samples.genome_info <- read_tsv("all_samples.genome_info.tsv", show_col_types = FALSE)

table(all_samples.genome_info$length == 0)
table(all_samples.genome_info$breadth_minCov == 0)
table(all_samples.genome_info$SNV_count == 0)

setDT(all_samples.genome_info)

all_samples.genome_info[, Genome1 := sub("_[0-9]+$", "", genome)]
all_samples.genome_info[, Genome2 := gsub("_k[0-9]+$", "", Genome1)]
all_samples.genome_info[, Genome3 := sub("(_genomic).*", "_genomic", Genome2)]
all_samples.genome_info[, Genome4 := sub("(_\\.\\d+.*)$", "", Genome3)]

all_samples.genome_info[, Genome5 :=
                          fifelse(grepl("_genomic", Genome4),
                                  sub("(_genomic).*", "\\1", Genome4),
                                  sub("_ERS.*$", "", Genome4))
]

table(filtered_SNVs_coverage5_top5$Genome1 %in% all_samples.genome_info$genome)
table(filtered_SNVs_coverage5_top5$genome %in% all_samples.genome_info$genome)
# TRUE - 1896337


#### Read genome species information

genome_species_info <- read_excel("D:/3_Projects/2_Children_jaundice/3_R_analysis/genome_species_info.xlsx")

genome_species_info2 <- data.frame(genome_species_info) %>%
  dplyr::select(Genome, Taxonomy.lineage..GTDB.)

strip_prefix <- function(x) sub("^[a-z]__", "", x)

genome_species_info3 <- genome_species_info2 %>%
  tidyr::separate(
    col = "Taxonomy.lineage..GTDB.",
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    Taxonomy.lineage..GTDB. = Taxonomy.lineage..GTDB.,
    Domain  = strip_prefix(Domain),
    Phylum  = strip_prefix(Phylum),
    Class   = strip_prefix(Class),
    Order   = strip_prefix(Order),
    Family  = strip_prefix(Family),
    Genus   = na_if(strip_prefix(Genus), ""),
    Species = na_if(strip_prefix(Species), "")
  )


#### Merge SNV table with species-level taxonomy

table(genome_species_info3$Genome %in% filtered_SNVs_coverage5_top5$Genome5)
# FALSE - 2022, TRUE - 150
table(filtered_SNVs_coverage5_top5$Genome5 %in% genome_species_info3$Genome)

setDT(filtered_SNVs_coverage5_top5)
setDT(genome_species_info3)

filtered_SNVs_coverage5_top5_species <- merge(
  filtered_SNVs_coverage5_top5,
  genome_species_info3,
  by.x = "Genome5",
  by.y = "Genome",
  all.x = FALSE,
  all.y = FALSE
)

print(dim(filtered_SNVs_coverage5_top5_species))

filtered_SNVs_coverage5_top5_species$Genus <- ifelse(
  is.na(filtered_SNVs_coverage5_top5_species$Genus),
  "Unknown",
  filtered_SNVs_coverage5_top5_species$Genus
)

filtered_SNVs_coverage5_top5_species$Species <- ifelse(
  is.na(filtered_SNVs_coverage5_top5_species$Species),
  "Unknown",
  filtered_SNVs_coverage5_top5_species$Species
)

table(!is.na(filtered_SNVs_coverage5_top5_species$Species))
table(filtered_SNVs_coverage5_top5_species$Species == "Unknown")

# Check unknown species labels.
test <- subset(filtered_SNVs_coverage5_top5_species,
               filtered_SNVs_coverage5_top5_species$Species == "Unknown")
dim(test) # 23726 * 30
table(test$Genus)
# Bifidobacterium: 23090
# Veillonella: 636
table(test$Genome1, test$Genus)

# Most unknown species were annotated only to Bifidobacterium or Veillonella genus level.
filtered_SNVs_coverage5_top5_species$Species <- ifelse(
  filtered_SNVs_coverage5_top5_species$Genus == "Veillonella" &
    filtered_SNVs_coverage5_top5_species$Species == "Unknown",
  paste0(filtered_SNVs_coverage5_top5_species$Genus, " sp."),
  filtered_SNVs_coverage5_top5_species$Species
)

filtered_SNVs_coverage5_top5_species$Species <- ifelse(
  filtered_SNVs_coverage5_top5_species$Genus == "Bifidobacterium" &
    filtered_SNVs_coverage5_top5_species$Species == "Unknown",
  paste0(filtered_SNVs_coverage5_top5_species$Genus, " spp."),
  filtered_SNVs_coverage5_top5_species$Species
)

table(is.na(filtered_SNVs_coverage5_top5_species$Species))
length(table(filtered_SNVs_coverage5_top5_species$Taxonomy.lineage..GTDB.))
# 61 taxa

table(filtered_SNVs_coverage5_top5_species$mutation_type)
# I: 661064
# M: 5844
# N: 410147
# S: 672753
# none_gene: 146529
# I = 34.9%
# M = 0.31%
# N = 21.6%
# S = 35.5%
# none_gene = 7.7%


#### Convert to data.table and retain selected mutation types

setDT(filtered_SNVs_coverage5_top5_species)

filtered_SNVs_coverage5_top5_species$mutation_type <- ifelse(
  is.na(filtered_SNVs_coverage5_top5_species$mutation_type),
  "NG",
  filtered_SNVs_coverage5_top5_species$mutation_type
)

dt <- filtered_SNVs_coverage5_top5_species[
  !is.na(sample2) &
    !is.na(Species) &
    !is.na(mutation_type) &
    mutation_type %in% c("S", "N", "I", "M", "NG")
]

dim(filtered_SNVs_coverage5_top5_species)
dim(dt)

# Build species-mutation-type SNV features.
dt[, Species_clean := make.names(Species)]
dt[, feature := paste(Species_clean, mutation_type, sep = "_")]

length(unique(dt$Taxonomy.lineage..GTDB.))
length(unique(dt$feature))
# 61 taxa and 174 species-mutation-type SNV features.


#### Count sample-level SNV feature records

snv_feature <- dt[, .N, by = .(sample2, feature)]

snv_matrix_dt <- dcast(
  snv_feature,
  sample2 ~ feature,
  value.var = "N",
  fill = 0
)


#### Convert to count matrix for DMM clustering

count_matrix <- as.matrix(snv_matrix_dt[, -1, with = FALSE])
rownames(count_matrix) <- snv_matrix_dt$sample2

dim(count_matrix)
head(rownames(count_matrix))
head(colnames(count_matrix))


#### Run DMM clustering

set.seed(123)

fit_list_k20 <- lapply(1:20, function(k) {
  dmn(count_matrix, k = k, verbose = TRUE)
})

save(fit_list_k20, file = "fit_list_k20.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/fit_list_k20.RData")


#### Compare models and select optimal cluster number

laplace_values <- sapply(fit_list_k20, laplace)

laplace_df <- data.frame(
  k = 1:20,
  Laplace = laplace_values
)

print(laplace_df)

p_dmm_laplace_k20 <- ggplot(laplace_df, aes(x = k, y = Laplace)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 2, alpha = 0.5, color = "red") +
  theme_bw() +
  labs(x = "K value", y = "Laplace") +
  geom_vline(xintercept = 7, size = 0.8, linetype = "dashed", color = "grey50")

p_dmm_laplace_k20
# Save as 4 * 3
# k = 7 was selected as an elbow solution based on the Laplace curve and biological interpretability.


#### Extract best model and assign samples to clusters

best_fit <- fit_list_k20[[7]]

posterior_prob <- mixture(best_fit)
cluster_assign <- apply(posterior_prob, 1, which.max)

cluster_result <- data.table(
  sample2 = rownames(count_matrix),
  DMM_cluster = paste0("Cluster_", cluster_assign)
)

head(cluster_result)

write_xlsx(cluster_result,
           "D:/3_Projects/2_Children_jaundice/3_R_analysis/DMM_cluster_result.xlsx")


#### Calculate cluster-level mean SNV feature counts

cluster_matrix_dt <- as.data.table(count_matrix)
cluster_matrix_dt[, sample2 := rownames(count_matrix)]
cluster_matrix_dt <- merge(cluster_result, cluster_matrix_dt, by = "sample2")

cluster_feature_mean <- cluster_matrix_dt[
  ,
  lapply(.SD, mean),
  by = DMM_cluster,
  .SDcols = setdiff(colnames(cluster_matrix_dt), c("sample2", "DMM_cluster"))
]

write_xlsx(cluster_feature_mean,
           "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 4_file/DMM_cluster_feature_means.xlsx")

cluster_feature_mean #

# Z-score transformation across clusters for each feature.
mat <- as.matrix(cluster_feature_mean[, -1, with = FALSE])
rownames(mat) <- cluster_feature_mean$DMM_cluster
mat_scaled <- scale(mat)


#### Fig. 4b Distinct SNV feature profiles across DMM clusters ####

cluster_order <- paste0("Cluster_", 1:7)
mat_scaled <- mat_scaled[cluster_order, ]

p_DMM_cluster_headmap <- pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  color = colorRampPalette(c("#0C7C59", "#F5F5F5", "#FF7F11"))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),
  border_color = NA,
  main = "DMM clusters driven by SNV features"
)

p_DMM_cluster_headmap
# Save as 9 * 3.3


#### Fig. 4c Cluster-specific dominant SNV features ####

cluster_feature_long <- melt(
  cluster_feature_mean,
  id.vars = "DMM_cluster",
  variable.name = "feature",
  value.name = "mean_feature_count"
)

top_features_per_cluster <- cluster_feature_long[
  order(-mean_feature_count),
  head(.SD, 10),
  by = DMM_cluster
]

top_features_per_cluster

p_DMM_cluster_top_feature <- ggplot(
  top_features_per_cluster,
  aes(x = reorder_within(feature, mean_feature_count, DMM_cluster),
      y = mean_feature_count,
      fill = DMM_cluster)
) +
  geom_bar(stat = "identity") +
  facet_wrap(~DMM_cluster, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  theme_bw() +
  labs(x = "Feature", y = "Mean SNV feature count") +
  scale_fill_manual(values = c(
    "Cluster_1" = "#E64B35",
    "Cluster_2" = "#4DBBD5",
    "Cluster_3" = "#00A087",
    "Cluster_4" = "#3C5488",
    "Cluster_5" = "#E69F00",
    "Cluster_6" = "#0072B2",
    "Cluster_7" = "#CC79A7"
  ))

p_DMM_cluster_top_feature
# Save as 11 * 9


#### Fig. 4d Temporal dynamics of DMM clusters during infant gut microbiome development ####

write.table(
  cluster_result,
  file = "dmm_cluster_result.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

dmm_cluster_result <- read_tsv(
  "dmm_cluster_result.tsv",
  show_col_types = FALSE
)

cluster_meta <- merge(dmm_cluster_result, meta_final, by = "sample2", all.x = TRUE)
dim(cluster_meta) # 5755

table(is.na(cluster_meta$DOL)) # 75 samples without DOL information
cluster_meta$DOL <- as.numeric(cluster_meta$DOL)

cluster_meta <- cluster_meta %>%
  mutate(Time_new3 = case_when(
    DOL >= 0 & DOL <= 1 ~ "0",
    DOL > 1 & DOL <= 14 ~ "0.5",
    DOL > 14 & DOL <= 30 ~ "1",
    DOL > 30 & DOL <= 90 ~ "3",
    DOL > 90 & DOL <= 180 ~ "6",
    DOL > 180 & DOL <= 270 ~ "9",
    DOL > 270 & DOL <= 360 ~ "12",
    DOL > 360 & DOL <= 540 ~ "18",
    DOL > 540 & DOL <= 720 ~ "24",
    DOL > 720 & DOL <= 900 ~ "30",
    DOL > 900 ~ "36"
  ))

cluster_meta$Time_new3 <- as.numeric(cluster_meta$Time_new3)

table(cluster_meta$Time_new2, cluster_meta$Time_new3)
table(cluster_meta$Time_new3)
# 0: 305
# 0.5: 1404
# 1: 1018
# 3: 822
# 6: 420
# 9: 432
# 12: 499
# 18: 375
# 24: 292
# 30: 72
# 36: 41

plot_dt <- cluster_meta %>%
  filter(!is.na(Time_new3)) %>%
  group_by(Time_new3, DMM_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Time_new3) %>%
  mutate(prop = n / sum(n))

plot_dt$Time_new3 <- factor(
  plot_dt$Time_new3,
  levels = c(0, 0.5, 1, 3, 6, 9, 12, 18, 24, 30, 36)
)

p_DMM_cluster_overtime <- ggplot(plot_dt, aes(x = Time_new3, y = prop, fill = DMM_cluster, group = DMM_cluster)) +
  geom_area(position = "fill", alpha = 0.9) +
  theme_bw() +
  ylab("Proportion of samples") +
  xlab("Age (months)") +
  scale_fill_manual(values = c(
    "Cluster_1" = "#E64B35",
    "Cluster_2" = "#4DBBD5",
    "Cluster_3" = "#00A087",
    "Cluster_4" = "#3C5488",
    "Cluster_5" = "#E69F00",
    "Cluster_6" = "#0072B2",
    "Cluster_7" = "#CC79A7"
  )) +
  geom_vline(xintercept = "3", size = 0.8, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = "12", size = 0.8, linetype = "dashed", color = "grey")

p_DMM_cluster_overtime
# Save as 4.8 * 3


#### Supplementary Fig. S3a Normalized SNV rate across DMM-defined SNV burden clusters ####

dmm_cluster_result <- read_tsv(
  "dmm_cluster_result.tsv",
  show_col_types = FALSE
)

cluster_meta <- merge(dmm_cluster_result, meta_final, by = "sample2", all.x = TRUE)
dim(cluster_meta) # 5755

range(cluster_meta$SNV_rate_mean)
# 5.644601e-05 1.343776e-01
summary(cluster_meta$SNV_rate_mean)
# Min. 5.645e-05
# 1st Qu. 1.497e-03
# Median 2.196e-03
# Mean 2.830e-03
# 3rd Qu. 3.233e-03
# Max. 1.344e-01

range(cluster_meta$nucl_diversity_mean)
# 0.0003611403 0.0769938624
summary(cluster_meta$nucl_diversity_mean)
# Min. 0.0003611
# 1st Qu. 0.0014622
# Median 0.0020269
# Mean 0.0028283
# 3rd Qu. 0.0030535
# Max. 0.0769939

# Order clusters by median sample-level mean SNV rate.
cluster_meta$DMM_cluster <- fct_reorder(
  cluster_meta$DMM_cluster,
  cluster_meta$SNV_rate_mean,
  .fun = median,
  .desc = TRUE
)

my_comparisons <- list(
  c("Cluster_6", "Cluster_2"),
  c("Cluster_2", "Cluster_5"),
  c("Cluster_5", "Cluster_3"),
  c("Cluster_3", "Cluster_7"),
  c("Cluster_7", "Cluster_1"),
  c("Cluster_1", "Cluster_4")
)

my_comparisons_df <- data.frame(
  group1 = sapply(my_comparisons, `[`, 1),
  group2 = sapply(my_comparisons, `[`, 2)
) %>%
  mutate(pair_id = paste(pmin(group1, group2), pmax(group1, group2), sep = "_vs_"))

stat_test <- compare_means(
  SNV_rate_mean ~ DMM_cluster,
  data = cluster_meta,
  method = "wilcox.test"
) %>%
  mutate(pair_id = paste(pmin(group1, group2), pmax(group1, group2), sep = "_vs_")) %>%
  inner_join(my_comparisons_df, by = "pair_id", suffix = c("", ".target")) %>%
  transmute(
    group1 = group1.target,
    group2 = group2.target,
    p = p
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    y.position = 10^c(-1.75, -1.60, -1.45, -1.30, -1.15, -1.00)
  )

stat_test #

p_DMM_cluster_snvrate <- ggplot(cluster_meta, aes(x = DMM_cluster, y = SNV_rate_mean, fill = DMM_cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.6, color = "grey50") +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "DMM SNV burden cluster",
    y = "Mean SNV rate",
    title = "Normalized SNV rate across DMM clusters"
  ) +
  scale_fill_manual(values = c(
    "Cluster_1" = "#E64B35",
    "Cluster_2" = "#4DBBD5",
    "Cluster_3" = "#00A087",
    "Cluster_4" = "#3C5488",
    "Cluster_5" = "#E69F00",
    "Cluster_6" = "#0072B2",
    "Cluster_7" = "#CC79A7"
  )) +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    size = 3,
    tip.length = 0.01,
    hide.ns = FALSE
  )

p_DMM_cluster_snvrate
# Save as 5 * 4

kruskal.test(SNV_rate_mean ~ DMM_cluster, data = cluster_meta)
# Kruskal-Wallis chi-squared = 567.38, df = 6, P < 2.2e-16

lm_snv_cluster <- lm(
  log10(SNV_rate_mean + 1e-6) ~ DMM_cluster + Time_new2,
  data = cluster_meta
)

anova(lm_snv_cluster)
# DMM_cluster: F = 83.305, P < 2.2e-16
# Time_new2: F = 101.484, P < 2.2e-16
# DMM burden clusters remain associated with normalized SNV rate after adjustment for age group.




#### Supplementary Fig. S3b MDS ordination of normalized species-level SNV-rate profiles colored by DMM cluster ####

sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")
sample_species_SNVs_rate #
head(sample_species_SNVs_rate)

snv_matrix <- sample_species_SNVs_rate %>%
  dplyr::select(
    sample2,
    Species,
    total_SNV_rate
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = total_SNV_rate,
    values_fill = 0
  )

head(snv_matrix)[, 1:6]

mat <- snv_matrix %>%
  column_to_rownames("sample2") %>%
  as.matrix()

dim(mat)

# Log-transform SNV rate because of its highly skewed distribution.
summary(as.vector(mat))
mat_log <- log10(mat + 1e-6)
summary(as.vector(mat_log))

# Scale each species to avoid dominance by high-SNV-rate species.
mat_scaled <- scale(mat_log)
summary(as.vector(mat_scaled))

dist_mat <- dist(
  mat_scaled,
  method = "euclidean"
)

mds_res <- cmdscale(
  dist_mat,
  k = 2,
  eig = TRUE
)

mds_df <- data.frame(
  sample2 = rownames(mat_scaled),
  Axis1 = mds_res$points[, 1],
  Axis2 = mds_res$points[, 2]
)

head(mds_df)

cluster_meta_mds <- cluster_meta %>%
  left_join(mds_df, by = "sample2")

head(cluster_meta_mds)

p_mds <- ggplot(
  cluster_meta_mds,
  aes(x = Axis1,
      y = Axis2,
      color = DMM_cluster)
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  theme_bw() +
  labs(
    x = "MDS1",
    y = "MDS2",
    color = "DMM cluster",
    title = "Euclidean distance of species-level SNV profiles"
  ) +
  scale_color_manual(values = c(
    "Cluster_1" = "#E64B35",
    "Cluster_2" = "#4DBBD5",
    "Cluster_3" = "#00A087",
    "Cluster_4" = "#3C5488",
    "Cluster_5" = "#E69F00",
    "Cluster_6" = "#0072B2",
    "Cluster_7" = "#CC79A7"
  ))

print(p_mds)
# Save as 5 * 4




#### PERMANOVA: DMM cluster explains normalized species-level SNV-rate profile

sample_species_SNVs_rate #
cluster_meta #

# Build sample × species SNV-rate matrix.
snv_rate_matrix <- sample_species_SNVs_rate %>%
  dplyr::select(sample2, Species, total_SNV_rate) %>%
  tidyr::pivot_wider(
    names_from = Species,
    values_from = total_SNV_rate,
    values_fill = 0
  )

head(snv_rate_matrix[, 1:6])

mat <- snv_rate_matrix %>%
  tibble::column_to_rownames("sample2") %>%
  as.matrix()

dim(mat)

# Log-transform species-level SNV rate.
mat_log <- log10(mat + 1e-6)
summary(as.vector(mat_log))

# Remove species with zero variance before scaling.
species_sd <- apply(mat_log, 2, sd, na.rm = TRUE)
table(species_sd == 0)

mat_log_filter <- mat_log[, species_sd > 0, drop = FALSE]
dim(mat_log_filter)

mat_scaled <- scale(mat_log_filter)

valid_cols <- apply(mat_scaled, 2, function(x) all(is.finite(x)))
mat_scaled <- mat_scaled[, valid_cols, drop = FALSE]

dim(mat_scaled)
sum(is.na(mat_scaled))
sum(is.nan(mat_scaled))
sum(is.infinite(mat_scaled))

dist_snv <- dist(
  mat_scaled,
  method = "euclidean"
)

dist_snv_mat <- as.matrix(dist_snv)
dim(dist_snv_mat)

meta_for_adonis <- cluster_meta %>%
  dplyr::filter(sample2 %in% rownames(mat_scaled)) %>%
  dplyr::arrange(match(sample2, rownames(mat_scaled)))

all(meta_for_adonis$sample2 == rownames(mat_scaled))

meta_for_adonis2 <- meta_for_adonis %>%
  dplyr::select(
    sample2,
    DMM_cluster,
    Time_new2,
    SNV_rate_mean
  ) %>%
  dplyr::filter(
    !is.na(DMM_cluster),
    !is.na(Time_new2),
    !is.na(SNV_rate_mean)
  )

dim(meta_for_adonis2)

common_samples <- meta_for_adonis2$sample2

dist_snv_sub_mat <- dist_snv_mat[
  common_samples,
  common_samples
]

dist_snv_sub <- as.dist(dist_snv_sub_mat)

all(meta_for_adonis2$sample2 == labels(dist_snv_sub))

meta_for_adonis2 <- meta_for_adonis2 %>%
  mutate(
    DMM_cluster = as.factor(DMM_cluster),
    Time_new2 = as.numeric(Time_new2),
    log_SNV_rate_mean = log10(SNV_rate_mean + 1e-6)
  )

summary(meta_for_adonis2)

# Test whether DMM cluster explains species-level SNV-rate profiles after adjustment for age group.
adonis_res_age <- adonis2(
  dist_snv_sub ~ Time_new2 + DMM_cluster,
  data = meta_for_adonis2,
  permutations = 999,
  by = "margin"
)

adonis_res_age
# Time_new2: R2 = 0.00712, F = 43.782, P = 0.001
# DMM_cluster: R2 = 0.04094, F = 41.954, P = 0.001
# This analysis tests whether DMM cluster explains species-level SNV-rate profile variation after adjustment for age group.

save(
  adonis_res_age,
  dist_snv_sub,
  meta_for_adonis2,
  file = "PERMANOVA_Euclidean_SNV_rate_results.RData"
)

