############### Gut microbial SNV ###################################

.libPaths(c("D:/2_Software/R/install/R-4.5.1/library_shuang",
            "D:/2_Software/R/install/R-4.5.1/library"))

# Set working directory
setwd("D:/3_Projects/2_Children_jaundice/3_R_analysis")


### Necessary packages
library(data.table)
library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
library(ggpubr)
library(tidyr)
library(readr)
library(ggrepel)
library(forcats)
library(scales)


### Load R workspace
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/Result5_R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Result5_R.RData")


############ Result: SNV rate is decoupled from species relative abundance ############

#### Metadata information ####

meta_final <- read_excel("meta_final.xlsx")


#### Read scaffold-to-genome mapping file ####

stb <- fread(
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/reps.concat.stb",
  col.names = c("scaffold", "Genome"),
  quote = "",
  header = FALSE,
  encoding = "UTF-8",
  verbose = TRUE
)

# Harmonize genome identifiers for bin-level mapping.
stb[, Genome2 := gsub("_k[0-9]+$", "", Genome)]
stb[, Genome3 := sub("(_genomic).*", "_genomic", Genome2)]
stb[, Genome4 := sub("(_\\.\\d+.*)$", "", Genome3)]
stb[, Genome5 :=
      fifelse(
        grepl("_genomic", Genome4),
        sub("(_genomic).*", "\\1", Genome4),
        sub("_ERS.*$", "", Genome4)
      )
]


#### Read genome taxonomy information ####

genome_species_info <- read_excel(
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/genome_species_info.xlsx"
)

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
    Domain  = strip_prefix(Domain),
    Phylum  = strip_prefix(Phylum),
    Class   = strip_prefix(Class),
    Order   = strip_prefix(Order),
    Family  = strip_prefix(Family),
    Genus   = na_if(strip_prefix(Genus), ""),
    Species = na_if(strip_prefix(Species), "")
  )


#### Merge scaffold-to-genome mapping with species-level taxonomy ####

setDT(stb)
setDT(genome_species_info3)

genome_species_df <- merge(
  stb,
  genome_species_info3,
  by.x = "Genome5",
  by.y = "Genome",
  all.x = FALSE,
  all.y = FALSE
) %>%
  dplyr::select(-Genome2, -Genome3, -Genome4)

# Replace missing genus and species labels with family/genus-based unknown labels.
genome_species_df$Genus <- ifelse(
  is.na(genome_species_df$Genus),
  paste0(genome_species_df$Family, "_unknown"),
  genome_species_df$Genus
)

genome_species_df$Species <- ifelse(
  is.na(genome_species_df$Species),
  paste0(genome_species_df$Genus, "_unknown"),
  genome_species_df$Species
)

# Remove duplicated genome mappings.
genome_species_df_unique <- genome_species_df %>%
  distinct(Genome5, Genome, .keep_all = TRUE)


#### Read genome-level SNV summary with mutation-type information ####

SNVs_with_geneandgenome_summary <- read_tsv(
  "SNVs_coverage5_with_mutationtype_genome_summary2.tsv",
  show_col_types = FALSE
)

# Merge SNV summary with species-level taxonomy.
genome_species_SNVs_merged <- merge(
  SNVs_with_geneandgenome_summary,
  genome_species_df_unique,
  by.x = "genome_sp",
  by.y = "Genome",
  all.x = TRUE
)


#### Read inStrain genome-level coverage information ####

all_samples.genome_info <- read_tsv(
  "all_samples.genome_info.tsv",
  show_col_types = FALSE
)

genome_species_SNVs_merged$genome <- genome_species_SNVs_merged$genome_sp

genome_species_SNVs_merged2 <- genome_species_SNVs_merged %>%
  left_join(
    all_samples.genome_info %>%
      dplyr::select(genome, length, breadth_minCov, sample2, coverage),
    by = c("sample2", "genome")
  ) %>%
  as.data.frame()


#### Calculate sample-level species SNV rate and coverage-based relative abundance ####

sample_species_SNVs_rate <- genome_species_SNVs_merged2 %>%
  mutate(
    eff_len = length * breadth_minCov
  ) %>%
  group_by(sample2, Species) %>%
  summarise(
    S_SNV_count = sum(S_SNV_richness, na.rm = TRUE),
    N_SNV_count = sum(N_SNV_richness, na.rm = TRUE),
    I_SNV_count = sum(I_SNV_richness, na.rm = TRUE),
    M_SNV_count = sum(M_SNV_richness, na.rm = TRUE),
    total_SNV_count = sum(total_SNV_richness, na.rm = TRUE),
    total_species_coverage = sum(coverage * eff_len, na.rm = TRUE),
    total_eff_len = sum(eff_len, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    S_SNV_rate     = S_SNV_count / total_eff_len,
    N_SNV_rate     = N_SNV_count / total_eff_len,
    I_SNV_rate     = I_SNV_count / total_eff_len,
    M_SNV_rate     = M_SNV_count / total_eff_len,
    total_SNV_rate = total_SNV_count / total_eff_len,
    species_coverage = total_species_coverage / total_eff_len
  ) %>%
  group_by(sample2) %>%
  mutate(
    species_coverage_pct = species_coverage / sum(species_coverage)
  ) %>%
  ungroup()

write_xlsx(sample_species_SNVs_rate, "sample_species_SNVs_rate.xlsx")
sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")


#### Optional filtering: retain species detected in more than 100 samples ####

species_filtered <- sample_species_SNVs_rate %>%
  group_by(Species) %>%
  mutate(n_samples = n_distinct(sample2)) %>%
  ungroup() %>%
  filter(n_samples > 100)


##### Supplementary Fig. S2a and Fig. 3a: Species-level SNV burden, SNV rate and abundance #####

species_summary <- sample_species_SNVs_rate %>%
  group_by(Species) %>%
  summarise(
    total_SNV_burden = sum(total_SNV_count, na.rm = TRUE),
    total_eff_len = sum(total_eff_len, na.rm = TRUE),
    mean_abundance = mean(species_coverage_pct, na.rm = TRUE),
    prevalence = n_distinct(sample2),
    .groups = "drop"
  ) %>%
  mutate(
    total_SNV_rate = total_SNV_burden / total_eff_len,
    log_rate = log10(total_SNV_rate + 1e-8),
    log_burden = log10(total_SNV_burden + 1),
    log_abundance = log10(mean_abundance + 1e-6),
    rate_level = ifelse(log_rate > median(log_rate), "High_rate", "Low_rate"),
    burden_level = ifelse(log_burden > median(log_burden), "High_burden", "Low_burden"),
    category = paste(rate_level, burden_level, sep = "_")
  )


##### Supplementary Fig. S2a: Relationship between abundance and total SNV burden #####

p_SNVburden_abundance <- ggplot(
  species_summary,
  aes(x = log_abundance, y = log_burden)
) +
  geom_point(size = 2, alpha = 0.7, color = "#1B9E77") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D95F0E",
    fill = "#D95F0E",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label.., sep = "")),
    label.x.npc = "left",
    label.y.npc = 0.85,
    size = 4
  ) +
  labs(
    x = "Log (abundance)",
    y = "Log (total SNV burden)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_SNVburden_abundance
# Save as 3 * 3


##### Fig. 3a: Relationship between abundance and SNV rate #####

p_SNVrate_abundance <- ggplot(
  species_summary,
  aes(x = log_abundance, y = log_rate)
) +
  geom_point(size = 2, alpha = 0.6, color = "#1B9E77") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D95F0E",
    fill = "#D95F0E",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label.., sep = "")),
    label.x.npc = "left",
    label.y.npc = 0.85,
    size = 4
  ) +
  labs(
    x = "Log (abundance)",
    y = "Log (SNV rate)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_SNVrate_abundance
# Save as 3 * 3


##### Fig. 3b and Fig. 3c: Observed N/S SNV ratio #####

# The observed N/S SNV ratio was calculated as nonsynonymous SNV count
# divided by synonymous SNV count with pseudocount correction.
# This metric is a descriptive coding-change burden measure and should not be
# interpreted as formal dN/dS or pN/pS.

sample_species_NSratio <- genome_species_SNVs_merged2 %>%
  mutate(
    eff_len = length * breadth_minCov
  ) %>%
  group_by(sample2, Species) %>%
  summarise(
    N_count = sum(N_SNV_richness, na.rm = TRUE),
    S_count = sum(S_SNV_richness, na.rm = TRUE),
    total_eff_len = sum(eff_len, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    NS_ratio = (N_count + 1) / (S_count + 1),
    log_NS_ratio = log10(NS_ratio)
  )

species_full <- species_summary %>%
  left_join(
    sample_species_NSratio %>%
      group_by(Species) %>%
      summarise(
        mean_NS_ratio = mean(NS_ratio, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Species"
  )


##### Fig. 3b: Relationship between abundance and observed N/S SNV ratio #####

p_abundance_NSratio <- ggplot(
  species_full,
  aes(x = log_abundance, y = mean_NS_ratio)
) +
  geom_point(size = 2.5, alpha = 0.6, color = "#1B9E77") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D95F0E",
    fill = "#D95F0E",
    alpha = 0.2,
    linewidth = 1
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85,
    size = 4
  ) +
  labs(
    x = "Log (abundance)",
    y = "Observed N/S SNV ratio"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_abundance_NSratio
# Save as 3 * 3


##### Fig. 3c: Relationship between SNV rate and observed N/S SNV ratio #####

p_snvrate_NSratio <- ggplot(
  species_full,
  aes(x = log_rate, y = mean_NS_ratio)
) +
  geom_point(size = 2.5, alpha = 0.6, color = "#1B9E77") +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#D95F0E",
    fill = "#D95F0E",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85,
    size = 4
  ) +
  labs(
    x = "Log (SNV rate)",
    y = "Observed N/S SNV ratio"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_snvrate_NSratio
# Save as 3 * 3


##### Fig. 3d: Relationship between SNV rate and total SNV burden #####

highlight_species_high_NS <- species_full %>%
  arrange(desc(mean_NS_ratio)) %>%
  slice(1:10)

highlight_species_low_NS <- species_full %>%
  arrange(mean_NS_ratio) %>%
  slice(1:10)

target_species <- c(
  "Escherichia coli",
  "Prevotella copri",
  "Klebsiella pneumoniae",
  "Enterobacter cloacae",
  "Faecalibacterium prausnitzii",
  "Bifidobacterium infantis",
  "Bifidobacterium breve",
  "Bifidobacterium bifidum",
  "Bifidobacterium longum",
  "Bifidobacterium pseudocatenulatum",
  "Bacteroides fragilis",
  "Phocaeicola vulgatus",
  "Bacteroides caccae",
  "Phocaeicola dorei",
  "Roseburia intestinalis",
  "Eubacterium rectale",
  "Streptococcus salivarius",
  "Streptococcus thermophilus",
  "Enterococcus faecalis",
  "Bacteroides thetaiotaomicron",
  "Bacteroides ovatus",
  "Akkermansia muciniphila",
  "Veillonella parvula",
  "Veillonella dispar",
  "Clostridium perfringens",
  "Collinsella aerofaciens",
  "Eggerthella lenta",
  "Citrobacter freundii"
)

species_full <- species_full %>%
  mutate(
    NS_ratio_group = ifelse(mean_NS_ratio > 1, ">1", "<=1")
  )

p_snvrate_burden_NSratio <- ggplot(
  species_full,
  aes(x = log_rate, y = log_burden)
) +
  geom_point(aes(color = NS_ratio_group), size = 2, alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "#D95F0E",
    fill = "grey",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.9,
    size = 4
  ) +
  scale_color_manual(
    values = c(">1" = "#8dd3c7", "<=1" = "#762a83"),
    name = "Observed N/S ratio"
  ) +
  geom_text_repel(
    data = highlight_species_high_NS,
    aes(label = Species),
    size = 3,
    color = "red",
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "red",
    segment.alpha = 0.8,
    segment.size = 0.5
  ) +
  geom_text_repel(
    data = highlight_species_low_NS,
    aes(label = Species),
    size = 3,
    color = "blue",
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "blue",
    segment.alpha = 0.8,
    segment.size = 0.5
  ) +
  geom_text_repel(
    data = species_full %>% filter(Species %in% target_species),
    aes(label = Species),
    size = 3,
    color = "purple",
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "purple",
    segment.alpha = 0.8,
    segment.size = 0.5
  ) +
  labs(
    x = "Log (SNV rate)",
    y = "Log (total SNV burden)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    vjust = 2,
    label = "Color: observed N/S ratio group\nRed: highest N/S ratio\nBlue: lowest N/S ratio\nPurple: selected infant gut taxa",
    size = 3
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

p_snvrate_burden_NSratio
# Save as 6.5 * 3.5 or 7.5 * 3.5


##### Supplementary Fig. S2b: Relationship between SNV rate and species prevalence #####

species_full$prevalence_pct <- species_full$prevalence / 5886

p_snvrate_prevalence_pct <- ggplot(
  species_full,
  aes(x = log_rate, y = prevalence_pct)
) +
  geom_point(alpha = 0.35, size = 2, color = "#2C7FB8") +
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "#D95F0E",
    linewidth = 1
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Log (SNV rate)",
    y = "Prevalence (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  stat_regline_equation(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  stat_cor(
    aes(label = paste(..p.label..)),
    label.x.npc = "left",
    label.y.npc = 0.85,
    size = 4
  )

p_snvrate_prevalence_pct
# Save as 3 * 3.5


##### Fig. 3e: Density distribution of log-transformed SNV rates across prevalence groups #####

species_full <- species_full %>%
  mutate(
    prev_group = case_when(
      prevalence > 1000 ~ "core",
      prevalence > 50 ~ "common",
      TRUE ~ "rare"
    ),
    prev_group = factor(prev_group, levels = c("rare", "common", "core"))
  )

p_snvrate_prevalence_group_density <- ggplot(
  species_full,
  aes(x = log_rate, fill = prev_group, color = prev_group)
) +
  geom_density(alpha = 0.6, linewidth = 1) +
  scale_fill_manual(
    values = c("rare" = "#c6dbef",
               "common" = "#6BAED6",
               "core" = "#2171B5")
  ) +
  scale_color_manual(
    values = c("rare" = "#c6dbef",
               "common" = "#6BAED6",
               "core" = "#2171B5")
  ) +
  labs(
    x = "Log (SNV rate)",
    y = "Density",
    fill = "Prevalence group",
    color = "Prevalence group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_snvrate_prevalence_group_density
# Save as 3.4 * 3.5


##### Fig. 3f: Mean and dispersion of log SNV rates across prevalence bins #####

summary_df <- species_full %>%
  mutate(
    prevalence_bin = cut(
      prevalence_pct,
      breaks = seq(0.0000, 0.56066, length.out = 6),
      include.lowest = TRUE
    )
  ) %>%
  group_by(prevalence_bin) %>%
  summarise(
    n = n(),
    mean_rate = mean(log_rate, na.rm = TRUE),
    sd_rate = sd(log_rate, na.rm = TRUE),
    .groups = "drop"
  )

p_prevalencebin_mean_var_rate <- ggplot(
  summary_df,
  aes(x = prevalence_bin, y = mean_rate, group = 1)
) +
  geom_line(color = "#756BB1", linewidth = 1) +
  geom_point(size = 3, color = "#54278F") +
  geom_errorbar(
    aes(ymin = mean_rate - sd_rate, ymax = mean_rate + sd_rate),
    width = 0.15,
    color = "#54278F"
  ) +
  geom_text(aes(label = n), vjust = -0.8, size = 3.5) +
  labs(
    x = "Prevalence bin",
    y = "Mean log SNV rate ± SD"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_prevalencebin_mean_var_rate
# Save as 3 * 4


##### Pairwise statistical tests across prevalence bins #####

# Use the same prevalence bins as Fig. 3f.
bin_test_df <- species_full %>%
  mutate(
    prevalence_bin = cut(
      prevalence_pct,
      breaks = seq(0.0000, 0.56066, length.out = 6),
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(prevalence_bin)) %>%
  filter(is.finite(log_rate))

# Summary statistics for each prevalence bin.
bin_summary <- bin_test_df %>%
  group_by(prevalence_bin) %>%
  summarise(
    n = n(),
    mean_log_rate = mean(log_rate, na.rm = TRUE),
    sd_log_rate = sd(log_rate, na.rm = TRUE),
    median_log_rate = median(log_rate, na.rm = TRUE),
    iqr_log_rate = IQR(log_rate, na.rm = TRUE),
    .groups = "drop"
  )

print(bin_summary)

# Global and pairwise tests for differences in log SNV rate.
kruskal_rate <- kruskal.test(log_rate ~ prevalence_bin, data = bin_test_df)

pairwise_rate <- compare_means(
  log_rate ~ prevalence_bin,
  data = bin_test_df,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

print(kruskal_rate)
print(pairwise_rate)

# Global and pairwise tests for differences in SNV-rate dispersion.
# Dispersion is measured as the absolute deviation from the bin-specific median.
dispersion_df <- bin_test_df %>%
  group_by(prevalence_bin) %>%
  mutate(
    bin_median_log_rate = median(log_rate, na.rm = TRUE),
    abs_dev_median = abs(log_rate - bin_median_log_rate)
  ) %>%
  ungroup()

fligner_dispersion <- fligner.test(log_rate ~ prevalence_bin, data = bin_test_df)

pairwise_dispersion <- compare_means(
  abs_dev_median ~ prevalence_bin,
  data = dispersion_df,
  method = "wilcox.test",
  p.adjust.method = "fdr"
)

print(fligner_dispersion)
print(pairwise_dispersion)

# Save statistical results.
write_xlsx(
  list(
    bin_summary = bin_summary,
    global_rate_test = data.frame(
      method = "Kruskal-Wallis test",
      statistic = unname(kruskal_rate$statistic),
      df = unname(kruskal_rate$parameter),
      p_value = kruskal_rate$p.value
    ),
    pairwise_rate_tests = pairwise_rate,
    global_dispersion_test = data.frame(
      method = "Fligner-Killeen test",
      statistic = unname(fligner_dispersion$statistic),
      df = unname(fligner_dispersion$parameter),
      p_value = fligner_dispersion$p.value
    ),
    pairwise_dispersion_tests = pairwise_dispersion
  ),
  "Fig3f_prevalence_bin_pairwise_tests.xlsx"
)



