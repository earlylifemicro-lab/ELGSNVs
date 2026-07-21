############### Gut microbial SNV ###################################

.libPaths()
.libPaths(c("D:/2_Software/R/install/R-4.5.1/library_shuang",
            "D:/2_Software/R/install/R-4.5.1/library"))
.libPaths()

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
library(tibble)
library(gratia)
library(mgcv)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(purrr)
library(forcats)
library(readr)
library(vegan)
library(cluster)
library(ggrepel)
library(nnet)
library(ggmosaic)


### Load data
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/Result6_R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Result6_R.RData")


############ Result: Prematurity is associated with lineage-specific microbial SNV remodeling and altered SNV-defined community states #####


### Shared colors ############

term_colors <- c(
  "full_term" = "#074cf5",
  "preterm" = "#f5b007"
)

cluster_colors <- c(
  "Cluster_1" = "#E64B35",
  "Cluster_2" = "#4DBBD5",
  "Cluster_3" = "#00A087",
  "Cluster_4" = "#3C5488",
  "Cluster_5" = "#E69F00",
  "Cluster_6" = "#0072B2",
  "Cluster_7" = "#CC79A7"
)

direction_colors <- c(
  "Not significant" = "grey70",
  "Higher in preterm" = "#f5b007",
  "Lower in preterm" = "#074cf5"
)


### Fig. 6a Distribution of sampling age by term status ############

meta_final <- read_excel("meta_final.xlsx")

table(meta_final$Term)

# Check sample structure across age, study, delivery mode, feeding pattern and continent.
# Term status can be confounded by cohort, age, delivery mode, feeding pattern and country.
meta_final %>% count(Term, Time_new2)
meta_final %>% count(Term, Study)
meta_final %>% count(Term, Delivery)
meta_final %>% count(Term, Feed)
meta_final %>% count(Term, Continent)

# Plot sampling-age distribution by term status
p_term_distribution_density2 <- meta_final %>%
  filter(!is.na(Term), !is.na(DOL)) %>%
  ggplot(aes(x = DOL + 1, fill = Term)) +
  geom_histogram(
    alpha = 0.6,
    position = "identity",
    bins = 40
  ) +
  scale_x_log10(
    breaks = c(1, 14, 30, 180, 540, 900, 1095),
    labels = c("0m", "0.5m", "1m", "6m", "12m", "24m", "36m")
  ) +
  scale_fill_manual(values = term_colors) +
  scale_color_manual(values = term_colors) +
  labs(
    x = "Chronologic age (months)",
    y = "Count",
    fill = "Term status",
    color = "Term status"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top")

p_term_distribution_density2
# save as 4 * 3





### Fig. 6b Overall microbial SNV rate in full-term and preterm infants ############

# Prepare data for overall and age-stratified comparisons
plot_term_df <- meta_final %>%
  filter(
    !is.na(Term),
    !is.na(Time_new2),
    !is.na(SNV_rate_mean),
    !is.na(pi_sample)
  ) %>%
  mutate(
    Term = factor(Term, levels = c("full_term", "preterm")),
    age_factor = factor(Time_new2, levels = c(0, 0.5, 1, 6, 12, 24, 36))
  )

p_snv_preterm_overall <- ggplot(
  plot_term_df,
  aes(x = Term, y = SNV_rate_mean, fill = Term)
) +
  geom_violin(alpha = 0.35, color = NA, trim = FALSE) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    alpha = 0.35,
    linewidth = 0.6
  ) +
  scale_y_log10() +
  stat_compare_means(
    method = "wilcox.test",
    size = 5
  ) +
  scale_fill_manual(values = term_colors) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Microbial SNV rate (log10)",
    title = "Overall difference in microbial SNV rate by term status"
  )

p_snv_preterm_overall
# save as 3.2 * 3





### Fig. 6c Overall microbial nucleotide diversity in full-term and preterm infants ############

p_pi_diversity_preterm_overall <- ggplot(
  plot_term_df,
  aes(x = Term, y = pi_sample, fill = Term)
) +
  geom_violin(alpha = 0.35, color = NA, trim = FALSE) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    alpha = 0.35,
    linewidth = 0.6
  ) +
  scale_y_log10() +
  stat_compare_means(
    method = "wilcox.test",
    size = 5
  ) +
  scale_fill_manual(values = term_colors) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_blank()
  ) +
  labs(
    x = "",
    y = "Nucleotide diversity (log10)",
    title = "Overall difference in nucleotide diversity by term status"
  )

p_pi_diversity_preterm_overall
# save as 3.2 * 3





### Supplementary Fig. S5a Age-stratified distributions of sample-level microbial SNV rate ############

plot_term_df <- plot_term_df %>%
  mutate(
    age_factor = factor(age_factor, levels = c("0", "0.5", "1", "6", "12", "24", "36")),
    Term = factor(Term, levels = c("full_term", "preterm"))
  )

# Wilcoxon tests comparing full-term and preterm infants within each age group
stat_snv_age_term <- ggpubr::compare_means(
  SNV_rate_mean ~ Term,
  data = plot_term_df,
  group.by = "age_factor",
  method = "wilcox.test"
) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p.adj.label = paste0("P = ", signif(p.adj, 3))
  )

# Set label positions
p_label_pos <- plot_term_df %>%
  group_by(age_factor) %>%
  summarise(
    y.position = max(SNV_rate_mean, na.rm = TRUE) * 1.8,
    .groups = "drop"
  )

stat_snv_age_term <- stat_snv_age_term %>%
  left_join(p_label_pos, by = "age_factor")

p_snv_age_term <- ggplot(
  plot_term_df,
  aes(x = age_factor, y = SNV_rate_mean, fill = Term)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = Term),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.75
    ),
    size = 0.45,
    alpha = 0.25
  ) +
  geom_text(
    data = stat_snv_age_term,
    aes(
      x = age_factor,
      y = y.position,
      label = p.adj.label
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_y_log10() +
  scale_color_manual(values = term_colors) +
  scale_fill_manual(values = term_colors) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(
    x = "Age (months)",
    y = "Microbial SNV rate (log10)",
    fill = "Term status",
    color = "Term status",
    title = "Prematurity-associated differences in microbial SNV rate across age"
  )

p_snv_age_term
# save as 8 * 4





### Supplementary Fig. S5b Age-stratified distributions of nucleotide diversity ############

# Wilcoxon tests comparing full-term and preterm infants within each age group
stat_pi_diversity_age_term <- ggpubr::compare_means(
  pi_sample ~ Term,
  data = plot_term_df,
  group.by = "age_factor",
  method = "wilcox.test"
) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p.adj.label = paste0("P = ", signif(p.adj, 3))
  )

# Set label positions
p_label_pos <- plot_term_df %>%
  group_by(age_factor) %>%
  summarise(
    y.position = max(pi_sample, na.rm = TRUE) * 1.8,
    .groups = "drop"
  )

stat_pi_diversity_age_term <- stat_pi_diversity_age_term %>%
  left_join(p_label_pos, by = "age_factor")

p_pi_diversity_age_term <- ggplot(
  plot_term_df,
  aes(x = age_factor, y = pi_sample, fill = Term)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = Term),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.75
    ),
    size = 0.45,
    alpha = 0.25
  ) +
  geom_text(
    data = stat_pi_diversity_age_term,
    aes(
      x = age_factor,
      y = y.position,
      label = p.adj.label
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_y_log10() +
  scale_color_manual(values = term_colors) +
  scale_fill_manual(values = term_colors) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(
    x = "Age (months)",
    y = "Weighted nucleotide diversity (log10)",
    fill = "Term status",
    color = "Term status",
    title = "Prematurity-associated differences in nucleotide diversity across age"
  )

p_pi_diversity_age_term
# save as 8 * 4





### Supplementary Fig. S6a+b Nonlinear age-dependent trajectories of microbial SNV rate ############

# Prepare data for GAM models
model_df <- meta_final %>%
  filter(
    !is.na(Term),
    !is.na(Time_new2),
    !is.na(SNV_rate_mean),
    !is.na(SubjectID),
    !is.na(Study)
  ) %>%
  mutate(
    log_SNV_rate = log10(SNV_rate_mean + 1e-6),
    log_pi = log10(pi_sample + 1e-6),
    Term = factor(Term),
    Time_new2 = as.numeric(Time_new2),
    Study = factor(Study),
    SubjectID = factor(SubjectID),
    Delivery = factor(Delivery, levels = c("vaginal", "CS"))
  )

dim(model_df)
colSums(is.na(model_df))
table(model_df$Time_new2)

# Feeding pattern is excluded because it is largely missing at 24 and 36 months,
# and its missingness is strongly age-dependent.

### Supplementary Fig. S6a: GAM without subject-level random effect

gam_term <- gam(
  log_SNV_rate ~ Term +
    s(Time_new2, by = Term, k = 5) +
    Delivery + Gender + Country +
    s(Study, bs = "re"),
  data = model_df,
  method = "REML"
)

gam.check(gam_term)
anova(gam_term)
summary(gam_term)

p_gam_term_woFeed_woSubjectID <- draw(
  gam_term,
  smooth_col = "#074cf5",
  ci_col = "#074cf5",
  ci_alpha = 0.2,
  residuals = FALSE
) &
  theme_bw(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_gam_term_woFeed_woSubjectID
# save as 6 * 4


### Supplementary Fig. S6b: GAM with subject-level random effect

gam_term2 <- gam(
  log_SNV_rate ~ Term +
    s(Time_new2, by = Term, k = 5) +
    Delivery + Gender + Country +
    s(Study, bs = "re") +
    s(SubjectID, bs = "re"),
  data = model_df,
  method = "REML"
)

anova(gam_term2)
summary(gam_term2)
# save(gam_term2, file = "gam_term2.RData")

p_gam_term_woFeed_wtSubjectID <- draw(
  gam_term2,
  smooth_col = "#074cf5",
  ci_col = "#074cf5",
  ci_alpha = 0.2,
  residuals = FALSE
) &
  theme_bw(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_gam_term_woFeed_wtSubjectID
# save as 6 * 4


### Difference smooths for microbial SNV rate

diff_term <- difference_smooths(
  gam_term,
  smooth = "s(Time_new2)",
  group_means = TRUE
)

draw(diff_term)

diff_df <- as.data.frame(diff_term) %>%
  mutate(
    sig_direction = case_when(
      .lower_ci > 0 ~ "full_term higher",
      .upper_ci < 0 ~ "preterm higher",
      TRUE ~ "not significant"
    )
  )

table(diff_df$sig_direction)

sig_age <- diff_df %>%
  filter(.lower_ci > 0 | .upper_ci < 0)

range(sig_age$Time_new2)





### Supplementary Fig. S6c+d Nonlinear age-dependent trajectories of weighted nucleotide diversity ############

### Supplementary Fig. S6c: GAM without subject-level random effect

gam_term_pi <- gam(
  log_pi ~ Term +
    s(Time_new2, by = Term, k = 5) +
    Delivery + Gender + Country +
    s(Study, bs = "re"),
  data = model_df,
  method = "REML"
)

gam.check(gam_term_pi)
anova(gam_term_pi)
summary(gam_term_pi)

p_gam_term_pi_woFeed_woSubjectID <- draw(
  gam_term_pi,
  smooth_col = "#074cf5",
  ci_col = "#074cf5",
  ci_alpha = 0.2,
  residuals = FALSE
) &
  theme_bw(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_gam_term_pi_woFeed_woSubjectID
# save as 6 * 4


### Supplementary Fig. S6d: GAM with subject-level random effect

gam_term2_pi <- gam(
  log_pi ~ Term +
    s(Time_new2, by = Term, k = 5) +
    Delivery + Gender + Country +
    s(Study, bs = "re") +
    s(SubjectID, bs = "re"),
  data = model_df,
  method = "REML"
)

anova(gam_term2_pi)
summary(gam_term2_pi)

save(gam_term2_pi, file = "gam_term2_pi.RData")

p_gam_term_pi_woFeed_wtSubjectID <- draw(
  gam_term2_pi,
  smooth_col = "#074cf5",
  ci_col = "#074cf5",
  ci_alpha = 0.2,
  residuals = FALSE
) &
  theme_bw(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_gam_term_pi_woFeed_wtSubjectID
# save as 6 * 4


### Difference smooths for weighted nucleotide diversity

diff_term_pi <- difference_smooths(
  gam_term2_pi,
  smooth = "s(Time_new2)",
  group_means = TRUE
)

draw(diff_term_pi)

diff_pi_df <- as.data.frame(diff_term_pi) %>%
  mutate(
    sig_direction = case_when(
      .lower_ci > 0 ~ "full_term higher",
      .upper_ci < 0 ~ "preterm higher",
      TRUE ~ "not significant"
    )
  )

table(diff_pi_df$sig_direction)

sig_age_pi <- diff_pi_df %>%
  filter(.lower_ci > 0 | .upper_ci < 0)

range(sig_age_pi$Time_new2)





### Term-associated overall species-level SNV-rate profile composition ############

# Build sample × species SNV-rate matrix
sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")

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

mat <- snv_matrix %>%
  column_to_rownames("sample2") %>%
  as.matrix()

mat_log <- log10(mat + 1e-6)

# Remove species with zero variance before scaling
species_sd <- apply(mat_log, 2, sd, na.rm = TRUE)
mat_log_filter <- mat_log[, species_sd > 0, drop = FALSE]

mat_scaled <- scale(mat_log_filter)

# Remove abnormal columns if present
valid_cols <- apply(mat_scaled, 2, function(x) all(is.finite(x)))
mat_scaled <- mat_scaled[, valid_cols, drop = FALSE]

# Calculate Euclidean distance
eucli_dist <- dist(
  mat_scaled,
  method = "euclidean"
)

eucli_mat <- as.matrix(eucli_dist)

write.csv(
  eucli_mat,
  "sample_species_SNVrate_Euclidean_distance_matrix.csv"
)

save(
  eucli_dist,
  eucli_mat,
  file = "sample_species_SNVrate_Euclidean_distance.RData"
)

load("sample_species_SNVrate_Euclidean_distance.RData")

# Prepare metadata for PERMANOVA
meta_beta <- meta_final %>%
  filter(sample2 %in% rownames(eucli_mat)) %>%
  dplyr::select(
    sample2,
    Term,
    Time_new2,
    Delivery,
    Feed,
    Gender,
    Study
  ) %>%
  na.omit()

common_samples <- meta_beta$sample2
eucli_mat_sub <- eucli_mat[common_samples, common_samples]
dist_snv <- as.dist(eucli_mat_sub)

all(meta_beta$sample2 == labels(dist_snv))

meta_beta <- meta_beta %>%
  mutate(
    Term = factor(Term),
    Delivery = factor(Delivery),
    Feed = factor(Feed),
    Gender = factor(Gender),
    Study = factor(Study),
    Time_new2 = as.numeric(Time_new2)
  )

adonis_term <- adonis2(
  dist_snv ~ Term + Time_new2 + Delivery + Feed + Gender + Study,
  data = meta_beta,
  method = "euclidean",
  permutations = 999,
  by = "margin"
)

adonis_term





### Fig. 6d Volcano plot showing species-level SNV rate differences associated with prematurity ############

meta_final <- read_excel("meta_final.xlsx")
sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")

# Merge species-level SNV-rate table with metadata
sample_species_SNVs_rate_metainfo <- merge(
  sample_species_SNVs_rate,
  meta_final,
  by = "sample2",
  all.x = TRUE
)

dim(sample_species_SNVs_rate_metainfo)

# Prepare species-level modelling data
species_term_df <- sample_species_SNVs_rate_metainfo %>%
  mutate(
    Time_new2 = as.numeric(as.character(Time_new2)),
    Term = factor(Term, levels = c("full_term", "preterm")),
    Delivery = factor(Delivery, levels = c("vaginal", "CS")),
    Gender = factor(Gender, levels = c("female", "male")),
    Study = factor(Study),
    SubjectID = factor(SubjectID),
    Country = factor(Country),
    Feed = factor(Feed),
    
    log_total_SNV_rate = log10(total_SNV_rate + 1e-6),
    log_S_SNV_rate = log10(S_SNV_rate + 1e-6),
    log_N_SNV_rate = log10(N_SNV_rate + 1e-6),
    log_species_coverage = log10(species_coverage + 1e-6),
    
    age_factor = factor(
      Time_new2,
      levels = c(0, 0.5, 1, 6, 12, 24, 36)
    )
  ) %>%
  filter(
    !is.na(Species),
    !is.na(sample2),
    !is.na(Term),
    !is.na(Time_new2),
    !is.na(total_SNV_rate),
    !is.na(species_coverage),
    !is.na(Study),
    !is.na(SubjectID)
  )

# Species-level filtering
species_term_df_filt0 <- species_term_df %>%
  filter(species_coverage >= 1)

species_summary <- species_term_df_filt0 %>%
  group_by(Species) %>%
  summarise(
    n_total = n_distinct(sample2),
    n_full_term = n_distinct(sample2[Term == "full_term"]),
    n_preterm = n_distinct(sample2[Term == "preterm"]),
    n_age = n_distinct(Time_new2),
    n_study = n_distinct(Study),
    median_coverage = median(species_coverage, na.rm = TRUE),
    mean_SNV_rate = mean(total_SNV_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_total))

species_keep <- species_summary %>%
  filter(
    n_total >= 20,
    n_full_term >= 10,
    n_preterm >= 10,
    n_age >= 1,
    n_study >= 1
  ) %>%
  pull(Species)

length(species_keep)

species_model_df <- species_term_df_filt0 %>%
  filter(Species %in% species_keep)


# Fit species-level mixed-effects model
fit_species_lmer <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(log_total_SNV_rate),
      !is.na(Term),
      !is.na(Time_new2),
      !is.na(log_species_coverage),
      !is.na(Delivery),
      !is.na(Gender),
      !is.na(Country),
      !is.na(Study),
      !is.na(SubjectID)
    ) %>%
    droplevels()
  
  if (
    nrow(df) < 80 ||
    n_distinct(df$Term) < 2 ||
    min(table(df$Term)) < 20 ||
    n_distinct(df$Time_new2) < 3 ||
    n_distinct(df$Study) < 2
  ) {
    return(tibble())
  }
  
  fit <- tryCatch(
    lmer(
      log_total_SNV_rate ~ Term + Time_new2 +
        Delivery + Gender + Country +
        (1 | Study) + (1 | SubjectID),
      data = df,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  res <- tryCatch(
    broom.mixed::tidy(fit, effects = "fixed"),
    error = function(e) tibble()
  )
  
  n_total <- nrow(df)
  n_full_term <- sum(df$Term == "full_term")
  n_preterm <- sum(df$Term == "preterm")
  n_age <- dplyr::n_distinct(df$Time_new2)
  n_study <- dplyr::n_distinct(df$Study)
  
  res %>%
    filter(term == "Termpreterm") %>%
    mutate(
      n = n_total,
      n_full_term = n_full_term,
      n_preterm = n_preterm,
      n_age = n_age,
      n_study = n_study
    )
}

species_term_results <- species_model_df %>%
  group_by(Species) %>%
  group_modify(~ fit_species_lmer(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    FC_preterm_vs_fullterm = 10^estimate,
    direction = case_when(
      p_adj < 0.05 & estimate > 0 ~ "Higher in preterm",
      p_adj < 0.05 & estimate < 0 ~ "Lower in preterm",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  species_term_results,
  "species_level_Term_differential_SNVrate_lmer.csv",
  row.names = FALSE
)

write_xlsx(
  species_term_results,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 6_file/species_level_Term_differential_SNVrate_lmer.xlsx"
)

# estimate > 0: higher SNV rate in preterm infants
# estimate < 0: lower SNV rate in preterm infants
# 10^estimate: preterm/full-term fold change

volcano_df <- species_term_results %>%
  mutate(
    neg_log10_FDR = -log10(p_adj),
    label = ifelse(p_adj < 0.05 & abs(estimate) > 0.05, Species, NA)
  )

p_volcano_species_term <- ggplot(
  volcano_df,
  aes(x = estimate, y = neg_log10_FDR)
) +
  geom_point(aes(color = direction), alpha = 0.8, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate(
    "text",
    x = -1.28,
    y = -log10(0.05),
    label = "FDR = 0.05",
    hjust = -0.1,
    vjust = -0.5,
    size = 3,
    color = "grey35"
  ) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  labs(
    x = "Effect size: preterm vs full-term (log10 SNV rate)",
    y = "-log10(FDR)",
    color = "",
    title = "Species-level SNV rate differences associated with prematurity"
  ) +
  theme(legend.position = "top")

p_volcano_species_term
# save as 7 * 5
# save as 3.9 * 3.2, p_volcano_species_term2





### Fig. 6f Coefficient plot of all significant species-level SNV rate differences associated with prematurity ############

top_species_effect <- species_term_results %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  mutate(Species = forcats::fct_reorder(Species, estimate))

p_species_allsig_coef <- ggplot(
  top_species_effect,
  aes(x = estimate, y = Species, color = direction)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(
    aes(
      xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error
    ),
    height = 0.2
  ) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  labs(
    x = "Adjusted effect of preterm birth on species-level SNV rate",
    y = "Species",
    color = "",
    title = "Top prematurity-associated species-level SNV rate differences"
  ) +
  theme(legend.position = "top")

p_species_allsig_coef
# save as 5 * 6
# save as 5.5 * 7





### Fig. 6e Volcano plot showing annotated gene-level SNV rate differences associated with prematurity ############

gene_sample_SNVrate_dNdS <- read_tsv(
  "gene_sample_SNVrate_dNdS.tsv",
  show_col_types = FALSE
)

table(!is.na(gene_sample_SNVrate_dNdS$mean_pi))
table(!is.na(gene_sample_SNVrate_dNdS$dNdS))
table(!is.na(gene_sample_SNVrate_dNdS$SNV_rate))

# Merge gene-level SNV-rate table with metadata
sample_gene_SNVs_rate_metainfo <- merge(
  gene_sample_SNVrate_dNdS,
  meta_final,
  by = "sample2",
  all.x = TRUE
)

dim(sample_gene_SNVs_rate_metainfo)

# Prepare gene-level modelling data
gene_term_df <- sample_gene_SNVs_rate_metainfo %>%
  mutate(
    Time_new2 = as.numeric(as.character(Time_new2)),
    Term = factor(Term, levels = c("full_term", "preterm")),
    Delivery = factor(Delivery, levels = c("vaginal", "CS")),
    Gender = factor(Gender, levels = c("female", "male")),
    Study = factor(Study),
    SubjectID = factor(SubjectID),
    Country = factor(Country),
    Feed = factor(Feed),
    
    log_SNV_rate = log10(SNV_rate + 1e-6),
    log_mean_pi = log10(mean_pi + 1e-6),
    log_dNdS = log10(dNdS + 1e-6),
    
    age_factor = factor(
      Time_new2,
      levels = c(0, 0.5, 1, 6, 12, 24, 36)
    )
  ) %>%
  filter(
    !is.na(gene_name),
    !is.na(sample2),
    !is.na(Term),
    !is.na(Time_new2),
    !is.na(SNV_rate),
    !is.na(mean_pi),
    !is.na(Study),
    !is.na(SubjectID)
  )

dim(gene_term_df)

# Gene-level filtering
gene_term_df_filt0 <- gene_term_df %>%
  filter(SNV_rate > 0)

gene_summary <- gene_term_df_filt0 %>%
  group_by(gene_name) %>%
  summarise(
    n_total = n_distinct(sample2),
    n_full_term = n_distinct(sample2[Term == "full_term"]),
    n_preterm = n_distinct(sample2[Term == "preterm"]),
    n_age = n_distinct(Time_new2),
    n_study = n_distinct(Study),
    median_SNV_rate = median(SNV_rate, na.rm = TRUE),
    mean_SNV_rate = mean(SNV_rate, na.rm = TRUE),
    median_pi = median(mean_pi, na.rm = TRUE),
    mean_dNdS = mean(dNdS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_total))

gene_keep <- gene_summary %>%
  filter(
    n_total >= 20,
    n_full_term >= 10,
    n_preterm >= 10,
    n_age >= 1,
    n_study >= 1
  ) %>%
  pull(gene_name)

length(gene_keep)

gene_model_df <- gene_term_df_filt0 %>%
  filter(gene_name %in% gene_keep)

dim(gene_model_df)


# Fit annotated gene-level mixed-effects model
fit_gene_lmer <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(log_SNV_rate),
      !is.na(Term),
      !is.na(Time_new2),
      !is.na(Delivery),
      !is.na(Gender),
      !is.na(Study),
      !is.na(SubjectID),
      !is.na(Country)
    ) %>%
    droplevels()
  
  if (
    nrow(df) < 80 ||
    dplyr::n_distinct(df$Term) < 2 ||
    min(table(df$Term)) < 20 ||
    dplyr::n_distinct(df$Time_new2) < 3 ||
    dplyr::n_distinct(df$Study) < 2
  ) {
    return(tibble())
  }
  
  n_total <- nrow(df)
  n_full_term <- sum(df$Term == "full_term")
  n_preterm <- sum(df$Term == "preterm")
  n_age <- dplyr::n_distinct(df$Time_new2)
  n_study <- dplyr::n_distinct(df$Study)
  
  fit <- tryCatch(
    lmer(
      log_SNV_rate ~ Term + Time_new2 +
        Delivery + Gender + Country +
        (1 | Study) + (1 | SubjectID),
      data = df,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  res <- tryCatch(
    broom.mixed::tidy(fit, effects = "fixed"),
    error = function(e) tibble()
  )
  
  res %>%
    filter(term == "Termpreterm") %>%
    mutate(
      n = n_total,
      n_full_term = n_full_term,
      n_preterm = n_preterm,
      n_age = n_age,
      n_study = n_study
    )
}

gene_term_results <- gene_model_df %>%
  group_by(gene_name) %>%
  group_modify(~ fit_gene_lmer(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    FC_preterm_vs_fullterm = 10^estimate,
    direction = case_when(
      p_adj < 0.05 & estimate > 0 ~ "Higher in preterm",
      p_adj < 0.05 & estimate < 0 ~ "Lower in preterm",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_adj)

write.csv(
  gene_term_results,
  "gene_level_Term_differential_SNVrate_lmer.csv",
  row.names = FALSE
)

write_xlsx(
  gene_term_results,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 6_file/gene_level_Term_differential_SNVrate_lmer.xlsx"
)

# estimate > 0: higher SNV rate in preterm infants
# estimate < 0: lower SNV rate in preterm infants
# 10^estimate: preterm/full-term fold change

volcano_gene_df <- gene_term_results %>%
  mutate(
    neg_log10_FDR = -log10(p_adj),
    label = ifelse(
      p_adj < 0.05 & abs(estimate) > 0.05,
      gene_name,
      NA
    )
  )

p_volcano_gene_term <- ggplot(
  volcano_gene_df,
  aes(x = estimate, y = neg_log10_FDR)
) +
  geom_point(
    aes(color = direction),
    alpha = 0.8,
    size = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = min(volcano_gene_df$estimate, na.rm = TRUE),
    y = -log10(0.05),
    label = "FDR = 0.05",
    hjust = -0.1,
    vjust = -0.5,
    size = 3,
    color = "grey35"
  ) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  labs(
    x = "Effect size: preterm vs full-term (log10 gene-level SNV rate)",
    y = "-log10(FDR)",
    color = "",
    title = "Gene-level SNV rate differences associated with prematurity"
  ) +
  theme(legend.position = "top")

p_volcano_gene_term
# save as 12 * 9
# save as 3.9 * 3.2, p_volcano_gene_term2





### Fig. 6g Coefficient plot of top 35 significant annotated gene-level SNV rate differences ############

top_gene_effect <- gene_term_results %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 35) %>%
  mutate(gene_name = as.character(gene_name))

# Shorten long gene annotation labels for plotting
top_gene_effect$gene_name <- ifelse(
  top_gene_effect$gene_name == "Essential cell division protein that forms a contractile ring structure (Z ring) at the future cell division site. The regulation of the ring assembly controls the timing and the location of cell division. One of the functions of the FtsZ ring is to recruit other cell division proteins to the septum to produce a new cell wall between the dividing cells. Binds GTP and shows GTPase activity",
  "Essential cell division protein",
  top_gene_effect$gene_name
)

top_gene_effect$gene_name <- ifelse(
  top_gene_effect$gene_name == "Forms oxaloacetate, a four-carbon dicarboxylic acid source for the tricarboxylic acid cycle",
  "A four-carbon dicarboxylic acid source for the tricarboxylic acid cycle",
  top_gene_effect$gene_name
)

top_gene_effect$gene_name <- ifelse(
  top_gene_effect$gene_name == "Participates in chromosomal partition during cell division. May act via the formation of a condensin-like complex containing Smc and ScpB that pull DNA away from mid-cell into both cell halves",
  "A condensin-like complex containing Smc and ScpB",
  top_gene_effect$gene_name
)

top_gene_effect$gene_name <- ifelse(
  top_gene_effect$gene_name == "Catalyzes the reversible interconversion of serine and glycine with tetrahydrofolate (THF) serving as the one-carbon carrier. This reaction serves as the major source of one-carbon groups required for the biosynthesis of purines, thymidylate, methionine, and other important biomolecules. Also exhibits THF- independent aldolase activity toward beta-hydroxyamino acids, producing glycine and aldehydes, via a retro-aldol mechanism",
  "Catalyzes the reversible interconversion of serine and glycine",
  top_gene_effect$gene_name
)

top_gene_effect <- top_gene_effect %>%
  arrange(p_adj) %>%
  mutate(gene_name = forcats::fct_reorder(gene_name, estimate))

p_gene_coef <- ggplot(
  top_gene_effect,
  aes(x = estimate, y = gene_name, color = direction)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(size = 2.5) +
  geom_errorbarh(
    aes(
      xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error
    ),
    height = 0.2
  ) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  labs(
    x = "Adjusted effect of preterm birth on gene-level SNV rate",
    y = "Gene",
    color = "",
    title = "Top prematurity-associated gene-level SNV rate differences"
  ) +
  theme(legend.position = "top")

p_gene_coef
# save as 8 * 7, p_gene_coef2





### Supplementary Fig. S7 Age-dependent distribution of SNV-defined DMM clusters by term status ############

# Test whether DMM cluster membership is associated with term status
table(meta_final$Term, meta_final$DMM_cluster)
chisq.test(table(meta_final$Term, meta_final$DMM_cluster))
# p-value < 2.2e-16

# Optional multinomial logistic regression
cluster_df <- meta_final %>%
  filter(!is.na(Term), !is.na(DMM_cluster), !is.na(Time_new2)) %>%
  mutate(
    DMM_cluster = factor(DMM_cluster),
    Term = factor(Term)
  )

m_cluster <- multinom(
  DMM_cluster ~ Term + Time_new2 + Delivery + Feed + Gender,
  data = cluster_df
)

summary(m_cluster)

# Plot age-dependent DMM cluster distribution by term status
cluster_term_df <- meta_final %>%
  filter(!is.na(Term), !is.na(DMM_cluster), !is.na(Time_new2)) %>%
  count(Time_new2, Term, DMM_cluster) %>%
  group_by(Time_new2, Term) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(age_factor = factor(Time_new2, levels = c(0, 0.5, 1, 6, 12, 24, 36)))

p_term_DMM_distribution <- ggplot(
  cluster_term_df,
  aes(x = age_factor, y = prop, fill = DMM_cluster)
) +
  geom_area(aes(group = DMM_cluster), position = "stack", alpha = 0.9) +
  facet_wrap(~ Term, ncol = 1) +
  theme_bw() +
  labs(
    x = "Age (months)",
    y = "Proportion of samples",
    fill = "DMM cluster",
    title = "Prematurity shifts SNV-defined microbial community states"
  ) +
  scale_fill_manual(values = cluster_colors)

p_term_DMM_distribution
# save as 7 * 6





### Prepare NEC-annotated metadata ############

meta_final <- read_excel("meta_final.xlsx")
NEC_meta_infor <- read_excel("NEC_meta_infor.xlsx")

# Check how many NEC-annotated samples are included in the global cohort
test <- NEC_meta_infor %>%
  subset(Study %in% c("BrooksB_2017", "RahmanSF_2018", "RavehSadkaT_2015", "RavehSadkaT_2016"))

dim(test)

table(NEC_meta_infor$sample2 %in% meta_final$sample2)
table(meta_final$sample2 %in% NEC_meta_infor$sample2)

# Add NEC information to metadata
meta_NEC_infor <- meta_final %>%
  left_join(
    NEC_meta_infor,
    by = c("sample2", "Study", "DOL", "Country", "SampleID")
  )

meta_NEC_infor$NEC_Shuang <- ifelse(
  meta_NEC_infor$PreNEC == "noNEC",
  "noNEC",
  "NEC"
)

# All NEC-annotated samples are preterm in this subset
fisher.test(table(meta_NEC_infor$PreNEC, meta_NEC_infor$Term))
fisher.test(table(meta_NEC_infor$NEC_Shuang, meta_NEC_infor$Term))

# Association tests between metadata variables and DMM clusters
chisq.test(table(meta_NEC_infor$Delivery, meta_NEC_infor$DMM_cluster))
chisq.test(table(meta_NEC_infor$Term, meta_NEC_infor$DMM_cluster))
fisher.test(table(meta_NEC_infor$NEC_Shuang, meta_NEC_infor$DMM_cluster))
# p-value = 0.02454





### Fig. 6h DMM cluster composition by term status ############

# Test association between term status and DMM cluster membership
table(meta_NEC_infor$Term)
table(meta_NEC_infor$DMM_cluster)
chisq.test(table(meta_NEC_infor$Term, meta_NEC_infor$DMM_cluster))
table(meta_NEC_infor$Term, meta_NEC_infor$DMM_cluster)

# Build contingency table
tab_term_DMM <- meta_NEC_infor %>%
  filter(!is.na(Term), !is.na(DMM_cluster)) %>%
  count(Term, DMM_cluster) %>%
  tidyr::pivot_wider(
    names_from = DMM_cluster,
    values_from = n,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Term") %>%
  as.matrix()

# Calculate within-term cluster proportions
prop_tab <- prop.table(tab_term_DMM, margin = 1)
round(prop_tab, 3)

# Plot mosaic
p_mosaic_term_DMM <- mosaicplot(
  tab_term_DMM,
  color = cluster_colors[colnames(tab_term_DMM)],
  border = "white",
  main = "Association between term status and DMM clusters",
  xlab = "Term status",
  ylab = "Proportion"
)

p_mosaic_term_DMM

# Add percentage labels manually
cum_prop <- t(apply(prop_tab, 1, function(x) cumsum(rev(x))))
mid_y <- 1 - (cum_prop - prop_tab[, ncol(prop_tab):1] / 2)
mid_y <- mid_y[, ncol(mid_y):1]

x_pos <- seq(
  1 / (2 * nrow(prop_tab)),
  1 - 1 / (2 * nrow(prop_tab)),
  length.out = nrow(prop_tab)
)

for (i in 1:nrow(prop_tab)) {
  for (j in 1:ncol(prop_tab)) {
    text(
      x = x_pos[i],
      y = mid_y[i, j],
      labels = paste0(round(rev(prop_tab[i, ])[j] * 100, 1), "%"),
      cex = 0.8,
      col = "white"
    )
  }
}

table(meta_NEC_infor$Term)
# full_term = 4015; preterm = 1651
# save as 4 * 4.3, p_mosaic_term_DMM





### Fig. 6i DMM cluster composition by NEC status in the NEC-annotated subset ############

# Test association between NEC status and DMM cluster membership
fisher.test(table(meta_NEC_infor$NEC_Shuang, meta_NEC_infor$DMM_cluster))
# p-value = 0.02454

table(meta_NEC_infor$NEC_Shuang, meta_NEC_infor$DMM_cluster)
table(meta_NEC_infor$NEC_Shuang, meta_NEC_infor$PreNEC)
fisher.test(table(meta_NEC_infor$PreNEC, meta_NEC_infor$DMM_cluster))

meta_NEC_infor$NEC_Shuang <- factor(
  meta_NEC_infor$NEC_Shuang,
  levels = c("noNEC", "NEC")
)

# Build contingency table
tab_NEC_DMM <- meta_NEC_infor %>%
  filter(
    !is.na(NEC_Shuang),
    !is.na(DMM_cluster)
  ) %>%
  count(NEC_Shuang, DMM_cluster) %>%
  tidyr::pivot_wider(
    names_from = DMM_cluster,
    values_from = n,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("NEC_Shuang") %>%
  as.matrix()

# Calculate within-NEC-group cluster proportions
prop_tab <- prop.table(tab_NEC_DMM, margin = 1)
round(prop_tab, 3)

# Plot mosaic
p_mosaic_NEC_DMM <- mosaicplot(
  tab_NEC_DMM,
  color = cluster_colors[colnames(tab_NEC_DMM)],
  border = "white",
  main = "Association between NEC status and DMM clusters",
  xlab = "NEC status",
  ylab = "Proportion"
)

p_mosaic_NEC_DMM

# Add percentage labels manually
cum_prop <- t(apply(prop_tab, 1, function(x) cumsum(rev(x))))
mid_y <- 1 - (cum_prop - prop_tab[, ncol(prop_tab):1] / 2)
mid_y <- mid_y[, ncol(mid_y):1]

x_pos <- seq(
  1 / (2 * nrow(prop_tab)),
  1 - 1 / (2 * nrow(prop_tab)),
  length.out = nrow(prop_tab)
)

for (i in 1:nrow(prop_tab)) {
  for (j in 1:ncol(prop_tab)) {
    text(
      x = x_pos[i],
      y = mid_y[i, j],
      labels = paste0(round(rev(prop_tab[i, ])[j] * 100, 1), "%"),
      cex = 0.8,
      col = "white"
    )
  }
}

prop_tab
table(meta_NEC_infor$NEC_Shuang, useNA = "ifany")
# NEC = 127; noNEC = 565; NA = 5063
# save as 4 * 4.3, p_mosaic_NEC_DMM

