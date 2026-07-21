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
library(vegan)
library(cluster)
library(readr)
library(cowplot)
library(rpart.plot)

### load phyloseq data
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/R.RData")


############ Result: Environmental and host factors shape the infant gut microbial SNV landscape #####
#### metadata information ####
load("meta_final.xlsx")
meta_final #
colnames(meta_final)



### Fig. 5a Mantel test  ############
############### Add Euclidean distance information
# Saved file
# sample_species_SNVs_rate.xlsx

# Mantel test based on Euclidean distance of species-level SNV rate

# 1. Read the original SNV-rate file
sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")

head(sample_species_SNVs_rate)
colnames(sample_species_SNVs_rate)

# 2. Build sample × species SNV-rate matrix
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
dim(snv_matrix)

# 3. Convert to matrix
mat <- snv_matrix %>%
  column_to_rownames("sample2") %>%
  as.matrix()
dim(mat)

# 4. Log transform
mat_log <- log10(mat + 1e-6)
summary(as.vector(mat_log))

# 5. Remove species with zero variance
#    Otherwise, scale() may generate NA/NaN values
species_sd <- apply(mat_log, 2, sd, na.rm = TRUE)
table(species_sd == 0)
mat_log_filter <- mat_log[, species_sd > 0, drop = FALSE]
dim(mat_log_filter)

# 6. Standardize values
#    Each species is scaled to mean 0 and variance 1
mat_scaled <- scale(mat_log_filter)

# Check for NA / NaN / Inf values
sum(is.na(mat_scaled))
sum(is.nan(mat_scaled))
sum(is.infinite(mat_scaled))

# If abnormal values still exist, remove abnormal columns
valid_cols <- apply(mat_scaled, 2, function(x) all(is.finite(x)))

mat_scaled <- mat_scaled[, valid_cols, drop = FALSE]

dim(mat_scaled)

# 7. Calculate Euclidean distance
eucli_dist <- dist(
  mat_scaled,
  method = "euclidean"
)
eucli_dist

# Convert to matrix for downstream subsetting
eucli_mat <- as.matrix(eucli_dist)

dim(eucli_mat)
head(rownames(eucli_mat))
head(colnames(eucli_mat))

# Save Euclidean distance matrix
write.csv(
  eucli_mat,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/sample_species_SNVrate_Euclidean_distance_matrix.csv"
)

save(
  eucli_dist,
  eucli_mat,
  file = "sample_species_SNVrate_Euclidean_distance.RData"
)


# 8. Prepare metadata

#    Here, the metadata table is assumed to be named meta_final
dmm_cluster_result <- read_tsv(
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/dmm_cluster_result.tsv",
  show_col_types = FALSE
)

cluster_meta <- merge(dmm_cluster_result, meta_final, by = "sample2", all.x = TRUE)
dim(cluster_meta) # 5755
head(cluster_meta)
range(cluster_meta$SNV_rate_mean) # 5.644601e-05 1.343776e-01
summary(cluster_meta$SNV_rate_mean) # 

colnames(cluster_meta)
dim(cluster_meta)


# Confirm the shared sample IDs
common_samples <- intersect(
  rownames(eucli_mat),
  cluster_meta$sample2
)
length(common_samples)

# Subset the Euclidean distance matrix
eucli_mat_sub <- eucli_mat[common_samples, common_samples]
eucli_dist_sub <- as.dist(eucli_mat_sub)

# Align metadata
meta_sub <- cluster_meta %>%
  dplyr::select(
    sample2,
    Time_new2,
    Delivery,
    Gender,
    Term,
    Feed,
    Continent,
    Country,
    DMM_cluster
  ) %>%
  filter(sample2 %in% common_samples) %>%
  arrange(match(sample2, common_samples))

# Check whether metadata and distance labels are fully aligned
all(meta_sub$sample2 == labels(eucli_dist_sub))

# This should return TRUE
head(meta_sub)


# 9. Check missing values
table(is.na(meta_sub$Time_new2)) # 75
table(is.na(meta_sub$Delivery)) # 9
table(is.na(meta_sub$Gender)) # 418
table(is.na(meta_sub$Term)) # 89
table(is.na(meta_sub$Feed)) # FALSE = 2829 | TRUE = 2926
table(meta_sub$Time_new2, meta_sub$Feed) # Feed metadata are available only within the first year, not at 2 or 3 years
table(is.na(meta_sub$Continent)) # 0





### Perform Mantel tests for each environmental factor separately
str(meta_sub)

# Principle:
# Continuous variables -> dist(method = "euclidean")
# Categorical variables -> Gower distance
# Missing values are handled separately for each factor

## Time_new2, continuous variable
idx <- !is.na(meta_sub$Time_new2)

# Subset at the matrix level
eucli_sub_time <- eucli_mat_sub[idx, idx]

# Convert to dist object
eucli_sub_dist_time <- as.dist(eucli_sub_time)

# Check whether dimensions are consistent
length(labels(eucli_sub_dist_time))
length(meta_sub$Time_new2[idx])

# Mantel test
mantel_time2 <- mantel(
  eucli_sub_dist_time,
  dist(meta_sub$Time_new2[idx], method = "euclidean"),
  method = "spearman",
  permutations = 999
)
save(mantel_time2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_time2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_time2.RData")

# This step takes a long time to run.
# Started at 20260509/10:14 and finished at 20260509/18:26.
View(mantel_time2)



## Feed, categorical variable
str(meta_sub)
table(meta_sub$Delivery)
table(meta_sub$Gender)
table(meta_sub$Term)
table(meta_sub$Feed)
table(meta_sub$Continent)
idx <- !is.na(meta_sub$Feed)
class(meta_sub)

eucli_sub_dist_feed <- as.dist(eucli_mat_sub[idx, idx])

env_dist_feed <- daisy(
  data.frame(Feed = as.factor(meta_sub$Feed[idx])),
  metric = "gower"
)
env_dist_feed <- as.dist(env_dist_feed)
View(env_dist_feed)

mantel_feed2 <- mantel(
  eucli_sub_dist_feed,
  env_dist_feed,
  method = "spearman",
  permutations = 999
)

length(labels(eucli_sub_dist_feed))
length(meta_sub$Feed[idx])

save(mantel_feed2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_feed2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_feed2.RData")





## Term
idx <- !is.na(meta_sub$Term)

eucli_sub_dist_term <- as.dist(eucli_mat_sub[idx, idx])

env_dist_term <- daisy(
  data.frame(Term = as.factor(meta_sub$Term[idx])),
  metric = "gower"
)
env_dist_term <- as.dist(env_dist_term)
View(env_dist_term)

length(labels(eucli_sub_dist_term))
length(meta_sub$Time_new2[idx])

mantel_term2 <- mantel(
  eucli_sub_dist_term,
  env_dist_term,
  method = "spearman",
  permutations = 999
)
save(mantel_term2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_term2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_term2.RData")


## Gender
idx <- !is.na(meta_sub$Gender)

eucli_sub_dist_gender <- as.dist(eucli_mat_sub[idx, idx])

env_dist_gender <- daisy(
  data.frame(Gender = as.factor(meta_sub$Gender[idx])),
  metric = "gower"
)
env_dist_gender <- as.dist(env_dist_gender)

length(labels(eucli_sub_dist_gender))
length(meta_sub$Time_new2[idx])

mantel_gender2 <- mantel(
  eucli_sub_dist_gender,
  env_dist_gender,
  method = "spearman",
  permutations = 999
)
save(mantel_gender2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_gender2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_gender2.RData")


## Delivery
idx <- !is.na(meta_sub$Delivery)

eucli_sub_dist_delivery <- as.dist(eucli_mat_sub[idx, idx])

env_dist_delivery <- daisy(
  data.frame(Delivery = as.factor(meta_sub$Delivery[idx])),
  metric = "gower"
)
env_dist_delivery <- as.dist(env_dist_delivery)

length(labels(eucli_sub_dist_delivery))
length(meta_sub$Time_new2[idx])

mantel_delivery2 <- mantel(
  eucli_sub_dist_delivery,
  env_dist_delivery,
  method = "spearman",
  permutations = 999
)
View(mantel_delivery2)
save(mantel_delivery2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_delivery2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_delivery2.RData")


### Continent
idx <- !is.na(meta_sub$Continent)

eucli_sub_dist_continent <- as.dist(eucli_mat_sub[idx, idx])

env_dist_continent <- daisy(
  data.frame(Continent = as.factor(meta_sub$Continent[idx])),
  metric = "gower"
)
env_dist_continent <- as.dist(env_dist_continent)

length(labels(eucli_sub_dist_continent))
length(meta_sub$Time_new2[idx])

mantel_continent2 <- mantel(
  eucli_sub_dist_continent,
  env_dist_continent,
  method = "spearman",
  permutations = 999
)
View(mantel_continent2)
save(mantel_continent2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_continent2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_continent2.RData")


### Country
idx <- !is.na(meta_sub$Country)

eucli_sub_dist_country <- as.dist(eucli_mat_sub[idx, idx])

env_dist_country <- daisy(
  data.frame(Country = as.factor(meta_sub$Country[idx])),
  metric = "gower"
)
env_dist_country <- as.dist(env_dist_country)

length(labels(eucli_sub_dist_country))
length(meta_sub$Time_new2[idx])

mantel_country2 <- mantel(
  eucli_sub_dist_country,
  env_dist_country,
  method = "spearman",
  permutations = 999
)
View(mantel_country2)
save(mantel_country2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_country2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_country2.RData")


### DMM_cluster
idx <- !is.na(meta_sub$DMM_cluster)

eucli_sub_dist_DMMcluster <- as.dist(eucli_mat_sub[idx, idx])

env_dist_DMMcluster <- daisy(
  data.frame(DMM_cluster = as.factor(meta_sub$DMM_cluster[idx])),
  metric = "gower"
)
env_dist_DMMcluster <- as.dist(env_dist_DMMcluster)

length(labels(eucli_sub_dist_DMMcluster))
length(meta_sub$Time_new2[idx])

mantel_DMMcluster2 <- mantel(
  eucli_sub_dist_DMMcluster,
  env_dist_DMMcluster,
  method = "spearman",
  permutations = 999
)
View(mantel_DMMcluster2)
save(mantel_DMMcluster2, file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_DMMcluster2.RData")
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/mantel_DMMcluster2.RData")








## Summarize Mantel results into a table
mantel_results2 <- data.frame(
  Factor = c("Time_new2", "Country", #"Continent", 
             "Feed", "Term", "Gender", "Delivery", "DMM_cluster"),
  r = c(
    mantel_time2$statistic,
    mantel_country2$statistic,
    # mantel_continent2$statistic,
    mantel_feed2$statistic,
    mantel_term2$statistic,
    mantel_gender2$statistic,
    mantel_delivery2$statistic,
    mantel_DMMcluster2$statistic
  ),
  p = c(
    mantel_time2$signif,
    mantel_country2$signif,
    # mantel_continent2$signif,
    mantel_feed2$signif,
    mantel_term2$signif,
    mantel_gender2$signif,
    mantel_delivery2$signif,
    mantel_DMMcluster2$signif
  )
)

# FDR correction
mantel_results2$p_adj <- p.adjust(mantel_results2$p, method = "BH")

mantel_results2

write.csv(mantel_results2, "Mantel_results_SNVrate_Eucli.csv", row.names = FALSE)


# Gower distance was used for categorical environmental variables, 
# whereas Euclidean distance was applied to continuous variables.


### Visualization
library(ggplot2)

mantel_results2$Significance <- ifelse(mantel_results2$p_adj < 0.05, "*", "")

p_mantel_envir2 <- ggplot(mantel_results2, aes(x = reorder(Factor, r), y = r)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Significance), vjust = -0.5, size = 6, color = "red") +
  labs(
    title = "Mantel Test: Environmental Factors vs. Microbial SNV Structure",
    x = "Environmental Factor",
    y = "Mantel r"
  ) +
  theme_minimal(base_size = 11) +
  coord_flip()
p_mantel_envir2
# save as 3.2 * 4






#### Supplementary Fig S4a. Univariate associations between environmental factors and microbial SNV features ######

### Associations between SNV rate, nucleotide diversity, Euclidean-distance MDS axes and environmental factors  
meta_df #
meta_final #

# dmm_cluster_result.tsv
dmm_cluster_result <- read_tsv(
  "dmm_cluster_result.tsv",
  show_col_types = FALSE
)

# Merge DMM cluster results with metadata
cluster_meta <- merge(dmm_cluster_result, meta_final, by = "sample2", all.x = TRUE)
dim(cluster_meta) # 5755
head(cluster_meta)
range(cluster_meta$SNV_rate_mean) # 5.644601e-05 1.343776e-01
summary(cluster_meta$SNV_rate_mean) # 
# Min.        1st Qu.     Median    Mean     3rd Qu.    Max. 
# 5.645e-05 1.497e-03 2.196e-03 2.830e-03 3.233e-03 1.344e-01 

range(cluster_meta$nucl_diversity_mean) # 0.0003611403 0.0769938624
summary(cluster_meta$nucl_diversity_mean) # 
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003611 0.0014622 0.0020269 0.0028283 0.0030535 0.0769939 



############### Add Euclidean distance information
# Saved files
# sample_species_SNVs_rate.xlsx
eucli_dist #
eucli_mat #
# sample_species_SNVrate_Euclidean_distance.RData



# Read file
sample_species_SNVs_rate <- read_excel("sample_species_SNVs_rate.xlsx")
sample_species_SNVs_rate #
head(sample_species_SNVs_rate)
View(sample_species_SNVs_rate)

load("D:/3_Projects/2_Children_jaundice/3_R_analysis/sample_species_SNVrate_Euclidean_distance.RData")

dist_mat <- eucli_dist


## 6. Classical MDS / PCoA
mds_res <- cmdscale(
  dist_mat,
  k = 2,
  eig = TRUE
)

## 7. Extract the first two axes
mds_df <- data.frame(
  sample2 = rownames(mat_scaled),
  Axis1 = mds_res$points[, 1],
  Axis2 = mds_res$points[, 2]
)
head(mds_df)

## 8. Merge with metadata
cluster_meta_mds <- cluster_meta %>%
  left_join(mds_df, by = "sample2")
head(cluster_meta_mds)














colnames(cluster_meta_mds)
table(cluster_meta_mds$DMM_cluster)

########______with direction
{
  ##### Single-factor linear model analysis
  cluster_meta_mds$Time_new2 <- factor(cluster_meta_mds$Time_new2, levels = c(0, 0.5, 1, 6, 12, 24, 36))
  # cluster_meta_mds$Time_new2 <- factor(cluster_meta_mds$Time_new2, levels = c("0", "0.5", "1", "6", "12", "24", "36"))
  cluster_meta_mds$Delivery <- factor(cluster_meta_mds$Delivery, levels = c("vaginal", "CS"))
  cluster_meta_mds$Gender <- factor(cluster_meta_mds$Gender, levels = c("male", "female"))
  cluster_meta_mds$Term <- factor(cluster_meta_mds$Term, levels = c("full_term", "preterm"))
  cluster_meta_mds$Feed <- factor(cluster_meta_mds$Feed, levels = c("breast", "combined", "formula"))
  
  # Set parameters
  group_vars <- c("Delivery", "Gender", "Term", "Feed", "Country" #, "Continent"
  )
  features <- c("SNV_rate_mean", "nucl_diversity_mean", # "shannon_entropy_median",
                "Axis1", "Axis2")
  # timepoints <- sort(unique(cluster_meta_mds$Time_new2))
  timepoints <- levels(cluster_meta_mds$Time_new2)
  
  # Convert the workflow into a function while retaining debug output
  run_lm_debug <- function(data, timepoint, group_var, feature) {
    message(">>> Time: ", timepoint, ", Group: ", group_var, ", Feature: ", feature)
    
    df <- data %>%
      filter(Time_new2 == timepoint) %>%
      select(all_of(group_var), all_of(feature)) %>%
      filter(!is.na(.data[[group_var]]), !is.na(.data[[feature]]))
    
    # Check whether the model can be fitted
    if (length(unique(df[[group_var]])) < 2) {
      message("  Only one group present, skipping.")
      return(NULL)
    }
    
    # Check whether the feature is numeric
    # If only one group is present or missing values prevent modelling, the function returns NULL.
    # Such cases will not be recorded in results_lm_debug_df.
    if (!is.numeric(df[[feature]])) {
      message("  Feature is not numeric, skipping.")
      return(NULL)
    }
    
    # Fit model
    model <- lm(as.formula(paste(feature, "~", group_var)), data = df) ### single-factor analysis
    anova_p <- tryCatch(anova(model)[1, "Pr(>F)"], error = function(e) NA)
    # r2 <- summary(model)$r.squared
    r2 <- summary(model)$adj.r.squared
    
    # R2 represents the proportion of variation in the response feature
    # explained by the grouping variable in this single-factor model.
    beta <- if (length(coef(model)) >= 2) coef(model)[2] else NA
    
    # Beta represents the difference in the feature value between the target group
    # and the reference group, including its direction.
    
    return(data.frame(
      Timepoint = timepoint,
      GroupVar = group_var,
      Feature = feature,
      p.value = anova_p,
      R2 = r2,
      Beta = beta,
      n = nrow(df)
    ))
  }
  
  ### Run one by one for debugging
  ## Feed
  # run_lm_debug(cluster_meta_mds, timepoint = 0, group_var = "Feed", feature = "BC_PCoA1")
  # run_lm_debug(cluster_meta_mds, timepoint = 0.5, group_var = "Feed", feature = "BC_PCoA1")
  
  # Execute all combinations and retain logs for debugging
  results_list <- list()
  row_counter <- 1
  
  for (tp in timepoints) {
    for (g in group_vars) {
      for (f in features) {
        res <- run_lm_debug(cluster_meta_mds, tp, g, f)
        if (!is.null(res)) {
          results_list[[row_counter]] <- res
          row_counter <- row_counter + 1
        }
      }
    }
  }
  
  results_lm_debug_df <- bind_rows(results_list)
  
  results_lm_debug_df %>%
    filter(is.na(p.value) | is.na(R2) | is.na(Beta)) %>%
    arrange(GroupVar, Feature)
  # View(results_lm_debug_df)
  
  # Generate all possible combinations
  all_combos <- expand.grid(
    Timepoint = unique(results_lm_debug_df$Timepoint),
    GroupVar = unique(results_lm_debug_df$GroupVar),
    Feature = unique(results_lm_debug_df$Feature),
    stringsAsFactors = FALSE
  )
  
  # Left join with original results and fill missing combinations with NA values
  results_lm_debug_df2 <- all_combos %>%
    left_join(results_lm_debug_df, by = c("Timepoint", "GroupVar", "Feature")) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("0", "0.5", "1", "6", "12", "24", "36")),
      
      sig_label = case_when(
        is.na(p.value) ~ "-",  # Missing values are shown as dashes
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        # p.value < 0.1 ~ "~",
        TRUE ~ ""
      )
    )
  
  # === Step 1: Add directional information ===
  results_lm_debug_df2 <- results_lm_debug_df2 %>%
    mutate(
      # Directional signs for binary variables
      signs = case_when(
        GroupVar %in% c("Delivery", "Gender", "Term") ~ sign(Beta),
        TRUE ~ NA_real_
      ),
      R2_signed = signs * R2  # Multiply R2 by direction
    )
  
  # === Step 2: Split data ===
  df_with_sign <- results_lm_debug_df2 %>% filter(!is.na(signs))
  df_no_sign   <- results_lm_debug_df2 %>% filter(is.na(signs))
  
  # Plot heatmap of features across time points and group variables
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  
  colnames(df_with_sign)
  table(df_with_sign$GroupVar)
  df_with_sign$GroupVar <- factor(df_with_sign$GroupVar,
                                  levels = c("Term", "Delivery", "Gender"))
  table(df_with_sign$Feature)
  df_with_sign$Feature <- factor(df_with_sign$Feature,
                                 levels = c("Axis2", "Axis1", "nucl_diversity_mean", "SNV_rate_mean"))
  
  
  df_no_sign$GroupVar <- factor(df_no_sign$GroupVar,
                                levels = c("Country", "Feed"))
  table(df_no_sign$Feature)
  df_no_sign$Feature <- factor(df_no_sign$Feature,
                               levels = c("Axis2", "Axis1", "nucl_diversity_mean", "SNV_rate_mean"))
  
  
  results_lm_debug_df2$GroupVar <- factor(results_lm_debug_df2$GroupVar,
                                          levels = c("Term", "Delivery", "Gender", "Country", "Feed"))
  table(results_lm_debug_df2$Feature)
  results_lm_debug_df2$Feature <- factor(results_lm_debug_df2$Feature,
                                         levels = c("Axis2", "Axis1", "nucl_diversity_mean", "SNV_rate_mean"))
  
  p_lm_single_envir <- ggplot() +
    # Layer with direction
    geom_tile(
      data = df_with_sign,
      aes(x = as.factor(Timepoint), y = Feature, fill = R2_signed),
      color = "white"
    ) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Signed Adjusted R²\n& Effect direction",
      na.value = "grey90",
      guide = guide_colorbar(order = 1)
    ) +
    
    # Reset fill scale
    ggnewscale::new_scale_fill() +
    
    # Layer without direction
    geom_tile(
      data = df_no_sign,
      aes(x = as.factor(Timepoint), y = Feature, fill = R2),
      color = "white"
    ) +
    scale_fill_gradient(
      low = "white", high = "darkgrey",
      name = "Adjusted R²\n(no direction)",
      na.value = "grey90",
      guide = guide_colorbar(order = 2)
    ) +
    
    # Significance labels
    geom_text(
      data = results_lm_debug_df2,
      aes(x = as.factor(Timepoint), y = Feature, label = sig_label),
      size = 4
    ) +
    facet_wrap(~ GroupVar, ncol = 1) +
    labs(
      x = "Infant age (months)",
      y = "Microbial Feature",
      title = "1.Group Differences in Microbial Features by Timepoint"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      plot.caption = element_text(size = 10, hjust = 0)
    )
  p_lm_single_envir
  
  
  # Add right-side annotation text
  caption_text <- "Significance codes:\n*** p < 0.001\n** p < 0.01\n* p < 0.05\n- not available"
  
  # Combine main plot and annotation text
  library("cowplot")
  p_lm_single_envir2 <- ggdraw() +
    draw_plot(p_lm_single_envir, x = 0, width = 0.85) +  # Main plot uses 85% of the width
    draw_label(caption_text, x = 0.655, y = 0.92, hjust = 0, size = 10)
  p_lm_single_envir2
  # save as 8.5 * 6
}

supp_s4a_pvalues <- results_lm_debug_df2 %>%
  dplyr::mutate(
    Figure = "Supplementary Fig. S4a",
    Model = "Univariate linear model",
    P_label = dplyr::case_when(
      is.na(p.value) ~ "not available",
      p.value < 0.001 ~ "P < 0.001",
      p.value < 0.01 ~ "P < 0.01",
      p.value < 0.05 ~ "P < 0.05",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    Figure,
    Model,
    Timepoint,
    GroupVar,
    Feature,
    p.value,
    P_label,
    R2,
    Beta,
    sig_label
  ) %>%
  dplyr::arrange(Timepoint, GroupVar, Feature)
supp_s4a_pvalues

write_xlsx(
  supp_s4a_pvalues,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/Supplementary_Fig_S4a_univariate_pvalues.xlsx"
)





### Fig. 5b Multivariable associations between environmental factors and microbial SNV features  ############

{
  ##### Multivariable linear model analysis
  
  ### Single time point + feature
  # Set variables to check
  ### 0 day
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Delivery)
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Gender)
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Term)
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Feed) # 24, 36
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Country) # 36
  
  cluster_meta_mds$Delivery <- factor(cluster_meta_mds$Delivery, levels = c("vaginal", "CS"))
  cluster_meta_mds$Gender <- factor(cluster_meta_mds$Gender, levels = c("male", "female"))
  cluster_meta_mds$Term <- factor(cluster_meta_mds$Term, levels = c("full_term", "preterm"))
  cluster_meta_mds$Feed <- factor(cluster_meta_mds$Feed, levels = c("breast", "combined", "formula"))
  
  my_time <- 0
  my_feature <- "SNV_rate_mean"
  my_group_vars <- c("Delivery", "Gender", "Term", "Feed", "Country")
  
  # Run one multivariable linear model
  ### Time = 0
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 0)
  lm_model_snvrate_0 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_snvrate_0
  summary(lm_model_snvrate_0)
  anova(lm_model_snvrate_0)
  
  lm_model_nucl_diversity_0 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_nucl_diversity_0
  summary(lm_model_nucl_diversity_0)
  anova(lm_model_nucl_diversity_0)
  
  lm_model_Eucli_MDS1_0 <- lm(Axis1 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS1_0
  summary(lm_model_Eucli_MDS1_0)
  anova(lm_model_Eucli_MDS1_0)
  
  lm_model_Eucli_MDS2_0 <- lm(Axis2 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS2_0
  summary(lm_model_Eucli_MDS2_0)
  anova(lm_model_Eucli_MDS2_0)
  
  library(car)
  # vif(lm_model_snvrate_0)
  
  
  ### Time = 0.5
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 0.5)
  lm_model_snvrate_0.5 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_snvrate_0.5
  summary(lm_model_snvrate_0.5)
  anova(lm_model_snvrate_0.5)
  
  lm_model_nucl_diversity_0.5 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_nucl_diversity_0.5
  summary(lm_model_nucl_diversity_0.5)
  anova(lm_model_nucl_diversity_0.5)
  
  lm_model_Eucli_MDS1_0.5 <- lm(Axis1 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS1_0.5
  summary(lm_model_Eucli_MDS1_0.5)
  anova(lm_model_Eucli_MDS1_0.5)
  
  lm_model_Eucli_MDS2_0.5 <- lm(Axis2 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS2_0.5
  summary(lm_model_Eucli_MDS2_0.5)
  anova(lm_model_Eucli_MDS2_0.5)
  
  
  ### Time = 1
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 1)
  lm_model_snvrate_1 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_snvrate_1
  summary(lm_model_snvrate_1)
  anova(lm_model_snvrate_1)
  
  lm_model_nucl_diversity_1 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_nucl_diversity_1
  summary(lm_model_nucl_diversity_1)
  anova(lm_model_nucl_diversity_1)
  
  lm_model_Eucli_MDS1_1 <- lm(Axis1 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS1_1
  summary(lm_model_Eucli_MDS1_1)
  anova(lm_model_Eucli_MDS1_1)
  
  lm_model_Eucli_MDS2_1 <- lm(Axis2 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS2_1
  summary(lm_model_Eucli_MDS2_1)
  anova(lm_model_Eucli_MDS2_1)
  
  
  ### Time = 6
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 6)
  lm_model_snvrate_6 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_snvrate_6
  summary(lm_model_snvrate_6)
  anova(lm_model_snvrate_6)
  
  lm_model_nucl_diversity_6 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_nucl_diversity_6
  summary(lm_model_nucl_diversity_6)
  anova(lm_model_nucl_diversity_6)
  
  lm_model_Eucli_MDS1_6 <- lm(Axis1 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS1_6
  summary(lm_model_Eucli_MDS1_6)
  anova(lm_model_Eucli_MDS1_6)
  
  lm_model_Eucli_MDS2_6 <- lm(Axis2 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS2_6
  summary(lm_model_Eucli_MDS2_6)
  anova(lm_model_Eucli_MDS2_6)
  
  
  ### Time = 12
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 12)
  lm_model_snvrate_12 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_snvrate_12
  summary(lm_model_snvrate_12)
  anova(lm_model_snvrate_12)
  
  lm_model_nucl_diversity_12 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_nucl_diversity_12
  summary(lm_model_nucl_diversity_12)
  anova(lm_model_nucl_diversity_12)
  
  lm_model_Eucli_MDS1_12 <- lm(Axis1 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS1_12
  summary(lm_model_Eucli_MDS1_12)
  anova(lm_model_Eucli_MDS1_12)
  
  lm_model_Eucli_MDS2_12 <- lm(Axis2 ~ Delivery + Gender + Term + Feed + Country, data = test)
  lm_model_Eucli_MDS2_12
  summary(lm_model_Eucli_MDS2_12)
  anova(lm_model_Eucli_MDS2_12)
  
  
  ### Time = 24
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Feed) # 24, 36
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Country) # 36
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 24)
  lm_model_snvrate_24 <- lm(SNV_rate_mean ~ Delivery + Gender + Term + Country, data = test)
  lm_model_snvrate_24
  summary(lm_model_snvrate_24)
  anova(lm_model_snvrate_24)
  
  lm_model_nucl_diversity_24 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term + Country, data = test)
  lm_model_nucl_diversity_24
  summary(lm_model_nucl_diversity_24)
  anova(lm_model_nucl_diversity_24)
  
  lm_model_Eucli_MDS1_24 <- lm(Axis1 ~ Delivery + Gender + Term + Country, data = test)
  lm_model_Eucli_MDS1_24
  summary(lm_model_Eucli_MDS1_24)
  anova(lm_model_Eucli_MDS1_24)
  
  lm_model_Eucli_MDS2_24 <- lm(Axis2 ~ Delivery + Gender + Term + Country, data = test)
  lm_model_Eucli_MDS2_24
  summary(lm_model_Eucli_MDS2_24)
  anova(lm_model_Eucli_MDS2_24)
  
  
  ### Time = 36
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Feed) # 24, 36
  table(cluster_meta_mds$Time_new2, cluster_meta_mds$Country) # 36
  test <- cluster_meta_mds %>%
    subset(Time_new2 == 36)
  lm_model_snvrate_36 <- lm(SNV_rate_mean ~ Delivery + Gender + Term, data = test)
  lm_model_snvrate_36
  summary(lm_model_snvrate_36)
  anova(lm_model_snvrate_36)
  
  lm_model_nucl_diversity_36 <- lm(nucl_diversity_mean ~ Delivery + Gender + Term, data = test)
  lm_model_nucl_diversity_36
  summary(lm_model_nucl_diversity_36)
  anova(lm_model_nucl_diversity_36)
  
  lm_model_Eucli_MDS1_36 <- lm(Axis1 ~ Delivery + Gender + Term, data = test)
  lm_model_Eucli_MDS1_36
  summary(lm_model_Eucli_MDS1_36)
  anova(lm_model_Eucli_MDS1_36)
  
  lm_model_Eucli_MDS2_36 <- lm(Axis2 ~ Delivery + Gender + Term, data = test)
  lm_model_Eucli_MDS2_36
  summary(lm_model_Eucli_MDS2_36)
  anova(lm_model_Eucli_MDS2_36)
  
  
  summary(lm_model_snvrate_0)$coefficients
  summary(lm_model_snvrate_0.5)$coefficients
  summary(lm_model_snvrate_1)$coefficients
  summary(lm_model_snvrate_12)$coefficients
  summary(lm_model_snvrate_24)$coefficients
  summary(lm_model_snvrate_36)$coefficients
  
  
  
  names(coef(lm_model_snvrate_36))
  
  # Update extraction function
  extract_anova_info <- function(model, time, feature_name) {
    if (is.null(model)) return(NULL)
    
    anova_tbl <- tryCatch({
      car::Anova(model, type = "II")
    }, error = function(e) return(NULL))
    
    if (is.null(anova_tbl)) return(NULL)
    
    ss_total <- sum(anova_tbl$`Sum Sq`, na.rm = TRUE)
    
    coefs <- coef(model)
    coef_names <- names(coefs)
    
    # Retain only the signs of main-effect variables, excluding intercept and interactions
    signs <- sapply(rownames(anova_tbl), function(var) {
      idx <- grep(paste0("^", var), coef_names)
      if (length(idx) > 0) {
        sign_val <- sign(mean(coefs[idx], na.rm = TRUE))
        if (is.na(sign_val)) return(0)
        return(sign_val)
      } else {
        return(0)
      }
    })
    
    result <- data.frame(
      Timepoint = time,
      GroupVar = rownames(anova_tbl),
      Feature = feature_name,
      p.value = anova_tbl$`Pr(>F)`,
      Partial_R2 = anova_tbl$`Sum Sq` / ss_total,
      signs = signs,
      R2_signed = (anova_tbl$`Sum Sq` / ss_total) * signs
    )
    
    return(result)
  }
  
  
  model_list <- list(
    
    # SNV rate
    list(model = lm_model_snvrate_0,   time = 0,   feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_0.5, time = 0.5, feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_1,   time = 1,   feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_6,   time = 6,   feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_12,  time = 12,  feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_24,  time = 24,  feature = "SNV_rate_mean"),
    list(model = lm_model_snvrate_36,  time = 36,  feature = "SNV_rate_mean"),
    
    # nucleotide diversity
    list(model = lm_model_nucl_diversity_0,   time = 0,   feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_0.5, time = 0.5, feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_1,   time = 1,   feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_6,   time = 6,   feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_12,  time = 12,  feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_24,  time = 24,  feature = "nucl_diversity_mean"),
    list(model = lm_model_nucl_diversity_36,  time = 36,  feature = "nucl_diversity_mean"),
    
    # Axis1
    list(model = lm_model_Eucli_MDS1_0,   time = 0,   feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_0.5, time = 0.5, feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_1,   time = 1,   feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_6,   time = 6,   feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_12,  time = 12,  feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_24,  time = 24,  feature = "Axis1"),
    list(model = lm_model_Eucli_MDS1_36,  time = 36,  feature = "Axis1"),
    
    # Axis2
    list(model = lm_model_Eucli_MDS2_0,   time = 0,   feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_0.5, time = 0.5, feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_1,   time = 1,   feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_6,   time = 6,   feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_12,  time = 12,  feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_24,  time = 24,  feature = "Axis2"),
    list(model = lm_model_Eucli_MDS2_36,  time = 36,  feature = "Axis2")
  )
  
  
  
  ## Extract and combine results
  results_list <- lapply(model_list, function(x) {
    extract_anova_info(x$model, x$time, x$feature)
  })
  
  results_df <- do.call(rbind, results_list)
  results_df$Timepoint <- as.character(results_df$Timepoint)
  
  ## Fill all combinations, including NA values
  # All combinations
  all_combos <- expand.grid(
    # Timepoint = sort(unique(cluster_meta_mds$Time_new2)),
    Timepoint = c("0", "0.5", "1", "6", "12", "24", "36"),
    GroupVar = c("Delivery", "Gender", "Term", "Feed", "Country"),
    Feature = c("SNV_rate_mean", "nucl_diversity_mean", "Axis1", "Axis2"),
    stringsAsFactors = FALSE
  )
  
  ## Update plotting data processing
  results_df_filled <- all_combos %>%
    left_join(results_df, by = c("Timepoint", "GroupVar", "Feature")) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("0", "0.5", "1", "6", "12", "24", "36")),
      
      sig_label = case_when(
        is.na(p.value) & is.na(Partial_R2) ~ "~",  # aliased
        is.na(p.value) ~ "-",                    # not available
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # For variables with more than two levels, remove directionality
  results_df_filled$signs <- ifelse(results_df_filled$GroupVar %in% c("Feed", "Country"),
                                    NA, results_df_filled$signs)
  
  results_df_filled$Partial_R2 <- ifelse(is.na(results_df_filled$signs),
                                         abs(results_df_filled$Partial_R2),
                                         results_df_filled$Partial_R2
  )
  results_df_filled #
  
  
  
  table(results_df_filled$GroupVar)
  results_df_filled$GroupVar <- factor(results_df_filled$GroupVar,
                                       levels = c("Term", "Delivery", "Gender", "Country", "Feed"))
  table(results_df_filled$Feature)
  results_df_filled$Feature <- factor(results_df_filled$Feature,
                                      levels = c("Axis2", "Axis1", "nucl_diversity_mean", "SNV_rate_mean", "DMM_Cluster"))
  
  ## Update heatmap fill to R2_signed
  # Split data
  df_with_sign  <- results_df_filled %>% filter(!is.na(signs))
  df_no_sign    <- results_df_filled %>% filter(is.na(signs))
  
  library(ggplot2)
  library(dplyr)
  library(ggnewscale)
  library(cowplot)
  
  ## Build main plot
  library("ggnewscale")
  p_lm_multi_envir <- ggplot() +
    ## Layer with direction, blue-white-red
    geom_tile(
      data = df_with_sign,
      aes(x = as.factor(Timepoint), y = Feature, fill = R2_signed),
      color = "white"
    ) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Signed Partial R²\n& Effect direction",
      na.value = "grey90",
      guide = guide_colorbar(order = 1)
    ) +
    
    ## Reset fill mapping
    ggnewscale::new_scale_fill() +
    
    ## Layer without direction, grey-scale gradient
    geom_tile(
      data = df_no_sign,
      aes(x = as.factor(Timepoint), y = Feature, fill = Partial_R2),
      color = "white"
    ) +
    scale_fill_gradient(
      low = "white", high = "darkgrey",
      name = "Partial R²\n(no direction)",
      na.value = "grey90",
      guide = guide_colorbar(order = 2)
    ) +
    
    ## Significance labels
    geom_text(
      data = results_df_filled,
      aes(x = as.factor(Timepoint), y = Feature, label = sig_label),
      size = 4
    ) +
    
    facet_wrap(~ GroupVar, ncol = 1) +
    labs(
      x = "Infant age (months)",
      y = "Microbial Feature",
      title = "2.Group Differences in Microbial Features by Timepoint"
    ) +
    theme_minimal(base_size = 11)
  
  ## Add right-side annotation text
  caption_text <- "Significance codes:\n*** p < 0.001\n** p < 0.01\n* p < 0.05\n~ aliased\n- not available"
  
  # Combine main plot and annotation text
  library("cowplot")
  p_lm_multi_envir2 <- ggdraw() +
    draw_plot(p_lm_multi_envir, x = 0, width = 0.85) +
    draw_label(caption_text, x = 0.685, y = 0.9, hjust = 0, size = 10)
  p_lm_multi_envir2
  # save as 8.5 * 6
  
  # Check specific directions
  summary(lm_model_snvrate_0)$coefficients
  summary(lm_model_snvrate_0.5)$coefficients
  summary(lm_model_snvrate_1)$coefficients
  summary(lm_model_snvrate_12)$coefficients
  summary(lm_model_snvrate_24)$coefficients
  summary(lm_model_snvrate_36)$coefficients
  
}

### Summarize Fig. 5b multivariable results
fig5b_pvalues <- results_df_filled %>%
  dplyr::mutate(
    Figure = "Fig. 5b",
    Model = "Multivariable linear model",
    Term_contribution = Partial_R2,
    Signed_term_contribution = R2_signed,
    P_label = dplyr::case_when(
      is.na(p.value) & is.na(Term_contribution) ~ "aliased",
      is.na(p.value) ~ "not available",
      p.value < 0.001 ~ "P < 0.001",
      p.value < 0.01 ~ "P < 0.01",
      p.value < 0.05 ~ "P < 0.05",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    Figure,
    Model,
    Timepoint,
    GroupVar,
    Feature,
    p.value,
    P_label,
    Term_contribution,
    Signed_term_contribution,
    signs,
    sig_label
  ) %>%
  dplyr::arrange(Timepoint, GroupVar, Feature)

write_xlsx(fig5b_pvalues, "Fig5b_multivariable_pvalues_term_contribution.xlsx")





### Fig. 5c Associations between environmental factors and SNV-defined DMM clusters across infant age  ############

### Associations between DMM clusters and environmental factors
library(dplyr)
library(tidyr)
library(ggplot2)
library(rcompanion)

## 1. Set factor order
table(cluster_meta_mds$Feed)
cluster_meta_mds <- cluster_meta_mds %>%
  mutate(
    Time_new2 = factor(
      as.character(Time_new2),
      levels = c("0", "0.5", "1", "6", "12", "24", "36")
    ),
    DMM_cluster = factor(DMM_cluster),
    Delivery = factor(Delivery, levels = c("vaginal", "CS")),
    Gender = factor(Gender, levels = c("male", "female")),
    Term = factor(Term, levels = c("full_term", "preterm")),
    Feed = factor(Feed, levels = c("breast", "combined", "formula")),
    Country = factor(Country)
  )

env_vars <- c("Delivery", "Gender", "Term", "Feed", "Country")
timepoints <- levels(cluster_meta_mds$Time_new2)

## 2. Define function: test DMM_cluster ~ environment at each time point
run_dmm_env_test <- function(data, timepoint, env_var) {
  
  df <- data %>%
    filter(Time_new2 == timepoint) %>%
    select(DMM_cluster, all_of(env_var)) %>%
    filter(!is.na(DMM_cluster), !is.na(.data[[env_var]])) %>%
    droplevels()
  
  if (nrow(df) < 10 ||
      length(unique(df$DMM_cluster)) < 2 ||
      length(unique(df[[env_var]])) < 2) {
    return(data.frame(
      Timepoint = timepoint,
      GroupVar = env_var,
      p.value = NA_real_,
      CramerV = NA_real_,
      n = nrow(df)
    ))
  }
  
  tab <- table(df[[env_var]], df$DMM_cluster)
  
  ## Chi-square test; use Fisher's exact test when expected counts are too small
  chi <- suppressWarnings(chisq.test(tab))
  
  pval <- if (any(chi$expected < 5)) {
    tryCatch(
      fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value,
      error = function(e) NA_real_
    )
  } else {
    chi$p.value
  }
  
  cv <- tryCatch(
    cramerV(tab),
    error = function(e) NA_real_
  )
  
  data.frame(
    Timepoint = timepoint,
    GroupVar = env_var,
    p.value = pval,
    CramerV = cv,
    n = nrow(df)
  )
}

## 3. Run tests in batch
dmm_env_results <- bind_rows(
  lapply(timepoints, function(tp) {
    bind_rows(
      lapply(env_vars, function(v) {
        run_dmm_env_test(cluster_meta_mds, tp, v)
      })
    )
  })
)

## 4. Add significance labels
dmm_env_results <- dmm_env_results %>%
  mutate(
    sig_label = case_when(
      is.na(p.value) ~ "-",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Timepoint = factor(Timepoint, levels = c("0", "0.5", "1", "6", "12", "24", "36")),
    GroupVar = factor(GroupVar, levels = c("Feed", "Country", "Gender", "Delivery", "Term"))
  )
# View(dmm_env_results)


# Heatmap
colnames(dmm_env_results)
table(dmm_env_results$GroupVar)

p_dmm_env2 <- ggplot(
  dmm_env_results,
  aes(x = Timepoint,
      y = GroupVar,
      fill = CramerV)
) +
  
  geom_tile(color = "white") +
  
  geom_text(aes(label = sig_label),
            size = 4) +
  
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "grey90",
    name = "Cramér's V"
  ) +
  
  labs(
    x = "Infant age (months)",
    y = "Environmental factor",
    title = "Associations between environmental factors and DMM clusters"
  ) +
  
  theme_minimal(base_size = 12)

print(p_dmm_env2)
# save as 4.5 * 3

fig5c_pvalues <- dmm_env_results %>%
  dplyr::mutate(
    Figure = "Fig. 5c",
    Model = "DMM cluster association test",
    P_label = dplyr::case_when(
      is.na(p.value) ~ "not available",
      p.value < 0.001 ~ "P < 0.001",
      p.value < 0.01 ~ "P < 0.01",
      p.value < 0.05 ~ "P < 0.05",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    Figure,
    Model,
    Timepoint,
    GroupVar,
    p.value,
    P_label,
    CramerV,
    n,
    sig_label
  ) %>%
  dplyr::arrange(Timepoint, GroupVar)

fig5c_pvalues

write_xlsx(
  fig5c_pvalues,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/Fig5c_DMM_cluster_environment_pvalues.xlsx"
)





### Fig. 5d Recursive partitioning tree of host and environmental determinants of DMM-defined SNV clusters  ############


### DMM cluster ~ environmental factors, partitioning tree
## First, use rpart to build the partitioning tree and identify which variables enter the tree and how each node is split.
## Then, calculate chi-square test P values for each rpart-derived split rule and summarize P values and branch proportions for the final figure.

cluster_meta_mds #
class(cluster_meta_mds)
colnames(cluster_meta_mds)
dim(cluster_meta_mds)
head(cluster_meta_mds)

### Data preparation
library(dplyr)
library(data.table)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(tidyr)
library(forcats)

set.seed(123)

# Safe function to convert NA values into an explicit factor level
na_to_level <- function(x, levels_order = NULL, na_level = "Unknown") {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- na_level
  
  if (!is.null(levels_order)) {
    factor(x, levels = c(levels_order, na_level))
  } else {
    factor(x)
  }
}

tree_df <- as.data.frame(cluster_meta_mds) %>%
  mutate(
    DMM_cluster = factor(DMM_cluster),
    
    # Use numeric age so the tree can automatically identify age split points
    age_month = as.numeric(as.character(Time_new2)),
    
    Time_new2 = factor(as.character(Time_new2),
                       levels = c("0", "0.5", "1", "6", "12", "24", "36")),
    
    Delivery = na_to_level(Delivery, levels_order = c("vaginal", "CS")),
    Gender = na_to_level(Gender, levels_order = c("male", "female")),
    Term = na_to_level(Term, levels_order = c("full_term", "preterm")),
    Feed = na_to_level(Feed, levels_order = c("breast", "combined", "formula")),
    Country = na_to_level(Country),
    Continent = na_to_level(Continent),
    Study = factor(Study),
    SubjectID = factor(SubjectID)
  ) %>%
  filter(!is.na(DMM_cluster), !is.na(age_month))

dim(tree_df) # 5680 * 20
table(tree_df$DMM_cluster)
table(tree_df$Time_new2, useNA = "ifany")
table(tree_df$age_month, useNA = "ifany")
table(tree_df$Term, useNA = "ifany")
table(tree_df$Feed, useNA = "ifany")
table(tree_df$Feed)



# Based on preliminary analysis, restrict tree complexity first
table(tree_df$age_month)
table(tree_df$Feed)
str(tree_df)

fit_tree_simple_raw <- rpart(
  DMM_cluster ~ age_month + Term + Delivery + Gender + Country + Feed,
  data = tree_df,
  method = "class",
  parms = list(split = "gini"),
  control = rpart.control(
    cp = 0.003,
    minsplit = 200,
    minbucket = 100,
    maxdepth = 4,
    xval = 10
  )
)

p_tree_all <- rpart.plot(
  fit_tree_simple_raw,
  type = 2,
  extra = 104,
  under = TRUE,
  fallen.leaves = TRUE,
  faclen = 0,
  varlen = 0,
  box.palette = "GnBu",
  cex = 0.85,
  main = "Determinants of SNV-defined DMM clusters"
)


##### Chi-square tests ####
dim(tree_df)

# Node 1: Time_new2 <= 6 months
class(tree_df$Time_new2)
tree_df$Time_new2_6m <- ifelse(
  as.numeric(as.character(tree_df$Time_new2)) > 6,
  "older",
  "younger"
)
chisq.test(table(tree_df$Time_new2_6m, tree_df$DMM_cluster))
# p-value < 2.2e-16
# 602/5680

### Branch 1
# table(test1$Term)
test1_1 <- tree_df %>% subset(Time_new2_6m == "younger")
test1_1$Term_group <- ifelse(test1_1$Term %in% c("preterm"), "G2", "G1")
chisq.test(table(test1_1$Term_group, test1_1$DMM_cluster))
# p-value < 2.2e-16

table(test1_1$Delivery)
test1_2 <- test1_1 %>%
  subset(Time_new2_6m == "younger" & Term_group %in% c("G1"))
test1_2$Delivery_group <- ifelse(test1_2$Delivery %in% c("CS"), "G2", "G1")
chisq.test(table(test1_2$Delivery_group, test1_2$DMM_cluster))
# p-value < 2.2e-16

table(test1_2$Time_new2)
test1_3 <- test1_2 %>%
  subset(Time_new2_6m == "younger" &
           Term_group %in% c("G1") &
           Delivery_group %in% c("G2"))
test1_3$Time_new2_1m <- ifelse(
  as.numeric(as.character(test1_3$Time_new2)) > 1,
  "older",
  "younger"
)
chisq.test(table(test1_3$Time_new2_1m, test1_3$DMM_cluster))
# p-value < 2.2e-16

### Branch 2
table(tree_df$Country)
test2_1 <- tree_df %>% subset(Time_new2_6m == "older")
table(test2_1$Country)
test2_1$Country_group <- ifelse(
  test2_1$Country %in%
    c("Estonia", "Finland", "Italy", "Russia", "Sweden", "UK", "US"),
  "G1",
  "G2"
)

tab <- table(test2_1$Country_group, test2_1$DMM_cluster)
tab2 <- tab[rowSums(tab) > 0, colSums(tab) > 0]
chisq.test(tab2, simulate.p.value = TRUE, B = 9999)
# p-value = 1e-04


table(test2_1$Country)
test2_2 <- test2_1 %>%
  subset(Country %in%
           c("Estonia", "Finland", "Italy", "Russia", "Sweden", "UK", "US") &
           Country_group == "G1")
table(test2_2$Country)
test2_2$Country_group <- ifelse(
  test2_2$Country %in%
    c("Estonia", "Finland", "Sweden"),
  "G1",
  "G2"
)

tab <- table(test2_2$Country_group, test2_2$DMM_cluster)
tab
tab2 <- tab[rowSums(tab) > 0, colSums(tab) > 0]
tab2
chisq.test(tab2, simulate.p.value = TRUE, B = 9999)
# p-value = 1e-04


table(test2_1$Country)
test2_3 <- test2_2 %>%
  subset(Country_group == "G2" &
           Country %in% c("Italy", "Russia", "UK", "US"))
table(test2_3$Country)
test2_3$Country_group <- ifelse(
  test2_3$Country %in%
    c("Russia", "UK"),
  "G1",
  "G2"
)

tab <- table(test2_3$Country_group, test2_3$DMM_cluster)
tab
tab2 <- tab[rowSums(tab) > 0, colSums(tab) > 0]
tab2
chisq.test(tab2, simulate.p.value = TRUE, B = 9999)
# p-value = 1e-04




#### Collect chi-square / Fisher test results for Fig. 5d splits

library(dplyr)
library(writexl)

# Function to calculate Cramer's V
calc_cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab))
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
}

# Function to run chi-square test or Fisher's exact test with simulated P value
run_split_test <- function(data, split_name, split_variable, split_group, cluster_var = "DMM_cluster") {
  
  df <- data %>%
    dplyr::select(
      all_of(split_group),
      all_of(cluster_var)
    ) %>%
    dplyr::filter(
      !is.na(.data[[split_group]]),
      !is.na(.data[[cluster_var]])
    ) %>%
    droplevels()
  
  tab <- table(df[[split_group]], df[[cluster_var]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
  
  chi <- suppressWarnings(chisq.test(tab))
  
  if (any(chi$expected < 5)) {
    test_res <- fisher.test(tab, simulate.p.value = TRUE, B = 9999)
    method_used <- "Fisher's exact test with simulated P value"
    statistic <- NA_real_
    df_value <- NA_real_
    p_value <- test_res$p.value
  } else {
    test_res <- chi
    method_used <- "Chi-square test"
    statistic <- as.numeric(test_res$statistic)
    df_value <- as.numeric(test_res$parameter)
    p_value <- test_res$p.value
  }
  
  group_n <- as.data.frame(table(df[[split_group]]))
  colnames(group_n) <- c("Group", "n")
  group_n_text <- paste0(group_n$Group, "=", group_n$n, collapse = "; ")
  
  cluster_n <- as.data.frame(table(df[[cluster_var]]))
  colnames(cluster_n) <- c("DMM_cluster", "n")
  cluster_n_text <- paste0(cluster_n$DMM_cluster, "=", cluster_n$n, collapse = "; ")
  
  data.frame(
    Split = split_name,
    Split_variable = split_variable,
    Test_method = method_used,
    n_total = nrow(df),
    group_n = group_n_text,
    cluster_n = cluster_n_text,
    statistic = statistic,
    df = df_value,
    p.value = p_value,
    CramerV = calc_cramers_v(tab),
    stringsAsFactors = FALSE
  )
}



#### Define split variables according to the tree structure

tree_df2 <- tree_df

# Split 1: age <= 6 months vs > 6 months
tree_df2$Split_age_6m <- ifelse(
  as.numeric(as.character(tree_df2$Time_new2)) > 6,
  ">6 months",
  "≤6 months"
)

# Younger branch: ≤6 months
test_young <- tree_df2 %>%
  dplyr::filter(Split_age_6m == "≤6 months")

# Split 2: term status among younger infants
test_young$Split_term_young <- ifelse(
  test_young$Term %in% c("preterm"),
  "Preterm",
  "Full-term/Unknown"
)

# Split 3: delivery mode among younger full-term/unknown infants
test_young_term_full <- test_young %>%
  dplyr::filter(Split_term_young == "Full-term/Unknown")

test_young_term_full$Split_delivery_young <- ifelse(
  test_young_term_full$Delivery %in% c("CS"),
  "C-section",
  "Vaginal/Unknown"
)

# Split 4: age <= 1 month vs 6 months among younger full-term/unknown C-section infants
test_young_term_full_cs <- test_young_term_full %>%
  dplyr::filter(Split_delivery_young == "C-section")

test_young_term_full_cs$Split_age_1m <- ifelse(
  as.numeric(as.character(test_young_term_full_cs$Time_new2)) > 1,
  "6 months",
  "≤1 month"
)

# Older branch: >6 months
test_old <- tree_df2 %>%
  dplyr::filter(Split_age_6m == ">6 months")

# Split 5: country group among older infants
test_old$Split_country_old_1 <- ifelse(
  test_old$Country %in% c("Estonia", "Finland", "Italy", "Russia", "Sweden", "UK", "US"),
  "Estonia/Finland/Italy/Russia/Sweden/UK/US",
  "Other countries"
)

# Split 6: country group among selected older countries
test_old_country_g1 <- test_old %>%
  dplyr::filter(Split_country_old_1 == "Estonia/Finland/Italy/Russia/Sweden/UK/US")

test_old_country_g1$Split_country_old_2 <- ifelse(
  test_old_country_g1$Country %in% c("Estonia", "Finland", "Sweden"),
  "Estonia/Finland/Sweden",
  "Italy/Russia/UK/US"
)

# Split 7: country group among Italy/Russia/UK/US
test_old_country_g2 <- test_old_country_g1 %>%
  dplyr::filter(Split_country_old_2 == "Italy/Russia/UK/US")

test_old_country_g2$Split_country_old_3 <- ifelse(
  test_old_country_g2$Country %in% c("Russia", "UK"),
  "Russia/UK",
  "Italy/US"
)


#### Run tests and combine into one table

tree_split_tests <- bind_rows(
  run_split_test(
    data = tree_df2,
    split_name = "Root split: ≤6 months vs >6 months",
    split_variable = "Age",
    split_group = "Split_age_6m"
  ),
  
  run_split_test(
    data = test_young,
    split_name = "Younger branch: preterm vs full-term/unknown",
    split_variable = "Term",
    split_group = "Split_term_young"
  ),
  
  run_split_test(
    data = test_young_term_full,
    split_name = "Younger full-term/unknown branch: C-section vs vaginal/unknown",
    split_variable = "Delivery",
    split_group = "Split_delivery_young"
  ),
  
  run_split_test(
    data = test_young_term_full_cs,
    split_name = "Younger full-term/unknown C-section branch: ≤1 month vs 6 months",
    split_variable = "Age",
    split_group = "Split_age_1m"
  ),
  
  run_split_test(
    data = test_old,
    split_name = "Older branch: selected countries vs other countries",
    split_variable = "Country",
    split_group = "Split_country_old_1"
  ),
  
  run_split_test(
    data = test_old_country_g1,
    split_name = "Older selected-country branch: Estonia/Finland/Sweden vs Italy/Russia/UK/US",
    split_variable = "Country",
    split_group = "Split_country_old_2"
  ),
  
  run_split_test(
    data = test_old_country_g2,
    split_name = "Older Italy/Russia/UK/US branch: Russia/UK vs Italy/US",
    split_variable = "Country",
    split_group = "Split_country_old_3"
  )
)

tree_split_tests <- tree_split_tests %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    p.label = case_when(
      p.value < 2.2e-16 ~ "P < 2.2e-16",
      p.value < 0.001 ~ "P < 0.001",
      p.value < 0.01 ~ "P < 0.01",
      p.value < 0.05 ~ "P < 0.05",
      TRUE ~ paste0("P = ", signif(p.value, 3))
    ),
    p.adj.label = case_when(
      p.adj < 2.2e-16 ~ "FDR < 2.2e-16",
      p.adj < 0.001 ~ "FDR < 0.001",
      p.adj < 0.01 ~ "FDR < 0.01",
      p.adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ paste0("FDR = ", signif(p.adj, 3))
    )
  )

tree_split_tests

write_xlsx(
  tree_split_tests,
  "D:/3_Projects/2_Children_jaundice/3_R_analysis/Fig5d_partitioning_tree_split_tests.xlsx"
)





#### Supplementary Fig S4b. Partitioning tree variable importance ######

### Variable importance: identify which variables the tree relies on most
var_imp <- data.frame(
  variable = names(fit_tree_simple_raw$variable.importance),
  importance = as.numeric(fit_tree_simple_raw$variable.importance)
) %>%
  arrange(desc(importance))

var_imp

p_var_imp <- ggplot(
  var_imp,
  aes(x = reorder(variable, importance), y = importance)
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = "",
    y = "Variable importance",
    title = "Variable importance in DMM cluster decision tree"
  )

p_var_imp
# save as 3.2 * 4


### Calculate tree explanatory performance
get_tree_explained <- function(fit) {
  frame <- fit$frame
  root_dev <- frame$dev[1]
  leaf_dev <- sum(frame$dev[frame$var == "<leaf>"])
  explained <- 1 - leaf_dev / root_dev
  explained
}

tree_explained_raw <- get_tree_explained(fit_tree_simple_raw)
tree_explained_raw
# 0.2820456
# The decision tree reduced classification impurity by 28.2%.
