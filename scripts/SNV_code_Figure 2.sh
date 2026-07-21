############### Gut microbial SNv ###################################
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

### Load phyloseq data
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/Result2_R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Result2_R.RData")


############ Result 2: The developmental trajectory of microbial SNVs in the gut from birth to 3 years old #####
#### Metadata information ####

meta_df <- read_excel("D:/3_Projects/2_Children_jaundice/3_R_analysis/41467_2022_32805_MOESM3_ESM.xlsx")
# SubjectID: child identifier
# DOL: days of life
dim(meta_df)# 6122
length(unique(meta_df$Run))# 3+7
length(unique(meta_df$SampleID))

colnames(meta_df)
str(meta_df)
summary(meta_df$DOL);mean(meta_df$DOL, na.rm = TRUE);sd(meta_df$DOL, na.rm = TRUE)
table(meta_df$Time_new)
meta_df$sample2 <- meta_df$Run
table(meta_df$Country);table(is.na(meta_df$Country))
unique(meta_df$Country)
# Italy, Sweden, US, Bangladesh, New_Zealand, Singapore, Finland, UK, Russia, Estonia, Luxembourg
# Russia spans Europe and Asia, but was assigned to Europe here because most samples are from the European region.
meta_df$Continent <- ifelse(meta_df$Country %in%c("Italy", "Sweden", "Finland", "UK","Estonia","Luxembourg","Russia"),"Europe",
                            ifelse(meta_df$Country %in%c("Bangladesh", "Singapore"),"Asia",
                                   ifelse(meta_df$Country %in%c("New_Zealand"),"Oceania",
                                          ifelse(meta_df$Country %in%c("US"),"North America",NA ))))
table(is.na(meta_df$Continent))
table(meta_df$Continent,meta_df$DOL)
colnames(meta_df)
meta_df$Study <- ifelse(meta_df$Study=="NA",NA,meta_df$Study)
meta_df$Run <- ifelse(meta_df$Run=="NA",NA,meta_df$Run)
meta_df$SampleID <- ifelse(meta_df$SampleID=="NA",NA,meta_df$SampleID)
meta_df$SubjectID <- ifelse(meta_df$SubjectID=="NA",NA,meta_df$SubjectID)
meta_df$DOL <- ifelse(meta_df$DOL=="NA",NA,meta_df$DOL)
meta_df$Delivery <- ifelse(meta_df$Delivery=="NA",NA,meta_df$Delivery)
meta_df$Country <- ifelse(meta_df$Country=="NA",NA,meta_df$Country)
meta_df$Gender <- ifelse(meta_df$Gender=="NA",NA,meta_df$Gender)
meta_df$Term <- ifelse(meta_df$Term=="NA",NA,meta_df$Term)
meta_df$Feed <- ifelse(meta_df$Feed=="NA",NA,meta_df$Feed)



##### Time column
meta_df$DOL <- as.numeric(meta_df$DOL)
meta_df <- meta_df %>%
  mutate(Time_new=case_when(DOL>=0 & DOL<=1 ~ "0",
                            DOL>1 & DOL<=14 ~ "0.5",
                            DOL>14 & DOL<=30 ~ "1",
                            DOL>30 & DOL<=90 ~ "3",
                            DOL>90 & DOL<=180 ~ "6",
                            DOL>180 & DOL<=360 ~ "12",
                            DOL>360 & DOL<=540 ~ "18",
                            DOL>540 & DOL<=720 ~ "24",
                            DOL>720 & DOL<=900 ~ "30",
                            DOL>900  ~ "36"),
         Time_new2=case_when(DOL>=0 & DOL<=1 ~ "0",
                             DOL>1 & DOL<=14 ~ "0.5",
                             DOL>14 & DOL<=30 ~ "1",
                             DOL>30 & DOL<=180 ~ "6",
                             DOL>180 & DOL<=540 ~ "12",
                             DOL>540 & DOL<=900 ~ "24",
                             DOL>900  ~ "36")
  )
meta_df$Time_new <- as.numeric(meta_df$Time_new)
meta_df$Time_new2 <- as.numeric(meta_df$Time_new2)
table(meta_df$Time_new,meta_df$Time_new2)
table(meta_df$Time_new)
#  0  0.5    1   12   18   24    3   30   36    6   months
# 326 1512 1128  942  390  309  896   76   44  424 
table(meta_df$Time_new2)
#  0  0.5    1    6   12   24   36 months
# 326 1512 1128 1320 1332  385   44 
meta_final <- meta_df
str(meta_final)
write_xlsx(meta_final, "meta_final.xlsx")




### Fig. 2a Microbial SNV rate across all time points ############
all_samples.genome_info.tsv #

library(tidyverse)

all_samples.genome_info <- read_tsv(
  "all_samples.genome_info.tsv",
  show_col_types = FALSE
)
View(all_samples.genome_info)

# Sample-level SNV rate 
# Count the total number of SNV sites in each sample across all genomes.
# Sum SNV counts across genomes for each sample.
colnames(all_samples.genome_info)
table(all_samples.genome_info$length==0)
table(all_samples.genome_info$breadth_minCov==0)# Contains zero values
table(all_samples.genome_info$SNV_count==0)# Contains zero values
snv_sample2 <- all_samples.genome_info %>%
  mutate(SNV_rate = SNV_count/(length*breadth_minCov+ 1e-10))

View(snv_sample2)

snv_sample3 <- snv_sample2 %>%
  group_by(sample2) %>%
  summarise(
    SNV_rate_mean = mean(SNV_rate),
    .groups = "drop"
  )
snv_sample3
table(snv_sample3$SNV_rate_mean==0)# T-10

# Merge data
meta_final <- meta_final %>%
  left_join(snv_sample3 %>% dplyr::select(sample2, SNV_rate_mean), by = "sample2")
str(meta_final)

## Plot: Time_new ~ SNV_rate_mean
str(meta_final)
table(meta_final$Time_new2)

meta_final$Time_new2 <- factor(meta_final$Time_new2)
table(meta_final$Time_new2)
library(ggpubr)

my_com <- list(c("0","0.5"), c("0.5","1"), c("1","6"), c("6","12"), c("12","24"), c("24","36"))

my_col <- c("0"="#8c510a", "0.5"="#bf812d", "1"="#dfc27d", "6"="#c7eae5", "12"="#80cdc1", 
            "24"="#35978f", "36"="#01665e")

test <- meta_final %>% subset(!is.na(Time_new2)) %>% as.data.frame()
stat.test <- compare_means(
  SNV_rate_mean ~ Time_new2,
  data = test,
  method = "wilcox.test",
  p.adjust.method = "fdr",
  comparisons = my_com
)

# Keep only the predefined pairwise comparisons in my_com.
stat.test <- stat.test %>%
  filter(interaction(group1, group2, sep = "") %in%
           interaction(map_chr(my_com, 1), map_chr(my_com, 2), sep = ""))

# Add y.position values for each comparison on the log10 scale.
#stat.test$y.position <- seq(0.005, 0.009, length.out = nrow(stat.test))
stat.test$y.position <- seq(0.00000001, 0.1, length.out = nrow(stat.test))

# Plot
p_snv_rate_overtime <- ggplot(test, aes(x = Time_new2, y = SNV_rate_mean)) +
  geom_jitter(aes(color = Time_new2), alpha = 0.5) +
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  scale_y_log10() +
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size = 2.5) +
  # ylim(0,0.1)+
  theme_bw() +
  scale_color_manual(values=my_col)+
  scale_fill_manual(values=my_col)+
  xlab("Chronologic age (months)")+
  ylab("Microbial SNV rate (log10)")+
  theme(panel.grid = element_blank(), text = element_text(size = 11))
p_snv_rate_overtime
# Save as 4.5 * 3.7





### Fig. 2b Microbial nucleotide diversity across all time points ############
# Sample-level nucleotide diversity
# Count the total number of SNV sites in each sample across all genomes.
# Sum SNV counts across genomes for each sample.
colnames(all_samples.genome_info)
table(all_samples.genome_info$length==0)
table(all_samples.genome_info$breadth_minCov==0)# Contains zero values
table(all_samples.genome_info$SNV_count==0)# Contains zero values
table(all_samples.genome_info$nucl_diversity==0)# Contains zero values
View(all_samples.genome_info)
range(all_samples.genome_info$coverage)
table(all_samples.genome_info$coverage==0)
table(all_samples.genome_info$coverage>1)

snv_sample4 <- all_samples.genome_info %>%
  subset(!is.na(nucl_diversity)) %>%
  mutate(eff_len = length * breadth_minCov * log10(coverage+1)) %>%
  group_by(sample2) %>%
  summarise(
    pi_sample = weighted.mean(nucl_diversity, w = eff_len, na.rm = TRUE),
    .groups = "drop"
  )

snv_sample4
table(snv_sample4$pi_sample==0)# T-2

# Merge data
meta_final <- meta_final %>%
  left_join(snv_sample4 %>% dplyr::select(sample2, pi_sample), by = "sample2")
str(meta_final)
write_xlsx(meta_final, "meta_final.xlsx")


## Plot: Time_new ~ nucl_diversity_mean
str(meta_final)
table(meta_final$Time_new2)

meta_final$Time_new2 <- factor(meta_final$Time_new2,levels = c(0,0.5,1,6,12,24,36))
table(meta_final$Time_new2)
library(ggpubr)

my_com <- list(c("0","0.5"), c("0.5","1"), c("1","6"), c("6","12"), c("12","24"), c("24","36"))
my_col <- c("0"="#8c510a", "0.5"="#bf812d", "1"="#dfc27d", "6"="#c7eae5", "12"="#80cdc1", 
            "24"="#35978f", "36"="#01665e")

test <- meta_final %>% subset(!is.na(Time_new2)) %>% as.data.frame()
stat.test <- compare_means(
  pi_sample ~ Time_new2,
  data = test,
  method = "wilcox.test",
  p.adjust.method = "fdr",
  comparisons = my_com
)

# Keep only the predefined pairwise comparisons in my_com.
stat.test <- stat.test %>%
  filter(interaction(group1, group2, sep = "") %in%
           interaction(map_chr(my_com, 1), map_chr(my_com, 2), sep = ""))

# Add y.position values for each comparison on the log10 scale.
#stat.test$y.position <- seq(0.2, 1, length.out = nrow(stat.test))
stat.test$y.position <- seq(0.00000001, 0.1, length.out = nrow(stat.test))

# Plot
p_pi_nucl_diversity_overtime <- ggplot(test, aes(x = Time_new2, y = pi_sample)) +
  geom_jitter(aes(color = Time_new2), alpha = 0.5) +
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  scale_y_log10() +
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size = 2.5) +
  theme_bw() +
  scale_color_manual(values=my_col)+
  scale_fill_manual(values=my_col)+
  xlab("Chronologic age (months)")+
  ylab("Weighted nucleotide diversity (log10)")+
  theme(panel.grid = element_blank(), text = element_text(size = 11))
p_pi_nucl_diversity_overtime
# Save as 4.5 * 3.7


p_snv_rate_overtime
p_nucl_diversity_overtime




### Fig. 2c Rate of change in SNV rate and identification of developmental turning points ############

meta_final#

colnames(meta_final)
test <- meta_final %>% subset(!is.na(Time_new2)) %>% as.data.frame()
str(test)
dim(test)# 6047

test$SNV_rate_mean
test$DOL

# Fit GAM: SNV rate ~ smooth(age)
gam_model <- mgcv::gam(SNV_rate_mean ~ s(DOL, bs = "cs"), data = test)
# snv_richness

# First derivative of the smooth with 95% CI
#remotes::install_github("gavinsimpson/gratia")
library("nanonext")
library("gratia")
library("mgcv")
derivs <- gratia::derivatives(gam_model, select = "s(DOL)", interval = "confidence")

# Compatibility shim: handle possible column name differences across gratia versions.
nm <- names(derivs)
if ("derivative" %in% nm)  names(derivs)[names(derivs) == "derivative"] <- ".derivative"
if ("lower" %in% nm)       names(derivs)[names(derivs) == "lower"]      <- ".lower_ci"
if ("upper" %in% nm)       names(derivs)[names(derivs) == "upper"]      <- ".upper_ci"

# Keep the derivative dataframe tidy.
deriv_df <- derivs %>%dplyr::select(DOL, .derivative, .lower_ci, .upper_ci)

# Turning points:
# 1) Maximum positive rate of change, representing the fastest increase.
peak1 <- deriv_df %>% slice_max(order_by = .derivative, n = 1)

# 2) Maximum negative rate of change, representing the fastest decrease.
peak3 <- deriv_df %>% slice_min(order_by = .derivative, n = 1)

# 3) Stabilization point: derivative close to 0 with confidence interval crossing 0;
# choose the point with the smallest absolute derivative.
deriv_df <- deriv_df %>%
  mutate(abs_deriv = abs(.derivative))

stable_zone <- deriv_df %>%
  filter(.lower_ci < 0 & .upper_ci > 0) %>% 
  slice_min(abs_deriv, n = 1)

# Combine turning points.
turning_points <- bind_rows(
  peak1 %>% mutate(type = "max increase"),
  stable_zone %>% mutate(type = "stable"),
  peak3 %>% mutate(type = "max decrease")
)
View(data.frame(turning_points))
print(turning_points)
# SNV_rate
#  DOL   .derivative   .lower_ci   .upper_ci type             abs_deriv
#  1  35.2  0.0000799     0.0000571   0.000103   max increase NA           
# 775.  -0.0000000371 -0.00000238  0.00000230 stable        0.0000000371
#   0   -0.000231     -0.000264   -0.000197   max decrease NA           


# Plot the first derivative with shaded 95% CI and vertical lines at turning points.
p_GAM_SNVrate_turningpoint_prediction <- ggplot(deriv_df, aes(x = DOL, y = .derivative)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "blue", alpha = 0.2) +
  geom_vline(data = turning_points, aes(xintercept = DOL, color = type),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("max increase" = "red",
                                "stable" = "orange",
                                "max decrease" = "green")) +
  scale_x_continuous(breaks = seq(0, max(test$DOL, na.rm = TRUE), by = 12)) +
  labs(#title = "First derivative of SNV rate: three turning points",
    x = "Age (days)", 
    y = "Rate of change in SNV rate\t
       (d(SNV rate)/d(Age))", 
    color = "Turning Point") +
  # d(SNV rate)/d(Age) quantifies the instantaneous rate of change in SNV rate with age.
  theme_minimal(base_size = 11)+
  theme(legend.position = "top")+
  scale_x_continuous(
    breaks = seq(0, max(test$DOL, na.rm = TRUE), by = 180),
    expand = expansion(mult = c(0.01, 0.01))) 
p_GAM_SNVrate_turningpoint_prediction
# Save as 4 * 3.5








### Fig. 2d Rate of change in nucleotide diversity and identification of developmental turning points ############

meta_final#

# Fit GAM: nucleotide diversity ~ smooth(age)
gam_model <- mgcv::gam(pi_sample ~ s(DOL, bs = "cs"), data = test)


# First derivative of the smooth with 95% CI
#remotes::install_github("gavinsimpson/gratia")
derivs <- gratia::derivatives(gam_model, select = "s(DOL)", interval = "confidence")

# Compatibility shim: handle possible column name differences across gratia versions.
nm <- names(derivs)
if ("derivative" %in% nm)  names(derivs)[names(derivs) == "derivative"] <- ".derivative"
if ("lower" %in% nm)       names(derivs)[names(derivs) == "lower"]      <- ".lower_ci"
if ("upper" %in% nm)       names(derivs)[names(derivs) == "upper"]      <- ".upper_ci"

# Keep the derivative dataframe tidy.
deriv_df <- derivs %>%dplyr::select(DOL, .derivative, .lower_ci, .upper_ci)

# Turning points:
# 1) Maximum positive rate of change, representing the fastest increase.
peak1 <- deriv_df %>% slice_max(order_by = .derivative, n = 1)

# 2) Maximum negative rate of change, representing the fastest decrease.
peak3 <- deriv_df %>% slice_min(order_by = .derivative, n = 1)

# 3) Stabilization point: derivative close to 0 with confidence interval crossing 0;
# choose the point with the smallest absolute derivative.
deriv_df <- deriv_df %>%
  mutate(abs_deriv = abs(.derivative))

stable_zone <- deriv_df %>%
  filter(.lower_ci < 0 & .upper_ci > 0) %>% 
  slice_min(abs_deriv, n = 1)

# Combine turning points.
turning_points <- bind_rows(
  peak1 %>% mutate(type = "max increase"),
  stable_zone %>% mutate(type = "stable"),
  peak3 %>% mutate(type = "max decrease")
)
View(data.frame(turning_points))
print(turning_points)
# pi----- nucleotide diversity
# DOL  .derivative   .lower_ci   .upper_ci type            abs_deriv
# 35.2  0.0000795    0.0000631   0.0000958  max increase     NA          
# 739.   0.000000147 -0.00000134  0.00000163 stable        0.000000147
#   0   -0.000158    -0.000182   -0.000134   max decrease    NA      


# Plot the first derivative with shaded 95% CI and vertical lines at turning points.
p_GAM_nucl_diversity_turningpoint_prediction <- ggplot(deriv_df, aes(x = DOL, y = .derivative)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "blue", alpha = 0.2) +
  geom_vline(data = turning_points, aes(xintercept = DOL, color = type),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("max increase" = "red",
                                "stable" = "orange",
                                "max decrease" = "green")) +
  scale_x_continuous(breaks = seq(0, max(test$DOL, na.rm = TRUE), by = 12)) +
  labs(#title = "First derivative of nucleotide diversity: three turning points",
    x = "Age (days)", 
    y = "Rate of change in nucleotide diversity\t
       (d(Nucleotide diversity)/d(Age))", 
    color = "Turning Point") +
  # d(Nucleotide diversity)/d(Age) quantifies the instantaneous rate of change in nucleotide diversity with age.
  theme_minimal(base_size = 11)+
  theme(legend.position = "top")+
  scale_x_continuous(
    breaks = seq(0, max(test$DOL, na.rm = TRUE), by = 180),
    expand = expansion(mult = c(0.01, 0.01))) 
p_GAM_nucl_diversity_turningpoint_prediction
# Save as 4 * 3.5

p_GAM_SNVrate_turningpoint_prediction
p_GAM_nucl_diversity_turningpoint_prediction





### Fig. 2e Bubble plots showing age-associated dynamics of SNV rate and relative abundance for the top 30 bacterial species ############


#### R code:

meta_final <- read_excel("meta_final.xlsx")
meta_final #
head(meta_final)
class(meta_final)

## Read sample-genome information
library(tidyverse)
library("readr")

all_samples.genome_info <- read_tsv("all_samples.genome_info.tsv", show_col_types = FALSE)
View(all_samples.genome_info)

colnames(genome_species_SNVs_merged)
table(genome_species_SNVs_merged$genome_sp %in% all_samples.genome_info$genome)
# TRUE -1790776
genome_species_SNVs_merged$genome <- genome_species_SNVs_merged$genome_sp
genome_species_SNVs_merged2 <- genome_species_SNVs_merged %>%
  left_join(all_samples.genome_info %>%
              dplyr::select(genome, length, breadth_minCov, sample2, coverage),
            by = c("sample2","genome")) %>% as.data.frame()
str(genome_species_SNVs_merged2)

# Calculate species-level SNV rate using callable-length weighting:
# total SNV counts divided by total effective callable genome length.
head(genome_species_SNVs_merged2)

# Calculate SNV rate for each sample × species pair.
library(dplyr)
sample_species_SNVs_rate  <- genome_species_SNVs_merged2 %>%
  
  # Calculate effective callable genome length.
  mutate(
    eff_len = length * breadth_minCov
  ) %>%
  
  # Aggregate by sample × species.
  group_by(sample2, Species) %>%
  summarise(
    
    # Sum SNV counts by mutation type.
    S_SNV_count = sum(S_SNV_richness, na.rm = TRUE),
    N_SNV_count = sum(N_SNV_richness, na.rm = TRUE),
    I_SNV_count = sum(I_SNV_richness, na.rm = TRUE),
    M_SNV_count = sum(M_SNV_richness, na.rm = TRUE),
    total_SNV_count = sum(total_SNV_richness, na.rm = TRUE),
    total_species_coverage = sum(coverage * eff_len, na.rm = TRUE),
    
    # Sum effective callable genome length.
    total_eff_len = sum(eff_len, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  # Calculate SNV rate for each mutation type and total SNV rate.
  mutate(
    S_SNV_rate     = S_SNV_count / total_eff_len,
    N_SNV_rate     = N_SNV_count / total_eff_len,
    I_SNV_rate     = I_SNV_count / total_eff_len,
    M_SNV_rate     = M_SNV_count / total_eff_len,
    total_SNV_rate = total_SNV_count / total_eff_len,
    species_coverage = total_species_coverage / total_eff_len
  )

View(sample_species_SNVs_rate)
table(!is.na(sample_species_SNVs_rate$total_SNV_rate))
head(sample_species_SNVs_rate)

# Calculate coverage-based relative abundance within each sample.
sample_species_SNVs_rate <- sample_species_SNVs_rate %>%
  group_by(sample2) %>%
  mutate(species_coverage_pct = species_coverage / sum(species_coverage)) %>%
  ungroup()


###
meta_final #
sample_species_SNVs_rate #
dim(sample_species_SNVs_rate) # 353261 * 16

# Merge species-level SNV table with metadata.
sample_species_SNVs_rate_metainfo <- merge(sample_species_SNVs_rate,
                                           meta_final,
                                           by="sample2",
                                           all.x = T)
dim(sample_species_SNVs_rate_metainfo) # 353261 * 32
View(sample_species_SNVs_rate_metainfo)
table(sample_species_SNVs_rate_metainfo$Time_new2)

### Data preprocessing
library(tidyverse)
library(mgcv)
library(patchwork)
library(viridis)

# Copy data.
df <- sample_species_SNVs_rate_metainfo
colnames(df)

# Set age order.
age_levels <- c(0, 0.5, 1, 6, 12, 24, 36)

df <- df %>%
  mutate(
    age = as.numeric(Time_new2),
    age_factor = factor(age, levels = age_levels),
    SNV_rate = total_SNV_rate,
    abundance = species_coverage_pct  # Coverage-based relative abundance.
  )

# Check zero values.
table(df$abundance==0) # FALSE -- 353261
table(df$total_SNV_rate==0) # FALSE -- 353261

# 1) Keep sample-level metadata; each sample corresponds to one age group.
sample_meta <- df %>%
  distinct(sample2, age, age_factor)

# 2) Extract all detected species.
all_species <- df %>%
  distinct(Species)

# 3) Construct a complete sample × species matrix.
# Undetected species in a given sample are assigned zero SNV rate and zero abundance.
df_complete <- sample_meta %>%
  crossing(all_species) %>%
  left_join(
    df %>%
      select(sample2, Species, SNV_rate, abundance),
    by = c("sample2", "Species")
  ) %>%
  mutate(
    SNV_rate = replace_na(SNV_rate, 0),
    abundance = replace_na(abundance, 0)
  )

# Check whether the complete matrix was generated correctly.
length(unique(df_complete$Species)) # 1127
test1 <- df_complete %>%
  subset(sample2=="ERR1600426")
length(unique(test1$Species)) # 1127
sum(test1$abundance)



## Calculate species prevalence across samples.
colnames(df_complete)
colnames(df)
dim(df_complete)

# Prevalence: proportion of samples in which each species has non-zero abundance.
length(unique(df_complete$sample2)) # 5886
prevalence_df <- df_complete %>%
  group_by(Species) %>%
  summarise(
    prevalence = mean(abundance > 0),
    .groups = "drop"
  )

## Calculate longitudinal consistency: whether each species was detected in at least three age groups.
timepoint_df <- df_complete %>%
  group_by(Species, age) %>%
  summarise(
    detected = any(abundance > 0),
    .groups = "drop"
  ) %>%
  group_by(Species) %>%
  summarise(
    n_timepoints = sum(detected),
    .groups = "drop"
  )


#### Identify species with significant GAM age effects.
# Calculate GAM results for each species.
colnames(df_complete)
length(unique(df_complete$sample2)) # 5886 
length(unique(df_complete$Species)) # 1127
str(df_complete)
table(df_complete$age)
range(df_complete$abundance)
range(df_complete$SNV_rate)

test <- df_complete %>% subset(abundance==1)
test
length(unique(test$Species)) # 5
length(unique(test$sample2)) # 16

# Check whether genome_info.tsv contains only genomes with detectable SNVs,
# and whether genomes without SNVs are absent from the table.
test2 <- df_complete %>% subset(sample2=="SRR7217609")
table(test2$abundance==0)
table(test2$abundance==1)
# Conclusion: some samples are dominated by only one detected species.

table(is.na(gam_results$age))
table(is.na(gam_results$SNV_rate))
table(is.na(gam_results$abundance))
range(gam_results$SNV_rate)
range(gam_results$abundance)


## GAM analysis
gam_results1 <- df_complete %>%
  group_by(Species) %>%
  group_modify(~{
    tryCatch({
      model <- gam(
        SNV_rate ~ s(age, k=6) + log10(abundance + 1e-6),
        # In this GAM, SNV_rate is modeled as a nonlinear smooth function of age.
        # The model does not assume a linear relationship; the smooth term allows the data to determine the shape of the age trajectory.
        # k is the basis dimension and controls the maximum complexity of the smooth function.
        # Smaller k values allow simpler curves, whereas larger k values allow more flexibility but may increase overfitting risk.
        # A GAM is appropriate here because SNV rate may change nonlinearly across early-life development.
        data = .x,
        method = "REML"
      )
      s_table <- summary(model)$s.table
      tibble(
        F_value = s_table[1, "F"],
        P_age = s_table[1, "p-value"],
        Adj_R2 = summary(model)$r.sq
      )
    }, error = function(e) {
      tibble(F_value = NA, P_age = NA, Adj_R2 = NA)
    })
  }) %>%
  ungroup() %>%
  mutate(
    P_label = case_when(
      P_age < 0.001 ~ "***",
      P_age < 0.01  ~ "**",
      P_age < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

gam_results1
table(gam_results1$P_age<0.05) # 432

sig_species <- gam_results1 %>%
  subset(P_age < 0.05) %>%
  pull(Species)



## Merge filtering information and select the top 30 species.
species_selection <- prevalence_df %>%
  left_join(timepoint_df, by = "Species") %>%
  left_join(dynamic_df, by = "Species") %>%
  filter(
    n_timepoints >= 3
  ) %>%
  arrange(desc(mean_abundance))

View(species_selection)
length(unique(species_selection$Species)) 

top30_species <- species_selection %>%
  subset(species_selection$Species %in% sig_species) %>%
  slice(1:30) %>%
  pull(Species)

top30_species


## Prepare data for the bubble plot.
plot_df <- df_complete %>%
  filter(Species %in% top30_species) %>%
  group_by(Species, age, age_factor) %>%
  summarise(mean_SNV_rate = mean(SNV_rate, na.rm = TRUE),
            mean_abundance = mean(abundance, na.rm = TRUE),
            .groups = "drop")

plot_df
range(plot_df$mean_abundance)
range(plot_df$mean_SNV_rate)


# Order species by mean abundance before clustering.
species_order <- species_selection %>%
  filter(Species %in% top30_species) %>%
  arrange(mean_abundance) %>%
  pull(Species)

plot_df$Species <- factor(plot_df$Species, levels = unique(species_order))

# Check duplicated species names.
species_order[duplicated(species_order)]


## Draw bubble plot.
library(RColorBrewer)
colnames(plot_df)
range(plot_df$mean_SNV_rate)



### 1. Prepare clustering matrix based on SNV rate and abundance.

# SNV-rate matrix in wide format.
snv_mat <- plot_df %>%
  select(Species, age, mean_SNV_rate) %>%
  pivot_wider(names_from = age, values_from = mean_SNV_rate, values_fill = 0) %>%
  as.data.frame()

# Relative-abundance matrix in wide format.
abun_mat <- plot_df %>%
  select(Species, age, mean_abundance) %>%
  pivot_wider(names_from = age, values_from = mean_abundance, values_fill = 0) %>%
  as.data.frame()

# Set Species as row names and remove the Species column.
rownames(snv_mat) <- snv_mat$Species
snv_mat$Species <- NULL

rownames(abun_mat) <- abun_mat$Species
abun_mat$Species <- NULL

# Row-scale each matrix by z-score to emphasize temporal patterns rather than absolute values.
snv_scaled <- t(scale(t(as.matrix(snv_mat))))
abun_scaled <- t(scale(t(as.matrix(abun_mat))))

# Replace NA values caused by constant rows with 0.
snv_scaled[is.na(snv_scaled)] <- 0
abun_scaled[is.na(abun_scaled)] <- 0

# Combine the two standardized matrices.
combined_mat <- cbind(snv_scaled, abun_scaled)

# Perform hierarchical clustering using Ward's minimum variance method.
hc_combined <- hclust(dist(combined_mat), method = "ward.D2")

# Extract clustered species order.
species_order_combined <- rownames(combined_mat)[hc_combined$order]


### 2. Set Species factor levels according to the clustered order.

plot_df_clustered <- plot_df %>%
  mutate(Species = factor(Species, levels = species_order_combined))


### 3. Draw the bubble plot.

p_bubble <- plot_df_clustered %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = age_factor, y = Species)) +
  geom_point(aes(size = mean_abundance,
                 color = mean_SNV_rate
  ),alpha = 0.9 ) +
  scale_size_continuous(name = "Relative abundance",
                        range = c(0, 10)
  ) +
  scale_color_gradient(low = "#d8daeb",
                       high = "#2d004b",
                       name = "SNV rate" ) +
  labs(x = "Age (months)",
       y = "Species (top 30)",
       title = "Age-associated dynamics of SNV rate across top 30 bacterial species") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Show plot.
print(p_bubble)


### Calculate the sign of correlation between SNV rate and age.
gam_results1 #

class(df_complete$age)
range(df_complete$SNV_rate)
range(df_complete$abundance)

cor_sign_df <- df_complete %>%
  subset(!is.na(age)) %>%
  filter(Species %in% top30_species) %>%
  group_by(Species) %>%
  summarise(rho = cor(SNV_rate, age, method = "spearman"),
            .groups = "drop") %>%
  mutate(cor_sign = ifelse(rho > 0, "+", ifelse(rho < 0, "-", "0")))

# Merge correlation direction into GAM results.
gam_results1 <- gam_results1 %>%
  left_join(cor_sign_df %>% dplyr::select(Species, cor_sign), by = "Species") 

# Check whether cor_sign was successfully added.
if(!"cor_sign" %in% names(gam_results1)) {
  gam_results1$cor_sign <- cor_sign_df$cor_sign[match(gam_results1$Species, cor_sign_df$Species)]
}


gam_results1$P_age2 <- ifelse(!is.na(gam_results1$P_age) & gam_results1$P_age < 2.2e-16,
                              2.2e-16,
                              formatC(signif(gam_results1$P_age, 2), format = "e", digits = 2))
table(gam_results1$P_age2==0)

gam_results1 <- gam_results1 %>%
  mutate(label = paste0(
    "P=", P_age2,
    ", R²=", round(Adj_R2, 2),
    ", ", cor_sign))


### Add GAM statistics.

### Ensure that the Species order in the GAM results matches the bubble plot.

# Convert Species to a factor with the clustered order.
gam_results <- gam_results1 %>%
  filter(Species %in% top30_species) %>%
  mutate(Species = factor(Species, levels = species_order_combined))

gam_results$Species

# Sort by factor level to align with the bubble plot.
gam_results <- gam_results %>% arrange(desc(Species))
View(gam_results)


### Right-side GAM statistics panel.

class(gam_results$P_age)

p_stats <- ggplot(gam_results,
                  aes(x = 1, y = Species)) +
  geom_text(aes(label = label), hjust = 0, size = 3) +
  xlim(1, 5) +
  theme_void() +
  labs(title = "GAM statistics") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


### Combine the bubble plot and the GAM statistics panel.
p_snvrate_relabun_top30species <- p_bubble + p_stats +
  plot_layout(widths = c(4, 1))

# Show plot.
print(p_snvrate_relabun_top30species)
# Save as 6.5 * 8.8





### Fig. 2f Bubble plots showing age-associated dynamics of SNV rate and observed N/S ratio for the top 30 genes ############

#### On the Linux server:


# Add gene_length, breadth_minCov and nucl_diversity columns from all_samples.gene_info.tsv
# to the filtered SNV table by matching the same gene and sample identifiers.

duckdb_cli

SET memory_limit='600GB';
SET threads=128;

COPY (
  SELECT
  snv.*,
  gi.gene_length,
  gi.breadth_minCov,
  gi.nucl_diversity
  FROM read_csv_auto(
    '/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/filtered_SNVs_coverage5_with_genome_genename_final2.tsv',
    delim = '\t',
    header = TRUE,
    sample_size = -1
  ) snv
  
  LEFT JOIN read_csv_auto(
    '/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/all_samples.gene_info.tsv',
    delim = '\t',
    header = TRUE,
    sample_size = -1
  ) gi
  
  ON snv.gene = gi.gene
  AND snv.sample2 = gi.sample2
  
)
TO '/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/filtered_SNVs_coverage5_with_genome_genename_final3.tsv'
WITH (HEADER, DELIMITER '\t');

# Check output file.
ls -alh /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/filtered_SNVs_coverage5_with_genome_genename_final3.tsv
# 91 G



# Calculate SNV rate and observed nonsynonymous-to-synonymous SNV ratio
# for each sample × gene_name pair.
duckdb_cli

COPY (
  WITH snv_filtered AS (
    SELECT *
      FROM read_csv_auto(
        '/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/filtered_SNVs_coverage5_with_genome_genename_final3.tsv',
        delim = '\t',
        header = TRUE,
        sample_size = -1
      )
    WHERE gene_name IS NOT NULL
    AND gene_name != ''
    AND gene IS NOT NULL
    AND gene != ''
  )
  
  SELECT
  gene_name,
  sample2,
  
  COUNT(*) AS SNV_count,
  
  SUM(CASE WHEN mutation_type = 'N' THEN 1 ELSE 0 END) AS SNV_N,
  SUM(CASE WHEN mutation_type = 'S' THEN 1 ELSE 0 END) AS SNV_S,
  
  CASE
  WHEN SUM(CASE WHEN gene_length IS NULL OR breadth_minCov IS NULL THEN 1 ELSE 0 END) > 0
  THEN NULL
  WHEN SUM(gene_length * breadth_minCov) <= 0
  THEN NULL
  ELSE COUNT(*) * 1.0 / (SUM(gene_length * breadth_minCov) + 1e-10)
  END AS SNV_rate,
  
  CASE
  WHEN SUM(CASE WHEN mutation_type = 'S' THEN 1 ELSE 0 END) = 0 THEN NULL
  ELSE SUM(CASE WHEN mutation_type = 'N' THEN 1 ELSE 0 END) * 1.0 /
    SUM(CASE WHEN mutation_type = 'S' THEN 1 ELSE 0 END)
  END AS NS_ratio,
  
  AVG(nucl_diversity) AS mean_pi
  
  FROM snv_filtered
  GROUP BY gene_name, sample2
)
TO '/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/gene_sample_SNVrate_NSratio.tsv'
WITH (HEADER, DELIMITER '\t');

# Check output file.
ls -alh /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/gene_sample_SNVrate_NSratio.tsv
# 1.4 G

less /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/gene_sample_SNVrate_NSratio.tsv




##### R code:

meta_final #
sample_species_SNVs_rate #
dim(sample_species_SNVs_rate) # 353261 * 16
gene_sample_SNVrate_NSratio.tsv #

# Read gene-level SNV summary table.
library("readr")
gene_sample_SNVrate_NSratio <- read_tsv("gene_sample_SNVrate_NSratio.tsv", show_col_types = FALSE)
View(gene_sample_SNVrate_NSratio)

colnames(gene_sample_SNVrate_NSratio)
table(!is.na(gene_sample_SNVrate_NSratio$mean_pi))# all T
table(!is.na(gene_sample_SNVrate_NSratio$NS_ratio))
table(!is.na(gene_sample_SNVrate_NSratio$SNV_rate))

# Set missing N/S ratio values to 0.
gene_sample_SNVrate_NSratio$NS_ratio <- ifelse(is.na(gene_sample_SNVrate_NSratio$NS_ratio)==TRUE,
                                               0,gene_sample_SNVrate_NSratio$NS_ratio)

# Check missing SNV-rate values.
test <- gene_sample_SNVrate_NSratio %>%
  subset(is.na(gene_sample_SNVrate_NSratio$SNV_rate)==TRUE)
dim(test)

# Check genes labeled as Unknown.
test <- gene_sample_SNVrate_NSratio %>%
  subset(gene_name=="Unknown")
dim(test)
View(test)



### Step 1: Build a clean GAM input table.
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mgcv)
library(broom)

gene_gam_df <- gene_sample_SNVrate_NSratio %>%
  left_join(
    meta_final %>%
      dplyr::select(sample2, age = Time_new2),
    by = "sample2"
  ) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(SNV_rate), !is.na(NS_ratio)) %>%
  filter(is.finite(SNV_rate), is.finite(NS_ratio)) %>%
  filter(SNV_rate > 0) %>%
  mutate(
    age = as.numeric(age),
    age_factor = factor(
      age,
      levels = c(0, 0.5, 1, 6, 12, 24, 36),
      labels = c("0", "0.5", "1", "6", "12", "24", "36")
    ),
    log_NS_ratio = log10(NS_ratio + 1e-6)
  )

View(gene_gam_df)


### Step 2: Fit GAM for each gene.
gam_gene_results <- gene_gam_df %>%
  group_by(gene_name) %>%
  group_modify(~{
    
    df <- .x
    
    # Set basis dimension according to the number of observed age groups, with a maximum of 6.
    k_use <- min(6, length(unique(df$age)) - 1)
    
    tryCatch({
      model <- gam(
        SNV_rate ~ s(age, k = k_use) + log_NS_ratio,
        data = df,
        method = "REML"
      )
      
      s_table <- summary(model)$s.table
      
      tibble(
        F_value = s_table[1, "F"],
        P_age = s_table[1, "p-value"],
        Adj_R2 = summary(model)$r.sq,
        n = nrow(df)
      )
      
    }, error = function(e) {
      tibble(
        F_value = NA_real_,
        P_age = NA_real_,
        Adj_R2 = NA_real_,
        n = nrow(df)
      )
    })
  }) %>%
  ungroup() %>%
  mutate(
    P_label = case_when(
      P_age < 0.001 ~ "***",
      P_age < 0.01  ~ "**",
      P_age < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

View(gam_gene_results)


### Step 3: Select significant genes.
sig_genes <- gam_gene_results %>%
  filter(!is.na(P_age)) %>%
  filter(P_age < 0.05) %>%         # Significant age effect.
  filter(Adj_R2 > 0.05) %>%        # Sufficient model explanatory power.
  filter(n > 100)                  # Sufficient sample size.

View(sig_genes)


### Step 4: Summarize gene × age table for bubble plot.

# Build a clean gene-age summary table.
gene_age_df_clean <- gene_gam_df %>%
  group_by(gene_name, age, age_factor) %>%
  summarise(
    n_samples = n_distinct(sample2),
    mean_SNV_rate = mean(SNV_rate, na.rm = TRUE),
    mean_NS_ratio = mean(NS_ratio, na.rm = TRUE),
    median_SNV_rate = median(SNV_rate, na.rm = TRUE),
    median_NS_ratio = median(NS_ratio, na.rm = TRUE),
    .groups = "drop"
  )


# Select the top 30 genes by mean SNV rate among significant genes.
top30_genes <- gene_age_df_clean %>%
  filter(gene_name %in% sig_genes$gene_name) %>%
  arrange(desc(mean_SNV_rate)) %>%
  slice_head(n = 30) %>%
  pull(gene_name)


# Extract plotting data for the top 30 genes.
plot_gene_df <- gene_age_df_clean %>%
  filter(gene_name %in% top30_genes)



### Step 5: Cluster genes using SNV-rate and N/S-ratio trajectories.
snv_gene_mat <- plot_gene_df %>%
  dplyr::select(gene_name, age, mean_SNV_rate) %>%
  pivot_wider(
    names_from = age,
    values_from = mean_SNV_rate,
    values_fill = 0
  ) %>%
  as.data.frame()

ns_gene_mat <- plot_gene_df %>%
  dplyr::select(gene_name, age, mean_NS_ratio) %>%
  pivot_wider(
    names_from = age,
    values_from = mean_NS_ratio,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(snv_gene_mat) <- snv_gene_mat$gene_name
snv_gene_mat$gene_name <- NULL

rownames(ns_gene_mat) <- ns_gene_mat$gene_name
ns_gene_mat$gene_name <- NULL

snv_gene_scaled <- t(scale(t(as.matrix(snv_gene_mat))))
ns_gene_scaled <- t(scale(t(as.matrix(ns_gene_mat))))

snv_gene_scaled[is.na(snv_gene_scaled)] <- 0
ns_gene_scaled[is.na(ns_gene_scaled)] <- 0

combined_gene_mat <- cbind(snv_gene_scaled, ns_gene_scaled)

hc_gene <- hclust(dist(combined_gene_mat), method = "ward.D2")
gene_order_combined <- rownames(combined_gene_mat)[hc_gene$order]


### Step 6: Left-side bubble plot.
colnames(plot_gene_df)
plot_gene_df_clustered <- plot_gene_df %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order_combined),
    
    # Cap extreme N/S ratios to avoid oversized bubbles in visualization.
    NS_ratio_plot = ifelse(mean_NS_ratio > 5, 5, mean_NS_ratio)
  )

## Left-side bubble plot.
range(plot_gene_df_clustered$NS_ratio_plot)
range(plot_gene_df_clustered$mean_SNV_rate)

dim(plot_gene_df_clustered)
length(top30_genes)

test <- plot_gene_df_clustered %>%
  subset(gene_name %in% top30_genes)
dim(test)
range(test$NS_ratio_plot)
range(test$mean_SNV_rate)
View(test)

p_gene_bubble <- plot_gene_df_clustered %>%
  filter(gene_name %in% top30_genes) %>%
  ggplot(aes(x = age_factor, y = gene_name)
  ) +
  geom_point(
    aes(
      size = NS_ratio_plot,
      color = mean_SNV_rate
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    name = "N/S ratio",
    range = c(0, 10)
  ) +
  scale_color_gradient(
    low = "#d8daeb",
    high = "#2d004b",
    limits = c(0, 0.26),
    name = "SNV rate"
  ) +
  labs(
    x = "Age (months)",
    y = "Gene",
    title = "Age-associated dynamics of gene-level SNV rate"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_gene_bubble


### Step 7: Add Spearman direction and draw the right-side GAM statistics panel.

## Calculate the Spearman correlation direction between gene-level SNV rate and age.
## GAMs model nonlinear associations and do not directly provide an overall positive or negative direction;
## therefore, Spearman rho is used here as an auxiliary trend-direction label.

cor_sign_df <- gene_gam_df %>%
  filter(gene_name %in% top30_genes) %>%
  group_by(gene_name) %>%
  summarise(
    rho = suppressWarnings(
      cor(SNV_rate, age, method = "spearman", use = "complete.obs")
    ),
    cor_sign = case_when(
      is.na(rho) ~ "NA",
      rho > 0 ~ "+",
      rho < 0 ~ "-",
      TRUE ~ "0"
    ),
    .groups = "drop"
  )

View(cor_sign_df)


## Prepare GAM results and trend-direction labels.
gam_results_plot <- gam_gene_results %>%
  filter(gene_name %in% top30_genes) %>%
  left_join(cor_sign_df, by = "gene_name") %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order_combined),
    
    P_age2 = case_when(
      is.na(P_age) ~ "NA",
      P_age < 2.2e-16 ~ "2.2e-16",
      TRUE ~ formatC(signif(P_age, 2), format = "e", digits = 2)
    ),
    
    label = paste0(
      "P=", P_age2,
      P_label,
      ", R²=", round(Adj_R2, 2),
      ", ", cor_sign
    )
  ) %>%
  arrange(gene_name)

View(gam_results_plot)


### Step 8: Right-side GAM statistics and trend-direction panel.
View(gam_gene_results)

p_gene_stats <- ggplot(
  gam_results_plot,
  aes(x = 1, y = gene_name)
) +
  geom_text(
    aes(label = label),
    hjust = 0,
    size = 3
  ) +
  xlim(1, 6.5) +
  theme_void() +
  labs(title = "GAM statistics") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_gene_stats)


### Step 9: Combine bubble plot and GAM statistics panel.
p_gene_GAM_SNPrate_top30 <- p_gene_bubble + p_gene_stats +
  plot_layout(widths = c(4, 1.5))

print(p_gene_GAM_SNPrate_top30)
# Save as 6.5 * 8.8

