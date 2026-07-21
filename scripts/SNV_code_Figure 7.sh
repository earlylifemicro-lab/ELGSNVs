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
library(patchwork)

### load phyloseq data
load("D:/3_Projects/2_Children_jaundice/3_R_analysis/Result7_R.RData")
save.image(file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Result7_R.RData")


############ Result: PreNEC-enriched bacterial SNV landscape and structural priorization of candidate nonsynonymous variants #####
#### metadata information ####

### prepare data 
meta_final <- read_excel("meta_final.xlsx")
meta_final #

# NEC_meta_infor.xlsx #
NEC_meta_infor <- read_excel("NEC_meta_infor.xlsx")
colnames(NEC_meta_infor)
dim(NEC_meta_infor)
table(NEC_meta_infor$Study)

# Check how many samples are included in the global cohort
test <- NEC_meta_infor %>%
  subset(Study%in% c("BrooksB_2017","RahmanSF_2018", "RavehSadkaT_2015", "RavehSadkaT_2016"))
dim(test) # 883 samples are included in this cohort, consistent with Chunlin's data

table(NEC_meta_infor$sample2 %in% meta_final$sample2)# T-692, F-3488
table(meta_final$sample2 %in% NEC_meta_infor$sample2)# T-692, F-5063

# Subset infants required for this analysis
NEC_meta_infor_sub <- NEC_meta_infor %>%
  subset(sample2 %in% meta_final$sample2)
dim(NEC_meta_infor_sub) # 692 samples
# 883 - 692 = 191 samples are missing, probably because inStrain profile failed for these samples

length(unique(NEC_meta_infor_sub$SubjectID)) # 94 
table(NEC_meta_infor_sub$Country) # 692
table(NEC_meta_infor_sub$Delivery_mode) #### CS-499, vaginal-193
table(NEC_meta_infor_sub$Sex) # female-326, male-366
table(NEC_meta_infor_sub$Feeding_pattern) # breast-72, combined-187, formula-78
table(NEC_meta_infor_sub$Infant_antibiotics) # 692
table(is.na(NEC_meta_infor_sub$Infant_antibiotics)) # yes-692
table(is.na(NEC_meta_infor_sub$PreNEC)) #0
table(NEC_meta_infor_sub$PreNEC)# noNEC-565,onsetNEC-7,postNEC-38,preNEC-82
table(NEC_meta_infor_sub$PreNEC, NEC_meta_infor_sub$Delivery_mode)
#           CS vaginal
# noNEC    416     149
# onsetNEC   1       6
# postNEC   25      13
# preNEC    57      25


range(NEC_meta_infor_sub$NEC_day_onset)
range(NEC_meta_infor_sub$Birth_weight_gram) # 550~2410
range(NEC_meta_infor_sub$Gestational_age_weeks) # 25~32
range(NEC_meta_infor_sub$DOL) # 0-82

### merge data
meta_final #
colnames(meta_final)

meta_sub <- meta_final %>%
  dplyr::select(sample2,DMM_cluster, Time_new2, Continent,SNV_rate_mean,pi_sample)

NEC_meta_infor_sub2 <- NEC_meta_infor_sub %>%
  left_join(meta_sub,by=c("sample2"))
NEC_meta_infor_sub2 #
dim(NEC_meta_infor_sub2 ) # 692 samples


# Define NEC infants
table(NEC_meta_infor_sub2$PreNEC)
NEC_meta_infor_sub2$NEC_children_Shuang <- ifelse(NEC_meta_infor_sub2$PreNEC=="noNEC","noNEC",
                                                  ifelse(NEC_meta_infor_sub2$PreNEC=="preNEC","preNEC",
                                                         "NEC"))
NEC_meta_infor_sub2$NEC_children_Shuang <- ifelse(NEC_meta_infor_sub2$PreNEC=="noNEC","noNEC","NEC")
table(NEC_meta_infor_sub2$NEC_children_Shuang) # NEC-127, noNEC-565
# Define whether NEC occurred

table(NEC_meta_infor_sub2$PreNEC, NEC_meta_infor_sub2$Time_new2)






### a. Baseline comparison of clinical variables between NEC and no NEC children ############
NEC_meta_infor_sub2 #
class(NEC_meta_infor_sub2)
colnames(NEC_meta_infor_sub2)

##### a1. Baseline comparison of continuous clinical variables ############
### Recheck birth weight and gestational age at the infant level
library(dplyr)
library(tidyr)

table(NEC_meta_infor_sub2$Infant_antibiotics)

child_df <- NEC_meta_infor_sub2 %>%
  mutate(NEC_group = ifelse(NEC_children_Shuang == "NEC", "NEC", "noNEC"),
         NEC_group = factor(NEC_group, levels = c("noNEC", "NEC"))) %>%
  group_by(SubjectID) %>%
  summarise(NEC_group = ifelse(any(NEC_group == "NEC", na.rm = TRUE), "NEC", "noNEC"),
            Birth_weight_gram = first(na.omit(Birth_weight_gram)),
            Gestational_age_weeks = first(na.omit(Gestational_age_weeks)),
            Sex = first(na.omit(Sex)),
            Delivery_mode = first(na.omit(Delivery_mode)),
            Country = first(na.omit(Country)),
            Study = first(na.omit(Study)),
            Infant_antibiotics = first(na.omit(Infant_antibiotics)),
            n_samples = n(),
            median_DOL = median(DOL, na.rm = TRUE),
            min_DOL = min(DOL, na.rm = TRUE),
            max_DOL = max(DOL, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(NEC_group = factor(NEC_group, levels = c("noNEC", "NEC")))

table(child_df$NEC_group)# NEC-17， noNEC-77
summary(child_df$n_samples)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   4.250   6.000   7.362   8.750  21.000



# First check whether fixed variables are consistent within the same infant
check_fixed_vars <- NEC_meta_infor_sub2 %>%
  group_by(SubjectID) %>%
  summarise(
    n_birth_weight = n_distinct(Birth_weight_gram, na.rm = TRUE),
    n_gest_age = n_distinct(Gestational_age_weeks, na.rm = TRUE),
    n_sex = n_distinct(Sex, na.rm = TRUE),
    n_delivery = n_distinct(Delivery_mode, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    n_birth_weight > 1 |
      n_gest_age > 1 |
      n_sex > 1 |
      n_delivery > 1
  )

check_fixed_vars
# The result shows that check_fixed_vars is empty, indicating that these variables are stable at the infant level.


### Infant-level comparison of continuous variables
calc_smd <- function(x, g) {
  g <- factor(g)
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE)
  s2 <- sd(x2, na.rm = TRUE)
  
  abs(m2 - m1) / sqrt((s1^2 + s2^2) / 2)
}

compare_child_cont <- function(df, var) {
  df2 <- df %>%
    filter(!is.na(.data[[var]]), !is.na(NEC_group))
  
  wt <- wilcox.test(
    as.formula(paste(var, "~ NEC_group")),
    data = df2
  )
  
  smd <- calc_smd(df2[[var]], df2$NEC_group)
  
  df2 %>%
    group_by(NEC_group) %>%
    summarise(
      n = n(),
      median = median(.data[[var]], na.rm = TRUE),
      Q1 = quantile(.data[[var]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[var]], 0.75, na.rm = TRUE),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd = sd(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      variable = var,
      p_value = wt$p.value,
      SMD = smd
    )
}

child_continuous_balance <- bind_rows(
  compare_child_cont(child_df, "Birth_weight_gram"),
  compare_child_cont(child_df, "Gestational_age_weeks"),
  compare_child_cont(child_df, "median_DOL")
)
child_continuous_balance

## Format as a table
child_continuous_table <- child_continuous_balance %>%
  mutate(
    value = paste0(
      round(median, 1),
      " [",
      round(Q1, 1),
      ", ",
      round(Q3, 1),
      "]"
    )
  ) %>%
  select(variable, NEC_group, value, p_value, SMD) %>%
  pivot_wider(
    names_from = NEC_group,
    values_from = value
  ) %>%
  distinct(variable, .keep_all = TRUE) %>%
  mutate(
    p_value = signif(p_value, 3),
    SMD = signif(SMD, 3)
  )
child_continuous_table
# Save table




##### a2. Baseline comparison of categorical clinical variables ############
### Infant-level comparison of categorical variables
compare_child_cat <- function(df, var) {
  df2 <- df %>%
    filter(!is.na(.data[[var]]), !is.na(NEC_group)) %>%
    droplevels()
  
  tab <- table(df2$NEC_group, df2[[var]])
  
  chi <- suppressWarnings(chisq.test(tab))
  
  if (any(chi$expected < 5)) {
    pval <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
    test_name <- "Fisher exact test, simulated"
  } else {
    pval <- chi$p.value
    test_name <- "Chi-square test"
  }
  
  data.frame(
    variable = var,
    test = test_name,
    p_value = pval,
    n = nrow(df2)
  )
}

child_categorical_balance <- bind_rows(
  compare_child_cat(child_df, "Sex"),
  compare_child_cat(child_df, "Delivery_mode"),
  compare_child_cat(child_df, "Study")#,
  #compare_child_cat(child_df, "Country"), # all are US
  #compare_child_cat(child_df, "Infant_antibiotics") # all are antibiotics---yes
)
child_categorical_balance #
# Save table

# Validation
table(child_df$Infant_antibiotics)#
table(is.na(child_df$Infant_antibiotics))#
table(child_df$Infant_antibiotics, child_df$NEC_group)

table(child_df$Country)#
table(is.na(child_df$Country))#
table(child_df$Country, child_df$NEC_group)

table(child_df$Sex)#
table(child_df$Delivery_mode)#
table(child_df$Study)#

table(NEC_meta_infor_sub2$Country)#
table(is.na(NEC_meta_infor_sub2$Country))#

table(NEC_meta_infor_sub2$Study)#
table(is.na(NEC_meta_infor_sub2$Study))#
table(NEC_meta_infor_sub2$NEC_children_Shuang, NEC_meta_infor_sub2$Study)#

table(child_df$Study)#
table(is.na(child_df$Study))#
chisq.test(table(child_df$Sex, child_df$NEC_group))
chisq.test(table(child_df$Delivery_mode, child_df$NEC_group))# chisq.test warning
chisq.test(table(child_df$Study, child_df$NEC_group))# chisq.test warning





##### Supplementary Fig. S8a Distribution of retained longitudinal samples per infant by NEC status #####

##### Distribution of retained longitudinal samples per infant
NEC_meta_infor_sub2 #
class(NEC_meta_infor_sub2)
colnames(NEC_meta_infor_sub2)


##### Number of samples collected per infant: faceted by NEC/noNEC and ordered from low to high
library(dplyr)
library(ggplot2)
library(forcats)

### 1. Build sample-level NEC stage
child_sample_count_stage <- NEC_meta_infor_sub2 %>%
  filter(!is.na(SubjectID),!is.na(sample2),!is.na(NEC_children_Shuang)) %>%
  mutate(NEC_group = ifelse(NEC_children_Shuang == "NEC", "NEC", "noNEC"),
         NEC_group = factor(NEC_group, levels = c("noNEC", "NEC")),
         
         # sample-level stage
         NEC_stage = case_when(
           NEC_group == "noNEC" ~ "noNEC",
           PreNEC %in% c("preNEC") ~ "preNEC",
           PreNEC %in% c("onsetNEC", "postNEC") ~ "postNEC",
           TRUE ~ NA_character_),
         NEC_stage = factor(NEC_stage,levels = c("noNEC",  "postNEC","preNEC"))) %>%
  filter(!is.na(NEC_stage)) %>%
  group_by(SubjectID, NEC_group, NEC_stage) %>%
  summarise(
    n_samples_stage = n_distinct(sample2),
    .groups = "drop")

### 2. Calculate total sample number per infant for ordering and annotation
child_sample_count_total <- child_sample_count_stage %>%
  group_by(SubjectID, NEC_group) %>%
  summarise(n_samples = sum(n_samples_stage),
            .groups = "drop") %>%
  group_by(NEC_group) %>%
  arrange(n_samples, SubjectID, .by_group = TRUE) %>%
  mutate(SubjectID_ordered = factor(
    paste(NEC_group, SubjectID, sep = "___"),
    levels = paste(NEC_group, SubjectID, sep = "___"))) %>%
  ungroup()

### 3. Merge ordering information
child_sample_count_stage_plot <- child_sample_count_stage %>%
  left_join(child_sample_count_total %>%
              select(SubjectID, NEC_group, n_samples, SubjectID_ordered),
            by = c("SubjectID", "NEC_group"))

### 4. Plot: stacked barplot
table(child_sample_count_stage_plot$SubjectID_ordered)
p_child_sample_count_stage <- ggplot(
  child_sample_count_stage_plot,
  aes(x = SubjectID_ordered,
      y = n_samples_stage,
      fill = NEC_stage )) +
  geom_col(width = 0.82,color = "white",linewidth = 0.08,alpha = 0.6) +
  
  # Annotate total sample number above each bar
  geom_text(data = child_sample_count_total,
            aes(x = SubjectID_ordered,y = n_samples,label = n_samples),
            inherit.aes = FALSE,vjust = -0.35,size = 2.6,color = "grey30") +
  scale_fill_manual(
    values = c("noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d"),
    
    labels = c("noNEC","postNEC","preNEC")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::pretty_breaks(n = 6)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90",linewidth = 0.4),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text( color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5,size = 13),
        plot.margin = margin(t = 8, r = 10, b = 8, l = 10)) +
  
  labs(x = "Infants",# Infants ordered by sample count
       y = "Number of retained fecal samples per infant",
       title = "Distribution of retained longitudinal samples per infant")
p_child_sample_count_stage
# save as 5 * 3.3






##### Fig. 7a Longitudinal fecal sampling across infants #####

##### Longitudinal fecal sampling across infants
NEC_meta_infor_sub2 #
class(NEC_meta_infor_sub2)
colnames(NEC_meta_infor_sub2)


library(dplyr)
library(ggplot2)
library(forcats)

table(NEC_meta_infor_sub2$PreNEC)
plot_subject_dol <- NEC_meta_infor_sub2 %>%
  filter(!is.na(SubjectID),
         !is.na(DOL),
         !is.na(PreNEC2)) %>%
  
  mutate(PreNEC = factor(PreNEC2,levels = c("noNEC","preNEC","postNEC"))) %>%
  group_by(SubjectID) %>%
  mutate(max_DOL = max(DOL, na.rm = TRUE),n_samples = n()) %>%
  ungroup() %>%
  mutate(SubjectID_ordered = fct_reorder(SubjectID,max_DOL))

p_subject_dol <- ggplot(
  plot_subject_dol,aes( x = DOL,y = SubjectID_ordered,
                        group = SubjectID,color = PreNEC)) +
  geom_line(linewidth = 0.4, alpha = 0.5,color = "grey70") +
  geom_point(size = 0.8,
             alpha = 0.95) +
  #scale_color_manual(values = c("noNEC" = "grey60",
  #                              "preNEC" = "#71E50C",#71E50C#6FDB0F
  #                              "postNEC" = "green")) +
  
  scale_color_manual(values = #c#("noNEC" = "grey60","preNEC" = "#fb9a99","postNEC" = "#e31a1c")
                       
                       c("noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Days of life (DOL)",
       y = "Infants ordered by latest sampling day",
       title = "Longitudinal fecal sampling across infants")

p_subject_dol +
  facet_grid(NEC_children_Shuang ~ .,
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text.y.left = element_text(angle = 90)) 
# save as 4 * 6, p_subject_dol





##### Supplementary Fig. S8b Comparison of microbial SNV rate across noNEC, preNEC and postNEC #####

###  NEC and no NEC children
NEC_meta_infor_sub2 #
class(NEC_meta_infor_sub2)
colnames(NEC_meta_infor_sub2)

table(NEC_meta_infor_sub2$NEC_children_Shuang)
table(NEC_meta_infor_sub2$DOL,NEC_meta_infor_sub2$PreNEC)
table(NEC_meta_infor_sub2$DMM_cluster,NEC_meta_infor_sub2$DOL,NEC_meta_infor_sub2$PreNEC)
fisher.test(table(NEC_meta_infor_sub2$DMM_cluster,NEC_meta_infor_sub2$NEC_children_Shuang))
NEC_meta_infor_sub2$NEC_children_Shuang <- factor(NEC_meta_infor_sub2$NEC_children_Shuang,
                                                  levels = c("noNEC","NEC"))

## noNEC vs. preNEC vs. postNEC
colnames(NEC_meta_infor_sub2)
table(NEC_meta_infor_sub2$PreNEC)

library(ggplot2)
library(ggpubr)
library(dplyr)

colnames(NEC_meta_infor_sub2)
NEC_meta_infor_sub2$PreNEC2 <- ifelse(NEC_meta_infor_sub2$PreNEC %in% c("noNEC", "preNEC"),
                                      NEC_meta_infor_sub2$PreNEC,"postNEC")

NEC_meta_infor_sub2$PreNEC2 <- factor(NEC_meta_infor_sub2$PreNEC2,
                                      levels = c("noNEC", "preNEC", "postNEC"))

my_comparisons <- combn(levels(NEC_meta_infor_sub2$PreNEC2),2,simplify = FALSE)

## SNV rate
p_PreNEC_SNVrate <- NEC_meta_infor_sub2 %>%
  filter(!is.na(PreNEC2),!is.na(SNV_rate_mean),SNV_rate_mean > 0 ) %>%
  ggplot(aes(x = PreNEC2,y = SNV_rate_mean,
             fill = PreNEC2,color = PreNEC2)) +
  geom_jitter(width = 0.12,alpha = 0.45,size = 1.2) +
  geom_violin(trim = FALSE,alpha = 0.45,linewidth = 0.4) +
  geom_boxplot( width = 0.14,outlier.shape = NA,alpha = 0.45,color = "black") +
  
  scale_y_log10() +
  scale_fill_manual(values = c( "noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d")) +
  scale_color_manual(values = c("noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 3,
                     tip.length = 0.01) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text( color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text( hjust = 0.5)) +
  labs(x = "",
       y = "Mean SNV rate",
       #title = "Microbial SNV rate across NEC stages"
  )
p_PreNEC_SNVrate
# save as 3 * 3.2





##### Supplementary Fig. S8c Comparison of nucleotide diversity across noNEC, preNEC and postNEC #####

## pi_sample
range(NEC_meta_infor_sub2$pi_sample)
p_PreNEC_pi <- NEC_meta_infor_sub2 %>%
  filter(!is.na(PreNEC2),!is.na(pi_sample)) %>%
  ggplot(aes(x = PreNEC2,y = pi_sample,
             fill = PreNEC2,color = PreNEC2)) +
  geom_jitter(width = 0.12,alpha = 0.45,size = 1.2) +
  geom_violin(trim = FALSE,alpha = 0.45,linewidth = 0.4) +
  geom_boxplot( width = 0.14,outlier.shape = NA,alpha = 0.45,color = "black") +
  
  scale_y_log10() +
  scale_fill_manual(values = c(#"noNEC" = "grey60","preNEC" = "#71E50C","postNEC" = "green"
    #"noNEC" = "grey60","preNEC" = "#fdbf6f","postNEC" = "#ff7f00"
    "noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d"
  )) +
  scale_color_manual(values = c(#"noNEC" = "grey60","preNEC" = "#71E50C","postNEC" = "green"
    #"noNEC" = "grey60","preNEC" = "#fb9a99","postNEC" = "#e31a1c"
    "noNEC" = "grey70","preNEC" = "#35978f","postNEC" = "#bf812d"
  )) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 3,
                     tip.length = 0.01) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text( color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "",
       y = "Nucleotide diversity",
       #title = "Microbial nucleotide diversity across NEC stages"
  )
p_PreNEC_pi
# save as 3 * 3.2





##### Supplementary Fig. S8d Shared mutated genomes among noNEC, preNEC and postNEC stages #####

### Genome difference caused by NEC
### data preparation
all_samples.genome_info.tsv #

library(tidyverse)

all_samples.genome_info <- read_tsv("D:/3_Projects/2_Children_jaundice/3_R_analysis/all_samples.genome_info.tsv",
                                    show_col_types = FALSE)


### Sample-level SNV rate 
# Total number of SNV sites in each sample across all scaffolds
# Sum SNV counts across all scaffolds for each sample
colnames(all_samples.genome_info)
table(all_samples.genome_info$length==0)
table(all_samples.genome_info$breadth_minCov==0)# contains zero values
table(all_samples.genome_info$SNV_count==0)# contains zero values
snv_sample2 <- all_samples.genome_info %>%
  mutate(SNV_rate_genome = SNV_count/(length*breadth_minCov+ 1e-10))


### Create Genome2 column for bin-level mapping:
library(data.table)
# Assume snv_sample2 has already been loaded
setDT(snv_sample2)
snv_sample2[, Genome1 := sub("_[0-9]+$", "", genome)]
snv_sample2[, Genome2 := gsub("_k[0-9]+$", "", Genome1)]
snv_sample2[, Genome3 := sub("(_genomic).*", "_genomic", Genome2)]
snv_sample2[, Genome4:= sub("(_\\.\\d+.*)$", "", Genome3)]
# Assume snv_sample2 is already a data.table and contains the Genome4 column
snv_sample2[, Genome5:=
              fifelse(grepl("_genomic", Genome4),
                      sub("(_genomic).*", "\\1", Genome4),  # keep up to _genomic
                      sub("_ERS.*$", "", Genome4)           # remove everything from _ERS onward
              )
]


### Read genome_species_info
genome_species_info <- read_excel("D:/3_Projects/2_Children_jaundice/3_R_analysis/genome_species_info.xlsx")
genome_species_info2 <- data.frame(genome_species_info )%>%
  dplyr::select(Genome, Taxonomy.lineage..GTDB.)
colnames(genome_species_info2)
# Check the first few rows

library(tidyverse)
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
    Taxonomy.lineage..GTDB.=Taxonomy.lineage..GTDB.,
    Domain  = strip_prefix(Domain),
    Phylum  = strip_prefix(Phylum),
    Class   = strip_prefix(Class),
    Order   = strip_prefix(Order),
    Family  = strip_prefix(Family),
    Genus   = na_if(strip_prefix(Genus), ""),
    Species = na_if(strip_prefix(Species), "")
  )


## Merge to the species level
table(genome_species_info3$Genome %in% snv_sample2$Genome5)
# FALSE-, TRUE-2170
table(snv_sample2$Genome5 %in% genome_species_info3$Genome)

# Ensure both tables are in data.table format
setDT(snv_sample2)
setDT(genome_species_info3)

# Merge using Genome5 and Genome columns
filtered_SNVs_coverage5_top5_species <- merge(
  snv_sample2,
  genome_species_info3,
  by.x = "Genome5",
  by.y = "Genome",
  all.x = FALSE,
  all.y = FALSE
)#%>% dplyr::select(-Genome1,-Genome2,-Genome3,-Genome4)

# Check the merged result
print(dim(filtered_SNVs_coverage5_top5_species))

# Handle NA values in Genus and Species columns
filtered_SNVs_coverage5_top5_species$Genus <- ifelse(is.na(filtered_SNVs_coverage5_top5_species$Genus),
                                                     paste0(filtered_SNVs_coverage5_top5_species$Family, "_Unknown"),
                                                     filtered_SNVs_coverage5_top5_species$Genus)
filtered_SNVs_coverage5_top5_species$Species <- ifelse(is.na(filtered_SNVs_coverage5_top5_species$Species),
                                                       paste0(filtered_SNVs_coverage5_top5_species$Genus, "_Unknown"),
                                                       filtered_SNVs_coverage5_top5_species$Species)

table(filtered_SNVs_coverage5_top5_species$Species)
table(!is.na(filtered_SNVs_coverage5_top5_species$Species))
table(filtered_SNVs_coverage5_top5_species$Species=="Unknown")
# Check which taxa are annotated as Unknown
test <- subset(filtered_SNVs_coverage5_top5_species, filtered_SNVs_coverage5_top5_species$Species=="Unknown")
dim(test) # 28263 * 46
table(test$Genus)
# Many species

table(test$Genome1,test$Genus)


#### Subset NEC samples
filtered_SNVs_coverage5_top5_species #
NEC_meta_infor_sub2 #

colnames(NEC_meta_infor_sub2)
dim(NEC_meta_infor_sub2)#  692
colnames(filtered_SNVs_coverage5_top5_species)

NEC_genome_SNVs_coverage5_top5_species <- filtered_SNVs_coverage5_top5_species %>%
  subset(sample2 %in% NEC_meta_infor_sub2$sample2)
dim(NEC_genome_SNVs_coverage5_top5_species)
length(unique(NEC_genome_SNVs_coverage5_top5_species$sample2)) # 692

## Merge PreNEC2, NEC_children_Shuang and PreNEC columns
library(dplyr)
NEC_genome_SNVs_coverage5_top5_species2 <-NEC_genome_SNVs_coverage5_top5_species %>%
  left_join(NEC_meta_infor_sub2 %>%
              dplyr::select(sample2,PreNEC2,NEC_children_Shuang,PreNEC,SubjectID),by = "sample2" )
dim(NEC_genome_SNVs_coverage5_top5_species2)

colnames(NEC_genome_SNVs_coverage5_top5_species2)

# check
table(NEC_genome_SNVs_coverage5_top5_species2$PreNEC2,useNA = "ifany")
table(NEC_genome_SNVs_coverage5_top5_species2$NEC_children_Shuang,useNA = "ifany")
table(NEC_genome_SNVs_coverage5_top5_species2$PreNEC,useNA = "ifany")
table(NEC_genome_SNVs_coverage5_top5_species2$SubjectID,useNA = "ifany")

test <- NEC_genome_SNVs_coverage5_top5_species2 %>%
  subset(genome %in%c("RavehSadkaT_2015_SRR1779110_bin.3_k141",
                      "ShaoY_2019_ERR3404912_bin.8_k141"))
table(test$Species)
# RavehSadkaT_2015_SRR1779110_bin.3_k141---Clostridioides difficile
# ShaoY_2019_ERR3404912_bin.8_k141---Clostridium_P perfringens







##### Bacterial SNV difference based on species-level 
NEC_genome_SNVs_coverage5_top5_species2 #

class(NEC_genome_SNVs_coverage5_top5_species2)
colnames(NEC_genome_SNVs_coverage5_top5_species2)



##### Genome-level Venn for noNEC / preNEC / postNEC groups

library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(tibble)

## Choose the genome ID column
## Use genome for the exact reference genome
## Use Genome5 for the curated genome label
#genome_id_col <- "Genome5"
genome_id_col <- "genome"

## Fix set order
sets_genome <- sets_genome[c("noNEC", "preNEC", "postNEC") ]

## Define colors for the three stages
stage_colors <- c("noNEC"   = "grey70","preNEC"  = "#35978f","postNEC" = "#bf812d")

## 1. Define mutated genomes
mut_genome_df <- NEC_genome_SNVs_coverage5_top5_species2 %>%
  mutate(genome_id = as.character(.data[[genome_id_col]]),
         PreNEC2 = factor(PreNEC2,levels = c("noNEC", "preNEC", "postNEC")) ) %>%
  filter(!is.na(PreNEC2),!is.na(genome_id),!is.na(SNV_count),SNV_count > 0) %>%
  distinct(PreNEC2, sample2, SubjectID, genome_id)

## Check
table(mut_genome_df$PreNEC2)
length(unique(mut_genome_df$genome_id))

## 2. Build genome sets for the three stages
sets_genome <- list(
  noNEC   = mut_genome_df %>% filter(PreNEC2 == "noNEC") %>% pull(genome_id) %>% unique(),
  preNEC  = mut_genome_df %>% filter(PreNEC2 == "preNEC") %>% pull(genome_id) %>% unique(),
  postNEC = mut_genome_df %>% filter(PreNEC2 == "postNEC") %>% pull(genome_id) %>% unique())

## Number of mutated genomes in each group
sapply(sets_genome, length)


## Build Venn object and extract plotting data
venn_object <- ggVennDiagram::Venn(sets_genome)
venn_data <- ggVennDiagram::process_data(venn_object)

## Mapping between set IDs and stage names
set_key <- tibble(id = as.character(seq_along(sets_genome)),
                  Stage = names(sets_genome))

## 1. Extract intersection regions
region_df <- ggVennDiagram::venn_regionedge(venn_data) %>%
  mutate(id = as.character(id),
         
         fill_group = case_when(id == "1" ~ "noNEC",
                                id == "2" ~ "preNEC",
                                id == "3" ~ "postNEC",
                                TRUE      ~ "Shared"))

## 2. Extract set borders
set_edge_df <- ggVennDiagram::venn_setedge(venn_data) %>%
  mutate(id = as.character(id)) %>%
  left_join(set_key,by = "id")

## 3. Extract set label positions
set_label_df <- ggVennDiagram::venn_setlabel(venn_data) %>%
  mutate(id = as.character(id) ) %>%
  left_join(set_key,by = "id")

## 4. Extract intersection-count label positions
region_label_df <- ggVennDiagram::venn_regionlabel(venn_data)

## 5. Redraw Venn plot
p_venn_genome <- ggplot() +
  
  ## Fill regions
  geom_polygon(data = region_df,
               aes(x = X,y = Y,group = id,fill = fill_group),color = NA,alpha = 0.75) +
  
  ## Borders of the three sets
  geom_path(data = set_edge_df,
            aes(x = X,y = Y,group = id,color = Stage),
            linewidth = 1.1,show.legend = FALSE) +
  
  ## Set labels
  geom_text(data = set_label_df,
            aes(x = X,y = Y,label = name,color = Stage),
            size = 5,fontface = "bold",show.legend = FALSE) +
  
  ## Genome counts in intersection regions
  geom_text(data = region_label_df,
            aes(x = X,y = Y,label = count),
            size = 4,color = "black") +
  
  ## Manually set region colors
  scale_fill_manual(values = c("noNEC"   = "grey70","preNEC"  = "#35978f",
                               "postNEC" = "#bf812d","Shared"  = "grey92"),
                    breaks = c("noNEC","preNEC","postNEC"),
                    name = "NEC stage") +
  
  ## Manually set border and set-label colors
  scale_color_manual(values = stage_colors) +
  coord_equal() +
  theme_void() +
  theme( legend.position = "right",
         plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Shared and stage-specific mutated genomes across NEC stages")
p_venn_genome
# save as 3 * 3


### Export genomes corresponding to each region
all_genomes <- sort(unique(unlist(sets_genome)))

venn_membership <- data.frame(genome_id = all_genomes) %>%
  mutate(noNEC   = genome_id %in% sets_genome$noNEC,
         preNEC  = genome_id %in% sets_genome$preNEC,
         postNEC = genome_id %in% sets_genome$postNEC,
         category = case_when(
           noNEC & !preNEC & !postNEC ~ "noNEC_only",
           !noNEC & preNEC & !postNEC ~ "preNEC_only",
           !noNEC & !preNEC & postNEC ~ "postNEC_only",
           noNEC & preNEC & !postNEC ~ "noNEC_preNEC_shared",
           noNEC & !preNEC & postNEC ~ "noNEC_postNEC_shared",
           !noNEC & preNEC & postNEC ~ "preNEC_postNEC_shared",
           noNEC & preNEC & postNEC ~ "shared_all_three",
           TRUE ~ "other"))
table(venn_membership$category)

#write.csv(venn_membership,"NEC_stage_mutated_genome_venn_membership.csv", row.names = FALSE)

## Export key regions separately
pre_post_shared <- venn_membership %>%
  filter(category %in% c("preNEC_postNEC_shared", "shared_all_three"))
colnames(pre_post_shared)

preNEC_only <- venn_membership %>%
  filter(category == "preNEC_only")
colnames(preNEC_only)

postNEC_only <- venn_membership %>%
  filter(category == "postNEC_only")
colnames(postNEC_only)

noNEC_only <- venn_membership %>%
  filter(category == "noNEC_only")
colnames(noNEC_only)

## Merge into one combined table
venn_membership_export <- bind_rows(
  pre_post_shared,
  preNEC_only,
  postNEC_only,
  noNEC_only)

## Save as an xlsx file
write_xlsx(
  venn_membership_export,
  path = "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 7_file/NEC_venn_membership_selected.xlsx")






##### Supplementary Fig. S8e Pairwise Jaccard similarity of mutated-genome sets across NEC stages #####


##### Jaccard index and hypergeometric enrichment test
### Test whether preNEC and postNEC share more mutated genomes
## Do not only examine overlap counts; use three metrics:
# Intersection count: number of shared elements;
# Jaccard index: intersection / union;
# Overlap coefficient: intersection / size of the smaller set.

pairwise_overlap <- function(setA, setB, nameA, nameB, universe) {
  
  inter_n <- length(intersect(setA, setB))
  union_n <- length(union(setA, setB))
  
  jaccard <- inter_n / union_n
  overlap_coef <- inter_n / min(length(setA), length(setB))
  
  ## hypergeometric enrichment test
  ## Test whether the overlap between A and B is greater than random expectation
  M <- length(universe) # total number of mutated genomes detected across the noNEC, preNEC and postNEC sets
  K <- length(setB) # number of mutated genomes in set B
  N <- length(setA) # number of mutated genomes in set A
  x <- inter_n # observed number of mutated genomes shared between set A and set B
  # Hypergeometric-test assumption: if N genomes are randomly sampled from M mutated genomes as set A,
  # and K genomes are randomly sampled as set B, what is the probability that the two sets share x or more genomes by chance?
  
  p_hyper <- phyper(q = x - 1,
                    m = K,
                    n = M - K,
                    k = N,
                    lower.tail = FALSE)
  
  data.frame(comparison = paste(nameA, nameB, sep = "_vs_"),
             n_A = length(setA),
             n_B = length(setB),
             intersection = inter_n,
             union = union_n,
             jaccard = jaccard,
             overlap_coefficient = overlap_coef,
             hypergeom_p = p_hyper)
}

universe_genomes <- sort(unique(unlist(sets_genome)))

overlap_summary <- bind_rows(
  pairwise_overlap(sets_genome$preNEC,sets_genome$postNEC,
                   "preNEC","postNEC",universe_genomes),
  pairwise_overlap(sets_genome$noNEC,sets_genome$preNEC,
                   "noNEC","preNEC",universe_genomes),
  pairwise_overlap(sets_genome$noNEC,sets_genome$postNEC,
                   "noNEC","postNEC",universe_genomes)) %>%
  mutate(hypergeom_FDR = p.adjust(hypergeom_p, method = "BH"))
overlap_summary

#write.csv(overlap_summary,"NEC_stage_pairwise_mutated_genome_overlap_summary.csv",row.names = FALSE )

overlap_summary %>% arrange(desc(jaccard))

## Add expected overlap and fold enrichment
M <- length(universe_genomes)

overlap_summary2 <- overlap_summary %>%
  mutate(expected_intersection = n_A * n_B / M,
         fold_enrichment = intersection / expected_intersection,
         hypergeom_p_label = case_when(
           hypergeom_p == 0 ~ "<2.2-16",
           TRUE ~ as.character(signif(hypergeom_p, 3))),
         hypergeom_FDR_label = case_when(hypergeom_FDR == 0 ~ "<2.2e-16",
                                         TRUE ~ as.character(signif(hypergeom_FDR, 3))) )
overlap_summary2



# plot
library(dplyr)
library(ggplot2)
library(scales)

## Format overlap_summary
overlap_plot_df <- overlap_summary %>%
  mutate(comparison_label = case_when(
    comparison == "preNEC_vs_postNEC" ~ "preNEC vs postNEC",
    comparison == "noNEC_vs_preNEC" ~ "noNEC vs preNEC",
    comparison == "noNEC_vs_postNEC" ~ "noNEC vs postNEC",
    TRUE ~ comparison),
    comparison_label = factor(comparison_label,
                              levels = c("noNEC vs preNEC","noNEC vs postNEC","preNEC vs postNEC")),
    jaccard_label = paste0("Jaccard = ", round(jaccard, 3),
                           "\nShared = ", intersection),
    FDR_label = case_when(hypergeom_FDR == 0 ~ "FDR < 1e-300",
                          TRUE ~ paste0("FDR = ", signif(hypergeom_FDR, 2))))

p_jaccard_genome <- ggplot(overlap_plot_df,
                           aes(x = comparison_label, y = jaccard, fill = comparison_label)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = jaccard_label),
            vjust = -0.25,size = 3.5) +
  scale_y_continuous(limits = c(0, max(overlap_plot_df$jaccard) * 1.25),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("noNEC vs preNEC" = "#9ecae1",
                               "noNEC vs postNEC" = "#c6dbef",
                               "preNEC vs postNEC" = "#756bb1")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank() ) +
  labs(y = "Jaccard index"#,
       #title = "Pairwise overlap of mutated genomes across NEC stages"
  )
p_jaccard_genome
# save as 3 * 3



##### Fig. 7b Overlap of mutated species across NEC stages #####

##### Species-level Venn for noNEC / preNEC / postNEC groups
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggVennDiagram)
library(data.table)

## 1. Prepare species-stage mutated table
genome_id_col <- "genome"

species_stage_df <- NEC_genome_SNVs_coverage5_top5_species2 %>%
  mutate(genome_id = as.character(.data[[genome_id_col]]),
         Species = as.character(Species),
         Genus = as.character(Genus),
         Family = as.character(Family),
         
         Species = ifelse(is.na(Species) | Species == "", "Unclassified", Species),
         Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
         Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family),
         
         NEC_stage = case_when(PreNEC2 == "noNEC" ~ "noNEC",
                               PreNEC2 == "preNEC" ~ "preNEC",
                               PreNEC2 %in% c("postNEC", "NEC", "onsetNEC") ~ "postNEC",
                               TRUE ~ NA_character_),
         NEC_stage = factor(NEC_stage, levels = c("noNEC", "preNEC", "postNEC"))) %>%
  filter(!is.na(NEC_stage),
         !is.na(Species),
         !is.na(sample2),
         !is.na(SubjectID),
         !is.na(genome_id),
         !is.na(SNV_count),
         SNV_count > 0) %>%
  distinct(NEC_stage,
           Species,
           Genus,
           Family,
           genome_id,
           sample2,
           SubjectID)

table(species_stage_df$NEC_stage)
length(unique(species_stage_df$Species))# 766 species



library(dplyr)
library(ggVennDiagram)
library(ggplot2)

## 1. Count species detection per NEC stage
species_stage_count <- species_stage_df %>%
  group_by(Species, Genus, Family, NEC_stage) %>%
  summarise(n_genomes = n_distinct(genome_id),
            n_samples = n_distinct(sample2),
            n_subjects = n_distinct(SubjectID),
            .groups = "drop")


## 2. Set minimum threshold for each group
##    Recommended to first use subject counts
min_subject_noNEC   <- 10
min_subject_preNEC  <- 8
min_subject_postNEC <- 8

qualified_species_stage <- species_stage_count %>%
  mutate(pass_stage_filter = case_when(
    NEC_stage == "noNEC"   & n_subjects >= min_subject_noNEC   ~ TRUE,
    NEC_stage == "preNEC"  & n_subjects >= min_subject_preNEC  ~ TRUE,
    NEC_stage == "postNEC" & n_subjects >= min_subject_postNEC ~ TRUE,
    TRUE ~ FALSE)) %>%
  filter(pass_stage_filter)

table(qualified_species_stage$NEC_stage)

## Check how many species remain after filtering in each stage
qualified_species_stage %>%
  count(NEC_stage, name = "n_species_pass")


## 3. Build filtered species-level Venn sets
sets_species_filtered <- list(noNEC = qualified_species_stage %>%
                                filter(NEC_stage == "noNEC") %>%
                                pull(Species) %>%
                                unique(),
                              
                              preNEC = qualified_species_stage %>%
                                filter(NEC_stage == "preNEC") %>%
                                pull(Species) %>%
                                unique(),
                              
                              postNEC = qualified_species_stage %>%
                                filter(NEC_stage == "postNEC") %>%
                                pull(Species) %>%
                                unique())

sapply(sets_species_filtered, length)



#### Visualization
## 1. Fix set order
## Must be run after sets_species_filtered has been created
sets_species_filtered <- sets_species_filtered[c("noNEC", "preNEC", "postNEC")]

## Check the number of species in each set
sapply(sets_species_filtered, length)

## 2. Define colors for the three stages
stage_colors <- c("noNEC"   = "grey70","preNEC"  = "#35978f","postNEC" = "#bf812d")

## 3. Build Venn object
venn_object_species <- ggVennDiagram::Venn(sets_species_filtered)
venn_data_species <- ggVennDiagram::process_data(venn_object_species)

## 4. Mapping between set IDs and stage names
## Because the order is fixed as noNEC, preNEC and postNEC,
## IDs 1, 2 and 3 correspond to the three stages, respectively
set_key_species <- tibble(id = as.character(seq_along(sets_species_filtered)),
                          Stage = names(sets_species_filtered))

## 5. Extract Venn intersection regions
species_region_df <- ggVennDiagram::venn_regionedge(venn_data_species) %>%
  mutate(id = as.character(id),
         
         ## Use the corresponding stage color for stage-specific regions;
         ## Use light grey for regions shared by two or three groups
         fill_group = case_when(id == "1" ~ "noNEC",
                                id == "2" ~ "preNEC",
                                id == "3" ~ "postNEC",
                                TRUE      ~ "Shared"))

## Optional check of region IDs
unique(species_region_df$id)

## 6. Extract set borders
species_set_edge_df <- ggVennDiagram::venn_setedge(venn_data_species) %>%
  mutate(id = as.character(id)) %>%
  left_join(set_key_species, by = "id")

## 7. Extract set label positions
species_set_label_df <- ggVennDiagram::venn_setlabel( venn_data_species) %>%
  mutate(id = as.character(id)) %>%
  left_join(set_key_species, by = "id")

## 8. Extract label positions for each intersection region
species_region_label_df <- ggVennDiagram::venn_regionlabel(venn_data_species)

## 9. Redraw the species-level Venn plot
p_species_venn_filtered_custom <- ggplot() +
  
  ## Fill regions
  geom_polygon(data = species_region_df,
               aes(x = X,y = Y,group = id,fill = fill_group ),
               color = NA,alpha = 0.75) +
  
  ## Borders of the three sets
  geom_path(data = species_set_edge_df,
            aes(x = X,y = Y,group = id,color = Stage),
            linewidth = 1.1,show.legend = FALSE) +
  
  ## Set labels
  geom_text(data = species_set_label_df,
            aes(x = X,y = Y,label = name,color = Stage),
            size = 5,fontface = "bold",show.legend = FALSE) +
  
  ## Species counts in each region
  geom_text(data = species_region_label_df,
            aes(x = X,y = Y,label = count),
            size = 4,color = "black") +
  
  ## Region colors
  scale_fill_manual(
    values = c("noNEC"   = "grey70","preNEC"  = "#35978f","postNEC" = "#bf812d","Shared"  = "grey92"),
    
    ## Do not show Shared separately in the legend
    breaks = c("noNEC", "preNEC", "postNEC"),name = "NEC stage" ) +
  
  ## Border and set-label colors
  scale_color_manual(values = stage_colors) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold",hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size = 10)) +
  
  labs(title = "Shared and stage-specific mutated species across NEC stages",
       subtitle = paste0("Stage-specific thresholds: noNEC ≥ ",min_subject_noNEC,
                         ", preNEC ≥ ",min_subject_preNEC,
                         ", postNEC ≥ ",min_subject_postNEC," infants"))
p_species_venn_filtered_custom
# save as 3 * 3


#
all_species_filtered <- sort(unique(unlist(sets_species_filtered)))

species_venn_membership_filtered <- data.frame(Species = all_species_filtered) %>%
  mutate(noNEC = Species %in% sets_species_filtered$noNEC,
         preNEC = Species %in% sets_species_filtered$preNEC,
         postNEC = Species %in% sets_species_filtered$postNEC,
         
         category = case_when(noNEC & !preNEC & !postNEC ~ "noNEC_only",
                              !noNEC & preNEC & !postNEC ~ "preNEC_only",
                              !noNEC & !preNEC & postNEC ~ "postNEC_only",
                              noNEC & preNEC & !postNEC ~ "noNEC_preNEC_shared",
                              noNEC & !preNEC & postNEC ~ "noNEC_postNEC_shared",
                              !noNEC & preNEC & postNEC ~ "preNEC_postNEC_shared_only",
                              noNEC & preNEC & postNEC ~ "shared_all_three",
                              TRUE ~ "other"))

table(species_venn_membership_filtered$category)
test <- species_venn_membership_filtered%>% subset(category=="preNEC_postNEC_shared_only")
test <- species_venn_membership_filtered%>% subset(category=="noNEC_postNEC_shared")
test <- species_venn_membership_filtered%>% subset(category=="noNEC_preNEC_shared")
test <- species_venn_membership_filtered%>% subset(category=="shared_all_three")
test$Species

#write.csv(species_venn_membership_filtered,
#  "filtered_species_level_mutated_species_venn_membership.csv",row.names = FALSE)
write_xlsx(
  species_venn_membership_filtered,
  path = "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 7_file/filtered_species_level_mutated_species_venn_membership.xlsx")





##### Fig. 7c. Candidate species-level SNV burden among 29 shared mutated species ####

##### panel 1. Shared genomes bumber ~ species number ####

### Extract only preNEC-NEC shared genomes
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

## 1. Extract the 1607 genomes shared by preNEC and NEC/postNEC
## Including:
## 1) preNEC_postNEC_shared: shared by preNEC and postNEC but absent from noNEC
## 2) shared_all_three: shared by noNEC, preNEC and postNEC
colnames(species_venn_membership_filtered)
table(species_venn_membership_filtered$category)

shared_29_species <- species_venn_membership_filtered %>%
  filter(category %in% c( "shared_all_three"#"preNEC_postNEC_shared",
  )) %>%
  distinct(Species, category)

dim(shared_29_species)
table(shared_29_species$category)

shared_species_ids <- shared_29_species$Species
length(shared_species_ids)


## Return to the original data and retain only these 1607 shared genomes
## 2. Retain only shared genomes and generate preNEC / NEC stages
shared_species_stage_df <- NEC_genome_SNVs_coverage5_top5_species2 %>%
  mutate(genome_id = as.character(.data[[genome_id_col]]),
         Species = as.character(Species),
         Species = ifelse(is.na(Species) | Species == "", "Unclassified", Species),
         NEC_stage = case_when(#PreNEC2 == "noNEC" ~ "noNEC",
           PreNEC2 == "preNEC" ~ "preNEC",
           PreNEC2 == "postNEC" ~ "postNEC",
           TRUE ~ NA_character_),
         NEC_stage = factor(NEC_stage, levels = c(#"noNEC",
           "preNEC", "postNEC")) ) %>%
  filter(Species %in% shared_species_ids,
         !is.na(NEC_stage),
         !is.na(sample2),
         !is.na(genome_id),
         # !is.na(Species),
         !is.na(SNV_count),
         SNV_count > 0) %>%
  ## sample-level detection: count the same genome only once within the same sample
  distinct(NEC_stage, Species, sample2,genome_id)

dim(shared_species_stage_df)
table(shared_species_stage_df$NEC_stage)
length(unique(shared_species_stage_df$genome_id))
length(unique(shared_species_stage_df$Species))


## Count shared mutated-genome detections in preNEC / NEC by species
## 3. Sample-level detection counts of shared genomes in preNEC / NEC for each species
shared_species_stage_count <- shared_species_stage_df %>%
  count(Species, NEC_stage, name = "n_genome_detections")

## Total detection count for each species, used for ordering
shared_species_total <- shared_species_stage_count %>%
  group_by(Species) %>%
  summarise(total_detections = sum(n_genome_detections),
            .groups = "drop" ) %>%
  arrange(desc(total_detections))



## 4. Long-tail stacked bar plot for all species
### Visualization method 1
species_order <- shared_species_total %>%
  arrange(desc(total_detections)) %>%
  pull(Species)

shared_species_stage_count_plot <- shared_species_stage_count %>%
  mutate(Species = factor(Species, levels = species_order),
         NEC_stage = factor(NEC_stage, levels = c("postNEC", "preNEC"#,"noNEC"
         )))






### Bacterial mutation sites differences 

## X species actually used in downstream analyses:
## Shared by preNEC and postNEC and passing the threshold in each group
X_species_df_filtered <- species_venn_membership_filtered %>%
  filter(preNEC == TRUE, postNEC == TRUE) %>%
  arrange(category, Species)

dim(X_species_df_filtered)
table(X_species_df_filtered$category)

#write.table(X_species_df_filtered,"filtered_X_preNEC_postNEC_shared_species.tsv",
#  sep = "\t",quote = FALSE,row.names = FALSE)

## Genome map for X species
genome_species_map_filtered <- NEC_genome_SNVs_coverage5_top5_species2 %>%
  mutate(genome_sp = as.character(.data[[genome_id_col]]),
         Species = as.character(Species),
         Genus = as.character(Genus),
         Family = as.character(Family),
         
         Species = ifelse(is.na(Species) | Species == "", "Unclassified", Species),
         Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus),
         Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)) %>%
  filter(Species %in% X_species_df_filtered$Species,!is.na(genome_sp)) %>%
  distinct(genome_sp,Species,Genus,Family)

dim(genome_species_map_filtered)
length(unique(genome_species_map_filtered$genome_sp))
length(unique(genome_species_map_filtered$Species))

#write.table(genome_species_map_filtered,"filtered_X_species_genome_map.tsv",
#  sep = "\t",quote = FALSE,row.names = FALSE)




### Run in Linux:
X_preNEC_postNEC_shared_species.tsv
X_species_genome_map.tsv
filtered_X_preNEC_postNEC_shared_species.tsv
filtered_X_species_genome_map.tsv

# Linux code:
duckdb_cli

.timer on

PRAGMA threads=20;
PRAGMA memory_limit='300GB';

CREATE OR REPLACE TEMP TABLE nec_samples AS
SELECT DISTINCT sample2
FROM read_csv(
  'NEC_sample_name.tsv',
  delim = '\t',
  header = false,
  columns = {'sample2': 'VARCHAR'}
);

CREATE OR REPLACE TEMP TABLE x_genome_map AS
SELECT *
  FROM read_csv_auto(
    'filtered_X_species_genome_map.tsv',
    delim = '\t',
    header = true,
    sample_size = -1
  );

SELECT COUNT(*) AS n_genomes FROM x_genome_map;
SELECT COUNT(DISTINCT Species) AS n_species FROM x_genome_map;

COPY (
  SELECT
  f.*,
  x.Species,
  x.Genus,
  x.Family
  FROM read_csv_auto(
    'all_samples_SNVs_with_genome.tsv',
    delim = '\t',
    header = true,
    sample_size = -1
  ) AS f
  INNER JOIN nec_samples AS n
  ON f.sample2 = n.sample2
  INNER JOIN x_genome_map AS x
  ON f.genome_sp = x.genome_sp
)
TO 'NEC_692samples_filtered_X_all_three_shared_species_SNVs.tsv'
WITH (
  DELIMITER '\t',
  HEADER true
);

## Check
CREATE OR REPLACE TEMP TABLE x_snv AS
SELECT *
  FROM read_csv_auto(
    'NEC_692samples_filtered_X_all_three_shared_species_SNVs.tsv',
    delim = '\t',
    header = true,
    sample_size = -1
  );

SELECT COUNT(*) AS n_rows FROM x_snv;

SELECT COUNT(DISTINCT sample2) AS n_samples FROM x_snv;

SELECT COUNT(DISTINCT genome_sp) AS n_genomes FROM x_snv;

SELECT COUNT(DISTINCT Species) AS n_species FROM x_snv;

SELECT Species, COUNT(*) AS n_snv_rows
FROM x_snv
GROUP BY Species
ORDER BY n_snv_rows DESC
LIMIT 30;


# linux code
ls -alh NEC_692samples_filtered_X_all_three_shared_species_SNVs.tsv
# 65 G
# 27G



##### Split each species into a separate file
.timer on

PRAGMA threads=40;
PRAGMA memory_limit='600GB';
PRAGMA temp_directory='/mnt/data/ShuangPeng/tmp_duckdb';
SET preserve_insertion_order=false;

COPY (
  SELECT *
    FROM read_csv_auto(
      'NEC_692samples_filtered_X_all_three_shared_species_SNVs.tsv',
      delim = '\t',
      header = true,
      sample_size = -1
    )
)
TO 'species_split_raw_SNVs_uncompressed'
WITH (
  FORMAT CSV,
  DELIMITER '\t',
  HEADER true,
  PARTITION_BY (Species),
  OVERWRITE_OR_IGNORE true
);


## Check
du -h --max-depth=2 species_split_raw_SNVs_uncompressed | sort -h
find species_split_raw_SNVs_uncompressed -type f | wc -l
find species_split_raw_SNVs_uncompressed -type f | head

less 'species_split_raw_SNVs_uncompressed/Species=Citrobacter_A amalonaticus/data_0.csv'


## Convert all outputs into independent species files in one directory
mkdir -p species_split_named_uncompressed

find species_split_raw_SNVs_uncompressed -name "data_0.csv" | while read f
do
spdir=$(basename "$(dirname "$f")")
sp=${spdir#Species=}
  safe=$(echo "$sp" | sed 's/[^A-Za-z0-9_]/_/g' | sed 's/_\+/_/g')
  cp "$f" "species_split_named_uncompressed/${safe}.tsv"
  done
  
  # Check
  ls -alh species_split_named_uncompressed | head
  du -h --max-depth=1 species_split_named_uncompressed | sort -h
  less species_split_named_uncompressed/Citrobacter_A_amalonaticus.tsv
  
  du -h --max-depth=2 species_split_raw_SNVs_uncompressed | sort -h
  rm -rf species_split_raw_SNVs_uncompressed
  
  
  
  
  
  
  ### In R
  Klebsiella_indica.tsv #
  
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  
  
  
  ######################## Batch run for panels 2-4: Direct SNV-site NEC association screen, parallel version ########
  ######################## In Linux
  options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))
  install.packages("ggplot2")
  install.packages("data.table")
  install.packages("dplyr")
  install.packages("tidyr")
  install.packages("vegan")
  install.packages("ggplot2")
  install.packages("ggrepel")
  install.packages("future")
  install.packages("future.apply")
  
  
  
  ######################## Calculation
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(future)
    library(future.apply)
  })
  
  ## 0. Run mode
  
  ## First test with TRUE; after successful testing, change to FALSE to run all
  #TEST_MODE <- TRUE
  TEST_MODE <- FALSE
  
  ## Test several specified species
  TEST_SPECIES <- c("Enterobacter hormaechei_A",
                    "Enterobacter quasihormaechei",
                    "Klebsiella michiganensis")
  
  N_WORKERS_TEST <- 3
  N_WORKERS_FULL <- 29
  
  ## site-level filtering
  MIN_COV <- 5
  MIN_ALT_COUNT <- 2
  MIN_ALT_FREQ <- 0.05
  
  MIN_SITE_SUBJECTS <- 5
  MIN_NONEC_SUBJECTS_SITE <- 5
  MIN_PRENEC_SUBJECTS_SITE <- 2
  MIN_ALT_SUBJECTS_SITE <- 2
  
  ## Maximum number of informative SNV sites retained for each species
  ## If a species is too large, this can be changed to 20000
  MAX_SITES_PER_SPECIES <- 50000
  
  ## candidate SNV-site criteria
  CANDIDATE_MIN_PRENEC_PREV <- 0.25
  CANDIDATE_MIN_DIFF_PREV <- 0.20
  CANDIDATE_MIN_ALT_PRENEC_SUBJECTS <- 3
  CANDIDATE_MAX_FISHER_P <- 0.05
  
  ## Whether to additionally calculate a Wilcoxon alt-frequency test for candidate sites
  DO_WILCOX_FOR_CANDIDATES <- TRUE
  
  
  ## 1. Parallel-safe settings
  data.table::setDTthreads(1)
  
  Sys.setenv(OMP_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1",
             NUMEXPR_NUM_THREADS = "1")
  options(future.globals.maxSize = 100 * 1024^3)
  
  
  ## 2. Load metadata
  NEC_meta_infor_sub2 <- readRDS("NEC_meta_infor_sub2.rds")
  
  meta_stage <- NEC_meta_infor_sub2 %>%
    mutate(NEC_stage = case_when(PreNEC2 == "noNEC" ~ "noNEC",
                                 PreNEC2 == "preNEC" ~ "preNEC",
                                 PreNEC2 %in% c("postNEC", "NEC", "onsetNEC") ~ "postNEC",
                                 TRUE ~ NA_character_),
           NEC_primary = case_when(PreNEC2 == "noNEC" ~ "noNEC",
                                   PreNEC2 == "preNEC" ~ "preNEC",
                                   TRUE ~ NA_character_)) %>%
    dplyr::select(sample2, SubjectID, NEC_stage, NEC_primary, DOL, Study) %>%
    distinct() %>%
    mutate(sample2 = as.character(sample2),
           SubjectID = as.character(SubjectID),
           NEC_stage = as.character(NEC_stage),
           NEC_primary = as.character(NEC_primary),
           DOL = as.numeric(DOL),
           Study = as.character(Study))
  cat("\nNEC_stage table:\n")
  print(table(meta_stage$NEC_stage, useNA = "ifany"))
  
  cat("\nNEC_primary table:\n")
  print(table(meta_stage$NEC_primary, useNA = "ifany"))
  
  
  ## 3. Species file map
  species_file_map <- fread("filtered_X_species_genome_map.tsv") %>%
    distinct(Species, Genus, Family) %>%
    mutate(Species = as.character(Species),
           Genus = as.character(Genus),
           Family = as.character(Family),
           safe_species = gsub("[^A-Za-z0-9_]+", "_", Species),
           safe_species = gsub("_+", "_", safe_species),
           file = file.path("species_split_named_uncompressed",
                            paste0(safe_species, ".tsv"))) %>%
    filter(file.exists(file)) %>%
    arrange(Species)
  
  cat("\nSpecies files found:\n")
  print(dim(species_file_map))
  print(species_file_map)
  
  
  ## 4. Helper functions
  safe_first_char <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) {
      return(NA_character_)
    } else {
      return(x[1])
    }
  }
  
  safe_first_num <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    } else {
      return(x[1])
    }
  }
  
  safe_fisher <- function(alt_pre, n_pre, alt_no, n_no) {
    
    tab <- matrix(
      c(
        alt_pre,
        n_pre - alt_pre,
        alt_no,
        n_no - alt_no
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    out <- tryCatch(
      fisher.test(tab),
      error = function(e) NULL
    )
    
    if (is.null(out)) {
      return(data.frame(fisher_p = NA_real_, fisher_OR = NA_real_))
    }
    
    data.frame(
      fisher_p = out$p.value,
      fisher_OR = unname(out$estimate)
    )
  }
  
  
  ## Function 1: build subject-level SNV-site table
  build_subject_site_table <- function(
    sp_file,
    sp_name,
    meta_stage,
    min_cov = 5,
    min_alt_count = 2,
    min_alt_freq = 0.05,
    min_site_subjects = 5,
    min_noNEC_subjects_site = 5,
    min_preNEC_subjects_site = 2,
    min_alt_subjects_site = 2,
    max_sites = 50000
  ) {
    
    message("Reading species file: ", sp_name)
    
    dt <- fread(sp_file)
    
    if (!"Species" %in% colnames(dt)) {
      dt[, Species := sp_name]
    }
    
    if (!"Genus" %in% colnames(dt)) {
      dt[, Genus := NA_character_]
    }
    
    if (!"Family" %in% colnames(dt)) {
      dt[, Family := NA_character_]
    }
    
    required_cols <- c(
      "sample2", "genome_sp", "scaffold", "position",
      "A", "C", "T", "G"
    )
    
    miss_cols <- setdiff(required_cols, colnames(dt))
    if (length(miss_cols) > 0) {
      stop("Missing columns: ", paste(miss_cols, collapse = ", "))
    }
    
    dt[, sample2 := as.character(sample2)]
    
    dt <- dt %>%
      left_join(meta_stage, by = "sample2")
    
    setDT(dt)
    
    dt <- dt[
      NEC_stage %in% c("noNEC", "preNEC", "postNEC") &
        !is.na(SubjectID) &
        !is.na(NEC_stage)
    ]
    
    if (nrow(dt) == 0) return(NULL)
    
    dt[, site_id := paste(genome_sp, scaffold, position, sep = "__")]
    
    dt[, A := as.numeric(A)]
    dt[, C := as.numeric(C)]
    dt[, T := as.numeric(T)]
    dt[, G := as.numeric(G)]
    
    ## 1. Global major allele per site
    site_major <- dt[
      ,
      .(
        A_total = sum(A, na.rm = TRUE),
        C_total = sum(C, na.rm = TRUE),
        T_total = sum(T, na.rm = TRUE),
        G_total = sum(G, na.rm = TRUE)
      ),
      by = site_id
    ]
    
    site_major[
      ,
      major_base := c("A", "C", "T", "G")[
        max.col(
          as.matrix(.SD),
          ties.method = "first"
        )
      ],
      .SDcols = c("A_total", "C_total", "T_total", "G_total")
    ]
    
    dt <- merge(
      dt,
      site_major[, .(site_id, major_base)],
      by = "site_id",
      all.x = TRUE
    )
    
    ## 2. alt count / alt freq / alt present
    dt[, total_count := A + C + T + G]
    
    dt[
      ,
      major_count := fcase(
        major_base == "A", A,
        major_base == "C", C,
        major_base == "T", T,
        major_base == "G", G,
        default = NA_real_
      )
    ]
    
    dt[, alt_count := total_count - major_count]
    dt[, alt_freq := alt_count / total_count]
    
    dt[
      ,
      alt_present := as.integer(
        alt_count >= min_alt_count &
          alt_freq >= min_alt_freq
      )
    ]
    
    dt <- dt[
      total_count >= min_cov &
        !is.na(alt_freq) &
        !is.na(site_id)
    ]
    
    if (nrow(dt) == 0) return(NULL)
    
    ## 3. Collapse to subject × NEC_stage × site
    subject_site <- dt[
      ,
      .(
        alt_present = as.integer(any(alt_present == 1, na.rm = TRUE)),
        alt_freq_subject = max(alt_freq, na.rm = TRUE),
        total_alt_count = sum(alt_count, na.rm = TRUE),
        total_coverage = sum(total_count, na.rm = TRUE),
        n_samples = uniqueN(sample2),
        
        NEC_primary = safe_first_char(NEC_primary),
        DOL = median(DOL, na.rm = TRUE),
        Study = safe_first_char(Study),
        
        Species = safe_first_char(Species),
        Genus = safe_first_char(Genus),
        Family = safe_first_char(Family),
        genome_sp = safe_first_char(genome_sp),
        scaffold = safe_first_char(scaffold),
        position = safe_first_num(position),
        major_base = safe_first_char(major_base),
        
        gene = if ("gene" %in% colnames(dt)) safe_first_char(gene) else NA_character_,
        mutation = if ("mutation" %in% colnames(dt)) safe_first_char(mutation) else NA_character_,
        mutation_type = if ("mutation_type" %in% colnames(dt)) safe_first_char(mutation_type) else NA_character_
      ),
      by = .(
        SubjectID,
        NEC_stage,
        site_id
      )
    ]
    
    subject_site[
      ,
      row_id := paste(SubjectID, NEC_stage, sep = "__")
    ]
    
    ## 4. Row universe per species
    ## All subject-stage records with SNV records in this species
    row_universe <- subject_site[
      ,
      .(
        SubjectID = safe_first_char(SubjectID),
        NEC_stage = safe_first_char(NEC_stage),
        NEC_primary = safe_first_char(NEC_primary),
        DOL = median(DOL, na.rm = TRUE),
        Study = safe_first_char(Study)
      ),
      by = row_id
    ]
    
    ## Primary comparison universe: noNEC vs preNEC
    row_primary <- row_universe[
      NEC_stage %in% c("noNEC", "preNEC")
    ]
    
    n_noNEC <- uniqueN(row_primary$row_id[row_primary$NEC_stage == "noNEC"])
    n_preNEC <- uniqueN(row_primary$row_id[row_primary$NEC_stage == "preNEC"])
    n_postNEC <- uniqueN(row_universe$row_id[row_universe$NEC_stage == "postNEC"])
    
    if (n_noNEC < 5 || n_preNEC < 2) {
      return(NULL)
    }
    
    ## 5. Site filtering based on primary comparison
    site_stats <- subject_site[
      NEC_stage %in% c("noNEC", "preNEC"),
      .(
        n_subject_stage = uniqueN(row_id),
        n_noNEC_subjects = uniqueN(row_id[NEC_stage == "noNEC"]),
        n_preNEC_subjects = uniqueN(row_id[NEC_stage == "preNEC"]),
        n_alt_subjects = sum(alt_present == 1, na.rm = TRUE),
        var_alt_freq = var(alt_freq_subject, na.rm = TRUE)
      ),
      by = site_id
    ]
    
    site_keep <- site_stats[
      n_subject_stage >= min_site_subjects &
        n_noNEC_subjects >= min_noNEC_subjects_site &
        n_preNEC_subjects >= min_preNEC_subjects_site &
        n_alt_subjects >= min_alt_subjects_site &
        !is.na(var_alt_freq) &
        var_alt_freq > 0
    ][
      order(-n_alt_subjects, -var_alt_freq)
    ]
    
    if (nrow(site_keep) == 0) {
      return(NULL)
    }
    
    if (nrow(site_keep) > max_sites) {
      site_keep <- site_keep[1:max_sites]
    }
    
    subject_site2 <- subject_site[
      site_id %in% site_keep$site_id
    ]
    
    site_anno <- subject_site2[
      ,
      .(
        Species = safe_first_char(Species),
        Genus = safe_first_char(Genus),
        Family = safe_first_char(Family),
        genome_sp = safe_first_char(genome_sp),
        scaffold = safe_first_char(scaffold),
        position = safe_first_num(position),
        major_base = safe_first_char(major_base),
        gene = safe_first_char(gene),
        mutation = safe_first_char(mutation),
        mutation_type = safe_first_char(mutation_type)
      ),
      by = site_id
    ]
    
    list(
      Species = sp_name,
      subject_site = subject_site2,
      row_universe = row_universe,
      row_primary = row_primary,
      site_stats = site_stats,
      site_keep = site_keep,
      site_anno = site_anno,
      n_noNEC = n_noNEC,
      n_preNEC = n_preNEC,
      n_postNEC = n_postNEC
    )
  }
  
  
  ## Function 2: per-site NEC association
  run_site_association <- function(obj) {
    
    ss <- as.data.table(obj$subject_site)
    site_anno <- as.data.table(obj$site_anno)
    row_universe <- as.data.table(obj$row_universe)
    row_primary <- as.data.table(obj$row_primary)
    
    n_noNEC <- obj$n_noNEC
    n_preNEC <- obj$n_preNEC
    n_postNEC <- obj$n_postNEC
    
    ## Primary noNEC vs preNEC counts
    ## Missing site in a subject-stage is treated as alt_present = 0
    
    primary_counts <- ss[
      NEC_stage %in% c("noNEC", "preNEC"),
      .(
        alt_noNEC_subjects = sum(alt_present[NEC_stage == "noNEC"] == 1, na.rm = TRUE),
        alt_preNEC_subjects = sum(alt_present[NEC_stage == "preNEC"] == 1, na.rm = TRUE),
        
        detected_noNEC_subjects = uniqueN(row_id[NEC_stage == "noNEC"]),
        detected_preNEC_subjects = uniqueN(row_id[NEC_stage == "preNEC"]),
        
        sum_altfreq_noNEC = sum(alt_freq_subject[NEC_stage == "noNEC"], na.rm = TRUE),
        sum_altfreq_preNEC = sum(alt_freq_subject[NEC_stage == "preNEC"], na.rm = TRUE),
        
        max_altfreq_noNEC = suppressWarnings(max(alt_freq_subject[NEC_stage == "noNEC"], na.rm = TRUE)),
        max_altfreq_preNEC = suppressWarnings(max(alt_freq_subject[NEC_stage == "preNEC"], na.rm = TRUE))
      ),
      by = site_id
    ]
    
    primary_counts[
      is.infinite(max_altfreq_noNEC),
      max_altfreq_noNEC := NA_real_
    ]
    primary_counts[
      is.infinite(max_altfreq_preNEC),
      max_altfreq_preNEC := NA_real_
    ]
    
    ## PostNEC retention counts
    
    post_counts <- ss[
      NEC_stage == "postNEC",
      .(
        alt_postNEC_subjects = sum(alt_present == 1, na.rm = TRUE),
        detected_postNEC_subjects = uniqueN(row_id),
        sum_altfreq_postNEC = sum(alt_freq_subject, na.rm = TRUE),
        max_altfreq_postNEC = suppressWarnings(max(alt_freq_subject, na.rm = TRUE))
      ),
      by = site_id
    ]
    
    if (nrow(post_counts) > 0) {
      post_counts[
        is.infinite(max_altfreq_postNEC),
        max_altfreq_postNEC := NA_real_
      ]
    }
    
    res <- merge(
      primary_counts,
      post_counts,
      by = "site_id",
      all.x = TRUE
    )
    
    for (cc in c(
      "alt_postNEC_subjects",
      "detected_postNEC_subjects",
      "sum_altfreq_postNEC"
    )) {
      if (cc %in% colnames(res)) {
        res[is.na(get(cc)), (cc) := 0]
      }
    }
    
    ## Denominators
    res[, n_noNEC := n_noNEC]
    res[, n_preNEC := n_preNEC]
    res[, n_postNEC := n_postNEC]
    
    ## Prevalence
    res[, noNEC_prev := alt_noNEC_subjects / n_noNEC]
    res[, preNEC_prev := alt_preNEC_subjects / n_preNEC]
    res[, postNEC_prev := fifelse(n_postNEC > 0, alt_postNEC_subjects / n_postNEC, NA_real_)]
    
    res[, diff_prev_preNEC_vs_noNEC := preNEC_prev - noNEC_prev]
    res[, diff_prev_postNEC_vs_noNEC := postNEC_prev - noNEC_prev]
    res[, diff_prev_postNEC_vs_preNEC := postNEC_prev - preNEC_prev]
    
    ## Mean alt frequency, missing site treated as 0
    res[, mean_altfreq_noNEC := sum_altfreq_noNEC / n_noNEC]
    res[, mean_altfreq_preNEC := sum_altfreq_preNEC / n_preNEC]
    res[, mean_altfreq_postNEC := fifelse(n_postNEC > 0, sum_altfreq_postNEC / n_postNEC, NA_real_)]
    
    res[, diff_mean_altfreq_preNEC_vs_noNEC := mean_altfreq_preNEC - mean_altfreq_noNEC]
    res[, diff_mean_altfreq_postNEC_vs_noNEC := mean_altfreq_postNEC - mean_altfreq_noNEC]
    
    ## Fisher exact test
    fisher_list <- lapply(seq_len(nrow(res)), function(i) {
      safe_fisher(
        alt_pre = res$alt_preNEC_subjects[i],
        n_pre = res$n_preNEC[i],
        alt_no = res$alt_noNEC_subjects[i],
        n_no = res$n_noNEC[i]
      )
    })
    
    fisher_df <- rbindlist(fisher_list)
    
    res[, fisher_p := fisher_df$fisher_p]
    res[, fisher_OR_preNEC_vs_noNEC := fisher_df$fisher_OR]
    
    res[, fisher_FDR_species := p.adjust(fisher_p, method = "BH")]
    
    res[, direction := case_when(
      diff_prev_preNEC_vs_noNEC > 0 ~ "preNEC-enriched",
      diff_prev_preNEC_vs_noNEC < 0 ~ "noNEC-enriched",
      TRUE ~ "no_difference"
    )]
    
    ## Merge annotation
    res <- merge(
      res,
      site_anno,
      by = "site_id",
      all.x = TRUE
    )
    
    ## Candidate sites
    
    candidate_nominal <- res[
      direction == "preNEC-enriched" &
        preNEC_prev >= CANDIDATE_MIN_PRENEC_PREV &
        diff_prev_preNEC_vs_noNEC >= CANDIDATE_MIN_DIFF_PREV &
        alt_preNEC_subjects >= CANDIDATE_MIN_ALT_PRENEC_SUBJECTS &
        fisher_p < CANDIDATE_MAX_FISHER_P
    ][
      order(fisher_p, -diff_prev_preNEC_vs_noNEC)
    ]
    
    candidate_FDR10 <- res[
      direction == "preNEC-enriched" &
        preNEC_prev >= CANDIDATE_MIN_PRENEC_PREV &
        diff_prev_preNEC_vs_noNEC >= CANDIDATE_MIN_DIFF_PREV &
        alt_preNEC_subjects >= CANDIDATE_MIN_ALT_PRENEC_SUBJECTS &
        fisher_FDR_species < 0.10
    ][
      order(fisher_FDR_species, fisher_p, -diff_prev_preNEC_vs_noNEC)
    ]
    
    ## Optional Wilcoxon only for nominal candidates
    
    if (DO_WILCOX_FOR_CANDIDATES && nrow(candidate_nominal) > 0) {
      
      candidate_ids <- candidate_nominal$site_id
      
      wilcox_res <- lapply(candidate_ids, function(sid) {
        
        tmp <- ss[site_id == sid & NEC_stage %in% c("noNEC", "preNEC")]
        
        vec <- data.table(
          row_id = row_primary$row_id,
          NEC_stage = row_primary$NEC_stage,
          alt_freq = 0
        )
        
        tmp2 <- tmp[, .(alt_freq = max(alt_freq_subject, na.rm = TRUE)), by = row_id]
        vec[tmp2, alt_freq := i.alt_freq, on = "row_id"]
        
        p <- tryCatch(
          wilcox.test(
            alt_freq ~ NEC_stage,
            data = vec,
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
        
        data.table(
          site_id = sid,
          wilcox_altfreq_p = p,
          median_altfreq_noNEC = median(vec$alt_freq[vec$NEC_stage == "noNEC"], na.rm = TRUE),
          median_altfreq_preNEC = median(vec$alt_freq[vec$NEC_stage == "preNEC"], na.rm = TRUE)
        )
      })
      
      wilcox_dt <- rbindlist(wilcox_res, fill = TRUE)
      wilcox_dt[, wilcox_altfreq_FDR_candidate := p.adjust(wilcox_altfreq_p, method = "BH")]
      
      res <- merge(
        res,
        wilcox_dt,
        by = "site_id",
        all.x = TRUE
      )
      
      candidate_nominal <- merge(
        candidate_nominal,
        wilcox_dt,
        by = "site_id",
        all.x = TRUE
      )
      
      if (nrow(candidate_FDR10) > 0) {
        candidate_FDR10 <- merge(
          candidate_FDR10,
          wilcox_dt,
          by = "site_id",
          all.x = TRUE
        )
      }
    } else {
      res[, wilcox_altfreq_p := NA_real_]
      res[, median_altfreq_noNEC := NA_real_]
      res[, median_altfreq_preNEC := NA_real_]
      res[, wilcox_altfreq_FDR_candidate := NA_real_]
    }
    
    list(
      all_sites = res,
      candidate_nominal = candidate_nominal,
      candidate_FDR10 = candidate_FDR10
    )
  }
  
  
  ## Function 3: run one species
  
  run_one_species_site_screen <- function(i) {
    
    suppressPackageStartupMessages({
      library(data.table)
      library(dplyr)
      library(tidyr)
    })
    
    data.table::setDTthreads(1)
    
    sp <- species_file_map$Species[i]
    sp_file <- species_file_map$file[i]
    safe_sp <- species_file_map$safe_species[i]
    
    log_file <- file.path(
      "site_assoc_logs",
      paste0(safe_sp, "_site_assoc_log.txt")
    )
    
    cat(
      "\n==============================\n",
      "Running site-level association for species ", i, "/", nrow(species_file_map), ": ", sp, "\n",
      "File: ", sp_file, "\n",
      "==============================\n",
      file = log_file,
      append = TRUE
    )
    
    obj <- tryCatch(
      build_subject_site_table(
        sp_file = sp_file,
        sp_name = sp,
        meta_stage = meta_stage,
        min_cov = MIN_COV,
        min_alt_count = MIN_ALT_COUNT,
        min_alt_freq = MIN_ALT_FREQ,
        min_site_subjects = MIN_SITE_SUBJECTS,
        min_noNEC_subjects_site = MIN_NONEC_SUBJECTS_SITE,
        min_preNEC_subjects_site = MIN_PRENEC_SUBJECTS_SITE,
        min_alt_subjects_site = MIN_ALT_SUBJECTS_SITE,
        max_sites = MAX_SITES_PER_SPECIES
      ),
      error = function(e) {
        cat("Build subject-site failed: ", e$message, "\n", file = log_file, append = TRUE)
        return(NULL)
      }
    )
    
    if (is.null(obj)) {
      cat("Object is NULL. Skip.\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    out <- tryCatch(
      run_site_association(obj),
      error = function(e) {
        cat("Site association failed: ", e$message, "\n", file = log_file, append = TRUE)
        return(NULL)
      }
    )
    
    if (is.null(out)) {
      cat("No valid site association result.\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    fwrite(
      out$all_sites,
      file.path(
        "site_assoc_results",
        paste0(safe_sp, "_SNV_site_NEC_assoc.tsv")
      ),
      sep = "\t"
    )
    
    fwrite(
      out$candidate_nominal,
      file.path(
        "site_assoc_results",
        paste0(safe_sp, "_preNEC_enriched_candidate_sites_nominal.tsv")
      ),
      sep = "\t"
    )
    
    fwrite(
      out$candidate_FDR10,
      file.path(
        "site_assoc_results",
        paste0(safe_sp, "_preNEC_enriched_candidate_sites_FDR10.tsv")
      ),
      sep = "\t"
    )
    
    summary_dt <- data.table(
      Species = sp,
      n_sites_tested = nrow(out$all_sites),
      n_nominal_preNEC_enriched_sites = nrow(out$candidate_nominal),
      n_FDR10_preNEC_enriched_sites = nrow(out$candidate_FDR10),
      n_noNEC = unique(out$all_sites$n_noNEC),
      n_preNEC = unique(out$all_sites$n_preNEC),
      n_postNEC = unique(out$all_sites$n_postNEC),
      min_fisher_p = suppressWarnings(min(out$all_sites$fisher_p, na.rm = TRUE)),
      min_fisher_FDR_species = suppressWarnings(min(out$all_sites$fisher_FDR_species, na.rm = TRUE))
    )
    
    cat("Finished successfully.\n", file = log_file, append = TRUE)
    
    list(
      summary = summary_dt,
      all_sites = out$all_sites,
      candidate_nominal = out$candidate_nominal,
      candidate_FDR10 = out$candidate_FDR10
    )
  }
  
  
  ## 5. Run TEST or FULL
  
  dir.create("site_assoc_results", showWarnings = FALSE)
  dir.create("site_assoc_logs", showWarnings = FALSE)
  
  if (TEST_MODE) {
    
    cat("\nRunning TEST mode for site-level association...\n")
    
    run_idx <- which(species_file_map$Species %in% TEST_SPECIES)
    
    if (length(run_idx) == 0) {
      stop("None of TEST_SPECIES found in species_file_map.")
    }
    
    n_workers <- min(N_WORKERS_TEST, length(run_idx))
    out_prefix <- "TEST_SNV_site_NEC_association"
    
  } else {
    
    cat("\nRunning FULL mode for site-level association...\n")
    
    run_idx <- seq_len(nrow(species_file_map))
    n_workers <- min(N_WORKERS_FULL, length(run_idx))
    out_prefix <- "ALL_species_SNV_site_NEC_association"
  }
  
  cat("\nWorkers: ", n_workers, "\n")
  cat("Species to run: ", length(run_idx), "\n")
  
  plan(multicore, workers = n_workers)
  
  set.seed(123)
  
  site_screen_list <- future_lapply(
    run_idx,
    run_one_species_site_screen,
    future.seed = TRUE
  )
  
  site_screen_list <- site_screen_list[!vapply(site_screen_list, is.null, logical(1))]
  
  if (length(site_screen_list) == 0) {
    stop("No valid species-level site association results. Check site_assoc_logs/")
  }
  
  summary_all <- rbindlist(
    lapply(site_screen_list, function(x) x$summary),
    fill = TRUE
  )
  
  all_sites <- rbindlist(
    lapply(site_screen_list, function(x) x$all_sites),
    fill = TRUE
  )
  
  candidate_nominal_all <- rbindlist(
    lapply(site_screen_list, function(x) x$candidate_nominal),
    fill = TRUE
  )
  
  candidate_FDR10_all <- rbindlist(
    lapply(site_screen_list, function(x) x$candidate_FDR10),
    fill = TRUE
  )
  
  ## Global FDR across all tested sites
  all_sites[, fisher_FDR_global := p.adjust(fisher_p, method = "BH")]
  
  ## Refresh global-level candidates
  candidate_global_FDR10 <- all_sites[
    direction == "preNEC-enriched" &
      preNEC_prev >= CANDIDATE_MIN_PRENEC_PREV &
      diff_prev_preNEC_vs_noNEC >= CANDIDATE_MIN_DIFF_PREV &
      alt_preNEC_subjects >= CANDIDATE_MIN_ALT_PRENEC_SUBJECTS &
      fisher_FDR_global < 0.10
  ][
    order(fisher_FDR_global, fisher_p, -diff_prev_preNEC_vs_noNEC)
  ]
  
  ## Write combined outputs
  fwrite(
    summary_all,
    paste0(out_prefix, "_species_summary.tsv"),
    sep = "\t"
  )
  
  fwrite(
    all_sites,
    paste0(out_prefix, "_all_sites.tsv"),
    sep = "\t"
  )
  
  fwrite(
    candidate_nominal_all,
    paste0(out_prefix, "_preNEC_enriched_candidate_sites_nominal.tsv"),
    sep = "\t"
  )
  
  fwrite(
    candidate_FDR10_all,
    paste0(out_prefix, "_preNEC_enriched_candidate_sites_species_FDR10.tsv"),
    sep = "\t"
  )
  
  fwrite(
    candidate_global_FDR10,
    paste0(out_prefix, "_preNEC_enriched_candidate_sites_global_FDR10.tsv"),
    sep = "\t"
  )
  
  ## 6. Gene-level burden from nominal candidate sites
  if (nrow(candidate_nominal_all) > 0) {
    
    gene_burden_nominal <- candidate_nominal_all %>%
      as.data.frame() %>%
      mutate(
        gene = ifelse(is.na(gene) | gene == "", "intergenic_or_unannotated", gene),
        mutation_type = ifelse(is.na(mutation_type) | mutation_type == "", "unknown", mutation_type)
      ) %>%
      group_by(Species, Genus, Family, genome_sp, gene) %>%
      summarise(
        n_candidate_sites = n_distinct(site_id),
        n_preNEC_alt_subjects_total = sum(alt_preNEC_subjects, na.rm = TRUE),
        mean_preNEC_prev = mean(preNEC_prev, na.rm = TRUE),
        mean_noNEC_prev = mean(noNEC_prev, na.rm = TRUE),
        mean_diff_prev = mean(diff_prev_preNEC_vs_noNEC, na.rm = TRUE),
        max_diff_prev = max(diff_prev_preNEC_vs_noNEC, na.rm = TRUE),
        min_fisher_p = min(fisher_p, na.rm = TRUE),
        min_fisher_FDR_species = min(fisher_FDR_species, na.rm = TRUE),
        n_postNEC_retained_sites = sum(postNEC_prev > 0, na.rm = TRUE),
        mean_postNEC_prev = mean(postNEC_prev, na.rm = TRUE),
        mutation_types = paste(sort(unique(mutation_type)), collapse = ";"),
        .groups = "drop"
      ) %>%
      arrange(desc(n_candidate_sites), min_fisher_p, desc(max_diff_prev))
    
  } else {
    
    gene_burden_nominal <- data.frame()
  }
  
  fwrite(
    gene_burden_nominal,
    paste0(out_prefix, "_preNEC_enriched_gene_burden_nominal.tsv"),
    sep = "\t"
  )
  
  ## 7. Print top outputs
  cat("\nSpecies summary:\n")
  print(summary_all[order(-n_nominal_preNEC_enriched_sites)])
  
  cat("\nTop nominal preNEC-enriched sites:\n")
  if (nrow(candidate_nominal_all) > 0) {
    print(
      candidate_nominal_all[
        order(fisher_p, -diff_prev_preNEC_vs_noNEC)
      ][
        ,
        .(
          Species,
          site_id,
          gene,
          mutation,
          mutation_type,
          n_noNEC,
          n_preNEC,
          alt_noNEC_subjects,
          alt_preNEC_subjects,
          noNEC_prev,
          preNEC_prev,
          diff_prev_preNEC_vs_noNEC,
          postNEC_prev,
          fisher_p,
          fisher_FDR_species
        )
      ][1:min(.N, 30)]
    )
  } else {
    cat("No nominal candidate sites found.\n")
  }
  
  cat("\nTop gene burden:\n")
  if (nrow(gene_burden_nominal) > 0) {
    print(head(gene_burden_nominal, 30))
  } else {
    cat("No gene burden result because no nominal candidate sites found.\n")
  }
  
  cat("\nDone.\n")
  
  ## Check test files
  tail -f NEC_direct_SNV_site_association_TEST.log
  
  less TEST_SNV_site_NEC_association_species_summary.tsv
  less TEST_SNV_site_NEC_association_all_sites.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_species_FDR10.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_global_FDR10.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_gene_burden_nominal.tsv
  
  
  
  less TEST_SNV_site_NEC_association_species_summary.tsv
  less TEST_SNV_site_NEC_association_all_sites.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_species_FDR10.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_candidate_sites_global_FDR10.tsv
  less TEST_SNV_site_NEC_association_preNEC_enriched_gene_burden_nominal.tsv
  
  
  
  # Real results:
  tail -f NEC_direct_SNV_site_association_FULL.log
  
  less ALL_species_SNV_site_NEC_association_all_sites.tsv
  less ALL_species_SNV_site_NEC_association_species_summary.tsv
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_species_FDR10.tsv
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_global_FDR10.tsv
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_gene_burden_nominal.tsv
  
  
  
  less ALL_species_SNV_site_NEC_association_species_summary.tsv
  
  less ALL_species_SNV_site_NEC_association_all_sites.tsv
  
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv
  
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_gene_burden_nominal.tsv
  
  
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_species_FDR10.tsv
  site_id alt_noNEC_subjects      alt_preNEC_subjects     detected_noNEC_subjects detected_preNEC_subjects        sum_altfreq_noNEC       sum_altfreq_preNEC        max_altfreq_noNEC       max_altfreq_preNEC      alt_postNEC_subjects    detected_postNEC_subjects       sum_altfreq_postNEC     max_altfreq_postNEC       n_noNEC n_preNEC        n_postNEC       noNEC_prev      preNEC_prev     postNEC_prev    diff_prev_preNEC_vs_noNEC       diff_prev_postNEC_vs_noNEC        diff_prev_postNEC_vs_preNEC     mean_altfreq_noNEC      mean_altfreq_preNEC     mean_altfreq_postNEC    diff_mean_altfreq_preNEC_vs_noNEC diff_mean_altfreq_postNEC_vs_noNEC      fisher_p        fisher_OR_preNEC_vs_noNEC       fisher_FDR_species      direction         Species Genus   Family  genome_sp       scaffold        position        major_base      gene    mutation        mutation_type
  
  
  less ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_global_FDR10.tsv
  site_id alt_noNEC_subjects      alt_preNEC_subjects     detected_noNEC_subjects detected_preNEC_subjects        sum_altfreq_noNEC       sum_altfreq_preNEC        max_altfreq_noNEC       max_altfreq_preNEC      alt_postNEC_subjects    detected_postNEC_subjects       sum_altfreq_postNEC     max_altfreq_postNEC       n_noNEC n_preNEC        n_postNEC       noNEC_prev      preNEC_prev     postNEC_prev    diff_prev_preNEC_vs_noNEC       diff_prev_postNEC_vs_noNEC        diff_prev_postNEC_vs_preNEC     mean_altfreq_noNEC      mean_altfreq_preNEC     mean_altfreq_postNEC    diff_mean_altfreq_preNEC_vs_noNEC diff_mean_altfreq_postNEC_vs_noNEC      fisher_p        fisher_OR_preNEC_vs_noNEC       fisher_FDR_species      direction         Species Genus   Family  genome_sp       scaffold        position        major_base      gene    mutation        mutation_type   wilcox_altfreq_p  median_altfreq_noNEC    median_altfreq_preNEC   wilcox_altfreq_FDR_candidate    fisher_FDR_global
  
  
  
  ALL_species_SNV_site_NEC_association_species_summary <- fread(
    "D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_species_summary.tsv",
    sep="\t",header=TRUE)
  
  
  ALL_species_SNV_site_NEC_association_all_sites <- fread(
    "D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_all_sites.tsv",
    sep="\t",header=TRUE)
  write_xlsx(ALL_species_SNV_site_NEC_association_all_sites, "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 7_file/ALL_species_SNV_site_NEC_association_all_sites.xlsx")
  
  
  test <- ALL_species_SNV_site_NEC_association_all_sites %>%
    subset(Species %in% c("Klebsiella pneumoniae"#,
                          #"Enterobacter hormaechei_A",
                          #"Escherichia coli",
                          #"Klebsiella michiganensis"
    ))
  
  
  
  
  ##### R code:
  ####### Merge all plotting data into one table
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(scales)
  library(openxlsx)
  
  ## 1. Read species-level summary
  sp <- fread(
    "D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_species_summary.tsv"
  )
  
  ## 2. Prepare species-level summary table
  sp2 <- sp %>%
    mutate(
      n_detected_subjects = n_noNEC + n_preNEC + n_postNEC,
      nominal_rate = n_nominal_preNEC_enriched_sites / n_sites_tested,
      FDR10_rate = ifelse(
        n_sites_tested > 0,
        n_FDR10_preNEC_enriched_sites / n_sites_tested,
        NA_real_
      )
    ) %>%
    arrange(desc(n_nominal_preNEC_enriched_sites))
  
  ## 3. Define species order
  species_order <- sp2 %>%
    arrange(n_nominal_preNEC_enriched_sites, n_sites_tested) %>%
    pull(Species)
  
  ## 4. Prepare shared mutated genome table
  shared_species_wide <- shared_species_stage_count_plot %>%
    group_by(Species, NEC_stage) %>%
    summarise(
      n_shared_mutated_genomes = sum(n_genome_detections, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = NEC_stage,
      values_from = n_shared_mutated_genomes,
      names_prefix = "shared_mutated_genomes_",
      values_fill = 0
    )
  colnames(shared_species_wide)
  
  ## 5. Merge all plotting data into one table
  plot_table <- sp2 %>%
    left_join(shared_species_wide, by = "Species") %>%
    mutate(
      # shared_mutated_genomes_noNEC = ifelse(
      #   is.na(shared_mutated_genomes_noNEC),
      #   0,
      #   shared_mutated_genomes_noNEC
      # ),
      shared_mutated_genomes_preNEC = ifelse(
        is.na(shared_mutated_genomes_preNEC),
        0,
        shared_mutated_genomes_preNEC
      ),
      shared_mutated_genomes_postNEC = ifelse(
        is.na(shared_mutated_genomes_postNEC),
        0,
        shared_mutated_genomes_postNEC
      ),
      Species = factor(Species, levels = species_order)
    )
  
  ## 6. Export the merged plotting table
  write.xlsx(plot_table,
             file = "D:/3_Projects/2_Children_jaundice/3_R_analysis/final_code/Figure 7_file/NEC_species_SNV_site_plot_table.xlsx",
             rowNames = FALSE)
  
  
  ####### Plot using the merged plot_table
  
  ## Plot 1: shared mutated genomes
  plot_shared <- plot_table %>%
    dplyr::select(
      Species,
      # shared_mutated_genomes_noNEC,
      shared_mutated_genomes_preNEC,
      shared_mutated_genomes_postNEC
    ) %>%
    pivot_longer(
      cols = c(#shared_mutated_genomes_noNEC,
        shared_mutated_genomes_preNEC, shared_mutated_genomes_postNEC),
      names_to = "NEC_stage",
      values_to = "n_genome_detections"
    ) %>%
    mutate(
      NEC_stage = factor(
        NEC_stage,
        levels = c(#"shared_mutated_genomes_postNEC",
          "shared_mutated_genomes_preNEC",
          "shared_mutated_genomes_noNEC"
          
        ),
        labels = c("postNEC", "preNEC"#,"noNEC"
        )
      )
    )
  
  p1_new <- ggplot(
    plot_shared,
    aes(x = Species, y = n_genome_detections, fill = NEC_stage)
  ) +
    geom_col(width = 0.8, color = NULL, alpha = 0.75) +
    coord_flip() +
    scale_fill_manual(
      values = c(#"noNEC"="red",
        "preNEC" = "#35978f",
        "postNEC" = "#bf812d"
      )
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks.y = element_line(),
      axis.title.y = element_blank()
    ) +
    labs(
      y = "Number of shared mutated genomes",
      x = NULL
    )
  p1_new
  
  ## Plot 2: tested SNV sites
  p2_new <- ggplot(
    plot_table,
    aes(x = Species, y = n_sites_tested)
  ) +
    geom_col(fill = "#77A88D", width = 0.8, alpha = 0.75) +
    geom_text(
      aes(label = comma(n_sites_tested)),
      hjust = -0.1,
      size = 3,
      color = "black"
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Number of tested SNV sites"
    )
  
  
  ## Plot 3: nominal preNEC-enriched SNV sites
  p3_new <- ggplot(
    plot_table,
    aes(x = Species, y = n_nominal_preNEC_enriched_sites)
  ) +
    geom_col(fill = "#E64B35", width = 0.8, alpha = 0.75) +
    geom_text(
      aes(label = n_nominal_preNEC_enriched_sites),
      hjust = -0.1,
      size = 3,
      color = "black"
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18))
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Number of nominal preNEC-enriched SNV sites"
    )
  
  
  ## Plot 4: nominal candidate rate
  p4_new <- ggplot(
    plot_table,
    aes(x = Species, y = nominal_rate)
  ) +
    geom_col(fill = "#3C5488", width = 0.8, alpha = 0.75) +
    geom_text(
      aes(label = percent(nominal_rate, accuracy = 0.01)),
      hjust = -0.1,
      size = 3,
      color = "black"
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.20))
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Nominal candidate rate"
    )
  
  ## Combine plots
  p_all <- (p1_new | p2_new | p3_new | p4_new) +
    plot_layout(widths = c(1, 1, 1, 1))
  
  p_all
  # save as 8.8 * 3.8, p_29species_mutatedgenome_snvsite_sigsnv_candidaterate 
  
  
  
  ##### Supplementary Fig. S9 Distribution of retained longitudinal samples per infant by NEC status #####
  
  ### Bacterial mutation sites ~ genes 
  
  ### Read candidate site and gene burden files
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
  
  ## Read files
  cand <- fread("D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv")
  gene_burden <- fread("D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_preNEC_enriched_gene_burden_nominal.tsv")
  
  dim(cand)
  dim(gene_burden)
  
  colnames(cand)
  colnames(gene_burden)
  
  ### Define priority species and format gene labels
  priority_species <- c("Klebsiella pneumoniae",
                        "Enterobacter hormaechei_A",
                        "Escherichia coli",
                        "Klebsiella michiganensis")
  
  make_gene_label <- function(gene, genome_sp) {
    gene <- as.character(gene)
    genome_sp <- as.character(genome_sp)
    
    if (is.na(gene) || gene == "" || gene == "intergenic_or_unannotated") {
      return("intergenic/unannotated")
    }
    
    prefix <- paste0(genome_sp, "_")
    
    if (!is.na(genome_sp) && startsWith(gene, prefix)) {
      return(paste0("gene_", substring(gene, nchar(prefix) + 1)))
    }
    
    return(gene)
  }
  
  make_gene_label_vec <- Vectorize(make_gene_label)
  
  
  
  ### Plot A: gene/region burden bubble plot for priority species
  # Demonstrate that candidate SNV sites are not randomly distributed, but concentrated in a few genes/regions of selected species
  gene_plot_df <- gene_burden %>%
    filter(Species %in% priority_species) %>%
    mutate(gene_short = make_gene_label_vec(gene, genome_sp),
           region_label = paste0(gene_short, " | ", genome_sp),
           has_N = grepl("(^|;)N($|;)", mutation_types),
           has_NS = ifelse(has_N, "Contains nonsynonymous", "No nonsynonymous")) %>%
    group_by(Species) %>%
    arrange(desc(n_candidate_sites), desc(mean_diff_prev), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    mutate(Species = factor(Species, levels = priority_species),
           region_label = fct_reorder(region_label, n_candidate_sites))
  
  ### Plot B: prevalence patterns of these genes/regions across noNEC, preNEC and postNEC
  # This plot is the most important. It supports the biological story:
  # These SNV clusters increase in the preNEC stage, and many are retained in postNEC samples.
  stage_prev_df <- gene_plot_df %>%
    dplyr::select(Species,
                  genome_sp,
                  gene,
                  gene_short,
                  region_label,
                  n_candidate_sites,
                  mean_noNEC_prev,
                  mean_preNEC_prev,
                  mean_postNEC_prev,
                  mean_diff_prev,
                  mutation_types) %>%
    pivot_longer(cols = c(mean_noNEC_prev, mean_preNEC_prev, mean_postNEC_prev),
                 names_to = "Stage",
                 values_to = "mean_prevalence" ) %>%
    mutate(Stage = case_when(
      Stage == "mean_noNEC_prev" ~ "noNEC",
      Stage == "mean_preNEC_prev" ~ "preNEC",
      Stage == "mean_postNEC_prev" ~ "postNEC"),
      Stage = factor(Stage, levels = c("noNEC", "preNEC", "postNEC")),
      region_label = factor(region_label, levels = levels(gene_plot_df$region_label)))
  
  
  
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
  
  ## 1. Align y-axis order
  gene_plot_df2 <- gene_plot_df %>%
    group_by(Species) %>%
    arrange(desc(n_candidate_sites), desc(mean_diff_prev), .by_group = TRUE) %>%
    ungroup()
  
  region_levels <- gene_plot_df2$region_label
  
  gene_plot_df2 <- gene_plot_df2 %>%
    mutate(region_label = factor(region_label, levels = region_levels),
           panel = "Summary",
           diff_label = paste0("Δ", percent(mean_diff_prev, accuracy = 1)))
  
  stage_prev_df2 <- stage_prev_df %>%
    mutate(region_label = factor(region_label, levels = region_levels) )
  
  ## Unified color scale: mean prevalence
  prev_limits <- c(0,max(c(gene_plot_df2$mean_preNEC_prev,
                           stage_prev_df2$mean_prevalence), 
                         na.rm = TRUE))
  
  shared_prev_color <- scale_color_gradient(
    low = "#FEE8C8",
    high = "#B2182B",
    limits = prev_limits,
    labels = percent_format(accuracy = 1),
    name = "Mean prevalence")
  
  ## Unified size scale: candidate SNV sites
  site_size_limits <- c(min(gene_plot_df2$n_candidate_sites, na.rm = TRUE),
                        max(gene_plot_df2$n_candidate_sites, na.rm = TRUE))
  
  shared_size_scale <- scale_size_continuous(
    range = c(2.5, 9),
    limits = site_size_limits,
    name = "Candidate\nSNV sites")
  
  
  ### Left plot: summary panel
  p_left <- ggplot(gene_plot_df2,aes(x = panel, y = region_label)) +
    geom_point(aes(size = n_candidate_sites,color = mean_preNEC_prev,
                   shape = has_NS),alpha = 0.9) +
    facet_grid(Species ~ .,scales = "free_y",space = "free_y") +
    shared_size_scale +
    shared_prev_color +
    scale_shape_manual(values = c("Contains nonsynonymous" = 16,"No nonsynonymous" = 17),
                       name = "Mutation type") +
    scale_x_discrete(expand = expansion(mult = c(0.25, 0.85))) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(size = 7, color = "black"),
          axis.title = element_text(color = "black"),
          strip.text.y = element_text(angle = 0, size = 8),
          strip.background = element_rect(fill = "grey95", color = "grey70"),
          legend.position = "right",
          plot.margin = margin(5.5, 20, 5.5, 5.5)) +
    labs(x = NULL,
         y = "Candidate gene / region")
  p_left
  
  
  ### Right plot: stage prevalence panel
  p_right <- ggplot(stage_prev_df2,aes(x = Stage, y = region_label)) +
    geom_line(aes(group = region_label),color = "grey75",
              linewidth = 0.35) +
    geom_point(aes(size = n_candidate_sites,
                   color = mean_prevalence,shape=has_NS),
               alpha = 0.95) +
    facet_grid(Species ~ .,
               scales = "free_y",
               space = "free_y") +
    shared_size_scale +
    shared_prev_color +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          strip.text.y = element_blank(),
          strip.background = element_blank(),
          legend.position = "right") +
    labs(x = NULL,
         y = NULL)
  p_right
  
  ### Combine plots
  p_combined <- p_left + p_right +
    plot_layout(widths = c(1.45, 2.2),guides = "collect") &
    theme(legend.position = "right")
  
  p_combined <- p_combined +
    plot_annotation(title = "Priority SNV-enriched genes/regions and their stage-specific prevalence patterns",
                    subtitle = "Color indicates mean prevalence; point size indicates candidate SNV-site burden; Δ indicates preNEC–noNEC prevalence difference")
  p_combined
  # save as 12 * 10, p_4bacs_allgenes_prevalence
  
  
  
  
  
  
  
  ##### Fig. 7d. Genome-level distribution of nominal preNEC-enriched candidate SNV sites in priority species ####
  
  ### Four bacterial genome-level candidate SNV burden plot
  # Klebsiella pneumoniae, Enterobacter hormaechei_A, Escherichia coli, Klebsiella michiganensis
  
  selected_genomes <- c("GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1",
                        "GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1",
                        "CokerMO_2019_SRR8692237_bin.9_k141",
                        "ShaoY_2019_ERR3405384_bin.5_k141")
  
  ### Format genome labels
  make_genome_label <- function(x) {
    x <- as.character(x)
    
    ## GCF reference genome: retain GCF accession + scaffold/accession
    if (grepl("^GCF_", x)) {
      acc <- sub("^(GCF_[0-9\\.]+).*", "\\1", x)
      scaffold <- sub("^.*_genomic_", "", x)
      scaffold <- sub("^NZ_", "", scaffold)
      return(paste0(acc, " | ", scaffold))
    }
    
    ## Metagenomic bin: remove k-mer suffix
    x2 <- sub("_k[0-9]+$", "", x)
    return(x2)
  }
  
  make_genome_label_vec <- Vectorize(make_genome_label)
  
  ### Count candidate SNV sites in each genome
  # This table summarizes:
  # Total number of candidate sites in each genome
  # Number of nonsynonymous / synonymous / intergenic sites
  # Number of involved genes
  # Number of postNEC retained sites
  # Mean preNEC prevalence
  # Mean preNEC-noNEC difference
  
  
  length(unique(cand$genome_sp))
  genome_summary <- cand %>%
    filter(Species %in% priority_species) %>%
    mutate(
      Species = factor(Species, levels = priority_species),
      mutation_type = ifelse(is.na(mutation_type) | mutation_type == "", "unknown", mutation_type),
      mutation_type = factor(mutation_type, levels = c("I", "S", "N", "unknown")),
      genome_label = make_genome_label_vec(genome_sp),
      selected_genome = genome_sp %in% selected_genomes
    ) %>%
    group_by(Species, genome_sp, genome_label) %>%
    summarise(
      n_candidate_sites = n_distinct(site_id),
      n_N = n_distinct(site_id[mutation_type == "N"]),
      n_S = n_distinct(site_id[mutation_type == "S"]),
      n_I = n_distinct(site_id[mutation_type == "I"]),
      n_unknown = n_distinct(site_id[mutation_type == "unknown"]),
      n_genes = n_distinct(gene[!is.na(gene) & gene != ""]),
      n_postNEC_retained_sites = n_distinct(site_id[postNEC_prev > 0]),
      postNEC_retained_fraction = mean(postNEC_prev > 0, na.rm = TRUE),
      mean_noNEC_prev = mean(noNEC_prev, na.rm = TRUE),
      mean_preNEC_prev = mean(preNEC_prev, na.rm = TRUE),
      mean_postNEC_prev = mean(postNEC_prev, na.rm = TRUE),
      mean_diff_prev = mean(diff_prev_preNEC_vs_noNEC, na.rm = TRUE),
      min_fisher_p = min(fisher_p, na.rm = TRUE),
      selected_genome = any(selected_genome),
      .groups = "drop"
    ) %>%
    arrange(Species, desc(n_N), desc(n_candidate_sites), desc(mean_diff_prev))
  genome_summary
  
  #fwrite(genome_summary,"Genome_level_candidate_SNV_site_summary_priority_species.tsv",
  #  sep = "\t")
  
  
  
  colnames(genome_summary)
  length(unique(genome_summary$genome_sp))
  ### Prepare mutation-type stacked barplot data
  top_n_genomes_per_species <- 20
  
  genomes_to_plot <- genome_summary %>%
    group_by(Species) %>%
    arrange(desc(n_candidate_sites), desc(n_N), .by_group = TRUE) %>%
    mutate(rank_by_sites = row_number()) %>%
    filter(rank_by_sites <= top_n_genomes_per_species | selected_genome) %>%
    ungroup()
  
  genome_type_count <- cand %>%
    filter(Species %in% priority_species) %>%
    mutate(Species = factor(Species, levels = priority_species),
           mutation_type = ifelse(is.na(mutation_type) | mutation_type == "", "unknown", mutation_type),
           mutation_type = factor(mutation_type, levels = c("I", "S", "N", "unknown")),
           genome_label = make_genome_label_vec(genome_sp),
           selected_genome = genome_sp %in% selected_genomes) %>%
    semi_join(genomes_to_plot %>% dplyr::select(Species, genome_sp),
              by = c("Species", "genome_sp")) %>%
    distinct(Species, genome_sp, genome_label, selected_genome, mutation_type, site_id) %>%
    count(Species, genome_sp, genome_label, selected_genome, mutation_type, name = "n_sites") %>%
    left_join(genome_summary %>%
                dplyr::select(Species, genome_sp, n_candidate_sites, n_N,
                              mean_preNEC_prev, mean_diff_prev, postNEC_retained_fraction),
              by = c("Species", "genome_sp"))
  
  ## Set y-axis order
  genome_order_df <- genome_type_count %>%
    distinct(Species, genome_sp, genome_label, selected_genome, n_candidate_sites, n_N, mean_diff_prev) %>%
    arrange(Species, n_candidate_sites, n_N, mean_diff_prev) %>%
    mutate(genome_plot_id = paste(Species, genome_sp, sep = "___"),
           genome_label_plot = ifelse(selected_genome, paste0(genome_sp, "  ★"), genome_sp))
  
  genome_levels <- genome_order_df$genome_plot_id
  
  genome_type_count_plot <- genome_type_count %>%
    mutate(genome_plot_id = paste(Species, genome_sp, sep = "___"),
           genome_plot_id = factor(genome_plot_id, levels = genome_levels),
           selected_label = ifelse(selected_genome, "Selected for downstream", "Other candidate genome"))
  table(genome_type_count_plot$genome_sp,genome_type_count_plot$Species)
  
  
  # Plot: candidate-site counts per genome, stratified by mutation type
  p_genome_candidate_bar <- ggplot(genome_type_count_plot,
                                   aes(x = genome_plot_id, y = n_sites, fill = mutation_type)) +
    geom_col(aes(color = selected_label),width = 0.78,linewidth = 0.25) +
    geom_text(data = genome_order_df %>%
                mutate(genome_plot_id = factor(genome_plot_id, levels = genome_levels)),
              aes(x = genome_plot_id,y = n_candidate_sites,label = n_candidate_sites),
              inherit.aes = FALSE,hjust = -0.15,size = 2.8,color = "black") +
    coord_flip() +
    facet_grid(Species ~ .,
               scales = "free_y",
               space = "free_y") +
    #scale_x_discrete(labels = function(x) {
    #    full_lab <- sub("^.*___", "", x)
    #    selected_full_labs <- genome_order_df$genome_sp[genome_order_df$selected_genome]
    #    ifelse(full_lab %in% selected_full_labs, paste0(full_lab, "  ★"), full_lab)
    #  }
    #) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(values = c("I" = "#b3de69",
                                 "S" = "#8dd3c7",
                                 "N" = "#fdb462",
                                 "unknown" = "grey70"),
                      name = "Mutation type") +
    scale_color_manual(values = c("Selected for downstream" = "black",
                                  "Other candidate genome" = "grey80"),
                       name = "") +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0, size = 8),
          strip.background = element_rect(fill = "grey95", color = "grey70"),
          axis.text.y = element_text(size = 6.5, color = "black"),
          axis.text.x = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          legend.position = "right") +
    labs(x = NULL,
         y = "Number of nominal preNEC-enriched candidate SNV sites",
         title = "Genome-level distribution of nominal preNEC-enriched candidate SNV sites",
         subtitle = "Bars show candidate SNV-site counts within each genome; stars indicate genomes selected for downstream nonsynonymous-site analysis")
  p_genome_candidate_bar
  # save as 6 * 11
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##### Supplementary Fig. S10. Linked SNV modules are enriched before NEC onset and identify nonsynonymous candidate loci ##########
  
  ### Bacterial mutation sites 
  # Four bacteria were selected from previous results: Klebsiella pneumoniae, Enterobacter hormaechei_A, Escherichia coli and Klebsiella michiganensis
  
  ##### Modularity of mutation sites in 4 selected bacteria and 8 genomes 
  # Priority 1: Klebsiella pneumoniae
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1
  # VatanenT_2016_SRR4408069_bin.6_k119
  # Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141
  # ShaoY_2019_ERR3404951_bin.2_k141
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020838.1
  # Priority 2: Enterobacter hormaechei_A
  # GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1
  # Priority 3: Escherichia coli
  # CokerMO_2019_SRR8692237_bin.9_k141
  # Priority 4: Klebsiella michiganensis
  # ShaoY_2019_ERR3405384_bin.5_k141_6236_14
  
  # Explore linked haplotypes of these 8 genomes across 692 samples
  
  ## Linked haplotype-like analysis for 8 priority genomes
  ## Purpose: Test whether candidate SNV sites in 8 selected genomes
  ##   form tightly co-occurring haplotype-like modules.
  
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(igraph)
  
  
  
  ## 1. Read candidate sites and metadata
  cand <- fread("D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed/ALL_species_SNV_site_NEC_association_preNEC_enriched_candidate_sites_nominal.tsv")
  
  NEC_meta_infor_sub2#
  
  cat("Candidate site table dimension:\n")
  print(dim(cand))
  
  cat("Metadata dimension:\n")
  print(dim(NEC_meta_infor_sub2))
  
  ## 3. Check required columns
  required_cand_cols <- c("Species",
                          "genome_sp",
                          "site_id",
                          "scaffold",
                          "position",
                          "gene",
                          "mutation",
                          "mutation_type",
                          "noNEC_prev",
                          "preNEC_prev",
                          "postNEC_prev",
                          "diff_prev_preNEC_vs_noNEC",
                          "fisher_p",
                          "fisher_FDR_species")
  
  missing_cand_cols <- setdiff(required_cand_cols, colnames(cand))
  
  if (length(missing_cand_cols) > 0) {
    stop(
      "cand is missing required columns: ",
      paste(missing_cand_cols, collapse = ", ")
    )
  }
  
  required_meta_cols <- c(
    "sample2",
    "SubjectID",
    "PreNEC2",
    "DOL",
    "Study"
  )
  
  missing_meta_cols <- setdiff(required_meta_cols, colnames(NEC_meta_infor_sub2))
  
  if (length(missing_meta_cols) > 0) {
    stop(
      "NEC_meta_infor_sub2 is missing required columns: ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
  
  
  ## 4. Define priority species and 8 target genomes
  priority_species <- c("Klebsiella pneumoniae",
                        "Enterobacter hormaechei_A",
                        "Escherichia coli",
                        "Klebsiella michiganensis")
  
  target_genomes_raw <- c(
    ## Priority 1: Klebsiella pneumoniae
    "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1",
    "VatanenT_2016_SRR4408069_bin.6_k119",
    "Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141",
    "ShaoY_2019_ERR3404951_bin.2_k141",
    "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020838.1",
    
    ## Priority 2: Enterobacter hormaechei_A
    "GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1",
    
    ## Priority 3: Escherichia coli
    "CokerMO_2019_SRR8692237_bin.9_k141",
    
    ## Priority 4: Klebsiella michiganensis
    "ShaoY_2019_ERR3405384_bin.5_k141"
  )
  
  ## Note:
  ## ShaoY_2019_ERR3405384_bin.5_k141_6236_14 may be a gene ID.
  ## If it is not found in cand$genome_sp, the code will try to automatically replace it with:
  ## ShaoY_2019_ERR3405384_bin.5_k141
  
  target_genomes <- target_genomes_raw
  
  for (i in seq_along(target_genomes)) {
    
    g <- target_genomes[i]
    
    if (!g %in% unique(cand$genome_sp)) {
      
      fallback_g <- sub("(_k[0-9]+).*", "\\1", g)
      
      if (fallback_g %in% unique(cand$genome_sp)) {
        message("Replace genome ID: ", g, " -> ", fallback_g)
        target_genomes[i] <- fallback_g
      }
    }
  }
  
  missing_target_genomes <- setdiff(target_genomes, unique(cand$genome_sp))
  
  if (length(missing_target_genomes) > 0) {
    
    cat("\nThese target genomes are NOT found in cand$genome_sp:\n")
    print(missing_target_genomes)
    
    cat("\nPossible matches:\n")
    for (g in missing_target_genomes) {
      pattern <- sub("_k[0-9]+.*$", "", g)
      print(grep(pattern, unique(cand$genome_sp), value = TRUE))
    }
    
    stop("Please fix target_genomes before continuing.")
  }
  
  cat("\nFinal target genomes:\n")
  print(target_genomes)
  
  
  ## 5. Prepare selected candidate sites
  ## Important:
  ## For linked haplotype discovery, it is recommended to use all I/S/N candidate sites.
  ## This is because haplotype / strain background is usually marked by a group of I/S/N sites.
  ## For downstream protein and 3D analyses, N sites can be selected later from strong linked modules.
  
  colnames(cand)
  #View(cand)
  selected_sites <- cand %>%
    filter(Species %in% priority_species,
           genome_sp %in% target_genomes) %>%
    mutate(Species = as.character(Species),
           genome_sp = as.character(genome_sp),
           scaffold = as.character(scaffold),
           position = as.numeric(position),
           site_id = as.character(site_id),
           gene = ifelse(is.na(gene) | gene == "","intergenic_or_unannotated",as.character(gene)),
           mutation = as.character(mutation),
           mutation_type = ifelse(is.na(mutation_type) | mutation_type == "","unknown",
                                  as.character(mutation_type)))
  
  cat("\nSelected candidate sites:\n")
  print(dim(selected_sites))
  
  cat("\nSelected sites by species:\n")
  print(table(selected_sites$Species))
  
  cat("\nSelected sites by genome:\n")
  print(table(selected_sites$genome_sp))
  
  cat("\nSelected sites by mutation type:\n")
  print(table(selected_sites$mutation_type))
  
  colnames(selected_sites)
  
  
  
  ## Local Windows path for the species file directory
  species_dir <- "D:/3_Projects/2_Children_jaundice/3_R_analysis/0_backup/0002_backup2_有修改的部分/species_split_named_uncompressed"
  
  ## Output results are also saved in this folder
  out_dir <- species_dir
  
  if (!dir.exists(species_dir)) {
    stop("species_dir does not exist: ", species_dir)
  }
  
  setwd(out_dir)
  
  cat("Species directory:\n")
  cat(species_dir, "\n\n")
  
  cat("Number of species files:\n")
  print(length(list.files(species_dir, pattern = "\\.tsv$", full.names = TRUE)))
  
  
  ## 6. Build species file map
  safe_species_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_]+", "_", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  species_files <- data.table(
    file = list.files(
      species_dir,
      pattern = "\\.tsv$",
      full.names = TRUE
    )) %>%
    mutate(
      safe_species = sub("\\.tsv$", "", basename(file))
    )
  
  target_species_file_map <- selected_sites %>%
    distinct(Species) %>%
    mutate(
      safe_species = safe_species_name(Species)
    ) %>%
    left_join(
      species_files,
      by = "safe_species"
    )
  
  cat("\nSpecies file map:\n")
  print(target_species_file_map)
  
  if (any(is.na(target_species_file_map$file))) {
    
    cat("\nMissing species files:\n")
    print(target_species_file_map %>% filter(is.na(file)))
    
    cat("\nAvailable species files:\n")
    print(head(species_files$safe_species, 100))
    
    stop("Some species files are missing. Check species_dir or file naming.")
  }
  
  
  
  ## 7. Prepare metadata
  meta_stage <- NEC_meta_infor_sub2 %>%
    mutate(
      NEC_stage = case_when(
        PreNEC2 == "noNEC" ~ "noNEC",
        PreNEC2 == "preNEC" ~ "preNEC",
        PreNEC2 %in% c("postNEC", "NEC", "onsetNEC") ~ "postNEC",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      sample2,
      SubjectID,
      NEC_stage,
      DOL,
      Study
    ) %>%
    distinct() %>%
    mutate(
      sample2 = as.character(sample2),
      SubjectID = as.character(SubjectID),
      NEC_stage = as.character(NEC_stage),
      DOL = as.numeric(DOL),
      Study = as.character(Study)
    )
  
  cat("\nNEC stage distribution:\n")
  print(table(meta_stage$NEC_stage, useNA = "ifany"))
  
  
  ## 8. Helper functions
  safe_first_char <- function(x) {
    
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    
    if (length(x) == 0) {
      return(NA_character_)
    } else {
      return(x[1])
    }
  }
  
  get_species_file <- function(sp, species_file_map) {
    
    hit <- species_file_map %>%
      filter(Species == sp) %>%
      pull(file)
    
    if (length(hit) == 0 || is.na(hit[1]) || !file.exists(hit[1])) {
      stop("Cannot find species file for: ", sp)
    }
    
    hit[1]
  }
  
  
  ## 9. Build one genome × candidate-site matrix
  build_selected_genome_site_matrix <- function(
    sp,
    genome_id,
    selected_sites,
    meta_stage,
    species_file_map,
    unit = c("sample", "subject_stage"),
    min_cov = 5,
    min_alt_count = 2,
    min_alt_freq = 0.05
  ) {
    
    unit <- match.arg(unit)
    
    sp_file <- get_species_file(
      sp = sp,
      species_file_map = species_file_map
    )
    
    message("Reading species file: ", sp_file)
    
    dt <- fread(sp_file)
    
    required_cols <- c("sample2",
                       "genome_sp",
                       "scaffold",
                       "position",
                       "A",
                       "C",
                       "T",
                       "G")
    
    miss <- setdiff(required_cols, colnames(dt))
    
    if (length(miss) > 0) {
      stop("Species file missing columns: ", paste(miss, collapse = ", "))
    }
    
    dt[, sample2 := as.character(sample2)]
    dt[, genome_sp := as.character(genome_sp)]
    dt[, scaffold := as.character(scaffold)]
    dt[, position := as.numeric(position)]
    
    ## Retain only the current genome
    dt <- dt[genome_sp == genome_id]
    
    if (nrow(dt) == 0) {
      message("No raw SNV records for genome: ", genome_id)
      return(NULL)
    }
    
    site_keep <- selected_sites %>%
      filter(
        Species == sp,
        genome_sp == genome_id
      ) %>%
      dplyr::select(
        site_id,
        genome_sp,
        scaffold,
        position,
        gene,
        mutation,
        mutation_type
      ) %>%
      distinct()
    
    if (nrow(site_keep) == 0) {
      message("No selected candidate sites for genome: ", genome_id)
      return(NULL)
    }
    
    site_keep_dt <- as.data.table(site_keep)
    site_keep_dt[, position := as.numeric(position)]
    
    ## Retain only selected candidate sites from the current genome
    dt <- merge(
      dt,
      site_keep_dt[, .(site_id, genome_sp, scaffold, position)],
      by = c("genome_sp", "scaffold", "position"),
      all = FALSE
    )
    
    if (nrow(dt) == 0) {
      message("Selected candidate sites not found in raw SNV file for genome: ", genome_id)
      return(NULL)
    }
    
    ## merge metadata
    dt <- dt %>%left_join(meta_stage,by = "sample2")
    setDT(dt)
    
    dt <- dt[NEC_stage %in% c("noNEC", "preNEC", "postNEC") &!is.na(SubjectID)]
    
    if (nrow(dt) == 0) {
      return(NULL)
    }
    
    ## Force numeric type
    dt[, A := as.numeric(A)]
    dt[, C := as.numeric(C)]
    dt[, T := as.numeric(T)]
    dt[, G := as.numeric(G)]
    
    ## Calculate the global major allele for each site
    site_major <- dt[,.(A_total = sum(A, na.rm = TRUE),
                        C_total = sum(C, na.rm = TRUE),
                        T_total = sum(T, na.rm = TRUE),
                        G_total = sum(G, na.rm = TRUE)),by = site_id]
    
    site_major[,major_base := c("A", "C", "T", "G")[
      max.col(as.matrix(.SD),ties.method = "first") ],
      .SDcols = c("A_total", "C_total", "T_total", "G_total")]
    
    dt <- merge(dt,site_major[, .(site_id, major_base)],
                by = "site_id",all.x = TRUE)
    
    dt[, total_count := A + C + T + G]
    
    dt[,major_count := fcase(major_base == "A", A,
                             major_base == "C", C,
                             major_base == "T", T,
                             major_base == "G", G,
                             default = NA_real_)]
    
    dt[, alt_count := total_count - major_count]
    dt[, alt_freq := alt_count / total_count]
    
    dt[,alt_present := as.integer(
      total_count >= min_cov &
        alt_count >= min_alt_count &
        alt_freq >= min_alt_freq)]
    
    if (unit == "sample") {
      
      dt[, row_id := sample2]
      
      row_meta <- dt[,.(sample2 = safe_first_char(sample2),
                        SubjectID = safe_first_char(SubjectID),
                        NEC_stage = safe_first_char(NEC_stage),
                        DOL = median(DOL, na.rm = TRUE),
                        Study = safe_first_char(Study)),by = row_id ]
      
    } else {
      
      dt[, row_id := paste(SubjectID, NEC_stage, sep = "__")]
      
      row_meta <- dt[,.(sample2 = paste(sort(unique(sample2)), collapse = ";"),
                        SubjectID = safe_first_char(SubjectID),
                        NEC_stage = safe_first_char(NEC_stage),
                        DOL = median(DOL, na.rm = TRUE),
                        Study = safe_first_char(Study)),by = row_id]
    }
    
    row_site <- dt[, .(alt_present = as.integer(any(alt_present == 1, na.rm = TRUE)),
                       alt_freq = max(alt_freq, na.rm = TRUE)),by = .(row_id, site_id)]
    
    wide_bin <- dcast(row_site,row_id ~ site_id,
                      value.var = "alt_present",fill = 0 )
    
    wide_af <- dcast(row_site,row_id ~ site_id,value.var = "alt_freq",fill = 0)
    
    row_meta <- row_meta[match(wide_bin$row_id, row_meta$row_id)]
    
    mat_bin <- as.matrix(wide_bin[, -1, with = FALSE])
    rownames(mat_bin) <- wide_bin$row_id
    
    mat_af <- as.matrix(wide_af[, -1, with = FALSE])
    rownames(mat_af) <- wide_af$row_id
    
    site_anno <- selected_sites %>%
      filter(Species == sp,genome_sp == genome_id) %>%
      distinct(Species,
               genome_sp,
               site_id,
               scaffold,
               position,
               gene,
               mutation,
               mutation_type,
               noNEC_prev,
               preNEC_prev,
               postNEC_prev,
               diff_prev_preNEC_vs_noNEC,
               fisher_p,
               fisher_FDR_species)
    
    list(Species = sp,
         genome_sp = genome_id,
         unit = unit,
         meta = as.data.frame(row_meta),
         mat_bin = mat_bin,
         mat_af = mat_af,
         site_anno = site_anno)
  }
  
  
  ## 10. Pairwise site co-occurrence
  calc_site_pair_cooccurrence <- function(obj) {
    
    mat <- obj$mat_bin
    
    keep <- colSums(mat, na.rm = TRUE) > 0
    mat <- mat[, keep, drop = FALSE]
    
    site_ids <- colnames(mat)
    
    if (ncol(mat) < 2) {
      return(NULL)
    }
    
    pairs <- combn(site_ids,2,simplify = FALSE)
    
    res <- lapply(pairs,
                  function(pp) {
                    
                    s1 <- pp[1]
                    s2 <- pp[2]
                    
                    x <- mat[, s1]
                    y <- mat[, s2]
                    
                    n11 <- sum(x == 1 & y == 1, na.rm = TRUE)
                    n10 <- sum(x == 1 & y == 0, na.rm = TRUE)
                    n01 <- sum(x == 0 & y == 1, na.rm = TRUE)
                    n00 <- sum(x == 0 & y == 0, na.rm = TRUE)
                    
                    jaccard <- ifelse((n11 + n10 + n01) > 0, n11 / (n11 + n10 + n01),NA_real_)
                    
                    phi <- suppressWarnings(cor(x, y, method = "pearson"))
                    
                    ft <- tryCatch(fisher.test(matrix(c(n11, n10, n01, n00),nrow = 2))$p.value,
                                   error = function(e) NA_real_)
                    
                    data.table(Species = obj$Species,
                               genome_sp = obj$genome_sp,
                               unit = obj$unit,
                               site_A = s1,
                               site_B = s2,
                               n_both = n11,
                               n_A_only = n10,
                               n_B_only = n01,
                               n_neither = n00,
                               jaccard = jaccard,
                               phi = phi,
                               fisher_cooccur_p = ft)
                  }
    )
    
    res <- rbindlist(res,fill = TRUE)
    
    res[,fisher_cooccur_FDR := p.adjust(fisher_cooccur_p,method = "BH")]
    
    anno <- as.data.table(obj$site_anno)
    
    anno_A <- anno %>%
      dplyr::select(site_A = site_id,
                    scaffold_A = scaffold,
                    position_A = position,
                    gene_A = gene,
                    mutation_A = mutation,
                    mutation_type_A = mutation_type)
    
    anno_B <- anno %>%
      dplyr::select(site_B = site_id,
                    scaffold_B = scaffold,
                    position_B = position,
                    gene_B = gene,
                    mutation_B = mutation,
                    mutation_type_B = mutation_type)
    
    res <- res %>%
      left_join(anno_A,by = "site_A" ) %>%
      left_join(anno_B,by = "site_B") %>%
      mutate(same_scaffold = scaffold_A == scaffold_B,
             same_gene = gene_A == gene_B,
             distance_bp = ifelse(same_scaffold,abs(as.numeric(position_A) - as.numeric(position_B)),
                                  NA_real_),
             both_N = mutation_type_A == "N" & mutation_type_B == "N",
             any_N = mutation_type_A == "N" | mutation_type_B == "N")
    
    as.data.table(res)
  }
  
  
  ## 11. Run 8 target genomes
  target_genome_map <- selected_sites %>%
    filter(genome_sp %in% target_genomes) %>%
    distinct(Species, genome_sp) %>%
    arrange(Species, genome_sp)
  
  cat("\nTarget genome map:\n")
  print(target_genome_map)
  
  matrix_objects <- list()
  pair_results <- list()
  
  for (i in seq_len(nrow(target_genome_map))) {
    
    sp <- target_genome_map$Species[i]
    genome_id <- target_genome_map$genome_sp[i]
    
    message("\n============================")
    message("Running: ", sp, " | ", genome_id)
    message("============================")
    
    obj <- tryCatch(
      build_selected_genome_site_matrix(
        sp = sp,
        genome_id = genome_id,
        selected_sites = selected_sites,
        meta_stage = meta_stage,
        species_file_map = target_species_file_map,
        unit = "sample",
        min_cov = 5,
        min_alt_count = 2,
        min_alt_freq = 0.05
      ),
      error = function(e) {
        message("Failed at matrix building: ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(obj)) {
      next
    }
    
    key <- paste(sp, genome_id, sep = " | ")
    matrix_objects[[key]] <- obj
    
    pair_res <- tryCatch(
      calc_site_pair_cooccurrence(obj),
      error = function(e) {
        message("Failed at pairwise co-occurrence: ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(pair_res) && nrow(pair_res) > 0) {
      pair_results[[key]] <- pair_res
    }
  }
  
  
  
  if (length(pair_results) == 0) {
    stop("No pairwise co-occurrence results generated.")
  }
  
  cooccur_all <- rbindlist(
    pair_results,
    fill = TRUE
  )
  
  #
  cat("\nFinished pairwise co-occurrence analysis.\n")
  cat("Rows in cooccur_all:\n")
  print(nrow(cooccur_all))
  
  cat("\nNumber of site pairs per genome:\n")
  cooccur_all %>%
    dplyr::count(Species, genome_sp, name = "n_site_pairs") %>%
    dplyr::arrange(desc(n_site_pairs)) %>%
    tibble::as_tibble() %>%
    head(n = 50)
  
  
  ## 12. Define tightly linked candidate pairs
  linked_pairs <- cooccur_all %>%
    filter(!is.na(phi),
           !is.na(jaccard),
           phi >= 0.70,
           jaccard >= 0.60,
           n_both >= 3,
           fisher_cooccur_FDR < 0.05) %>%
    arrange(Species,
            genome_sp,
            desc(phi),
            desc(jaccard),
            distance_bp)
  
  #fwrite(linked_pairs, "Tightly_linked_candidate_SNV_pairs_phi0.7_jaccard0.6.tsv", sep = "\t")
  
  cat("\nTightly linked pairs summary:\n")
  linked_pairs %>%
    dplyr::count(Species, genome_sp, name = "n_linked_pairs") %>%
    dplyr::arrange(desc(n_linked_pairs)) %>%
    head(n = 50)
  
  
  ## 13. Relaxed linked pairs
  linked_pairs_relaxed <- cooccur_all %>%
    filter(!is.na(phi),
           !is.na(jaccard),
           phi >= 0.60,
           jaccard >= 0.50,
           n_both >= 3,
           fisher_cooccur_FDR < 0.05) %>%
    arrange( Species,
             genome_sp,
             desc(phi),
             desc(jaccard),
             distance_bp)
  
  #fwrite(linked_pairs_relaxed,"Tightly_linked_candidate_SNV_pairs_phi0.6_jaccard0.5_relaxed.tsv",sep = "\t")
  
  cat("\nRelaxed linked pairs summary:\n")
  linked_pairs_relaxed %>%
    dplyr::count(Species, genome_sp, name = "n_linked_pairs") %>%
    dplyr::arrange(desc(n_linked_pairs)) %>%
    head(n = 50)
  
  
  ## 14. Choose pair set for module detection
  if (nrow(linked_pairs) > 0) {
    
    linked_pairs_for_module <- linked_pairs
    module_threshold_label <- "strict_phi0.7_jaccard0.6"
    
  } else {
    
    linked_pairs_for_module <- linked_pairs_relaxed
    module_threshold_label <- "relaxed_phi0.6_jaccard0.5"
  }
  
  if (nrow(linked_pairs_for_module) == 0) {
    stop("No linked pairs found even under relaxed threshold.")
  }
  
  cat("\nModule detection uses threshold:\n")
  print(module_threshold_label)
  
  
  ## 15. Identify haplotype-like modules
  identify_linked_modules <- function(linked_pairs_one) {
    
    if (nrow(linked_pairs_one) == 0) {
      return(NULL)
    }
    
    g <- graph_from_data_frame(linked_pairs_one %>%
                                 dplyr::select(site_A, site_B),directed = FALSE)
    
    comp <- components(g)
    
    data.frame(site_id = names(comp$membership),
               module_id = paste0("Module_", comp$membership),
               stringsAsFactors = FALSE)
  }
  
  module_list <- linked_pairs_for_module %>%
    group_by(Species, genome_sp) %>%
    group_split()
  
  module_sites_all <- lapply(
    module_list,
    function(df) {
      
      tmp <- identify_linked_modules(df)
      
      if (is.null(tmp)) {
        return(NULL)
      }
      
      tmp$Species <- unique(df$Species)
      tmp$genome_sp <- unique(df$genome_sp)
      
      tmp
    }
  ) %>%
    bind_rows() %>%
    dplyr::select(
      Species,
      genome_sp,
      module_id,
      site_id)
  
  #fwrite(module_sites_all,"Tightly_linked_haplotype_like_modules_sites.tsv",sep = "\t")
  
  cat("\nModule site table:\n")
  
  
  ## 16. Module annotation and summary
  module_site_anno <- module_sites_all %>%
    left_join(selected_sites %>%
                dplyr::select(Species,
                              genome_sp,
                              site_id,
                              scaffold,
                              position,
                              gene,
                              mutation,
                              mutation_type,
                              noNEC_prev,
                              preNEC_prev,
                              postNEC_prev,
                              diff_prev_preNEC_vs_noNEC,
                              fisher_p,
                              fisher_FDR_species),
              by = c("Species", "genome_sp", "site_id"))
  
  module_summary <- module_site_anno %>%
    group_by(Species, genome_sp, module_id) %>%
    summarise( n_sites = n_distinct(site_id),
               n_N = n_distinct(site_id[mutation_type == "N"]),
               n_S = n_distinct(site_id[mutation_type == "S"]),
               n_I = n_distinct(site_id[mutation_type == "I"]),
               n_unknown = n_distinct(site_id[mutation_type == "unknown"]),
               n_genes = n_distinct(gene),
               genes = paste(sort(unique(gene)), collapse = ";"),
               mutations = paste(sort(unique(mutation)), collapse = ";"),
               mutation_types = paste(sort(unique(mutation_type)), collapse = ";"),
               mean_noNEC_prev = mean(noNEC_prev, na.rm = TRUE),
               mean_preNEC_prev = mean(preNEC_prev, na.rm = TRUE),
               mean_postNEC_prev = mean(postNEC_prev, na.rm = TRUE),
               mean_diff_prev = mean(diff_prev_preNEC_vs_noNEC, na.rm = TRUE),
               n_postNEC_retained_sites = sum(postNEC_prev > 0, na.rm = TRUE),
               postNEC_retained_fraction = mean(postNEC_prev > 0, na.rm = TRUE),
               min_fisher_p = min(fisher_p, na.rm = TRUE),
               .groups = "drop" ) %>%
    arrange(Species,
            genome_sp,
            desc(n_N),
            desc(n_sites),
            desc(mean_diff_prev))
  
  fwrite(module_site_anno,"Tightly_linked_haplotype_like_module_site_annotation.tsv",sep = "\t")
  fwrite(module_summary,"D:/3_Projects/2_Children_jaundice/3_R_analysis/Tightly_linked_haplotype_like_module_summary.tsv",sep = "\t")
  write_xlsx(module_summary, path = "D:/3_Projects/2_Children_jaundice/3_R_analysis/Tightly_linked_haplotype_like_module_summary.xlsx")
  
  
  cat("\nModule summary:\n")
  print(module_summary)
  
  
  ## 17. Calculate module score by sample
  calc_module_score_for_obj <- function(obj, module_sites_all) {
    
    mat <- obj$mat_bin
    meta <- obj$meta
    
    modules <- module_sites_all %>%
      filter(Species == obj$Species,
             genome_sp == obj$genome_sp)
    
    if (nrow(modules) == 0) {
      return(NULL)
    }
    
    out <- list()
    
    for (mid in unique(modules$module_id)) {
      
      sites <- modules %>%
        filter(module_id == mid) %>%
        pull(site_id)
      
      sites <- intersect(sites, colnames(mat))
      
      if (length(sites) == 0) {
        next
      }
      
      score <- rowMeans(mat[, sites, drop = FALSE],na.rm = TRUE)
      
      tmp <- data.frame(Species = obj$Species,
                        genome_sp = obj$genome_sp,
                        module_id = mid,
                        row_id = rownames(mat),
                        module_score = score,
                        n_module_sites_in_matrix = length(sites)) %>%
        left_join(meta,by = "row_id")
      
      out[[mid]] <- tmp
    }
    
    bind_rows(out)
  }
  
  module_score_all <- bind_rows(lapply(matrix_objects,
                                       calc_module_score_for_obj,
                                       module_sites_all = module_sites_all))
  
  module_score_all$NEC_stage <- factor(module_score_all$NEC_stage,
                                       levels = c("noNEC", "preNEC", "postNEC"))
  
  fwrite(module_score_all,"Tightly_linked_haplotype_like_module_score_by_sample.tsv",sep = "\t")
  
  
  ## 18. Module score statistics: noNEC vs preNEC
  module_score_stats <- module_score_all %>%
    filter(NEC_stage %in% c("noNEC", "preNEC")) %>%
    group_by(Species,
             genome_sp,
             module_id) %>%
    summarise(n_noNEC = sum(NEC_stage == "noNEC"),
              n_preNEC = sum(NEC_stage == "preNEC"),
              median_noNEC = median(
                module_score[NEC_stage == "noNEC"],
                na.rm = TRUE),
              median_preNEC = median(
                module_score[NEC_stage == "preNEC"],
                na.rm = TRUE),
              diff_median_preNEC_vs_noNEC = median_preNEC - median_noNEC,
              p_wilcox = tryCatch(wilcox.test(module_score ~ NEC_stage)$p.value,
                                  error = function(e) NA_real_),.groups = "drop") %>%
    mutate(FDR_wilcox = p.adjust(p_wilcox,method = "BH")) %>%
    arrange(p_wilcox)
  
  fwrite(module_score_stats,"Tightly_linked_haplotype_like_module_score_stats.tsv",sep = "\t")
  
  cat("\nModule score statistics:\n")
  print(module_score_stats)
  
  
  ## 19. Select strong linked modules
  strong_modules <- module_summary %>%
    left_join(module_score_stats %>%
                dplyr::select(Species,
                              genome_sp,
                              module_id,
                              median_noNEC,
                              median_preNEC,
                              diff_median_preNEC_vs_noNEC,
                              p_wilcox,
                              FDR_wilcox),
              by = c("Species", "genome_sp", "module_id")) %>%
    filter(n_sites >= 3,
           n_N >= 1,
           mean_diff_prev >= 0.20,# 0.2
           postNEC_retained_fraction >= 0.50, # 0.5
           diff_median_preNEC_vs_noNEC > 0) %>%
    arrange(FDR_wilcox,
            p_wilcox,
            desc(n_N),
            desc(n_sites),
            desc(mean_diff_prev))
  
  colnames(module_summary$postNEC_retained_fraction)
  
  fwrite(strong_modules,"Strong_linked_haplotype_like_modules_for_downstream.tsv", sep = "\t")
  
  cat("\nStrong modules for downstream analysis:\n")
  print(strong_modules)
  
  
  
  
  ## 20. Select nonsynonymous sites for downstream protein/3D
  selected_N_sites_for_3D <- module_site_anno %>%
    semi_join(strong_modules %>%
                dplyr::select(Species,
                              genome_sp,
                              module_id),
              by = c("Species", "genome_sp", "module_id")) %>%
    filter(mutation_type == "N") %>%
    arrange(Species,
            genome_sp,
            module_id,
            gene,
            position)
  
  fwrite(selected_N_sites_for_3D,
         "Selected_nonsynonymous_sites_in_strong_linked_modules_for_3D.tsv",sep = "\t")
  
  cat("\nNonsynonymous sites selected for downstream protein/3D:\n")
  print(selected_N_sites_for_3D)
  
  cat("\nN sites by gene:\n")
  selected_N_sites_for_3D %>%
    dplyr::count(Species,
                 genome_sp,
                 module_id,
                 gene,
                 name = "n_N_sites") %>%
    dplyr::arrange(desc(n_N_sites)) %>%
    head(n = 100)
  
  
  ## 21. Final output files
  cat("\nMain output files generated:\n")
  cat("1. Selected_8genomes_candidate_sites_for_haplotype_analysis.tsv\n")
  cat("2. Eight_genomes_candidate_sites_pairwise_cooccurrence_sample_level.tsv\n")
  cat("3. Tightly_linked_candidate_SNV_pairs_phi0.7_jaccard0.6.tsv\n")
  cat("4. Tightly_linked_candidate_SNV_pairs_phi0.6_jaccard0.5_relaxed.tsv\n")
  cat("5. Tightly_linked_haplotype_like_modules_sites.tsv\n")
  cat("6. Tightly_linked_haplotype_like_module_site_annotation.tsv\n")
  cat("7. Tightly_linked_haplotype_like_module_summary.tsv\n")
  cat("8. Tightly_linked_haplotype_like_module_score_by_sample.tsv\n")
  cat("9. Tightly_linked_haplotype_like_module_score_stats.tsv\n")
  cat("10. Strong_linked_haplotype_like_modules_for_downstream.tsv\n")
  cat("11. Selected_nonsynonymous_sites_in_strong_linked_modules_for_3D.tsv\n")
  
  
  
  
  ###### Visualization: Method 1 - Supplementary figure
  ## Candidate NEC-associated linked SNV modules
  ## 4 species + 8 priority genomes
  ##
  ## Input objects or files:
  ##   module_summary
  ##   module_score_all
  ##   strong_modules
  ##   module_site_anno
  ##
  ## If objects are not in memory, this script reads:
  ##   Tightly_linked_haplotype_like_module_summary.tsv
  ##   Tightly_linked_haplotype_like_module_score_by_sample.tsv
  ##   Strong_linked_haplotype_like_modules_for_downstream.tsv
  ##   Tightly_linked_haplotype_like_module_site_annotation.tsv
  ##
  ## Main outputs:
  ##   Candidate_8genomes_all_module_summary.tsv
  ##   Candidate_8genomes_strong_module_summary.tsv
  ##   Candidate_8genomes_module_gene_summary.tsv
  ##   Candidate_8genomes_module_site_annotation.tsv
  ##
  ## Figures:
  ##   Candidate_8genomes_Fig1_genome_module_burden.pdf/png
  ##   Candidate_8genomes_Fig2_module_prevalence_heatmap.pdf/png
  ##   Candidate_8genomes_Fig3_module_score_by_NEC_stage.pdf/png
  ##   Candidate_8genomes_Fig4_gene_region_prevalence_bubble.pdf/png
  
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
  library(stringr)
  
  module_summary #
  module_score_all#
  strong_modules#
  module_site_anno#
  
  ## 1. Define 4 candidate species and 8 candidate genomes
  target_info <- data.frame(
    Species = c(rep("Klebsiella pneumoniae", 5),
                "Enterobacter hormaechei_A",
                "Escherichia coli",
                "Klebsiella michiganensis"),
    genome_sp = c("GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1",
                  "VatanenT_2016_SRR4408069_bin.6_k119",
                  "Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141",
                  "ShaoY_2019_ERR3404951_bin.2_k141",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020838.1",
                  "GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1",
                  "CokerMO_2019_SRR8692237_bin.9_k141",
                  "ShaoY_2019_ERR3405384_bin.5_k141"),
    Priority = c(rep("Priority 1", 5),
                 "Priority 2",
                 "Priority 3",
                 "Priority 4"),
    genome_order = 1:8, stringsAsFactors = FALSE)
  priority_species <- unique(target_info$Species)
  target_genomes <- target_info$genome_sp
  
  ## If the plot is too crowded, only plot the top modules for each genome
  ## To plot all modules, set this value to Inf
  top_modules_per_genome <- Inf
  
  ## Define thresholds for preNEC-enriched modules and persistent modules
  delta_cutoff <- 0.05
  retained_ratio_cutoff <- 0.70
  
  
  ## 2. Helper functions
  make_genome_label <- function(x) {
    x <- as.character(x)
    
    if (grepl("^GCF_", x)) {
      acc <- sub("^(GCF_[0-9\\.]+).*", "\\1", x)
      scaffold <- sub("^.*_genomic_", "", x)
      scaffold <- sub("^NZ_", "", scaffold)
      return(paste0(acc, " | ", scaffold))
    }
    
    return(x)
  }
  
  
  make_species_short <- function(x) {
    x <- as.character(x)
    
    dplyr::case_when(x == "Klebsiella pneumoniae" ~ "K. pneumoniae",
                     x == "Enterobacter hormaechei_A" ~ "E. hormaechei_A",
                     x == "Escherichia coli" ~ "E. coli",
                     x == "Klebsiella michiganensis" ~ "K. michiganensis",
                     TRUE ~ x)
  }
  
  
  make_gene_label <- function(gene, genome_sp) {
    gene <- as.character(gene)
    genome_sp <- as.character(genome_sp)
    
    if (is.na(gene) || gene == "" || gene == "intergenic_or_unannotated") {
      return("intergenic/unannotated")
    }
    
    prefix <- paste0(genome_sp, "_")
    
    if (!is.na(genome_sp) && startsWith(gene, prefix)) {
      return(paste0("gene_", substring(gene, nchar(prefix) + 1)))
    }
    
    return(gene)
  }
  
  make_gene_label_vec <- Vectorize(make_gene_label)
  
  
  check_required <- function(df, required_cols, df_name) {
    missing_cols <- setdiff(required_cols, colnames(df))
    
    if (length(missing_cols) > 0) {
      stop(
        paste0(
          df_name,
          " is missing required columns:\n",
          paste(missing_cols, collapse = ", ")
        )
      )
    }
  }
  
  
  ## 3. Standardize input tables
  module_summary2 <- module_summary %>%
    as.data.frame() %>%
    mutate(Species = as.character(Species),
           genome_sp = as.character(genome_sp),
           module_id = as.character(module_id))
  
  module_score_all2 <- module_score_all %>%
    as.data.frame() %>%
    mutate(Species = as.character(Species),
           genome_sp = as.character(genome_sp),
           module_id = as.character(module_id))
  
  strong_modules2 <- strong_modules %>%
    as.data.frame() %>%
    mutate(Species = as.character(Species),
           genome_sp = as.character(genome_sp),
           module_id = as.character(module_id))
  
  module_site_anno2 <- module_site_anno %>%
    as.data.frame() %>%
    mutate(Species = as.character(Species),
           genome_sp = as.character(genome_sp),
           module_id = as.character(module_id))
  
  check_required(module_summary2, c("Species", "genome_sp", "module_id"), "module_summary")
  
  check_required(module_site_anno2,c("Species", "genome_sp", "module_id"),
                 "module_site_anno")
  
  if (!"position" %in% colnames(module_site_anno2)) {
    module_site_anno2$position <- NA
  }
  
  if (!"site_id" %in% colnames(module_site_anno2)) {
    module_site_anno2 <- module_site_anno2 %>%
      mutate(site_id = paste(genome_sp, position, sep = "_pos_"))
  }
  
  if (!"mutation_type" %in% colnames(module_site_anno2)) {
    module_site_anno2$mutation_type <- "unknown"
  }
  
  if (!"gene" %in% colnames(module_site_anno2)) {
    module_site_anno2$gene <- NA_character_
  }
  
  if (!"mutation" %in% colnames(module_site_anno2)) {
    module_site_anno2$mutation <- NA_character_
  }
  
  module_site_anno2 <- module_site_anno2 %>%
    mutate(site_id = as.character(site_id),
           position = as.character(position),
           gene = as.character(gene),
           mutation = as.character(mutation),
           mutation_type = ifelse(is.na(mutation_type) | mutation_type == "",
                                  "unknown",as.character(mutation_type)) )
  
  site_count_by_module <- module_site_anno2 %>%
    group_by(Species, genome_sp, module_id) %>%
    summarise(n_sites_from_site = n_distinct(site_id),
              n_N_from_site = sum(mutation_type == "N", na.rm = TRUE),
              n_S_from_site = sum(mutation_type == "S", na.rm = TRUE),
              n_I_from_site = sum(mutation_type == "I", na.rm = TRUE),
              n_unknown_from_site = sum(!mutation_type %in% c("N", "S", "I"),
                                        na.rm = TRUE),
              n_genes_from_site = n_distinct(gene[!is.na(gene) & gene != ""]),
              .groups = "drop")
  
  module_summary2 <- module_summary2 %>%
    left_join(site_count_by_module,by = c("Species", "genome_sp", "module_id"))
  
  numeric_cols <- c( "n_sites",
                     "n_N",
                     "n_S",
                     "n_I",
                     "n_genes",
                     "mean_noNEC_prev",
                     "mean_preNEC_prev",
                     "mean_postNEC_prev",
                     "mean_diff_prev",
                     "postNEC_retained_fraction")
  
  for (cc in numeric_cols) {
    if (!cc %in% colnames(module_summary2)) {
      module_summary2[[cc]] <- NA_real_
    }
    
    module_summary2[[cc]] <- as.numeric(module_summary2[[cc]])
  }
  
  module_summary2 <- module_summary2 %>%
    mutate(n_sites = coalesce(n_sites, n_sites_from_site),
           n_N = coalesce(n_N, n_N_from_site, 0),
           n_S = coalesce(n_S, n_S_from_site, 0),
           n_I = coalesce(n_I, n_I_from_site, 0),
           n_genes = coalesce(n_genes, n_genes_from_site, 0),
           n_unknown = pmax(n_sites - n_N - n_S - n_I, 0),
           mean_diff_prev = coalesce(mean_diff_prev,
                                     mean_preNEC_prev - mean_noNEC_prev),
           postNEC_retained_fraction = coalesce(postNEC_retained_fraction,
                                                ifelse(mean_preNEC_prev > 0,mean_postNEC_prev / mean_preNEC_prev,NA_real_)))
  
  
  ## 4. Extract modules from 4 species and 8 genomes
  strong_key <- strong_modules2 %>%
    mutate(module_key = paste(Species, genome_sp, module_id, sep = "||")) %>%
    pull(module_key)
  
  target_modules_all <- module_summary2 %>%
    inner_join(target_info,by = c("Species", "genome_sp")) %>%
    mutate(module_key = paste(Species, genome_sp, module_id, sep = "||"),
           is_strong_module = module_key %in% strong_key,
           species_short = make_species_short(Species),
           genome_short = sapply(genome_sp, make_genome_label),
           genome_label = paste0(species_short,"\n",str_trunc(genome_short, width = 55)),
           has_N = ifelse(n_N > 0, "Contains nonsynonymous", "No nonsynonymous"),
           delta_pre_no = mean_preNEC_prev - mean_noNEC_prev,
           retained_ratio = ifelse(mean_preNEC_prev > 0,
                                   mean_postNEC_prev / mean_preNEC_prev,NA_real_ ),
           module_pattern = case_when(delta_pre_no >= delta_cutoff &!is.na(retained_ratio) &
                                        retained_ratio >= retained_ratio_cutoff ~ "Persistent into postNEC",
                                      
                                      delta_pre_no >= delta_cutoff &!is.na(retained_ratio) &
                                        retained_ratio < retained_ratio_cutoff ~ "Transient preNEC module",
                                      
                                      delta_pre_no >= delta_cutoff ~ "preNEC-enriched",
                                      
                                      TRUE ~ "Other / weak pattern")) %>%
    group_by(Species, genome_sp) %>%
    arrange(desc(is_strong_module),desc(delta_pre_no),desc(n_N),desc(n_sites),
            .by_group = TRUE) %>%
    mutate(module_rank_within_genome = row_number()) %>%
    ungroup() %>%
    arrange(genome_order, module_rank_within_genome)
  
  if (nrow(target_modules_all) == 0) {
    cat("\nNo modules were found for the 4 candidate species and 8 genomes.\n")
    cat("\nAvailable Species in module_summary:\n")
    print(sort(unique(module_summary2$Species)))
    
    cat("\nAvailable genome_sp in module_summary for target species:\n")
    print(module_summary2 %>%filter(Species %in% priority_species) %>%
            distinct(Species, genome_sp) %>%arrange(Species, genome_sp))
    
    stop("Please check whether Species and genome_sp names match your module_summary table.")
  }
  
  cat("\nNumber of modules found in 4 species / 8 genomes:\n")
  print(nrow(target_modules_all))
  
  cat("\nModule count per genome:\n")
  print(target_modules_all %>%
          count(Species, genome_sp, name = "n_modules") %>%
          arrange(match(genome_sp, target_genomes)))
  
  missing_target_genomes <- target_info %>%
    anti_join(target_modules_all %>%
                distinct(Species, genome_sp),
              by = c("Species", "genome_sp"))
  
  if (nrow(missing_target_genomes) > 0) {
    cat("\nWarning: these target genomes have no modules in module_summary:\n")
    print(missing_target_genomes)
  }
  
  
  ## 5. Prepare module site and gene summaries
  target_sites_all <- module_site_anno2 %>%
    inner_join(target_info,by = c("Species", "genome_sp")) %>%
    semi_join(target_modules_all %>%select(Species, genome_sp, module_id),
              by = c("Species", "genome_sp", "module_id")) %>%
    mutate(species_short = make_species_short(Species),
           genome_short = sapply(genome_sp, make_genome_label),
           genome_label = paste0(species_short,"\n",
                                 str_trunc(genome_short, width = 55)),
           gene_short = make_gene_label_vec(gene, genome_sp),
           mutation_label = ifelse(is.na(mutation) | mutation == "",
                                   mutation_type,paste0(mutation_type, ":", mutation)),
           site_label = paste0(gene_short," | pos ",position," | ",mutation_label))
  
  gene_prev_cols <- intersect(c("noNEC_prev", "preNEC_prev", "postNEC_prev"),
                              colnames(target_sites_all))
  
  if (length(gene_prev_cols) < 3) {
    target_sites_all$noNEC_prev <- NA_real_
    target_sites_all$preNEC_prev <- NA_real_
    target_sites_all$postNEC_prev <- NA_real_
  }
  
  target_sites_all <- target_sites_all %>%
    mutate(noNEC_prev = as.numeric(noNEC_prev),
           preNEC_prev = as.numeric(preNEC_prev),
           postNEC_prev = as.numeric(postNEC_prev))
  
  target_gene_summary <- target_sites_all %>%
    group_by(Species,
             genome_sp,
             Priority,
             genome_order,
             module_id,
             gene_short) %>%
    summarise(n_sites = n_distinct(site_id),
              n_N = sum(mutation_type == "N", na.rm = TRUE),
              n_S = sum(mutation_type == "S", na.rm = TRUE),
              n_I = sum(mutation_type == "I", na.rm = TRUE),
              n_unknown = sum(!mutation_type %in% c("N", "S", "I"), na.rm = TRUE),
              mean_noNEC_prev = mean(noNEC_prev, na.rm = TRUE),
              mean_preNEC_prev = mean(preNEC_prev, na.rm = TRUE),
              mean_postNEC_prev = mean(postNEC_prev, na.rm = TRUE),
              mutations = paste(sort(unique(na.omit(mutation_label))),collapse = ";"),
              .groups = "drop") %>%
    mutate(mean_noNEC_prev = ifelse(is.nan(mean_noNEC_prev), NA, mean_noNEC_prev),
           mean_preNEC_prev = ifelse(is.nan(mean_preNEC_prev), NA, mean_preNEC_prev),
           mean_postNEC_prev = ifelse(is.nan(mean_postNEC_prev), NA, mean_postNEC_prev),
           delta_pre_no = mean_preNEC_prev - mean_noNEC_prev,
           retained_ratio = ifelse(mean_preNEC_prev > 0, mean_postNEC_prev / mean_preNEC_prev,NA_real_),
           contains_N = ifelse(n_N > 0, "Contains nonsynonymous", "No nonsynonymous")) %>%
    arrange(genome_order,
            module_id,
            desc(delta_pre_no),
            desc(n_N),
            desc(n_sites))
  
  
  ## 6. Save module and gene tables
  target_strong_modules <- target_modules_all %>%
    filter(is_strong_module)
  
  
  ## 7. Select modules for plotting
  if (is.finite(top_modules_per_genome)) {
    target_modules_plot <- target_modules_all %>%
      group_by(Species, genome_sp) %>%
      arrange(desc(is_strong_module),
              desc(delta_pre_no),
              desc(n_N),
              desc(n_sites),
              .by_group = TRUE) %>%
      slice_head(n = top_modules_per_genome) %>%
      ungroup()
  } else {
    target_modules_plot <- target_modules_all
  }
  
  target_modules_plot <- target_modules_plot %>%
    mutate(module_plot_id = paste(
      Species,
      genome_sp,
      module_id,
      sep = "||"),
      #module_label = paste0(
      #  module_id," | ",species_short," | ",str_trunc(genome_short, width = 35),
      #  " | ",n_sites," sites"," | N=",n_N),
      
      module_label = paste0(
        module_id," | ",species_short," | ",str_trunc(genome_short, width = 35)),
      
      module_pattern = factor(
        module_pattern,
        levels = c("Persistent into postNEC",
                   "Transient preNEC module",
                   "preNEC-enriched",
                   "Other / weak pattern")),
      has_N = factor(has_N,levels = c("Contains nonsynonymous", "No nonsynonymous"))) %>%
    arrange(genome_order,
            desc(is_strong_module),
            desc(delta_pre_no),
            desc(n_N),
            desc(n_sites))
  
  module_levels <- rev(target_modules_plot$module_plot_id)
  
  module_label_map <- target_modules_plot$module_label
  names(module_label_map) <- target_modules_plot$module_plot_id
  
  target_modules_plot <- target_modules_plot %>%
    mutate(module_plot_id = factor(module_plot_id, levels = module_levels))
  
  
  ## 8. Figure 2: module prevalence heatmap across NEC stages
  sample_col <- intersect(c("sample2", "sample", "Run", "SampleID"),
                          colnames(module_score_all2))
  
  if (length(sample_col) == 0) {
    sample_col <- NA_character_
  } else {
    sample_col <- sample_col[1]
  }
  
  stage_labels <- c("noNEC" = "noNEC","preNEC" = "preNEC","postNEC" = "postNEC")
  
  if (!is.na(sample_col) && "NEC_stage" %in% colnames(module_score_all2)) {
    stage_n <- module_score_all2 %>%
      mutate(NEC_stage = as.character(NEC_stage)) %>%
      distinct(.data[[sample_col]], NEC_stage) %>%
      count(NEC_stage, name = "n")
    
    for (ss in names(stage_labels)) {
      nn <- stage_n$n[stage_n$NEC_stage == ss]
      
      if (length(nn) == 1) {
        stage_labels[ss] <- paste0(ss, "\n(n=", nn, ")")
      }
    }
  }
  
  module_prev_long <- target_modules_plot %>%
    select(Species,
           genome_sp,
           module_id,
           module_plot_id,
           module_label,
           n_sites,
           n_N,
           n_S,
           n_I,
           n_unknown,
           has_N,
           is_strong_module,
           module_pattern,
           mean_noNEC_prev,
           mean_preNEC_prev,
           mean_postNEC_prev) %>%
    pivot_longer(cols = c(mean_noNEC_prev,mean_preNEC_prev,mean_postNEC_prev),
                 names_to = "Stage",values_to = "prevalence") %>%
    mutate(Stage = case_when(Stage == "mean_noNEC_prev" ~ "noNEC",
                             Stage == "mean_preNEC_prev" ~ "preNEC",
                             Stage == "mean_postNEC_prev" ~ "postNEC"),
           Stage = factor(Stage, levels = c("noNEC", "preNEC", "postNEC")),
           module_plot_id = factor(module_plot_id, levels = module_levels))
  
  max_prev <- max(module_prev_long$prevalence, na.rm = TRUE)
  
  if (!is.finite(max_prev)) {
    max_prev <- 1
  }
  
  p_heat <- ggplot(module_prev_long,
                   aes(x = Stage, y = module_plot_id, fill = prevalence)) +
    geom_tile(color = "white",linewidth = 0.25,width = 0.95,height = 0.95) +
    scale_x_discrete(labels = stage_labels) +
    scale_y_discrete(labels = module_label_map) +
    scale_fill_gradient(low = "white",high = "#B2182B",
                        limits = c(0, max_prev),
                        labels = percent_format(accuracy = 1),
                        name = "Mean prevalence") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_text(color = "black", size = 5.8),
          axis.title = element_text(color = "black"),
          legend.position = "top") +
    labs(x = NULL,
         y = "Linked SNV module")
  
  p_sites <- ggplot(target_modules_plot,
                    aes(x = "SNV sites", y = module_plot_id)) +
    geom_tile(fill = "grey98",color = "grey75",linewidth = 0.25,
              width = 0.95,height = 0.95) +
    geom_text(aes(label = n_sites),size = 2.8, color = "black") +
    scale_y_discrete(labels = module_label_map) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank() )
  
  p_N <- ggplot(target_modules_plot,aes(x = "Nonsynonymous", y = module_plot_id)) +
    geom_tile(fill = "grey98",color = "grey75",linewidth = 0.25,width = 0.95,height = 0.95) +
    geom_point(aes(shape = has_N),size = 2.6,color = "black") +
    scale_shape_manual(values = c("Contains nonsynonymous" = 16,
                                  "No nonsynonymous" = 17),name = "") +
    scale_y_discrete(labels = module_label_map) +
    geom_label(aes(label=n_N,hjust = -1))+
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          legend.position = "top")
  
  p_pattern <- ggplot(target_modules_plot,
                      aes(x = "Pattern", y = module_plot_id, fill = module_pattern)) +
    geom_tile(color = "white",linewidth = 0.25,width = 0.95,height = 0.95) +
    scale_fill_manual(values = c("Persistent into postNEC" = "#B2182B",
                                 "Transient preNEC module" = "#2166AC",
                                 "preNEC-enriched" = "#F4A261",
                                 "Other / weak pattern" = "grey85"),name = "Module pattern" ) +
    scale_y_discrete(labels = module_label_map) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          legend.position = "top")
  
  p_module_heatmap <- p_heat + p_sites + p_N + p_pattern +
    plot_layout(widths = c(3.8, 0.55, 0.55, 1.0))
  
  p_module_heatmap <- p_module_heatmap +
    plot_annotation( title = "Module prevalence across NEC stages in 4 candidate species and 8 genomes",
                     subtitle = "Tiles show mean module prevalence across noNEC, preNEC and postNEC. N indicates whether the module contains nonsynonymous SNV sites.")
  p_module_heatmap
  # save as 10 * 10
  
  


  ##### Supplementary Fig. S11. Co-occurrence network representation of NEC-associated linked SNV modules across disease stages ##########
  
  
  ######### Visualization: Method 2 - Supplementary figure 
  ## Module prevalence network plot for 4 species / 8 genomes
  ##
  ## Input objects or files:
  ##   module_summary
  ##   strong_modules
  ##   module_site_anno
  ##
  ## If objects are not in memory, this script reads:
  ##   Tightly_linked_haplotype_like_module_summary.tsv
  ##   Strong_linked_haplotype_like_modules_for_downstream.tsv
  ##   Tightly_linked_haplotype_like_module_site_annotation.tsv
  ##
  ## Main figure:
  ##   Candidate_8genomes_module_prevalence_network_all_stages.pdf/png
  ##
  ## Optional figures:
  ##   Candidate_8genomes_module_prevalence_network_noNEC.pdf/png
  ##   Candidate_8genomes_module_prevalence_network_preNEC.pdf/png
  ##   Candidate_8genomes_module_prevalence_network_postNEC.pdf/png
  
  
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(igraph)
  library(purrr)
  
  
  ## 1. Read data if objects are not in memory
  if (!exists("module_summary")) {
    module_summary <- fread("Tightly_linked_haplotype_like_module_summary.tsv")
  }
  
  if (!exists("strong_modules")) {
    strong_modules <- fread("Strong_linked_haplotype_like_modules_for_downstream.tsv")
  }
  
  if (!exists("module_site_anno")) {
    module_site_anno <- fread("Tightly_linked_haplotype_like_module_site_annotation.tsv")
  }
  
  cat("module_summary: ")
  print(dim(module_summary))
  
  cat("strong_modules: ")
  print(dim(strong_modules))
  
  cat("module_site_anno: ")
  print(dim(module_site_anno))
  
  
  
  table(module_summary$Species,module_summary$module_id)
  
  ## 2. Define 4 species and 8 genomes
  target_info <- data.frame(
    Species = c(rep("Klebsiella pneumoniae", 2),
                "Enterobacter hormaechei_A",
                #"Escherichia coli",
                "Klebsiella michiganensis"),
    genome_sp = c("GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1",
                  #"VatanenT_2016_SRR4408069_bin.6_k119",
                  "Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141",
                  #"ShaoY_2019_ERR3404951_bin.2_k141",
                  #### "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020838.1",
                  "GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1",
                  #"CokerMO_2019_SRR8692237_bin.9_k141",
                  "ShaoY_2019_ERR3405384_bin.5_k141" ),
    Priority = c( rep("Priority 1", 2),
                  "Priority 2",
                  "Priority 3"#,#"Priority 4"
    ),
    genome_order = 1:4,
    stringsAsFactors = FALSE)
  
  
  priority_species <- unique(target_info$Species)
  target_genomes <- target_info$genome_sp
  
  
  ## 3. User options
  ## Whether to plot only strong modules
  only_strong_modules <- FALSE
  
  ## Maximum number of modules to plot per genome
  ## Set to Inf if you want to plot all modules
  top_modules_per_genome <- Inf
  
  ## Maximum number of edges retained per module to avoid an overly dense network
  max_edges_per_module <- 120
  
  ## Random seed for network layout
  layout_seed <- 123
  
  
  ## 4. Helper functions
  make_genome_label <- function(x) {
    x <- as.character(x)
    
    if (grepl("^GCF_", x)) {
      acc <- sub("^(GCF_[0-9\\.]+).*", "\\1", x)
      scaffold <- sub("^.*_genomic_", "", x)
      scaffold <- sub("^NZ_", "", scaffold)
      return(paste0(acc, " | ", scaffold))
    }
    
    return(x)
  }
  
  
  make_species_short <- function(x) {
    x <- as.character(x)
    
    dplyr::case_when(
      x == "Klebsiella pneumoniae" ~ "K. pneumoniae",
      x == "Enterobacter hormaechei_A" ~ "E. hormaechei_A",
      x == "Escherichia coli" ~ "E. coli",
      x == "Klebsiella michiganensis" ~ "K. michiganensis",
      TRUE ~ x
    )
  }
  
  
  make_gene_label <- function(gene, genome_sp) {
    gene <- as.character(gene)
    genome_sp <- as.character(genome_sp)
    
    if (is.na(gene) || gene == "" || gene == "intergenic_or_unannotated") {
      return("intergenic/unannotated")
    }
    
    prefix <- paste0(genome_sp, "_")
    
    if (!is.na(genome_sp) && startsWith(gene, prefix)) {
      return(paste0("gene_", substring(gene, nchar(prefix) + 1)))
    }
    
    return(gene)
  }
  
  make_gene_label_vec <- Vectorize(make_gene_label)
  
  
  make_module_edges <- function(site_keys, max_edges = 120) {
    site_keys <- unique(site_keys)
    
    if (length(site_keys) < 2) {
      return(data.frame(
        from_key = character(),
        to_key = character(),
        stringsAsFactors = FALSE
      ))
    }
    
    pair_mat <- t(combn(site_keys, 2))
    
    edge_df <- data.frame(
      from_key = pair_mat[, 1],
      to_key = pair_mat[, 2],
      stringsAsFactors = FALSE
    )
    
    if (nrow(edge_df) > max_edges) {
      set.seed(123)
      edge_df <- edge_df[sample(seq_len(nrow(edge_df)), max_edges), ]
    }
    
    return(edge_df)
  }
  
  
  make_hull <- function(df) {
    if (nrow(df) < 3) {
      return(df[0, ])
    }
    
    hull_idx <- chull(df$x, df$y)
    df[hull_idx, , drop = FALSE]
  }
  
  
  ## 5. Standardize module_summary and module_site_anno
  module_summary2 <- module_summary %>%
    as.data.frame() %>%
    mutate(
      Species = as.character(Species),
      genome_sp = as.character(genome_sp),
      module_id = as.character(module_id)
    )
  
  strong_modules2 <- strong_modules %>%
    as.data.frame() %>%
    mutate(
      Species = as.character(Species),
      genome_sp = as.character(genome_sp),
      module_id = as.character(module_id)
    )
  
  module_site_anno2 <- module_site_anno %>%
    as.data.frame() %>%
    mutate(
      Species = as.character(Species),
      genome_sp = as.character(genome_sp),
      module_id = as.character(module_id)
    )
  
  if (!"site_id" %in% colnames(module_site_anno2)) {
    if ("position" %in% colnames(module_site_anno2)) {
      module_site_anno2 <- module_site_anno2 %>%
        mutate(site_id = paste(genome_sp, position, sep = "_pos_"))
    } else {
      stop("module_site_anno needs either site_id or position.")
    }
  }
  
  if (!"mutation_type" %in% colnames(module_site_anno2)) {
    module_site_anno2$mutation_type <- "unknown"
  }
  
  if (!"gene" %in% colnames(module_site_anno2)) {
    module_site_anno2$gene <- NA_character_
  }
  
  if (!"mutation" %in% colnames(module_site_anno2)) {
    module_site_anno2$mutation <- NA_character_
  }
  
  required_prev_cols <- c("noNEC_prev", "preNEC_prev", "postNEC_prev")
  missing_prev_cols <- setdiff(required_prev_cols, colnames(module_site_anno2))
  
  if (length(missing_prev_cols) > 0) {
    stop(
      paste0(
        "module_site_anno is missing prevalence columns:\n",
        paste(missing_prev_cols, collapse = ", "),
        "\nPlease make sure your site annotation table already contains noNEC_prev / preNEC_prev / postNEC_prev."
      )
    )
  }
  
  module_site_anno2 <- module_site_anno2 %>%
    mutate(
      site_id = as.character(site_id),
      mutation_type = ifelse(
        is.na(mutation_type) | mutation_type == "",
        "unknown",
        as.character(mutation_type)
      ),
      gene = as.character(gene),
      mutation = as.character(mutation),
      noNEC_prev = as.numeric(noNEC_prev),
      preNEC_prev = as.numeric(preNEC_prev),
      postNEC_prev = as.numeric(postNEC_prev)
    )
  
  if (!"mean_diff_prev" %in% colnames(module_summary2)) {
    if (all(c("mean_preNEC_prev", "mean_noNEC_prev") %in% colnames(module_summary2))) {
      module_summary2 <- module_summary2 %>%
        mutate(mean_diff_prev = as.numeric(mean_preNEC_prev) - as.numeric(mean_noNEC_prev))
    } else {
      module_summary2$mean_diff_prev <- NA_real_
    }
  }
  
  
  ## 6. Select modules from 8 genomes
  strong_key <- strong_modules2 %>%
    mutate(module_key = paste(Species, genome_sp, module_id, sep = "||")) %>%
    pull(module_key)
  
  target_modules_all <- module_summary2 %>%
    inner_join(target_info, by = c("Species", "genome_sp")) %>%
    mutate(module_key = paste(Species, genome_sp, module_id, sep = "||"),
           is_strong_module = module_key %in% strong_key,
           species_short = make_species_short(Species),
           genome_short = sapply(genome_sp, make_genome_label),
           genome_label = paste0(species_short,
                                 "\n",
                                 str_trunc(genome_short, width = 55))) %>%
    group_by(Species, genome_sp) %>%
    arrange(desc(is_strong_module),
            desc(mean_diff_prev),
            desc(n_N),
            desc(n_sites),
            .by_group = TRUE) %>%
    mutate(module_rank_within_genome = row_number()) %>%
    ungroup() %>%
    arrange(genome_order, module_rank_within_genome)
  
  if (only_strong_modules) {
    target_modules_plot <- target_modules_all %>%
      filter(is_strong_module)
  } else {
    target_modules_plot <- target_modules_all
  }
  
  if (is.finite(top_modules_per_genome)) {
    target_modules_plot <- target_modules_plot %>%
      group_by(Species, genome_sp) %>%
      arrange(
        desc(is_strong_module),
        desc(mean_diff_prev),
        desc(n_N),
        desc(n_sites),
        .by_group = TRUE
      ) %>%
      slice_head(n = top_modules_per_genome) %>%
      ungroup()
  }
  
  cat("\nModules selected for network plotting:\n")
  print(target_modules_plot %>%
          count(Species, genome_sp, name = "n_modules") %>%
          arrange(match(genome_sp, target_genomes)))
  
  
  ## 7. Build node table
  nodes_base <- module_site_anno2 %>%
    inner_join(target_info,
               by = c("Species", "genome_sp")) %>%
    semi_join(target_modules_plot %>%
                select(Species, genome_sp, module_id),
              by = c("Species", "genome_sp", "module_id") ) %>%
    mutate(species_short = make_species_short(Species),
           genome_short = sapply(genome_sp, make_genome_label),
           genome_panel = paste0(
             species_short,
             "\n",
             str_trunc(genome_short, width = 50)),
           gene_short = make_gene_label_vec(gene, genome_sp),
           mutation_label = ifelse(is.na(mutation) | mutation == "",
                                   mutation_type,
                                   paste0(mutation_type, ":", mutation)),
           site_label = paste0(gene_short,
                               " | pos ",
                               ifelse("position" %in% colnames(.), as.character(position), site_id),
                               " | ",),
           node_key = paste(Species, genome_sp, module_id, site_id, sep = "||")) %>%
    distinct(Species, genome_sp, module_id, site_id, .keep_all = TRUE)
  
  if (nrow(nodes_base) == 0) {
    stop("No module sites found for the selected 8 genomes.")
  }
  
  
  ## 8. Build edge table within each module
  edges_base <- nodes_base %>%
    group_by(Species, genome_sp, module_id) %>%
    summarise(node_keys = list(unique(node_key)),
              .groups = "drop") %>%
    mutate(edge_df = map(node_keys, ~ make_module_edges(.x, max_edges = max_edges_per_module))) %>%
    select(-node_keys) %>%
    tidyr::unnest(edge_df)
  
  if (nrow(edges_base) == 0) {
    cat("\nWarning: No edges generated. This usually means most modules contain only 1 site.\n")
  }
  
  
  ## 9. Compute a stable layout for each genome
  layout_one_genome <- function(df_nodes, df_edges, seed = 123) {
    
    vertex_df <- df_nodes %>%
      distinct(node_key, .keep_all = TRUE) %>%
      mutate(name = node_key)
    
    if (nrow(df_edges) > 0) {
      edge_df2 <- df_edges %>%
        filter(
          from_key %in% vertex_df$node_key,
          to_key %in% vertex_df$node_key
        ) %>%
        transmute(
          from = from_key,
          to = to_key
        )
    } else {
      edge_df2 <- data.frame(
        from = character(),
        to = character(),
        stringsAsFactors = FALSE
      )
    }
    
    if (nrow(edge_df2) == 0) {
      g <- make_empty_graph(n = nrow(vertex_df), directed = FALSE)
      V(g)$name <- vertex_df$node_key
    } else {
      g <- graph_from_data_frame(
        d = edge_df2,
        vertices = vertex_df %>% select(name),
        directed = FALSE
      )
    }
    
    set.seed(seed)
    xy <- layout_with_fr(g)
    
    out <- vertex_df %>%
      mutate(
        x = xy[, 1],
        y = xy[, 2]
      )
    
    return(out)
  }
  
  nodes_layout <- nodes_base %>%
    group_by(genome_panel) %>%
    group_split() %>%
    map_dfr(function(df_one) {
      
      current_genome <- unique(df_one$genome_panel)
      
      df_edges_one <- edges_base %>%
        inner_join(
          df_one %>%
            distinct(Species, genome_sp, module_id),
          by = c("Species", "genome_sp", "module_id")
        )
      
      layout_one_genome(
        df_nodes = df_one,
        df_edges = df_edges_one,
        seed = layout_seed
      )
    })
  
  
  ## 10. Prepare edges with coordinates
  edges_xy <- edges_base %>%
    left_join(nodes_layout %>%
                select(from_key = node_key,x_from = x,y_from = y,
                       genome_panel_from = genome_panel),
              by = "from_key" ) %>%
    left_join(nodes_layout %>%
                select(to_key = node_key,x_to = x,y_to = y,
                       genome_panel_to = genome_panel ), by = "to_key" ) %>%
    filter(!is.na(x_from),!is.na(y_from),!is.na(x_to),!is.na(y_to),
           genome_panel_from == genome_panel_to) %>%
    mutate(genome_panel = genome_panel_from)
  
  
  ## 11. Prepare module hull polygons
  hull_df <- nodes_layout %>%
    group_by(genome_panel, module_id) %>%
    group_modify(~ make_hull(.x)) %>%
    ungroup()
  
  module_centers <- nodes_layout %>%
    group_by(genome_panel, module_id) %>%
    summarise(x = mean(x, na.rm = TRUE),
              y = mean(y, na.rm = TRUE),
              n_sites = n(),
              .groups = "drop")
  
  
  ## 12. Build stage-specific node table
  nodes_prev_long <- nodes_layout %>%
    select( Species, genome_sp, module_id, site_id, node_key,
            genome_panel, Priority, genome_order,
            mutation_type, gene_short, site_label,
            noNEC_prev, preNEC_prev, postNEC_prev,
            x, y) %>%
    pivot_longer(cols = c(noNEC_prev, preNEC_prev, postNEC_prev),
                 names_to = "Stage",
                 values_to = "prevalence") %>%
    mutate(Stage = case_when(Stage == "noNEC_prev" ~ "noNEC",
                             Stage == "preNEC_prev" ~ "preNEC",
                             Stage == "postNEC_prev" ~ "postNEC"),
           Stage = factor(Stage, levels = c("noNEC", "preNEC", "postNEC")),
           mutation_type = factor(mutation_type,
                                  levels = c("N", "S", "I", "unknown")))
  
  max_prev <- max(nodes_prev_long$prevalence, na.rm = TRUE)
  if (!is.finite(max_prev)) {
    max_prev <- 1
  }
  
  
  ## 13. Plot combined network: 8 genomes x 3 stages
  max_prev <- max(nodes_prev_long$prevalence, na.rm = TRUE)
  
  p_network_strong <- ggplot() +
    geom_polygon(data = hull_df,
                 aes(x = x,y = y, group = interaction(genome_panel, module_id)),
                 fill = "grey85",alpha = 0.20, color = "grey55",
                 linewidth = 0.30,linetype = "dashed",show.legend = FALSE ) +
    geom_segment(data = edges_xy,
                 aes(x = x_from,y = y_from,xend = x_to,yend = y_to),
                 color = "grey70",linewidth = 0.22,alpha = 0.65) +
    geom_point(data = nodes_prev_long,
               aes(x = x, y = y,fill = prevalence,shape = mutation_type),
               size = 2.5,color = "grey20",stroke = 0.30,alpha = 0.95) +
    geom_text(data = module_centers,
              aes(x = x,y = y,label = module_id),
              size = 2.8,fontface = "bold",color = "black",vjust = -1.0) +
    facet_grid(genome_panel ~ Stage,scales = "free",switch = "y") +
    
    scale_fill_gradient(low = "white",high = "#B2182B",
                        limits = c(0, max_prev),labels = percent_format(accuracy = 1),name = "Site prevalence") +
    scale_shape_manual(values = c("N" = 21,"S" = 24,"I" = 22,"unknown" = 23),
                       name = "Mutation type") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text.x = element_text(size = 9, face = "bold"),
          strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
          legend.position = "right") +
    labs(title = "Module prevalence network across 4 candidate species and 8 genomes",
         subtitle = "Nodes are SNV sites; edges connect sites within the same module; node color shows site prevalence across NEC stages.")
  p_network_strong
  # save as 10 * 16
  
  
  
  ##### Fig. 7e+f. The 3D structure of (e) genes GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432 and (f) GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436 in K. pneumoniae Module_9 ####
  
  ##### e+f1. 3D structure plot ~ 4 target modules ############
  # Extract directly from .faa:
  # Klebsiella pneumoniae
  # Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141
  # Module_4
  # Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141_13392_1
  # N:G381R;N:Q702E;N:R705*
  # Indicates a premature stop codon, which may cause protein truncation. This may have greater biological interpretability than a regular missense mutation.
  
  # Klebsiella pneumoniae
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1
  # Module_9
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2431;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2433;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2434;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2435;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436;intergenic_or_unannotated
  # ;N:I171L;N:K1178N;N:V472A;S:128;S:131;S:1352;S:149;S:161;S:188;S:2060;S:35;S:479;S:512;S:59;S:95
  
  # Klebsiella michiganensis
  # ShaoY_2019_ERR3405384_bin.5_k141
  # Module_1
  # ShaoY_2019_ERR3405384_bin.5_k141_6236_14;
  # ShaoY_2019_ERR3405384_bin.5_k141_6236_15
  # N:P88Q;N:P96S;N:R94H;N:S103N;N:T105S;S:143;S:152;S:182;S:86;S:98
  # Gene structure plot
  
  # Enterobacter hormaechei_A
  # GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1
  # Module_2
  # GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1_91
  # N:A280V;N:K105N;N:T670I;S:42;S:654
  # 3D structure, because all mutations are located in one gene
  
  # .gff is used to draw gene neighborhood / gene structure plots, such as the upper arrow plot in Supplementary figure 2.
  # .faa + structure prediction / homologous structures are used to draw protein 3D structures, such as the protein 3D structure on the right side of Supplementary figure 1.
  
  
  ### Enter Linux; Linux code
  # Enter the gene folder
  cd /mnt/data/ShuangPeng/Database/ELGG_representatives_2172/gene
  
  mkdir -p NEC_3D_faa_download
  cd NEC_3D_faa_download
  
  # Write exact protein IDs that have already been confirmed
  # These protein IDs already include the final ORF number:
  cat > exact_protein_ids.txt << 'EOF_INNER'
  Baumann-DudenhoefferAM_2018_SRR7217883_bin.2_k141_13392_1
  ShaoY_2019_ERR3405384_bin.5_k141_6236_14
  ShaoY_2019_ERR3405384_bin.5_k141_6236_15
  GCF_013743755.1_ASM1374375v1_genomic_NZ_CP056748.1_91
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2431
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2433
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2434
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2435
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436
  EOF_INNER
  
  # Write protein IDs that need to be extracted by prefix
  # These are currently scaffold/gene prefixes and may not be complete protein IDs, so extract them by prefix:
  cat > prefix_protein_ids.txt << 'EOF_INNER'
  EOF_INNER
  
  
  ### Extract these proteins from the combined .faa file
  #### If the format becomes incorrect after pasting, paste the code directly from GPT
  cat > extract_candidate_faa.py << 'EOF_INNER'
  import os
  
  faa_file = "../reps.concat.genes.fixed.faa"
  
  exact_file = "exact_protein_ids.txt"
  prefix_file = "prefix_protein_ids.txt"
  
  out_raw = "Candidate_4modules_WT.raw.faa"
  out_clean = "Candidate_4modules_WT.clean_for_AlphaFold.faa"
  out_table = "Candidate_4modules_extracted_protein_ids.tsv"
  out_missing = "Candidate_4modules_missing_ids.txt"
  
  with open(exact_file) as f:
    exact_ids = set(x.strip() for x in f if x.strip())
  
  with open(prefix_file) as f:
    prefix_ids = [x.strip() for x in f if x.strip()]
  
  found_exact = set()
  found_prefix = {p: [] for p in prefix_ids}
  
  records = []
  
  current_header = None
  current_seq = []
  
  
  def flush_record(header, seq_lines):
    if header is None:
    return
  
  seq = "".join(seq_lines)
  protein_id = header.split()[0]
  
  matched_by = []
  
  if protein_id in exact_ids:
    matched_by.append("exact")
  found_exact.add(protein_id)
  
  for p in prefix_ids:
    if protein_id.startswith(p + "_"):
    matched_by.append(p)
  found_prefix[p].append(protein_id)
  
  if len(matched_by) > 0:
    records.append((protein_id, header, seq, ";".join(matched_by)))
  
  
  with open(faa_file) as f:
    for line in f:
    line = line.rstrip("\n")
  
  if line.startswith(">"):
    flush_record(current_header, current_seq)
  current_header = line[1:]
  current_seq = []
  else:
    current_seq.append(line.strip())
  
  flush_record(current_header, current_seq)
  
  
  with open(out_raw, "w") as o_raw, \
  open(out_clean, "w") as o_clean, \
  open(out_table, "w") as o_tab:
    
    o_tab.write("protein_id\tmatched_by\tlength_raw\tlength_clean\theader\n")
  
  for protein_id, header, seq, matched_by in records:
    clean_seq = seq.replace("*", "")
  
  o_raw.write(">" + header + "\n")
  for i in range(0, len(seq), 60):
    o_raw.write(seq[i:i+60] + "\n")
  
  o_clean.write(">" + protein_id + "\n")
  for i in range(0, len(clean_seq), 60):
    o_clean.write(clean_seq[i:i+60] + "\n")
  
  o_tab.write(
    protein_id + "\t" +
      matched_by + "\t" +
      str(len(seq)) + "\t" +
      str(len(clean_seq)) + "\t" +
      header + "\n"
  )
  
  
  with open(out_missing, "w") as o:
    for x in sorted(exact_ids - found_exact):
    o.write("missing_exact\t" + x + "\n")
  
  for p in prefix_ids:
    if len(found_prefix[p]) == 0:
    o.write("missing_prefix\t" + p + "\n")
  
  
  print("Finished.")
  print("Extracted proteins:", len(records))
  print("Raw fasta:", out_raw)
  print("Clean fasta for AlphaFold/ColabFold:", out_clean)
  print("Extracted table:", out_table)
  print("Missing ID report:", out_missing)
  EOF_INNER
  
  ## Rerun
  python extract_candidate_faa.py
  
  ## Check whether extraction was successful
  cat Candidate_4modules_extracted_protein_ids.tsv
  
  # Check whether there are missing IDs:
  cat Candidate_4modules_missing_ids.txt
  
  # Count how many proteins were extracted:
  grep -c "^>" Candidate_4modules_WT.raw.faa
  grep -c "^>" Candidate_4modules_WT.clean_for_AlphaFold.faa
  
  
  # If successful, split into individual .faa files
  mkdir -p individual_WT_faa
  
  seqkit split \
  -i \
  Candidate_4modules_WT.clean_for_AlphaFold.faa \
  -O individual_WT_faa
  
  # Package for download
  tar -czvf Candidate_4modules_WT_faa_for_3D.tar.gz \
  Candidate_4modules_WT.raw.faa \
  Candidate_4modules_WT.clean_for_AlphaFold.faa \
  Candidate_4modules_extracted_protein_ids.tsv \
  Candidate_4modules_missing_ids.txt \
  exact_protein_ids.txt \
  prefix_protein_ids.txt \
  individual_WT_faa
  
  # Finally download this file: Candidate_4modules_WT_faa_for_3D.tar.gz
  
  
  #### Remaining steps
  # Path: sequences of the 10 genes extracted under Linux
  less /mnt/data/ShuangPeng/Database/ELGG_representatives_2172/gene/3D_input/structure_candidate_proteins.clean.renamed.faa
  # Save as species 3_ 4 genomes_10 genes.tsv
  
  
  
  
  ##### e+f2. Comparison of 3D structures before and after mutation ############
  
  #### First identify the selected sites in Klebsiella pneumoniae
  ### In R
  selected_sites #
  colnames(selected_sites)
  
  # "Klebsiella pneumoniae"
  # "Klebsiella_pneumoniae"
  
  
  # bacteria
  # Klebsiella pneumoniae.tsv
  
  # gene
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2431;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2433;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2434;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2435;
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436;
  # intergenic_or_unannotated
  
  # mutation
  # ;
  # N:I171L;
  # N:K1178N;
  # N:V472A;
  # S:128;
  # S:131;
  # S:1352;
  # S:149;
  # S:161;
  # S:188;
  # S:2060;
  # S:35;
  # S:479;
  # S:512;
  # S:59;
  # S:95
  
  # subset
  test2_1 <- selected_sites %>%
    subset(#gene == "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436" &
      gene %in% c("GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2431",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2433",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2434",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2435",
                  "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436")&
        mutation %in% c("N:I171L","N:K1178N","N:V472A","S:128","S:131","S:1352",
                        "S:149","S:161","S:188","S:2060","S:35","S:479","S:512",
                        "S:59","S:95"))%>%
    distinct(gene, mutation, .keep_all = TRUE)
  
  
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2431 # all S
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432 # "N:I171L","N:V472A"
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2433 # all S
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2434 # all S
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2435 # all S
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436 # "N:K1178N"
  
  ## Focus on two genes: GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432,
  ## GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436
  
  
  #### In Linux
  
  ## First extract the 6 WT proteins in Module_9
  cd /mnt/data/ShuangPeng/Database/ELGG_representatives_2172/gene
  
  # Build nonsynonymous mutation map
  cat > module9_nonsyn_mutation_map.tsv << 'EOF_INNER'
  gene_id	gene_nt_pos	aa_pos	ref_aa	alt_aa	protein_mutation	genome_position	original_label
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432	171	57	I	L	I57L	2538128	N:I171L
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432	472	158	V	A	V158A	2538429	N:V472A
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436	1178	393	K	N	K393N	2541448	N:K1178N
  EOF_INNER
  
  # The conversion formula is: aa_pos = floor((gene_nt_pos - 1) / 3) + 1
  
  
  ## Extract WT .faa files for gene2432 and gene2436
  cat > module9_nonsyn_gene_ids.txt << 'EOF_INNER'
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436
  EOF_INNER
  
  seqkit grep -f module9_nonsyn_gene_ids.txt reps.concat.genes.fixed.faa > module9_nonsyn_WT.faa
  
  grep -c "^>" module9_nonsyn_WT.faa
  seqkit fx2tab -n -l module9_nonsyn_WT.faa
  
  # Expected result: GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432    176
  # GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436    805
  
  ## Generate mutant.faa
  cat > make_module9_mutant_faa.py << 'EOF_INNER'
  #!/usr/bin/env python3
  
  import os
  import pandas as pd
  
  wt_faa = "module9_nonsyn_WT.faa"
  mutation_map = "module9_nonsyn_mutation_map.tsv"
  
  out_dir = "module9_mutant_faa"
  os.makedirs(out_dir, exist_ok=True)
  
  
  def read_fasta(path):
    seqs = {}
  order = []
  header = None
  chunks = []
  
  with open(path) as f:
    for line in f:
    line = line.strip()
  if not line:
    continue
  
  if line.startswith(">"):
    if header is not None:
    seq = "".join(chunks).replace("*", "")
  seqs[header] = seq
  order.append(header)
  
  header = line[1:].split()[0]
  chunks = []
  else:
    chunks.append(line)
  
  if header is not None:
    seq = "".join(chunks).replace("*", "")
  seqs[header] = seq
  order.append(header)
  
  return seqs, order
  
  
  def write_fasta(name, seq, path):
    with open(path, "w") as out:
    out.write(f">{name}\n")
  for i in range(0, len(seq), 80):
    out.write(seq[i:i+80] + "\n")
  
  
  seqs, order = read_fasta(wt_faa)
  mut_df = pd.read_csv(mutation_map, sep="\t")
  
  report = []
  
  for gene_id in order:
    wt_seq = seqs[gene_id]
  
  short_id = gene_id.replace(
    "GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_",
    "Kpn_CP020837_gene"
  )
  
  # Output WT, with stop codon * removed
  write_fasta(
    f"{short_id}|WT",
    wt_seq,
    os.path.join(out_dir, f"{short_id}_WT.faa")
  )
  
  sub = mut_df[mut_df["gene_id"] == gene_id].copy()
  
  if sub.shape[0] == 0:
    continue
  
  mut_seq = list(wt_seq)
  applied = []
  
  for _, row in sub.iterrows():
    aa_pos = int(row["aa_pos"])
  ref_aa = str(row["ref_aa"])
  alt_aa = str(row["alt_aa"])
  protein_mutation = str(row["protein_mutation"])
  original_label = str(row["original_label"])
  
  if aa_pos < 1 or aa_pos > len(wt_seq):
    observed = "OUT_OF_RANGE"
  status = "FAILED_OUT_OF_RANGE"
  else:
    observed = wt_seq[aa_pos - 1]
  if observed != ref_aa:
    status = "FAILED_REF_AA_MISMATCH"
  else:
    mut_seq[aa_pos - 1] = alt_aa
  applied.append(protein_mutation)
  status = "APPLIED"
  
  report.append({
    "gene_id": gene_id,
    "short_id": short_id,
    "original_label": original_label,
    "protein_mutation": protein_mutation,
    "aa_pos": aa_pos,
    "expected_ref_aa": ref_aa,
    "observed_aa": observed,
    "alt_aa": alt_aa,
    "status": status
  })
  
  mut_label = "_".join(applied) if applied else "no_valid_mutation"
  
  write_fasta(
    f"{short_id}|MUT|{mut_label}",
    "".join(mut_seq),
    os.path.join(out_dir, f"{short_id}_MUT_{mut_label}.faa")
  )
  
  report_df = pd.DataFrame(report)
  report_df.to_csv(
    os.path.join(out_dir, "module9_mutant_faa_validation_report.tsv"),
    sep="\t",
    index=False
  )
  
  print(report_df)
  print("[Done] Output directory:", out_dir)
  EOF_INNER
  
  python make_module9_mutant_faa.py
  
  ## Check whether it succeeded
  ls -lh module9_mutant_faa
  
  cat module9_mutant_faa/module9_mutant_faa_validation_report.tsv
  
  ## Check FASTA length
  seqkit fx2tab -n -l module9_mutant_faa/*.faa
  
  ## Next upload files to AlphaFold
  # module9_mutant_faa/Kpn_CP020837_gene2432_MUT_I57L_V158A.faa
  # module9_mutant_faa/Kpn_CP020837_gene2436_MUT_K393N.faa
  
  
  
  ##### e+f3. Codon-level validation ############
  
  ### First build a table for the 3 nonsynonymous sites
  cd /mnt/data/ShuangPeng/Database/ELGG_representatives_2172/gene
  
  cat > module9_nonsyn_sites.tsv << 'EOF_INNER'
  site_id	scaffold	position	gene_id	major_base	mutation
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__2538128	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2538128	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432	A	N:I171L
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__2538429	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2538429	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432	T	N:V472A
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1__2541448	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2541448	GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436	A	N:K1178N
  EOF_INNER
  
  ### Extract .fna files for these two genes
  # Use nucleotide CDS rather than .faa:
  cat > module9_nonsyn_gene_ids.txt << 'EOF_INNER'
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2432
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1_2436
  EOF_INNER
  
  seqkit grep -f module9_nonsyn_gene_ids.txt reps.concat.genes.fixed.fna > module9_nonsyn_WT.fna
  
  grep -c "^>" module9_nonsyn_WT.fna
  seqkit fx2tab -n -l module9_nonsyn_WT.fna
  
  # Expected length: gene2432    531 nt; gene2436    2418 nt
  
  ### Automatically infer codon changes using CDS
  cat > infer_module9_codon_changes.py << 'EOF_INNER'
  #!/usr/bin/env python3
  
  import re
  import pandas as pd
  
  sites_file = "module9_nonsyn_sites.tsv"
  fna_file = "module9_nonsyn_WT.fna"
  out_file = "module9_codon_changes.tsv"
  
  genetic_code = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
  }
  
  
  def read_fasta(path):
    seqs = {}
  header = None
  chunks = []
  with open(path) as f:
    for line in f:
    line = line.strip()
  if not line:
    continue
  if line.startswith(">"):
    if header is not None:
    seqs[header] = "".join(chunks).upper().replace("*", "")
  header = line[1:].split()[0]
  chunks = []
  else:
    chunks.append(line)
  if header is not None:
    seqs[header] = "".join(chunks).upper().replace("*", "")
  return seqs
  
  
  def translate(codon):
    return genetic_code.get(codon.upper(), "X")
  
  
  seqs = read_fasta(fna_file)
  sites = pd.read_csv(sites_file, sep="\t")
  
  rows = []
  
  for _, r in sites.iterrows():
    gene_id = r["gene_id"]
  mutation = r["mutation"]
  major_base = str(r["major_base"]).upper()
  
  m = re.match(r"N:([A-Z\*])(\d+)([A-Z\*])", mutation)
  if not m:
    raise ValueError(f"Cannot parse mutation: {mutation}")
  
  ref_aa_label = m.group(1)
  nt_label = int(m.group(2))
  alt_aa_label = m.group(3)
  
  seq = seqs[gene_id]
  
  candidates = []
  
  # Test nt_label-1, nt_label and nt_label+1 simultaneously
  # because 0-based / 1-based differences are common between inStrain/SNV positions and GFF/CDS coordinates
  for cds_nt_pos in [nt_label - 1, nt_label, nt_label + 1]:
    if cds_nt_pos < 1 or cds_nt_pos > len(seq):
    continue
  
  codon_start = ((cds_nt_pos - 1) // 3) * 3
  codon_end = codon_start + 3
  ref_codon = seq[codon_start:codon_end]
  codon_pos = (cds_nt_pos - 1) % 3 + 1
  aa_pos = codon_start // 3 + 1
  
  ref_aa = translate(ref_codon)
  
  for alt_base in ["A", "C", "G", "T"]:
    if alt_base == ref_codon[codon_pos - 1]:
    continue
  
  alt_codon_list = list(ref_codon)
  alt_codon_list[codon_pos - 1] = alt_base
  alt_codon = "".join(alt_codon_list)
  alt_aa = translate(alt_codon)
  
  match = (ref_aa == ref_aa_label and alt_aa == alt_aa_label)
  
  candidates.append({
    "site_id": r["site_id"],
    "gene_id": gene_id,
    "genomic_position_from_site_table": int(r["position"]),
    "mutation_label": mutation,
    "label_nt_number": nt_label,
    "tested_cds_nt_pos": cds_nt_pos,
    "nt_shift_vs_label": cds_nt_pos - nt_label,
    "aa_pos": aa_pos,
    "codon_pos": codon_pos,
    "ref_codon": ref_codon,
    "alt_codon": alt_codon,
    "ref_base_from_CDS": ref_codon[codon_pos - 1],
    "alt_base_inferred": alt_base,
    "major_base_from_site_table": major_base,
    "ref_aa_from_CDS": ref_aa,
    "alt_aa_if_mutated": alt_aa,
    "ref_aa_label": ref_aa_label,
    "alt_aa_label": alt_aa_label,
    "match_label": match
  })
  
  matched = [x for x in candidates if x["match_label"]]
  
  if matched:
    rows.extend(matched)
  else:
    # If there is no complete match, output all candidates for troubleshooting
    rows.extend(candidates)
  
  out = pd.DataFrame(rows)
  out.to_csv(out_file, sep="\t", index=False)
  
  print("[Done]", out_file)
  print(out[out["match_label"] == True].to_string(index=False))
  EOF_INNER
  
  python infer_module9_codon_changes.py
  
  ## Check results
  column -t -s $'\t' module9_codon_changes.tsv | less -S
  
  
  ##### For the SNVs.tsv file
  ### First extract these 3 sites from all_samples.SNVs.tsv
  find /mnt/data/ShuangPeng -name "/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/all_samples.SNVs.tsv" 2>/dev/null
  
  ###
  python - << 'PY'
  import pandas as pd
  
  snv_file = "/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/all_samples.SNVs.tsv"
  
  df = pd.read_csv(snv_file, sep="\t", nrows=3)
  for i, c in enumerate(df.columns, 1):
    print(i, c)
  PY
  
  ### Extract 3 target genomic positions
  # First build the target site file:
  cat > module9_three_nonsyn_positions.tsv << 'EOF_INNER'
  scaffold	position	mutation_label
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2538128	N:I171L
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2538429	N:V472A
  GCF_002156725.1_ASM215672v1_genomic_NZ_CP020837.1	2541448	N:K1178N
  EOF_INNER
  
  # Then extract
  cat > extract_module9_three_sites_from_SNVS.py << 'EOF_INNER'
  #!/usr/bin/env python3
  
  import pandas as pd
  
  snv_file = "/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/all_samples.SNVs.tsv"
  target_file = "module9_three_nonsyn_positions.tsv"
  out_file = "module9_three_nonsyn_raw_SNVS.tsv"
  
  targets = pd.read_csv(target_file, sep="\t")
  target_keys = set(zip(targets["scaffold"].astype(str), targets["position"].astype(int)))
  
  chunks = []
  
  for chunk in pd.read_csv(snv_file, sep="\t", chunksize=500000):
    # Automatically adapt to the real column names in SNVs.tsv
    scaffold_col = None
  position_col = None
  
  for c in ["scaffold", "scaffold_id", "genome", "scaffold_name"]:
    if c in chunk.columns:
    scaffold_col = c
  break
  
  for c in ["position", "pos", "location"]:
    if c in chunk.columns:
    position_col = c
  break
  
  if scaffold_col is None or position_col is None:
    print("Cannot find scaffold/position columns.")
  print(chunk.columns.tolist())
  raise SystemExit
  
  chunk[scaffold_col] = chunk[scaffold_col].astype(str)
  chunk[position_col] = chunk[position_col].astype(int)
  
  mask = [
    (s, p) in target_keys
    for s, p in zip(chunk[scaffold_col], chunk[position_col])
  ]
  
  sub = chunk.loc[mask].copy()
  if sub.shape[0] > 0:
    chunks.append(sub)
  
  if chunks:
    out = pd.concat(chunks, ignore_index=True)
  out.to_csv(out_file, sep="\t", index=False)
  print("[Done]", out_file, out.shape)
  else:
    print("No records found.")
  EOF_INNER
  
  ## Then run
  python extract_module9_three_sites_from_SNVS.py
  
  # Check which base-count columns are present in SNVs.tsv
  head -1 module9_three_nonsyn_raw_SNVS.tsv | tr '\t' '\n' | nl | less
  
  ### Summarize alleles for each site
  cat > summarize_module9_three_nonsyn_alleles.py << 'EOF_INNER'
  #!/usr/bin/env python3
  
  import pandas as pd
  
  infile = "module9_three_nonsyn_raw_SNVS.tsv"
  
  df = pd.read_csv(infile, sep="\t")
  
  for b in ["A", "C", "T", "G"]:
    df[b] = pd.to_numeric(df[b], errors="coerce").fillna(0)
  
  genetic_code = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
  }
  
  site_info = {
    2538128: {
      "mutation_label": "N:I171L",
      "protein_mutation": "I58L",
      "ref_codon": "ATA",
      "codon_pos": 1,
      "target_aa": "L"
    },
    2538429: {
      "mutation_label": "N:V472A",
      "protein_mutation": "V158A",
      "ref_codon": "GTG",
      "codon_pos": 2,
      "target_aa": "A"
    },
    2541448: {
      "mutation_label": "N:K1178N",
      "protein_mutation": "K393N",
      "ref_codon": "AAA",
      "codon_pos": 3,
      "target_aa": "N"
    }
  }
  
  rows = []
  chosen = []
  
  for pos, info in site_info.items():
    sub = df[df["position"].astype(int) == pos].copy()
  
  if sub.shape[0] == 0:
    print(f"[WARN] position {pos} not found")
  continue
  
  ref_codon = info["ref_codon"]
  codon_pos = info["codon_pos"]
  ref_aa = genetic_code[ref_codon]
  target_aa = info["target_aa"]
  
  base_totals = sub[["A", "C", "T", "G"]].sum()
  total_depth = base_totals.sum()
  
  for base in ["A", "C", "T", "G"]:
    alt_codon_list = list(ref_codon)
  alt_codon_list[codon_pos - 1] = base
  alt_codon = "".join(alt_codon_list)
  aa = genetic_code.get(alt_codon, "X")
  
  if base == ref_codon[codon_pos - 1]:
    role = "reference_base"
  elif aa == target_aa:
    role = "target_nonsyn_base"
  elif aa == ref_aa:
    role = "synonymous_alt_base"
  else:
    role = "other_alt_base"
  
  rows.append({
    "position": pos,
    "mutation_label": info["mutation_label"],
    "protein_mutation": info["protein_mutation"],
    "ref_codon": ref_codon,
    "ref_aa": ref_aa,
    "codon_pos": codon_pos,
    "base": base,
    "alt_codon_if_this_base": alt_codon,
    "aa_if_this_base": aa,
    "role": role,
    "total_read_count_across_extracted_records": int(base_totals[base]),
    "fraction_across_extracted_records": float(base_totals[base] / total_depth) if total_depth > 0 else 0
  })
  
  target_rows = [
    r for r in rows
    if r["position"] == pos and r["role"] == "target_nonsyn_base"
  ]
  
  if len(target_rows) > 0:
    best = max(target_rows, key=lambda x: x["total_read_count_across_extracted_records"])
  chosen.append(best)
  
  out = pd.DataFrame(rows)
  out.to_csv("module9_three_nonsyn_allele_codon_summary.tsv", sep="\t", index=False)
  
  chosen_df = pd.DataFrame(chosen)
  chosen_df.to_csv("module9_three_nonsyn_chosen_codon_change.tsv", sep="\t", index=False)
  
  print("\n[Allele/codon summary]")
  print(out.to_string(index=False))
  
  print("\n[Chosen target nonsyn codon by highest total read count]")
  print(chosen_df.to_string(index=False))
  
  print("\n[Mutation column counts by position]")
  print(df.groupby(["position", "mutation"]).size().reset_index(name="n_records").to_string(index=False))
  
  print("\n[con_base counts by position]")
  print(df.groupby(["position", "con_base"]).size().reset_index(name="n_records").to_string(index=False))
  
  print("\n[var_base counts by position]")
  print(df.groupby(["position", "var_base"]).size().reset_index(name="n_records").to_string(index=False))
  EOF_INNER
  
  python summarize_module9_three_nonsyn_alleles.py
  
  ## Check these two output files: they will show the most likely codon change for each site based on read count.
  column -t -s $'\t' module9_three_nonsyn_chosen_codon_change.tsv
  
  # Check the complete allele table
  column -t -s $'\t' module9_three_nonsyn_allele_codon_summary.tsv | less -S


