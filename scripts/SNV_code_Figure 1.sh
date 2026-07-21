
############### Gut microbial SNV ###################################

.libPaths()
.libPaths(c("D:/2_Software/R/install/R-4.5.1/library_shuang",
            "D:/2_Software/R/install/R-4.5.1/library"))
.libPaths()

# Set working directory
setwd("D:/3_Projects/2_Children_jaundice/3_R_analysis")


### Necessary packages 
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
library(ggpubr)
library(tidyr)



############ Result 1: Children gut bacterial SNV development over time #####


###### On the Linux server:
######## Merge inStrain profile output files into single combined files #######################


##### Generate the combined SNVs.tsv file
# Set output file path
OUTFILE="all_samples.SNVs.tsv"
> "$OUTFILE"

# Write unified header
echo -e "scaffold\tposition\tposition_coverage\tallele_count\tref_base\tcon_base\tvar_base\tref_freq\tcon_freq\tvar_freq\tA\tC\tT\tG\tgene\tmutation\tmutation_type\tcryptic\tclass\tsample\tsample2" >> "$OUTFILE"

# Iterate through all *_SNVs.tsv files
find /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/ \
-type f \( -name "*_SNVs.tsv" -o -name "*_SNVs.tsv.gz" \) | while read file; do

# Extract sample name without path and suffix
sample=$(basename "$file" | sed -E 's/_SNVs\.tsv(\.gz)?$//')

# Skip empty files
[ -s "$file" ] || continue

# Process each file: skip the header and retain only rows with the expected number of fields
# tail -n +2 "$file" | awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" 'NF==19 {OFS="\t"; print $0, s1, s2}' >> "$OUTFILE"

# Use cat or gunzip depending on the file type
if [[ "$file" == *.gz ]]; then
gunzip -c "$file"
else
  cat "$file"
fi | tail -n +2 | \
awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" '
    NF==19 {
      OFS="\t";
      print $0, s1, s2
    }
  ' >> "$OUTFILE"
done

# tail -n +2: skip the original header
# NF==19: retain only correctly formatted rows with 19 fields
# s1="$sample": add the full sample name
# s2="${sample%%_*}": extract the simplified sample name as sample2
# echo -e: add a unified header to ensure DuckDB can correctly parse the table

## Check the merged file:
# Check the number of fields; all rows should have 21 columns
awk -F'\t' '{print NF}' all_samples.SNVs.tsv | sort | uniq -c

# Inspect the first few rows
head -n 3 all_samples.SNVs.tsv | column -t
cut -f20 all_samples.SNVs.tsv | tail -n +2 | sort | uniq | wc -l
# 5897 sample names



##### Generate the combined genome_info.tsv file
# Set output file path
OUTFILE="all_samples.genome_info.tsv"
> "$OUTFILE"

# Write unified header, including the original genome_info columns plus sample and sample2
echo -e "genome\tcoverage\tbreadth\tnucl_diversity\tlength\ttrue_scaffolds\tdetected_scaffolds\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_rarefied\tconANI_reference\tpopANI_reference\tiRep\tiRep_GC_corrected\tlinked_SNV_count\tSNV_distance_mean\tr2_mean\td_prime_mean\tconsensus_divergent_sites\tpopulation_divergent_sites\tSNS_count\tSNV_count\tfiltered_read_pair_count\treads_unfiltered_pairs\treads_mean_PID\treads_unfiltered_reads\tdivergent_site_count\tsample\tsample2" >> "$OUTFILE"

# Iterate through all *_genome_info.tsv files
find /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/ \
-type f \( -name "*_genome_info.tsv" -o -name "*_genome_info.tsv.gz" \) | while read file; do

# Extract sample name
sample=$(basename "$file" | sed -E 's/_genome_info\.tsv(\.gz)?$//')

# Skip empty files
[ -s "$file" ] || continue

# Use cat or gunzip depending on the file type
if [[ "$file" == *.gz ]]; then
gunzip -c "$file"
else
  cat "$file"
fi | tail -n +2 | \
awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" '
    NF==30 {
      OFS="\t";
      print $0, s1, s2
    }
  ' >> "$OUTFILE"

done

less all_samples.genome_info.tsv
awk -F'\t' '{print NF}' all_samples.genome_info.tsv | sort | uniq -c
head -n 3 all_samples.genome_info.tsv | column -t
cut -f31 all_samples.genome_info.tsv | tail -n +2 | sort | uniq | wc -l
# 5897 sample names



##### Generate the combined scaffold_info.tsv file
# Set output file path
OUTFILE="all_samples.scaffold_info.tsv"
> "$OUTFILE"

# Write unified header, including the original scaffold_info columns plus sample and sample2
echo -e "scaffold\tlength\tcoverage\tbreadth\tnucl_diversity\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_median\tnucl_diversity_rarefied\tnucl_diversity_rarefied_median\tbreadth_rarefied\tconANI_reference\tpopANI_reference\tSNS_count\tSNV_count\tdivergent_site_count\tpopulation_divergent_sites\tconsensus_divergent_sites\tsample\tsample2" >> "$OUTFILE"

# Iterate through all *_scaffold_info.tsv files
find /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/ \
-type f \( -name "*_scaffold_info.tsv" -o -name "*_scaffold_info.tsv.gz" \) | while read file; do

# Extract sample name
sample=$(basename "$file" | sed -E 's/_scaffold_info\.tsv(\.gz)?$//')

# Skip empty files
[ -s "$file" ] || continue

# Use cat or gunzip depending on the file type
if [[ "$file" == *.gz ]]; then
gunzip -c "$file"
else
  cat "$file"
fi | tail -n +2 | \
awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" '
    NF==21 {
      OFS="\t";
      print $0, s1, s2
    }
  ' >> "$OUTFILE"

done

# Expected number of columns: 23 columns, including 21 original columns plus sample and sample2
less all_samples.scaffold_info.tsv
awk -F'\t' '{print NF}' all_samples.scaffold_info.tsv | sort | uniq -c
head -n 3 all_samples.scaffold_info.tsv | column -t
cut -f22 all_samples.scaffold_info.tsv | tail -n +2 | sort | uniq | wc -l
# 5897 sample names



##### Generate the combined gene_info.tsv file
# Set output file path
OUTFILE="all_samples.gene_info.tsv"
> "$OUTFILE"

# Write unified header, including the original gene_info columns plus sample and sample2
echo -e "scaffold\tgene\tgene_length\tcoverage\tbreadth\tbreadth_minCov\tnucl_diversity\tstart\tend\tdirection\tpartial\tdNdS_substitutions\tpNpS_variants\tSNV_count\tSNV_S_count\tSNV_N_count\tSNS_count\tSNS_S_count\tSNS_N_count\tdivergent_site_count\tsample\tsample2" >> "$OUTFILE"

# Iterate through all *_gene_info.tsv files
find /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/ \
-type f \( -name "*_gene_info.tsv" -o -name "*_gene_info.tsv.gz" \) | while read file; do

# Extract sample name
sample=$(basename "$file" | sed -E 's/_gene_info\.tsv(\.gz)?$//')

# Skip empty files
[ -s "$file" ] || continue

# Use cat or gunzip depending on the file type
if [[ "$file" == *.gz ]]; then
gunzip -c "$file"
else
  cat "$file"
fi | tail -n +2 | \
awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" '
    NF==20 {
      OFS="\t";
      print $0, s1, s2
    }
  ' >> "$OUTFILE"

done

# Expected number of columns: 22 columns, including 20 original columns plus sample and sample2
less all_samples.gene_info.tsv
awk -F'\t' '{print NF}' all_samples.gene_info.tsv | sort | uniq -c
head -n 3 all_samples.gene_info.tsv | column -t
cut -f21 all_samples.gene_info.tsv | tail -n +2 | sort | uniq | wc -l
# 5897 sample names



##### Generate the combined mapping_info.tsv file
# Set output file path
OUTFILE="all_samples.mapping_info.tsv"
> "$OUTFILE"

# Write unified header, including the original mapping_info columns plus sample and sample2
echo -e "scaffold\tpass_pairing_filter\tfiltered_pairs\tmean_mapq_score\tunfiltered_priority_reads\tmean_pair_length\tmean_mistmaches\tpass_min_insert\tpass_max_insert\tmean_insert_distance\tmean_PID\tunfiltered_singletons\tpass_min_mapq\tfiltered_singletons\tunfiltered_pairs\tunfiltered_reads\tmedian_insert\tfiltered_priority_reads\tpass_min_read_ani\tsample\tsample2" >> "$OUTFILE"

# Iterate through all *_mapping_info.tsv files
find /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/ \
-type f \( -name "*_mapping_info.tsv" -o -name "*_mapping_info.tsv.gz" \) | while read file; do

# Extract sample name
sample=$(basename "$file" | sed -E 's/_mapping_info\.tsv(\.gz)?$//')

# Skip empty files
[ -s "$file" ] || continue

# Use cat or gunzip depending on the file type
if [[ "$file" == *.gz ]]; then
gunzip -c "$file"
else
  cat "$file"
fi | tail -n +3 | \
awk -F'\t' -v s1="$sample" -v s2="${sample%%_*}" '
    NF==19 {
      OFS="\t";
      print $0, s1, s2
    }
  ' >> "$OUTFILE"

done

# Expected number of columns: 21 columns, including 19 original columns plus sample and sample2
less all_samples.mapping_info.tsv
awk -F'\t' '{print NF}' all_samples.mapping_info.tsv | sort | uniq -c
head -n 3 all_samples.mapping_info.tsv | column -t
cut -f20 all_samples.mapping_info.tsv | tail -n +2 | sort | uniq | wc -l
# 5897 sample names



### Combined inStrain output files
all_samples.gene_info.tsv
all_samples.genome_info.tsv
all_samples.mapping_info.tsv
all_samples.scaffold_info.tsv
all_samples.SNVs.tsv

mkdir -p final_result_20260105
mv all_samples.*.tsv final_result_20260105/
  cd /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/
  
  
  
  ### Generate small test files by randomly selecting 1000 rows
  
  # all_samples.gene_info.tsv -> small
  (head -n 1 all_samples.gene_info.tsv && \
   tail -n +2 all_samples.gene_info.tsv | shuf -n 1000) \
> all_samples.gene_info_small.tsv
less all_samples.gene_info_small.tsv
cat all_samples.gene_info_small.tsv | wc -l
# done

# all_samples.genome_info.tsv -> small
(head -n 1 all_samples.genome_info.tsv && \
  tail -n +2 all_samples.genome_info.tsv | shuf -n 1000) \
> all_samples.genome_info_small.tsv
less all_samples.genome_info_small.tsv
cat all_samples.genome_info_small.tsv | wc -l
# done

# all_samples.mapping_info.tsv -> small
(head -n 1 all_samples.mapping_info.tsv && \
  tail -n +2 all_samples.mapping_info.tsv | shuf -n 1000) \
> all_samples.mapping_info_small.tsv
less all_samples.mapping_info_small.tsv
cat all_samples.mapping_info_small.tsv | wc -l
# done

# all_samples.scaffold_info.tsv -> small
(head -n 1 all_samples.scaffold_info.tsv && \
  tail -n +2 all_samples.scaffold_info.tsv | shuf -n 1000) \
> all_samples.scaffold_info_small.tsv
less all_samples.scaffold_info_small.tsv
cat all_samples.scaffold_info_small.tsv | wc -l
# done

# all_samples.SNVs.tsv -> small
(head -n 1 all_samples.SNVs.tsv && \
  tail -n +2 all_samples.SNVs.tsv | shuf -n 1000) \
> all_samples.SNVs_small.tsv
less all_samples.SNVs_small.tsv
cat all_samples.SNVs_small.tsv | wc -l
# done



#### Check the data
# Count the total number of unique SNV sites detected across all samples in all_samples.SNVs.tsv
# A unique SNV site is defined by scaffold + position, regardless of how many samples it appears in
# Method: awk + sort + uniq on Linux
awk -F'\t' 'NR>1 {print $1"\t"$2}' all_samples.SNVs.tsv \
| sort -u \
| wc -l
# Explanation: $1 = scaffold; $2 = position; sort -u removes duplicates; wc -l counts unique SNV sites
# Result: 226,801,474 SNV sites

## Count the number of scaffolds with detected SNVs
awk -F'\t' 'NR>1 {print $1}' all_samples.SNVs.tsv \
| sort -u \
| wc -l
# Result: 827,157 scaffolds





####### calculation by duckDB 

mkdir -p /mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/duckdb_tmp

duckdb_cli

### Note: when running commands in DuckDB, paste only the commands and do not copy comment lines

PRAGMA threads=160;
PRAGMA temp_directory='/mnt/data/ShuangPeng/Project/SNP/workplace/ELGG_results/instrain/final_result_20260105/duckdb_tmp';
PRAGMA memory_limit='400GB';


############ Filter SNVs with coverage >= 5
{
  PRAGMA threads=64;
  PRAGMA memory_limit='600GB';
  
  CREATE TABLE snvs_filtered AS
  SELECT *
    FROM read_csv(
      'all_samples.SNVs.tsv',
      delim = '\t',
      header = TRUE,
      sample_size = -1
    )
  WHERE
  allele_count >= 2
  AND class IN ('SNV', 'con_SNV', 'pop_SNV')
  AND cryptic = 'False'
  AND con_freq >= 0.05;
  
  COPY snvs_filtered
  TO 'filtered_SNVs3.tsv'
  (HEADER, DELIMITER '\t');
  
  # Linux command
  mv filtered_SNVs3.tsv filtered_SNVs_coverage5.tsv
}


# Check whether the command is still running
watch -n 30 'ls -lh filtered_SNVs_coverage5.tsv'

# Check whether the file exists and whether the file size is stable
ls -lh filtered_SNVs_coverage5.tsv

## Check file size
du -sh filtered_SNVs_coverage5.tsv
# 64G

cut -f20 all_samples.SNVs.tsv | tail -n +2 | sort | uniq | wc -l
cut -f20 filtered_SNVs_coverage5.tsv | tail -n +2 | sort | uniq | wc -l
# 5886 samples

# Count the number of unique SNV sites after filtering
# A unique SNV site is defined by scaffold + position, regardless of how many samples it appears in
awk -F'\t' 'NR>1 {print $1"\t"$2}' all_samples.SNVs.tsv \
| sort -u \
| wc -l

awk -F'\t' 'NR>1 {print $1"\t"$2}' filtered_SNVs_coverage5.tsv \
| sort -u \
| wc -l
# Explanation: $1 = scaffold; $2 = position; sort -u removes duplicates; wc -l counts unique SNV sites
# Result: 121,518,501 SNV sites


## Count the number of scaffolds with detected SNVs
awk -F'\t' 'NR>1 {print $1}' all_samples.SNVs.tsv \
| sort -u \
| wc -l

awk -F'\t' 'NR>1 {print $1}' filtered_SNVs_coverage5.tsv \
| sort -u \
| wc -l
# Result: 753,261 scaffolds 




### Fig. 1b. Genome-wide SNV density, sequencing coverage, and nucleotide diversity across a representative bacterial genome ############

### Select one representative genome from one sample
single_SNVs_df <- fread("SRR1779122_1_kneaddata_paired_sorted_SNVs.tsv", sep = "\t", header = TRUE, fill = TRUE)
head(single_SNVs_df,5)

names(which(table(single_SNVs_df$scaffold) == 12632))
# "GCF_005145085.1_ASM514508v1_genomic_NZ_CP039705.1"

head(single_SNVs_df,5)

single_genome_SNVs_df <- single_SNVs_df %>%
  subset(scaffold%in%c("GCF_005145085.1_ASM514508v1_genomic_NZ_CP039705.1")&
           position<100000)

table(!is.na(single_genome_SNVs_df$position))
# No NA values


### Mean sequencing depth and nucleotide diversity
# mean_coverage: average sequencing depth
# nuc_diversity: mean nucleotide diversity, calculated as 1 - sum(p_i^2)
# snv_density: SNVs per bp

# Set window size
window_size <- 500

genome_summary <- single_genome_SNVs_df %>%
  mutate(window = floor(position / window_size) * window_size) %>%
  group_by(window) %>%
  summarise(mean_coverage = mean(position_coverage, na.rm = TRUE),
            snv_density = sum(class != "SNS") / window_size, 
            nuc_diversity = mean(1 - ((A/position_coverage)^2 + 
                                        (C/position_coverage)^2 +
                                        (T/position_coverage)^2 +
                                        (G/position_coverage)^2), na.rm = TRUE))%>%
  ungroup()

all_windows <- data.frame(window = seq(0, max(single_genome_SNVs_df$position), by = window_size))

# Merge with all genomic windows and replace missing values with 0
genome_summary_full <- all_windows %>%
  left_join(genome_summary, by = "window") %>%
  replace_na(list(mean_coverage = 0, snv_density = 0, nuc_diversity = 0))

# Convert to long format
plot_df <- genome_summary_full %>%
  pivot_longer(cols = c(snv_density, mean_coverage, nuc_diversity),
               names_to = "metric",
               values_to = "value" )

# Plot
plot_df$metric <- factor(plot_df$metric,
                         levels = c("snv_density", "mean_coverage", "nuc_diversity"))

p_SNV_density_coverage_nuc_diversity <- ggplot(plot_df, aes(x = window, y = value)) +
  geom_line(color = "blue", size = 0.6) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(
               snv_density = "SNV density (SNVs/bp)",
               mean_coverage = "Coverage depth",
               nuc_diversity = "Nucleotide diversity"))) +
  theme_bw(base_size = 12) +
  labs(x = "Genome position (bp)", y = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        plot.margin = margin(5, 5, 5, 5))

p_SNV_density_coverage_nuc_diversity 
# save as 6 * 5





### Fig. 1c. Position coverage vs. number of SNVs ############

### Filtering based on coverage > 5

## Before filtering: 827,157 scaffolds with 226,801,474 SNV sites
## After filtering: 753,261 scaffolds with 121,518,501 SNV sites

## Calculate the fraction of retained and filtered SNVs:
121518501/226801474
(226801474-121518501)/226801474
# 0.5357924 = 53.6%
# 0.4642076 = 46.4%
# Result: 53.6% of SNVs were retained, while 46.4% were filtered out

## 3D pie plot
library(plotrix)

# SNV counts
num_data <- c(
  Retained = 121518501,
  Filtered = 226801474 - 121518501
)

# Percentages for labels
percent <- round(num_data / sum(num_data) * 100, 1)

labels <- paste0(
  names(num_data), "\n",
  percent, "%"
)

p_SNV_filtering_summary  <- pie3D(
  num_data,
  labels = labels,
  explode = 0.25,
  col = c("#762a83", "grey60")
  # main = "Impact of SNV filtering",
  # labelcex = 1.2,
  # start = 60,
  # edgecol = "black"
)

p_SNV_filtering_summary
# save as 3 * 3







### Fig. 1d. Relationship between median SNV coverage and genome-wide sequencing depth ############

###### On the Linux server:
### Add genome information to filtered_SNVs.tsv and name the new column as genome_sp

duckdb_cli

SET memory_limit='400GB';
SET threads=64;

COPY (
  SELECT
  snv.*,
  stb.genome AS genome_sp
  FROM read_csv(
    'filtered_SNVs_coverage5.tsv',
    delim='\t',
    header=true
  ) AS snv
  LEFT JOIN read_csv(
    'stb_scaffold_genome_name.tsv',
    delim='\t',
    header=true
  ) AS stb
  ON snv.scaffold = stb.scaffold
)
TO 'filtered_SNVs_with_genome_coverage5.tsv'
(HEADER, DELIMITER '\t');

# All rows from filtered_SNVs_coverage5.tsv are retained
# If scaffold can be matched, genome_sp is assigned as genome
# If scaffold cannot be matched, genome_sp is set as NULL

# Check whether the new column exists
head -n 3 filtered_SNVs_with_genome_coverage5.tsv

# Check whether NA values exist, which is expected for scaffolds not present in the scaffold-to-genome table
cut -f1,17 filtered_SNVs_with_genome_coverage5.tsv | head



### Calculate median position coverage for each genome in each sample
SET memory_limit='300GB';
SET threads=64;

COPY (
  SELECT
  sample2,
  genome_sp,
  median(position_coverage) AS median_coverage
  FROM read_csv(
    'filtered_SNVs_with_genome_coverage5.tsv',
    delim = E'\t',
    header = true
  )
  WHERE genome_sp IS NOT NULL
  GROUP BY
  sample2,
  genome_sp
)
TO 'sample_genome_snp_median_coverage_coverage5.tsv'
(HEADER, DELIMITER '\t');



### Add genome-wide coverage from all_samples.genome_info.tsv to sample_genome_snp_median_coverage.tsv
## The join is based on matching genome in all_samples.genome_info.tsv and genome_sp in sample_genome_snp_median_coverage.tsv

COPY (
  SELECT
  s.sample2,
  s.genome_sp,
  s.median_coverage,
  g.coverage AS genome_coverage
  FROM read_csv(
    'sample_genome_snp_median_coverage_coverage5.tsv',
    delim = E'\t',
    header = true
  ) AS s
  LEFT JOIN (
    SELECT
    genome,
    sample2,
    coverage
    FROM read_csv_auto(
      'all_samples.genome_info.tsv',
      delim = E'\t',
      header = true
    )
  ) AS g
  ON s.genome_sp = g.genome
  AND s.sample2   = g.sample2
)
TO 'sample_genome_snp_median_coverage_coverage5_with_genome_cov.tsv'
(HEADER, DELIMITER '\t');



#### R code:

### Use coverage-filtered SNVs
snp_median_coverage_with_genome <- read_tsv(
  "sample_genome_snp_median_coverage_coverage5_with_genome_cov.tsv",
  show_col_types = FALSE
)

View(snp_median_coverage_with_genome)
colnames(snp_median_coverage_with_genome)
dim(snp_median_coverage_with_genome)

cor.test(
  log10(snp_median_coverage_with_genome$genome_coverage),
  log10(snp_median_coverage_with_genome$median_coverage),
  method = "spearman"
) 
# p-value < 2.2e-16

cor(
  log10(snp_median_coverage_with_genome$genome_coverage),
  log10(snp_median_coverage_with_genome$median_coverage),
  method = "spearman"
) 
# R = 0.5255964

# Median SNV coverage was strongly correlated with genome-wide sequencing depth 
# Spearman's R = 0.53

# Plot
## Randomly sample 30,000 rows for plotting
set.seed(123)

plot_df <- snp_median_coverage_with_genome %>%
  slice_sample(n = 30000)

p_snp_median_coverage_with_genome <- ggplot(plot_df, aes(
  x = genome_coverage,
  y = median_coverage
)) +
  geom_point(alpha = 0.3, size = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Genome-wide coverage (log10)",
    y = "Median SNV position coverage (log10)"
  ) +
  theme_bw()

p_snp_median_coverage_with_genome
# Each dot represents one genome in one sample
# save as 3.5 * 3






### Fig. 1e. Distribution of sequencing coverage across all detected SNVs ############

##### On the Linux server:

### Calculate the coverage distribution of SNV sites
# Main QC question: Are SNVs supported by sufficient sequencing depth rather than random noise?
# Dot plot: x = position_coverage, y = number of SNVs on a log scale
# Expected pattern: most SNVs have coverage above the filtering threshold and are not concentrated only at the threshold edge

## First calculate the frequency of position_coverage directly without considering sample2/scaffold/position
{duckdb <<'EOF'
COPY (
  SELECT
  position_coverage,
  COUNT(*) AS snp_count
  FROM read_csv(
    'SNVs_ATCG_subset_coverage5.tsv',
    delim = '\t',
    header = TRUE
  )
  GROUP BY
  position_coverage
)
TO 'SNP_coverage_distribution_coverage5.tsv'
(HEADER, DELIMITER '\t');
EOF

du -sh *
  # 100K
} 
  
##### R code:
  
### Use coverage-filtered SNVs
  cov_dist <- read_tsv(
    "SNP_coverage_distribution_coverage5.tsv",
    show_col_types = FALSE
  )

# Summary statistics
summary(cov_dist$position_coverage)

# coverage >= 10
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10    3207    6444   10075   11138   87402 

# coverage >= 5
# Min. 1st Qu.  Median    Mean  3rd Qu.   Max. 
# 5     1237     2472     3125    3989    87402

# coverage >= 5
# Min. 1st Qu.  Median    Mean  3rd Qu.   Max. 
# 5     3203     6442     10071   11133   87402 

median(cov_dist$position_coverage)
mean(cov_dist$position_coverage)
min(cov_dist$position_coverage)
max(cov_dist$position_coverage)

quantile(cov_dist$position_coverage, c(0.25, 0.5, 0.75))

summary(cov_dist$snp_count)
median(cov_dist$snp_count)
mean(cov_dist$snp_count)
min(cov_dist$snp_count)
max(cov_dist$snp_count)


p_SNP_coverage_distribution <- ggplot(cov_dist, aes(x = position_coverage, y = snp_count)) +
  geom_point(size = 1, alpha = 0.2, color="#252525") +
  geom_vline(xintercept = 5,
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  scale_y_log10() +
  # scale_x_log10() +
  # coord_cartesian(ylim = c(1, 1000000)) +
  labs(x = "SNV coverage (reads)",
       y = "Number of SNVs (log10)"
       # title = "Coverage distribution of SNVs",
       # subtitle = "Dashed line indicates coverage threshold (>=5)"
  ) +
  theme_bw()

p_SNP_coverage_distribution
# save as 3.5 * 3






### Fig. 1f. Distribution of SNV richness across 5,886 samples ############

#### On the Linux server:

##### Calculate SNV richness for each sample and scaffold
duckdb <<'EOF'
COPY (
  SELECT
  sample2,
  scaffold,
  COUNT(DISTINCT position) AS snv_richness
  FROM read_csv(
    'SNVs_ATCG_subset_coverage5.tsv',
    delim = '\t',
    header = TRUE,
    sample_size = -1
  )
  GROUP BY
  sample2,
  scaffold
)
TO 'SNV_richness_coverage5_sample_scaffold.tsv'
(HEADER, DELIMITER '\t');
EOF

## Explanation:
# read_csv(...): read the filtered SNV file
# GROUP BY sample2, scaffold: group rows by sample and scaffold
# COUNT(DISTINCT position): count the number of distinct SNV positions for each sample-scaffold pair
# COPY ... TO: export the result as a TSV file named SNV_richness_coverage5_sample_scaffold.tsv

# Linux command
du -sh *
  # 1.1 G
  
  
  
  ##### R code:
  
  ### Alpha diversity: observed SNV richness
  library(tidyverse)

snv_scaffold <- read_tsv(
  "SNV_richness_coverage5_sample_scaffold.tsv",
  show_col_types = FALSE
)

# Sample-level SNV richness
# This is the most commonly used metric in the main analysis
# It represents the total number of SNV sites in one sample across all scaffolds
# Sum SNV richness across all scaffolds for each sample
snv_sample <- snv_scaffold %>%
  group_by(sample2) %>%
  summarise(
    snv_richness = sum(snv_richness),
    .groups = "drop"
  )

dim(snv_sample)
# 5886

range(snv_sample$snv_richness)
# 1-1470043

sum(snv_sample$snv_richness)

summary(snv_sample$snv_richness)

# coverage >= 5
# Min. 1st Qu.  Median    Mean  3rd Qu.     Max. 
# 1     11890    26771     52203   62449     1470043 

# Plot 1: distribution of SNV richness across samples
p_SNV_richness_per_samples <- ggplot(snv_sample, aes(x = snv_richness)) +
  geom_histogram(
    bins = 100,
    fill = "#40004b",
    color = "white",
    linewidth=0.05,
    alpha=0.75
  ) +
  scale_x_log10() +
  labs(
    x = "SNV richness per sample (log10)",
    y = "Number of samples"
    # title = "Distribution of SNV richness across samples"
  ) +
  # xlim(1,1470043) +
  theme_bw()

p_SNV_richness_per_samples
# save as 3.5 * 3





### Fig. 1g. Distribution of detected SNVs across scaffolds ############

# Scaffold-level SNV richness, used for supplementary analysis
snv_scaffold #

# Calculate how SNV-rich each scaffold is
snv_scaffold_level <- snv_scaffold %>%
  group_by(scaffold) %>%
  summarise(
    mean_snv_richness = mean(snv_richness),
    .groups = "drop"
  )

head(snv_scaffold_level, 10)
dim(snv_scaffold_level)
# 753261 scaffolds

range(snv_scaffold_level$mean_snv_richness)
# 1.00-32515.28

summary(snv_scaffold_level$mean_snv_richness)

# coverage >= 10
# Min.   1st Qu.   Median   Mean   3rd Qu.   Max. 
# 1.000  2.000     3.333    6.799  5.714     32355.412 

# coverage >= 5
# Min.   1st Qu.   Median   Mean   3rd Qu.   Max. 
# 1.000  2.000     3.378    6.876  5.750     32515.278 

# test <- snv_scaffold_level %>% subset(mean_snv_richness > 16.353)
# summary(test$mean_snv_richness)

## Plot 2: distribution of SNV richness across scaffolds

# snv_scaffold_level <- snv_scaffold_level %>%
#   mutate(log10_mean_snv = log10(mean_snv_richness))

# ggplot(snv_scaffold_level, aes(x = log10_mean_snv)) +
#   geom_histogram(bins = 100, fill = "#40004b", color = "white", alpha = 0.75) +
#   labs(x = "Mean SNV richness per scaffold (log10)",
#        y = "Number of scaffolds",
#        title = "SNV richness across scaffolds") +
#   theme_bw()

# ggplot(snv_scaffold_level, aes(x = mean_snv_richness)) +
#   geom_density(fill = "#40004b", alpha = 0.5) +
#   scale_x_log10() +
#   labs(x = "Mean SNV richness per scaffold (log10)",
#        y = "Density",
#        title = "Distribution of SNV richness across scaffolds") +
#   theme_bw()

dim(snv_scaffold_level)

p_SNV_richness_per_scaffold <- ggplot(snv_scaffold_level, aes(x = mean_snv_richness)) +
  geom_histogram(bins = 100, fill = "#40004b", color = "white",
                 linewidth=0.05, alpha=0.75) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Mean SNV richness per scaffold (log10)",
    y = "Number of scaffolds (log10)"
    # title = "SNV richness across scaffolds"
  ) +
  theme_bw()

p_SNV_richness_per_scaffold
# save as 3.5 * 3





### Fig. 1h. Proportions of SNVs by mutation type ############

###### On the Linux server:

##### Calculate proportions of mutation types
## Biological rationale:
## This analysis evaluates whether detected SNVs show biologically plausible mutation-type patterns rather than random sequencing errors.
## Bar/pie plot: synonymous, nonsynonymous, intergenic and multiple mutation types
## Expected pattern: synonymous SNVs should generally be more frequent than nonsynonymous SNVs.
## If these were random errors, mutation types would not show such biological constraint.

# Check whether mutation_type has NA or empty values
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="mutation_type") col=i} 
NR>1 && ($col=="" || $col=="NA") {print}' filtered_SNVs_coverage5.tsv | head
# Output means NA or empty values exist; no output means all values are valid
# Result: output exists

# Replace NA or empty mutation_type values with "non_gene"
# Export key columns into a smaller file, with one row per mutation
duckdb_cli

COPY (
  SELECT 
  CASE 
  WHEN mutation_type IS NULL OR mutation_type = '' THEN 'non_gene'
  ELSE mutation_type
  END AS mutation_type,
  sample2,
  class,
  scaffold,
  position,
  position_coverage
  FROM read_csv(
    'filtered_SNVs_coverage5.tsv',
    delim = '\t',
    header = TRUE,
    sample_size = -1
  )
) 
TO 'SNP_coverage5_mutation_types.tsv' 
(HEADER, DELIMITER '\t');

# COPY (
#     SELECT mutation_type, sample2, class, scaffold, position, position_coverage
#     FROM read_csv(
#     'filtered_SNVs_coverage5.tsv',
#     delim = '\t',
#     header = TRUE,
#     sample_size = -1
#   )
#     WHERE mutation_type IS NOT NULL
# ) 
# TO 'SNP_coverage5_mutation_types.tsv' 
# (HEADER, DELIMITER '\t');

# The file SNP_coverage5_mutation_types.tsv has been generated from this command.
# Next, use DuckDB to calculate the proportions of mutation_type and class across all samples,
# and then generate pie plots in R.

# Check whether mutation_type still has NA or empty values
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="mutation_type") col=i} 
NR>1 && ($col=="" || $col=="NA") {print}' SNP_coverage5_mutation_types.tsv | head
# Output means NA or empty values exist; no output means all values are valid
# Result: no output

### Use DuckDB to calculate the overall proportion of mutation types
duckdb_cli

COPY (
  SELECT
  mutation_type,
  COUNT(*) AS snv_count,
  COUNT(*) * 1.0 / SUM(COUNT(*)) OVER () AS proportion
  FROM read_csv(
    'SNP_coverage5_mutation_types.tsv',
    delim = '\t',
    header = true
  )
  WHERE mutation_type IS NOT NULL
  GROUP BY mutation_type
)
TO 'mutation_type_proportion_coverage5.tsv'
(HEADER, DELIMITER '\t');



# Check whether class has NA or empty values
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="mutation_type") col=i} 
NR>1 && ($col=="" || $col=="NA") {print}' SNP_coverage5_mutation_types.tsv | head
# Output means NA or empty values exist; no output means all values are valid
# Result: no output



###### R code:

### Percentage of synonymous, nonsynonymous, intergenic and multiple SNVs
mutation_type_proportion <- read_tsv(
  "mutation_type_proportion_coverage5.tsv",
  show_col_types = FALSE
)

mutation_type_proportion
colnames(mutation_type_proportion)
dim(mutation_type_proportion)

mutation_type_proportion$mutation_type <- ifelse(mutation_type_proportion$mutation_type=="N","Nonsynonymous",
                                                 ifelse(mutation_type_proportion$mutation_type=="S","Synonymous",
                                                        ifelse(mutation_type_proportion$mutation_type=="I","Intergenic",  
                                                               ifelse(mutation_type_proportion$mutation_type=="M","Multiple",NA
                                                               ))))

# Prepare labels and order
mutation_type_proportion_df <- mutation_type_proportion %>%
  arrange(desc(proportion)) %>%
  mutate(
    mutation_type = factor(mutation_type, levels = mutation_type),
    label = paste0(
      mutation_type, "\n",
      scales::percent(proportion, accuracy = 0.1)
    )
  )

# Donut chart
p_mutation_type_proportion <- ggplot(mutation_type_proportion_df, aes(x = 2, y = proportion, fill = mutation_type)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "white"
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Synonymous" = "#8dd3c7",
      "Nonsynonymous" = "#fdb462",
      "Intergenic" = "#b3de69",
      "Multiple" = "#8C5BAA"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

p_mutation_type_proportion
# save as 3 * 3





### Fig. 1i. Proportions of SNVs by variant class ############

###### On the Linux server:

# Check whether class has NA or empty values
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="mutation_type") col=i} 
NR>1 && ($col=="" || $col=="NA") {print}' SNP_coverage5_mutation_types.tsv | head
# Output means NA or empty values exist; no output means all values are valid
# Result: no output

### Use DuckDB to calculate the overall proportion of SNV classes
duckdb_cli

COPY (
  SELECT
  class,
  COUNT(*) AS snv_count,
  COUNT(*) * 1.0 / SUM(COUNT(*)) OVER () AS proportion
  FROM read_csv(
    'SNP_coverage5_mutation_types.tsv',
    delim = '\t',
    header = true
  )
  WHERE class IS NOT NULL
  GROUP BY class
)
TO 'snp_class_proportion_coverage5.tsv'
(HEADER, DELIMITER '\t');



###### R code:

snp_class_proportion <- read_tsv(
  "snp_class_proportion_coverage5.tsv",
  show_col_types = FALSE
)

snp_class_proportion
colnames(snp_class_proportion)
dim(snp_class_proportion)

# Prepare labels and order
snp_class_proportion_df <- snp_class_proportion %>%
  arrange(desc(proportion)) %>%
  mutate(
    classe = factor(class, levels = class),
    label = paste0(
      class, "\n",
      scales::percent(proportion, accuracy = 0.1)
    )
  )

# Donut chart
p_snp_class_proportion <- ggplot(snp_class_proportion_df, aes(x = 2, y = proportion, fill = class)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "white"
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "SNV" = "#ffffb3",
      "con_SNV" = "#bebada",
      "pop_SNV" = "#fb8072"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

p_snp_class_proportion
# save as 3 * 3





### Fig. 1j. Distribution of consensus allele frequencies across SNV categories ############

###### On the Linux server:

duckdb_cli

COPY (
  SELECT
  class,
  FLOOR(con_freq * 100) / 100.0 AS con_freq_bin,
  COUNT(*) AS snv_count
  FROM read_csv(
    'filtered_SNVs_coverage5.tsv',
    delim = '\t',
    header = true
  )
  WHERE con_freq IS NOT NULL
  AND class IS NOT NULL
  GROUP BY
  class,
  con_freq_bin
)
TO 'con_freq_distribution_by_class_coverage5.tsv'
(HEADER, DELIMITER E'\t');

# FLOOR rounds values down
# con_freq_bin represents the binned consensus allele frequency
# In other words, it indicates which frequency interval each SNV falls into
# Directly using continuous con_freq values is not suitable because nearly every site may have a unique decimal value,
# making it difficult to visualize the distribution shape



###### R code:

### Allele-frequency structure, a core feature of microdiversity
# Consensus allele frequency distribution
# Main QC question: Do SNVs show a reasonable population-frequency structure rather than sequencing-error-like patterns?
# x = con_freq
# y = SNV count
# facet = SNV/con_SNV/pop_SNV
# Expected pattern: sequencing errors would cluster at extremely low frequencies, whereas true SNVs should show structured allele-frequency distributions.

con_freq_distribution_by_class <- read_tsv(
  "con_freq_distribution_by_class_coverage5.tsv", show_col_types = FALSE)

con_freq_distribution_by_class$class <- factor(con_freq_distribution_by_class$class,
                                               levels = c("SNV","con_SNV","pop_SNV")) 

con_freq_distribution_by_class

# Bar plot 
p_con_freq_distribution_by_class <- ggplot(
  con_freq_distribution_by_class,
  aes(x = con_freq_bin, y = snv_count, group=class,fill=class)) +
  geom_col(width = 0.01 ) +
  facet_wrap(~ class, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)) +
  labs(x = "Consensus allele frequency",
       y = "SNV count" ) +
  scale_fill_manual(values = c(
    "SNV" = "#ffffb3",
    "con_SNV" = "#bebada",
    "pop_SNV" = "#fb8072"))+
  theme_bw()+
  theme( legend.position = "none")

p_con_freq_distribution_by_class
# save as 9.5 * 4     





### Supplementary Fig. S1. Distribution and summary of assembled bacterial genomes based on genomic length and quality ############


#### Distribution of assembled microbial genomes based on genomic length and quality 

# genome_species_info.xlsx
library(readxl)
library(dplyr)

# Read Excel file; the first sheet is used by default
genome_species_info <- read_excel("D:/3_Projects/2_Children_jaundice/3_R_analysis/genome_species_info.xlsx")
class(genome_species_info)
genome_species_info <- data.frame(genome_species_info)

# Inspect the first few rows
head(genome_species_info,10)

## Step 1: classify each genome by genome quality
library(data.table)
library(ggplot2)
library(dplyr)

colnames(genome_species_info)
dt <- as.data.table(genome_species_info)

range(dt$Completeness)
range(dt$Contamination)
range(dt$quality)

# Genome quality classification
dt[, quality := case_when(
  Completeness >= 90 & Contamination < 5 ~ "Complete",
  Completeness >= 75 ~ "High-quality",
  TRUE ~ "Medium-quality"
)]

## Step 2: convert genome length to log10 scale for the x axis
dt[, genome_length_log10 := log10(Genome.size)]

## Step 3: plot stacked histogram of genome length by quality category
summary(dt$Genome.size)

p_genome_length_hist <- ggplot(dt, 
                               aes(x = genome_length_log10, fill = quality)) +
  geom_histogram(binwidth = 0.1, color = "white", size = 0.2) +
  scale_fill_manual(
    values = c(
      "Complete" = "#D55E00",
      "High-quality" = "#56B4E9",
      "Medium-quality" = "#009E73"
    )
  ) +
  labs(
    x = "Genome length (log10, bp)",
    y = "Number of genomes",
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "top"
  )

p_genome_length_hist
# save as 3.5 * 3


### Pie plot
library(ggplot2)
library(dplyr)

# Count the number of genomes in each quality category
df_pie <- dt %>%
  count(quality) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(sprintf("%.1f", percent), "%")
  )

# Define colors
quality_colors <- c(
  "Complete" = "#D55E00",
  "High-quality" = "#56B4E9",
  "Medium-quality" = "#009E73"
)

# Pie plot
p_completeness_pie <- ggplot(df_pie, aes(x = "", y = n, fill = quality)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = quality_colors) +
  labs(
    title = "Summary statistics for genome completeness",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

p_completeness_pie
# save as 4 * 3

