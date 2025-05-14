library(tidyverse)
library(ggplot2)
library(gridExtra)

count_lines <- function(file_path) {
  if(!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NA)
  }

  result <- system(paste("wc -l", file_path), intern = TRUE)
  as.numeric(strsplit(result, " +")[[1]][1])
}

initial_num_variants <- count_lines("initial/gwas_data.bim")

qc_steps <- tribble(
  ~step, ~description, ~samples, ~variants,
  0, "Initial", count_lines("initial/gwas_data.fam"), initial_num_variants,
  1, "After sex check", count_lines("qc/QC_output_sex_filtered.fam"), count_lines("qc/QC_output_sex_filtered.bim"),
  2, "After heterozygosity filter", count_lines("qc/QC_output_het_filtered.fam"), count_lines("qc/QC_output_het_filtered.bim"),
  3, "After relatedness check", count_lines("qc/QC_output_sample_final.fam"), count_lines("qc/QC_output_sample_final.bim"),
  4, "Final", count_lines("qc_final/QC_output_final.fam"), count_lines("qc_final/QC_output_final.bim")
)

group_data <- tribble(
  ~group, ~samples, ~variants_before, ~variants_after,
  "Group 1", count_lines("qc/QC_output_geno_filtered_group_1.fam"), initial_num_variants, count_lines("qc/QC_output_geno_filtered_group_1.bim"),
  "Group 2", count_lines("qc/QC_output_geno_filtered_group_2.fam"), initial_num_variants, count_lines("qc/QC_output_geno_filtered_group_2.bim"),
  "Group 3", count_lines("qc/QC_output_geno_filtered_group_3.fam"), initial_num_variants, count_lines("qc/QC_output_geno_filtered_group_3.bim"),
  "Group 4", count_lines("qc/QC_output_geno_filtered_group_4.fam"), initial_num_variants, count_lines("qc/QC_output_geno_filtered_group_4.bim"),
  "Group 5", count_lines("qc/QC_output_geno_filtered_group_5.fam"), initial_num_variants, count_lines("qc/QC_output_geno_filtered_group_5.bim"),
)


p1 <- ggplot(qc_steps, aes(x = factor(step), y = samples)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = samples), vjust = -0.5, size = 3.5) +
  labs(x = "QC Step", y = "Number of Samples",
       title = "Sample Count Through QC Steps") +
  theme_minimal() +
  scale_x_discrete(labels = c("Initial", "Sex Check", "Heterozygosity", "Relatedness", "Final"))

p2 <- ggplot(group_data, aes(x = group, y = samples)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = samples), vjust = -0.5, size = 3.5) +
  labs(x = "Group", y = "Number of Samples",
       title = "Sample Size per Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

group_data_long <- group_data %>%
  pivot_longer(cols = starts_with("variants"),
               names_to = "filter_status",
               values_to = "variant_count") %>%
  mutate(filter_status = ifelse(filter_status == "variants_before", "Before Filtering", "After Filtering"))

p3 <- ggplot(group_data_long, aes(x = group, y = variant_count, fill = filter_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = variant_count), position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(x = "Group", y = "Number of Variants",
       title = "Variants Before and After Filtering by Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Before Filtering" = "lightblue", "After Filtering" = "coral"))

group_data <- group_data %>%
  mutate(percent_retained = variants_after / variants_before * 100)


grid.arrange(p1, p2, p3, ncol = 2)

ggsave("plots/qc_visualization.png", arrangeGrob(p1, p2, p3, ncol = 2),
       width = 12, height = 14, dpi = 300)
