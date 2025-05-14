library(qqman)
library(tidyverse)

linear <- read.table(snakemake@input[['gwas_results']], head = TRUE)
linear_clean <- linear %>%
  filter(TEST == 'ADD') %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P)) %>%
  filter(is.finite(P)) %>%
  filter(P > 0)

write.table(linear_clean,
            file = snakemake@output[['linear_clean']],
            sep = '\t',
            quote = FALSE,
            row.names = FALSE)



mlma <- read.table(snakemake@input[['mlma_results']], head = TRUE)
colnames(mlma) <- c('CHR', 'SNP', 'BP', 'A1', 'A2', 'FREQ', 'B', 'SE', 'P')
mlma_clean <- mlma %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P)) %>%
  filter(is.finite(P)) %>%
  filter(P > 0)

write.table(mlma_clean,
            file = snakemake@output[['mlma_clean']],
            sep = '\t',
            quote = FALSE,
            row.names = FALSE)
