library(tidyverse)

heights <- read.table(snakemake@input[['height']], head=T)
covariates <- read.table(snakemake@input[['covariates']], head=T)

colnames(heights) <- c('FID', 'IID', 'HEIGHT')
colnames(covariates) <- c('FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

covariates <- covariates %>% semi_join(heights, by = 'FID')
heights <- heights %>% semi_join(covariates, by = 'FID')
heights <- heights %>% unique()
all_samples <- heights %>% select('FID', 'IID')

heights <- heights %>% arrange(FID, IID)
covariates <- covariates %>% arrange(FID, IID)
all_samples <- all_samples %>% arrange(FID, IID)

set.seed(1)
train_indices <- sample(1:nrow(heights), 0.7 * nrow(heights))

write_split_data <- function(data, train_indices, name, output_dir = 'pgs') {
  train_data <- data[train_indices,]
  test_data <- data[-train_indices,]

  write.table(train_data, file.path(output_dir, paste0('train_', name, '.txt')),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

  write.table(test_data, file.path(output_dir, paste0('test_', name, '.txt')),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_split_data(all_samples, train_indices, 'samples')
write_split_data(heights, train_indices, 'heights')
write_split_data(covariates, train_indices, 'covariates')
