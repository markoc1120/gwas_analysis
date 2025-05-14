library(tidyverse)
library(ggplot2)

test_heights <- read.table(snakemake@input[['test_heights']], col.names = c('FID', 'IID', 'HEIGHT'))
test_covariates <- read.table(snakemake@input[['test_covariates']], col.names = c('FID', 'IID', paste0('PC', 1:10)))

pgs_scores <- list(
  '0.01' = read.table(snakemake@input[['pgs_profiles']][1], header = TRUE),
  '0.05' = read.table(snakemake@input[['pgs_profiles']][2], header = TRUE),
  '0.1' = read.table(snakemake@input[['pgs_profiles']][3], header = TRUE)
)

merged <- merge(test_heights, test_covariates)
merged_data <- list()
model_results <- data.frame(
  threshold = c('0.01', '0.05', '0.1'),
  pgs_r2 = numeric(3)
)

for (i in 1:3) {
  threshold <- model_results$threshold[i]
  merged_data[[threshold]] <- merge(pgs_scores[[threshold]], merged)

  full_model <- lm(
    HEIGHT ~ SCORE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = merged_data[[threshold]]
  )
  full_r2 <- summary(full_model)$r.squared

  baseline_model <- lm(
    HEIGHT ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = merged_data[[threshold]]
  )
  baseline_r2 <- summary(baseline_model)$r.squared

  model_results$pgs_r2[i] <- full_r2 - baseline_r2
}

best_idx <- which.max(model_results$pgs_r2)
best_threshold <- model_results$threshold[best_idx]
best_merged <- merged_data[[best_threshold]]

height_prediction_plot <- ggplot(best_merged, aes(x = SCORE, y = HEIGHT)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm', se = FALSE, color = 'red') +
  labs(
    title = paste0('Height ~ PGS (p < ', best_threshold, ')'),
    subtitle = paste0('PGS contribution to R² = ', round(model_results$pgs_r2[best_idx], 4)),
    x = 'Polygenic Score (PGS)',
    y = 'Measured Height (cm)'
  ) +
  theme_minimal()

ggsave(
  filename = snakemake@output[['evaluation_plot']],
  plot = height_prediction_plot,
  device = 'png',
  width = 7,
  height = 5,
  units = 'in',
  dpi = 300
)
