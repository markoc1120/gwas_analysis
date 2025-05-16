library(tidyverse)
library(gridExtra)

pca <- read_table2(snakemake@input[['eigenvec']], col_names = FALSE)
eigenval <- scan(snakemake@input[['eigenval']])

pca <- pca[, -1]
names(pca)[1] <- 'ind'
names(pca)[2:ncol(pca)] <- paste0('PC', 1:(ncol(pca) - 1))
pve <- data.frame(PC = 1:20, pve = eigenval / sum(eigenval) * 100)

pca_plot <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(size = 3) +
  theme_light() +
  xlab(paste0('PC1 (', signif(pve$pve[1], 3), '%)')) +
  ylab(paste0('PC2 (', signif(pve$pve[2], 3), '%)'))

scree_plot <- ggplot(pve, aes(PC, pve)) +
  geom_bar(stat = 'identity') +
  ylab('Percentage variance explained') +
  theme_light()


grid.arrange(scree_plot, pca_plot, ncol = 2)
ggsave(snakemake@output[['pca_plot']], arrangeGrob(scree_plot, pca_plot, ncol = 2),
       width = 12, height = 7, dpi = 300)
