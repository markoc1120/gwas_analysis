library(tidyverse)

ibd <- read.table(snakemake@input[['genome']], header = TRUE)
members <- ibd$FID1
members <- unique(members)
write.table(cbind(members,members), file = snakemake@output[['ibd_fails']], col.names = F, row.names = F)