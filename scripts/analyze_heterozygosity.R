library(tidyverse)
library(dplyr)

d_miss <- read.table(snakemake@input[['imiss']],header=T)
d_het <- read.table(snakemake@input[['het']],header=T)

d_het <- d_het %>% 
  mutate(HET = (N.NM. - O.HOM.)/ N.NM.)
df <- dplyr::inner_join(d_miss,d_het)

na_rows <- df[!complete.cases(df), ]
df <- drop_na(df)

mean_het <- mean(df$HET)
sd_het <- sd(df$HET)
r_tail <- mean_het + 3*sd_het
l_tail <- mean_het - 3*sd_het

ggplot(data=df) +
  geom_point(mapping=aes(x=HET, y=F_MISS)) +
  geom_vline(mapping = aes(xintercept = r_tail))+
  geom_vline(mapping = aes(xintercept = l_tail))

filtered_df <- df %>% filter(
  HET < r_tail, 
  HET > l_tail,
  # F_MISS <= 0.4
)

filtered_out_entries <- anti_join(df, filtered_df)
filtered_out_entries <- bind_rows(filtered_out_entries, na_rows)
write.table(filtered_out_entries[,1:2], file = snakemake@output[['outliers']], col.names = F, row.names = F)