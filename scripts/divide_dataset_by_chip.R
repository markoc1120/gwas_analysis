library(tidyverse)

metadata <- read_tsv(snakemake@input[['metadata']])

output_files <- snakemake@output[['group_files']]
output_dir <- dirname(output_files[1])

metadata <- metadata %>%
  mutate(
    group = if_else(
      chip == '' & chip_version == '',
      paste0('source: ', source),
      paste0('chip: ', chip)
    )
  )

user_groups <- metadata %>%
  group_by(group) %>%
  summarise(
    users = list(user),
    n_users = n(),
    .groups = 'drop'
  ) %>%
  filter(n_users >= 20)


group_numbers <- c(1, 2, 3, 4, 5)
file_map <- data.frame(
  number = group_numbers,
  filepath = output_files,
  stringsAsFactors = FALSE
)

for (i in seq_along(group_numbers)) {
  output_file <- file_map$filepath[i]
  if (i <= nrow(user_groups)) {
    users <- user_groups$users[[i]]
    lines <- paste(users, users, sep = '\t')
    writeLines(lines, con = output_file)
    message(paste('Created group file', output_file, 'with', length(users), 'users'))
  }
}