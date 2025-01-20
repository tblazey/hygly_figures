library(tidyverse)
library(readxl)

col_order = c(
  1, 2, 12, 3, 13, 7, 17, 5, 15, 4, 14, 6, 16, 19, 11, 8, 10, 9, 18
)

regional_data = read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')
regional_bool = regional_data  %>% 
  filter(Region == 'Deep White Matter' & Modality != 'Task.Retrieval') %>% 
  mutate_at('Value', .funs=is.numeric) %>%
  dplyr::select(-Region, -Visit.Order) %>% 
  group_by(ID, Condition, Modality) %>% 
  summarise_all(any) %>%
  pivot_wider(
    names_from=c(Modality, Condition), values_from=Value, values_fill=FALSE
  ) %>% ungroup()
demo_table = regional_bool[, col_order]

# Get number of subjects with data at each condition
demo_sums = colSums(demo_table[, -1])

# Add totals to text and replace bools with checkmarks/empties
demo_table[, -1] = sapply(demo_table[, -1], as.character)
demo_table = rbind(demo_table, c('Totals', as.character(demo_sums)))
demo_table[demo_table == 'TRUE'] = '\U2713'
demo_table[demo_table == 'FALSE'] = ''

# Write out table
write.csv(demo_table, 'table_s3.csv', quote=FALSE, row.names=FALSE)

