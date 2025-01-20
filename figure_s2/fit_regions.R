library(tidyverse)
library(readxl)
library(nlme)
rm(list=ls())

fit_model = function(data) {
  fit = lme(Value ~ Condition, random=~1|ID, data = data[[1]])
  fixed = as_tibble(intervals(fit, which='fixed')$fixed) %>% rename(est='est.')
  fixed$Parameter = c('Basal', 'Delta')
  fixed = fixed %>% 
    pivot_wider(names_from=Parameter, values_from=c(lower, 'est', 'upper'))
  return(fixed)
}

# Load in data
df = read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')

# Run mixed model on each region for relative PET
df_sum = df %>% 
  filter(Modality %in% c('cmrglc')) %>%
  group_by(Modality, Region) %>% 
  nest() %>%
  summarise(coefs=fit_model(data)) %>%
  unnest(coefs)

write_csv(df_sum, 'regional_cmrglc_deltas.csv')