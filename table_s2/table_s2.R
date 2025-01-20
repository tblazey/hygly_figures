library(tidyverse)
library(nlme)
library(gt)
library(readxl)
rm(list=ls())

# Function to get mixed model estimates at each region
fit_model = function(data) {
  fit = lme(Value ~ Condition, random=~1|ID, data = data[[1]])
  fixed = as_tibble(intervals(fit, which='fixed')$fixed) %>% rename(est='est.')
  fixed$Parameter = c('Basal', 'Delta')
  fixed$ci = fixed$upper - fixed$est
  fixed = fixed %>% 
    pivot_wider(
      names_from=Parameter, values_from=c(lower, 'est', 'upper', 'ci')
    )
  return(fixed)
}

# Create a function to map values outside the domain to the boundary colors
# Thanks chatgpt
custom_pal = function(palette, domain){
  function(vector) {
    palette_function = scales:::col_numeric(palette=palette, domain=domain)
    n_points = length(vector)
    colors = rep(0, n_points)
    for (i in 1:n_points){
      value = vector[i]
      if (value < domain[1]) {
        colors[i] = palette[1]
      } else if (value > domain[2]) {
        colors[i] = palette[length(palette)]
      } else {
        colors[i] = palette_function(value)
      }
    }
    return(colors)
  }
}

# Load in data
df = read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')

# Run mixed model on each region for relative PET
df_sum = df %>% 
  filter(Modality %in% c('fdg', 'ho', 'om', 'oc', 'oef', 'ogi')) %>%
  group_by(Modality, Region) %>% 
  nest() %>%
  summarise(coefs=fit_model(data)) %>%
  unnest(coefs)

# Make a nice table
tab = df_sum %>%
  dplyr:::select(starts_with("ci") | starts_with('est'), Region, Modality) %>%
  pivot_wider(
    id_cols=Region, names_from=Modality, values_from=-c(Modality, Region)
  ) %>%
  gt(rowname_col = 'Region') %>%
  cols_merge(
    columns=c('est_Basal_fdg', 'ci_Basal_fdg'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Basal_ho', 'ci_Basal_ho'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Basal_om', 'ci_Basal_om'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Basal_oc', 'ci_Basal_oc'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Basal_oef', 'ci_Basal_oef'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Basal_ogi', 'ci_Basal_ogi'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_fdg', 'ci_Delta_fdg'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_ho', 'ci_Delta_ho'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_om', 'ci_Delta_om'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_oc', 'ci_Delta_oc'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_oef', 'ci_Delta_oef'),
    pattern = "{1} ± {2}"
  ) %>%
  cols_merge(
    columns=c('est_Delta_ogi', 'ci_Delta_ogi'),
    pattern = "{1} ± {2}"
  ) %>%
  tab_stubhead(label = md('**Region**')) %>%
  fmt_number(decimals=2) %>%
  cols_label(
    est_Basal_fdg=md('**[<sup>18</sup>F]FDG**'),
    est_Basal_ho=md('**[<sup>15</sup>O]H<sub>2</sub>O**'),
    est_Basal_om=md('**[<sup>15</sup>O]O<sub>2</sub>**'),
    est_Basal_oc=md('**[<sup>15</sup>O]CO**'),
    est_Basal_oef=md('**rOEF**'),
    est_Basal_ogi=md('**rOGI**'),
    est_Delta_fdg=md('**[<sup>18</sup>F]FDG**'),
    est_Delta_ho=md('**[<sup>15</sup>O]H<sub>2</sub>O**'),
    est_Delta_om=md('**[<sup>15</sup>O]O<sub>2</sub>**'),
    est_Delta_oc=md('**[<sup>15</sup>O]CO**'),
    est_Delta_oef=md('**rOEF**'),
    est_Delta_ogi=md('**rOGI**')
  ) %>%
  tab_spanner(
    columns=c('est_Basal_fdg',
              'est_Basal_ho',
              'est_Basal_om',
              'est_Basal_oc',
              'est_Basal_oef',
              'est_Basal_ogi'),
    label=md('**Basal (Mean [95% CI])**'),
    id='Basal'
  ) %>%
  tab_spanner(
    columns=c('est_Delta_fdg',
              'est_Delta_ho',
              'est_Delta_om',
              'est_Delta_oc',
              'est_Delta_oef', 
              'est_Delta_ogi'),
    label=md('**Hygly. - Basal  (Mean [95% CI])**'),
    id='Delta'
  ) %>% 
  data_color(
    columns = starts_with("est_Delta"),
    method = "numeric",
    fn = custom_pal(rev(RColorBrewer::brewer.pal(7, 'RdBu')), c(-0.15, 0.15)),
  ) %>%
  data_color(
    columns = starts_with("est_Basal"),
    method = "numeric",
    fn =custom_pal(viridis::plasma(7), c(0.7, 1.4))
  ) %>%
  gtsave(
    str_glue('table_s2.docx')
  )




