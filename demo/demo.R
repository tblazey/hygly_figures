library(tidyverse)
library(readxl)
rm(list=ls())

# Load in demographic data
data = read_excel('../common/data_values_file.xlsx', 'Demographics', na='NA')

# Separate out first and followup visits
base_data = data[data$Visit.Order == 1, ]
follow_data = data[data$Visit.Order != 1, ]

mean_age = round(mean(base_data$Age.Visit.1), 1)
sd_age = round(sd(base_data$Age.Visit.1), 1)
print(paste('Age:', mean_age, '+/-', sd_age))

mean_bmi = round(mean(base_data$BMI.Visit.1, na.rm = TRUE), 1)
sd_bmi = round(sd(base_data$BMI.Visit.1, na.rm = TRUE), 1)
print(paste('BMI:', mean_bmi, '+/-', sd_bmi))

median_delta = median(follow_data$Delta.Visit.1)
min_delta = min(follow_data$Delta.Visit.1)
max_delta = max(follow_data$Delta.Visit.1)
print('Time between visits:')
print(paste('Median =', median_delta))
print(paste('Min =', min_delta))
print(paste('Max =', max_delta))

n_both =sum(
    (data %>% group_by(ID) %>% summarize(u_cond=length(unique(Condition)) > 1))$u_cond
)
print(paste('Number of subjects with both conditions:', n_both))

n_first = data %>% filter(Visit.Order == 1, Visit.Count > 1) %>% count(Condition)
print(n_first)
