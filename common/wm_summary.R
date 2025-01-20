#Load in libraries
library(tidyverse)
library(readxl)
library(nlme)
library(multcomp)
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

#Load in data file
if (args[1] == 'cmrglc'){
  mod = 'cmrglc'
} else if (args[1] == 'oxy'){
  mod = 'om'
} else if (args[1] == 'ho'){
  mod = 'ho'
} else if (args[1] == 'cbf'){
  mod = 'cbf'
}

all_data = read_excel('../common/data_values_file.xlsx', "Freesurfer.Rois")
mod_data = all_data %>% filter(Modality==mod & Region == 'Deep White Matter')

#Make factors
mod_data$Subject = as.factor(mod_data$ID)
mod_data$Group = factor(mod_data$Condition,
                         levels=c('basal', 'hypergly'),
                         labels=c('Eugly.', 'Hyper.')) 
mod_data$Rev = factor(mod_data$Condition, levels = c('hypergly', 'basal'))

#Run fit
fit = lme(Value~Group, random=~1|Subject, data=mod_data)
fit2 = lme(Value~Rev, random=~1|Subject, data=mod_data)
fit_sum = summary(fit)
print(fit_sum)

# Print percent change
print(fit_sum$tTable[2, 1] / fit_sum$tTable[1, 1] * 100.0)

print(intervals(fit, which='fixed'))
print(intervals(fit2, which='fixed'))

# Add visit time to data frame
demo_data = read_excel('../common/data_values_file.xlsx', "Demographics") %>% 
  dplyr:::select(c(ID, Visit.Order, Delta.Visit.1))
mod_data = mod_data %>% left_join(demo_data, by=c('ID', 'Visit.Order'))

# Fit with time from baseline
fit3 = lme(Value~Group + Delta.Visit.1, random=~1|Subject, data=mod_data)
print(summary(fit3))

# Get an equivalence test for 5%
eq_val = abs(fixef(fit)[1] * 0.05)
eq_less = summary(
  glht(fit, linfct = matrix(c(0, 1), nrow=1), rhs=eq_val, alternative='less')
)
eq_greater = summary(
  glht(fit, linfct = matrix(c(0, 1), nrow=1), rhs=-eq_val, alternative='greater')
)
print(eq_less)
print(eq_greater)
print(eq_less$test$pvalues[1] < 0.05 & eq_greater$test$pvalues[1] < 0.05)


