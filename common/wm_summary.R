#Load in libraries
library(tidyverse)
library(readxl)
library(nlme)
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

#Load in data file
if (args[1] == 'cmrglc'){
  fig = '2'
  mod = 'cmrglc'
} else if (args[1] == 'oxy'){
  fig = '3'
  mod = 'om'
} else if (args[1] == 'ho'){
  fig = '5'
  mod = 'ho'
} else if (args[1] == 'cbf'){
  fig = 's6'
  mod = 'cbf'
}

all_data = read_excel('../common/data_values_file.xlsx', "Freesurfer.Rois")
mod_data = all_data %>% filter(Modality==mod & Region == 'Deep White Matter')

#Make factors
mod_data$Subject = as.factor(mod_data$ID)
mod_data$Group = factor(mod_data$Condition,
                         levels=c('basal', 'hypergly'),
                         labels=c('Eugly.', 'Hyper.')) 

#Run fit
fit = lme(Value~Group, random=~1|Subject, data=mod_data)
fit_sum = summary(fit)
print(fit_sum)

# Print coefficients and standard errors
print(fit_sum$tTable[1, 1])
print(fit_sum$tTable[1, 2] * 1.96)
print(fit_sum$tTable[1, 1] + fit_sum$tTable[2, 1])
print(sqrt(as.numeric(c(1,1) %*% vcov(fit) %*% c(1,1))) * 1.96)
print(fit_sum$tTable[2, 1])
print(fit_sum$tTable[2, 2] * 1.96)
print(fit_sum$tTable[2, 5])

# Print percent change
print(fit_sum$tTable[2, 1] / fit_sum$tTable[1, 1] * 100.0)
