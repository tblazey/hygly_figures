#Load in libraries
library(readxl)
library(nlme)
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

#Load in data file
sess_info = read_excel('../common/data_values_file.xlsx', "Ref.CMRglc")

#Make factors
sess_info$Subject = as.factor(sess_info$ID)
sess_info$Group = factor(sess_info$Condition,
                         levels=c('basal', 'hypergly'),
                         labels=c('Eugly.', 'Hyper.')) 

#Run fit
fit = lme(Value~Group, random=~1|Subject, data=sess_info)
fit_sum = summary(fit)
print(fit_sum)
print(intervals(fit_sum, which='fixed'))