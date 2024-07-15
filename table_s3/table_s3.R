#Load libraries
library(tidyverse)
library(readxl)
library(nlme)
rm(list=ls())

#Load in data
wb_data = read_excel('../common/data_values_file.xlsx', "Whole.Brain.FDG")

#Get values for table
pars = unique(wb_data$Parameter)
n_par = length(pars)
wb_table = array(0, dim=c(n_par, 7))
rownames(wb_table) = pars
colnames(wb_table) = c(
  'Eugly', 'Eugly.CI', 'Hygly', 'Hygly.CI', 'Diff', 'Diff.CI', 'Diff.p.value'
)

for (i in 1:n_par){

    #Make frame for fit
    df = wb_data %>% filter(Parameter == pars[i])

    #Run fit
    wb_fit = lme(Value~Condition, random=~1|ID, data=df)
    wb_sum = summary(wb_fit)

    #Extract means and CI for baseline
    wb_table[i, 1] = wb_sum$tTable[1, 1]
    wb_table[i, 3] = wb_sum$tTable[1, 1] + wb_sum$tTable[2, 1]
    wb_table[i, 2] = wb_sum$tTable[1, 2] * 1.96

    #Compute 95% confidence interval for hypeglycemia
    wb_table[i, 4] = sqrt(as.numeric(c(1,1) %*% vcov(wb_fit) %*% c(1,1))) * 1.96

    #Get mean difference and CI
    wb_table[i, 5] = wb_sum$tTable[2, 1]
    wb_table[i, 6] = wb_sum$tTable[2, 2] * 1.96

    #Save p-value
    wb_table[i, 7] = wb_sum$tTable[2, 5]

}
write.csv(wb_table, 'table_s3.csv', quote=FALSE)
print(wb_table)



