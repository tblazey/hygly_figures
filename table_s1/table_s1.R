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

    # Make frame for fit. Get a reversed factor so we can get confidence 
    # intervals for hyperglycemia factor. A total hack :(
    df = wb_data %>% filter(Parameter == pars[i])
    df$Rev = factor(df$Condition, levels = c('hypergly', 'basal'))
    
    #Run fit
    wb_fit = lme(Value~Condition, random=~1|ID, data=df)
    wb_fit2 = lme(Value~Rev, random=~1|ID, data=df)
    wb_sum = summary(wb_fit)
    
    # Extract intervals
    wb_int = intervals(wb_fit, which='fixed')
    wb_int2 = intervals(wb_fit2, which='fixed')

    # Extract means and CI for baseline
    wb_table[i, 1] = wb_int$fixed[1, 2]
    wb_table[i, 2] = wb_int$fixed[1, 2] - wb_int$fixed[1, 1]
    
    # Means and CI for hyperglycemia level
    wb_table[i, 3] = wb_int2$fixed[1, 2]
    wb_table[i, 4] = wb_int2$fixed[1, 2] - wb_int2$fixed[1, 1]

    #Get mean difference and CI
    wb_table[i, 5] = wb_int$fixed[2, 2]
    wb_table[i, 6] = wb_int$fixed[2, 2] - wb_int$fixed[2, 1]

    #Save p-value
    wb_table[i, 7] = wb_sum$tTable[2, 5]

}
write.csv(wb_table, 'table_s1.csv', quote=FALSE)
print(wb_table)
print(p.adjust(wb_table[,7], method='fdr'))


