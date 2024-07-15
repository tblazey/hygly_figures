#Load libraries
library(tidyverse)
library(readxl)
library(nlme)
library(cowplot)
library(multcomp)
rm(list=ls())

#Options for ggplot2
theme_set(theme_light())
gg_options = theme(axis.text.x=element_text(size=14),
                   axis.text.y=element_text(size=14),
                   axis.title.x=element_text(size=14,
                                             face="bold",
                                             margin = margin(t = 10,
                                                             r = 0,
                                                             b = 0,
                                                             l = 0)),
                   axis.title.y=element_text(size=16,
                                             angle=90,
                                             margin = margin(t = 0,
                                                             r = 10,
                                                             b = 0,
                                                             l = 0),
                                             face="bold"),
                   legend.text=element_text(size=12, face="bold"),
                   plot.title=element_text(size=16, face="bold", hjust=0.5),
                   legend.title=element_text(size=14, face="bold"),
                   legend.title.align=0.5, legend.position="None")

#Load in blood data
blood_long = read_excel('../common/data_values_file.xlsx', "Blood", na="NA")
blood_long$Group.Var = interaction(blood_long$ID, blood_long$Visit.Order)

#Make a regressor for piecewise regression
knot = 55
blood_long$Piece = pmax(0, blood_long$Time - knot)

# Get individual data frames
glu_long = blood_long %>% filter(Compound == 'Glucose') %>% drop_na(Value)
ins_long = blood_long %>% filter(Compound == 'Insulin') %>% drop_na(Value)

#Run piecewise fits
g_fit = lme(
  Value~Time+Piece+Condition+Time*Condition+Piece*Condition,
  random=~1|ID/Visit.Order,
  data=glu_long
)
i_fit = lme(Value~Time+Piece+Condition+Time*Condition+Piece*Condition,
            random=~1|ID/Visit.Order,
            data=ins_long
)
print(summary(g_fit))
print(summary(i_fit))
g_int = intervals(g_fit, which='fixed')
i_int = intervals(i_fit, which='fixed')
print(g_int)
print(i_int)

#Get baseline values
g_eu_base = glht(g_fit, linfct = matrix(c(1, 0, 0, 0, 0, 0), nrow=1))
g_hy_base = glht(g_fit, linfct = matrix(c(1, 0, 0, 1, 0, 0), nrow=1))
i_eu_base = glht(i_fit, linfct = matrix(c(1, 0, 0, 0, 0, 0), nrow=1))
i_hy_base = glht(i_fit, linfct = matrix(c(1, 0, 0, 1, 0, 0), nrow=1))

#Compute slopes before and after knot
g_eu_b = summary(glht(g_fit, linfct = matrix(c(0, 1, 0, 0, 0, 0), nrow=1)))
g_eu_a = summary(glht(g_fit, linfct = matrix(c(0, 1, 1, 0, 0, 0), nrow=1)))
g_hy_b = summary(glht(g_fit, linfct = matrix(c(0, 1, 0, 0, 1, 0), nrow=1))) 
g_hy_a = summary(glht(g_fit, linfct = matrix(c(0, 1, 1, 0, 1, 1), nrow=1)))
i_eu_b = summary(glht(i_fit, linfct = matrix(c(0, 1, 0, 0, 0, 0), nrow=1)))
i_eu_a = summary(glht(i_fit, linfct = matrix(c(0, 1, 1, 0, 0, 0), nrow=1)))
i_hy_b = summary(glht(i_fit, linfct = matrix(c(0, 1, 0, 0, 1, 0), nrow=1))) 
i_hy_a = summary(glht(i_fit, linfct = matrix(c(0, 1, 1, 0, 1, 1), nrow=1)))

#Show me the coefficients and normal CIs
print(paste(signif(g_eu_b$test$coefficients, 2), "±", signif(g_eu_b$test$sigma * 1.96, 2)))
print(paste(signif(g_eu_a$test$coefficients, 2), "±", signif(g_eu_a$test$sigma * 1.96, 2)))
print(paste(signif(g_hy_b$test$coefficients, 2), "±", signif(g_hy_b$test$sigma * 1.96, 2)))
print(paste(signif(g_hy_a$test$coefficients, 2), "±", signif(g_hy_a$test$sigma * 1.96, 2)))
print(paste(signif(i_eu_b$test$coefficients, 2), "±", signif(i_eu_b$test$sigma * 1.96, 2)))
print(paste(signif(i_eu_a$test$coefficients, 2), "±", signif(i_eu_a$test$sigma * 1.96, 2)))
print(paste(signif(i_hy_b$test$coefficients, 2), "±", signif(i_hy_b$test$sigma * 1.96, 2)))
print(paste(signif(i_hy_a$test$coefficients, 2), "±", signif(i_hy_a$test$sigma * 1.96, 2)))
             
#Make design matrix for predictions
t_pred = seq(from=0, to=300, length.out=100)
p_pred = pmax(0, t_pred - knot)
design = cbind(rep(1, 200),
                rep(t_pred, 2),
                rep(p_pred, 2),
                rep(c(0, 1), each=100),
                rep(c(0, 1), each=100) * t_pred,
                rep(c(0, 1), each=100) * p_pred)

#Get predictions
g_pred = data.frame(Time=design[, 2],
                    Value=design %*% fixef(g_fit),
                    Condition=rep(c('basal', 'hypergly'), each=100))
i_pred = data.frame(Time=design[, 2],
                    Value=design %*% fixef(i_fit),
                    Condition=rep(c('basal', 'hypergly'), each=100))

#Get limits
g_cov = vcov(g_fit)
i_cov = vcov(i_fit)
g_lower = rep(0, 200)
g_upper = rep(0, 200)
i_lower = rep(0, 200)
i_upper = rep(0, 200)
for (i in 1:200){
  
  #Compute standard error
  grad = c(0, design[i, 2:6])
  g_err = grad %*% g_cov %*% grad
  g_se = sqrt(g_err[1, 1])
  i_err = grad %*% i_cov %*% grad
  i_se = sqrt(i_err[1, 1])
  
  #Get upper and lower limits
  g_lower[i] = g_pred[i, 2] - (1.96 * g_se)
  g_upper[i] = g_pred[i, 2] + (1.96 * g_se)
  i_lower[i] = i_pred[i, 2] - (1.96 * i_se)
  i_upper[i] = i_pred[i, 2] + (1.96 * i_se)
  
}
g_pred$Lower = g_lower
g_pred$Upper = g_upper
i_pred$Lower = i_lower
i_pred$Upper = i_upper

#Get predictions at knot
g_eu_knot = c(1, knot, 0, 0, 0, 0) %*% fixef(g_fit)
g_hy_knot = c(1, knot, 0, 1, knot, 0) %*% fixef(g_fit)
g_knot_df = data.frame(Time=c(knot, knot),
                       Value=c(g_eu_knot, g_hy_knot),
                       Condition=c('basal', 'hypergly'))
i_eu_knot = c(1, knot, 0, 0, 0, 0) %*% fixef(i_fit)
i_hy_knot = c(1, knot, 0, 1, knot, 0) %*% fixef(i_fit)
i_knot_df = data.frame(Time=c(knot, knot),
                       Value=c(i_eu_knot, i_hy_knot),
                       Condition=c('basal', 'hypergly'))

g_eu_knot_se =  sqrt(
  c(0, knot, 0, 0, 0, 0) %*% g_cov %*% c(0, knot, 0, 0, 0, 0)
  )[1, 1]
g_hy_knot_se =  sqrt(
  c(0, knot, 0, 1, knot, 0) %*% g_cov %*% c(0, knot, 0, 1, knot, 0)
)[1, 1]
i_eu_knot_se =  sqrt(
  c(0, knot, 0, 0, 0, 0) %*% i_cov %*% c(0, knot, 0, 0, 0, 0)
)[1, 1]
i_hy_knot_se =  sqrt(
  c(0, knot, 0, 1, knot, 0) %*% i_cov %*% c(0, knot, 0, 1, knot, 0)
)[1, 1]

print(paste(signif(g_eu_knot, 4), "±", signif(g_eu_knot_se * 1.96, 4)))
print(paste(signif(g_hy_knot, 4), "±", signif(g_hy_knot_se * 1.96, 4)))
print(paste(signif(i_eu_knot, 4), "±", signif(i_eu_knot_se * 1.96, 4)))
print(paste(signif(i_hy_knot, 4), "±", signif(i_hy_knot_se * 1.96, 4)))

#Make a glucose plot
g_plot = ggplot(glu_long, aes(x=Time, y=Value)) 
g_point = geom_point(aes(group=Group.Var, color=Condition), alpha=0.3)
g_line = geom_line(aes(group=Group.Var, color=Condition), alpha=0.3)
g_color = scale_color_manual('Cond.',
                             labels=c('Eugly.', 'Hyper.'),
                             values=c('#006bb6', '#b6006b'))
g_x = scale_x_continuous('Time (Minutes)')
g_y = scale_y_continuous('Concentration (mg/dL)')
g_title = ggtitle('Glucose')
g_hat_line = geom_line(data=g_pred, aes(color=Condition), linewidth=1.75)
g_ribbon = geom_ribbon(data=g_pred,
                       aes(ymin=Lower, ymax=Upper, group=Condition),
                       alpha=0.4,
                       fill='gray60')
g_knot = geom_point(data=g_knot_df, color='black', size=3.25)
g_fig = g_plot + g_line + g_point + g_color + g_x +
        g_y + g_title + g_ribbon + g_hat_line + g_knot + gg_options

#Make a insulin plot
i_plot = ggplot(ins_long, aes(x=Time, y=Value)) 
i_point = geom_point(aes(group=Group.Var, color=Condition), alpha=0.3)
i_line = geom_line(aes(group=Group.Var, color=Condition), alpha=0.3)
i_color = scale_color_manual('Cond.',
                             labels=c('Eugly.', 'Hyper.'),
                             values=c('#006bb6', '#b6006b'))
i_x = scale_x_continuous('Time (Minutes)')
i_y = scale_y_continuous('Conc. (pmol/L)', limits=c(0, 500))
i_title = ggtitle('Insulin')
i_hat_line = geom_line(data=i_pred, aes(color=Condition), linewidth=1.75)
i_ribbon = geom_ribbon(data=i_pred,
                       aes(ymin=Lower, ymax=Upper, group=Condition),
                       alpha=0.4,
                       fill='gray60')
i_knot = geom_point(data=i_knot_df, color='black', size=3.25)
i_ref = geom_hline(yintercept=425, linetype="dashed", color = "black", linewidth=.9)
i_fast = geom_hline(yintercept=140, linetype="dashed", color = "black", linewidth=.9)
i_text = geom_text(x=165, y=450, label='Glucose Load', size=4.5)
i_text_two = geom_text(x=30, y=166, label='Fasting', size=4.5)
i_fig = i_plot + i_fast + i_line + i_point + i_color + i_x +
        i_y + i_title + i_ribbon + i_hat_line + 
        i_knot + i_ref + i_text + i_text_two + gg_options

#Join up all the plots
c_fig = plot_grid(NULL, g_fig, NULL, i_fig,
                  labels=c('', 'A)', '',  'B)'),
                  label_size=16, nrow=1, ncol=4,
                  rel_widths=c(0.015, 1, 0.025, 1))

#Add legend
c_legend = get_legend(g_fig + theme(legend.position="right"))
c_fig = plot_grid(c_fig, c_legend, rel_widths=c(3, 0.3))

#Save figure
ggsave('figure_1.tiff', plot=c_fig, width=12, height=5, units='in')

