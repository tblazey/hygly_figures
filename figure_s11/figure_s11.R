#Load libraries
library(tidyverse)
library(readxl)
library(nlme)
library(cowplot)
rm(list=ls())

#Options for ggplot2
theme_set(theme_light())
gg_options = theme(axis.text.x=element_text(size=16),
                   axis.text.y=element_text(size=16),
                   axis.title.x=element_text(size=20, face="bold",
                                             margin = margin(t = 20,
                                                             r = 0,
                                                             b = 0,
                                                             l = 0)),
                   axis.title.y=element_text(size=20, angle=90, 
                                             margin = margin(t = 0,
                                                             r = 20,
                                                             b = 0,
                                                             l = 0), 
                                             face="bold"),
                   legend.text=element_text(size=20, face="bold"),
                   plot.title=element_text(size=22, face="bold", hjust=0.5,
                                           margin = margin(t = 0,
                                                           r = 0,
                                                           b = 10,
                                                           l = 0)),
                   legend.title=element_text(size=20, face="bold"),
                   legend.title.align=0.5,
                   plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"),
                   legend.position="none",
                   strip.text.x = element_text(size = 16))

# Read in data
behav_df = read_excel('../common/data_values_file.xlsx', 'Task.Behav')

#Run model fits with individual trials
enc_rt_fit = lme(
  RT~Condition + Resp + Condition*Resp, random=~Resp|ID/Visit.Order/Run,
  data=behav_df %>% filter(Task == 'Encoding')
)
ret_rt_fit = lme(
  RT~Condition + Resp + Condition*Resp, random=~Resp|ID/Visit.Order/Run,
  data=behav_df %>% filter(Task == 'Retrieval')
)
print(summary(enc_rt_fit)$tTable)
print(summary(ret_rt_fit)$tTable)

#Get group level fixed effects for plots
enc_fix = fixef(enc_rt_fit)
enc_fit_avg = data.frame(RT=c(enc_fix[1] + enc_fix[3],
                              enc_fix[1],
                              sum(enc_fix),
                              enc_fix[1] + enc_fix[2]),
                         Resp=rep(c("Fits", "Doesn't Fit"), 2),
                         Condition=rep(c("basal", "hypergly"), each=2))
ret_fix = fixef(ret_rt_fit)
ret_fit_avg = data.frame(RT=c(ret_fix[1] + ret_fix[3],
                              ret_fix[1],
                              sum(ret_fix),
                              ret_fix[1] + ret_fix[2]),
                         Resp=rep(c("Remembered", "Forgotten"), 2),
                         Condition=rep(c("basal", "hypergly"), each=2))

#Compute individual means for plots
behav_means = behav_df %>% 
  group_by(ID, Condition, Resp, Task) %>%
  summarise(RT=mean(RT))

#Run model fits with average data
enc_rt_means = behav_means %>% filter(Task == 'Encoding')
enc_rt_avg_fit = lme(
  RT~Condition + Resp + Condition*Resp, random=~1|ID, data=enc_rt_means
)
ret_rt_means = behav_means %>% filter(Task == 'Retrieval')
ret_rt_avg_fit = lme(
  RT~Condition + Resp + Condition*Resp, random=~1|ID, data=ret_rt_means
)
print(summary(enc_rt_avg_fit)$tTable)
print(summary(ret_rt_avg_fit)$tTable)

# Make sure there isn't an overall condition effect when combining all times
rt_avg_fit = lme(
  RT~Condition + Task + Condition*Task, random=~1|ID, data=behav_means
)
print(summary(rt_avg_fit)$tTable)

#Make encoding rt plot
enc_rt_plot = ggplot(enc_rt_means, aes(x=Condition, y=RT, color=Resp))
enc_rt_line = geom_line(aes(group=ID), alpha=0.2) 
enc_rt_point = geom_point(alpha=0.2)
enc_rt_y = scale_y_continuous('Reaction Time (ms)', limits=c(600, 2500))
enc_rt_xlabel = scale_x_discrete('Condition', labels= c('Eugly.', 'Hygyl.'))
enc_rt_title = ggtitle('Encoding')
enc_rt_color = scale_color_manual("", values=c('#006bb6', '#b6006b'), )
enc_rt_facet = facet_grid(. ~ Resp)
enc_rt_sum = geom_line(
  data=enc_fit_avg, aes(x=Condition, y=RT, group=1), linewidth=3
)
enc_rt_fig = enc_rt_plot + enc_rt_line + enc_rt_point + enc_rt_y +
             enc_rt_xlabel + enc_rt_title + enc_rt_facet + enc_rt_sum +
             enc_rt_color + gg_options

#Make retrieval rt plot
ret_rt_plot = ggplot(ret_rt_means, aes(x=Condition, y=RT, color=Resp))
ret_rt_line = geom_line(aes(group=ID), alpha=0.2) 
ret_rt_point = geom_point(alpha=0.2)
ret_rt_y = scale_y_continuous('', limits=c(600, 2500))
ret_rt_xlabel = scale_x_discrete('Condition', labels= c('Eugly.', 'Hygyl.'))
ret_rt_title = ggtitle('Retrieval')
ret_rt_color = scale_color_manual("", values=c('#006bb6', '#b6006b'), )
ret_rt_facet = facet_grid(. ~ Resp)
ret_rt_sum = geom_line(
  data=ret_fit_avg, aes(x=Condition, y=RT, group=1), linewidth=3
)
ret_rt_fig = ret_rt_plot + ret_rt_line + ret_rt_point + ret_rt_y +
  ret_rt_xlabel + ret_rt_title + ret_rt_facet + ret_rt_sum +
  ret_rt_color + gg_options

#Join up plots
empty_y_theme = theme(axis.text.y=element_blank(),
                     axis.title.y = element_blank())
grid_plot = plot_grid(
  enc_rt_fig,
  ret_rt_fig + empty_y_theme,
  rel_widths=c(1.1, 1),
  labels=c('A)', 'B)'),
  label_size=20)
ggsave('figure_s11.tiff', plot=grid_plot,width=12, height=6)
