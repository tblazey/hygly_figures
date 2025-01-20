#Load libs
library(tidyverse)
library(readxl)
library(nlme)
library(cowplot)
library(ggbrain)
library(multcomp)
rm(list=ls())

# Which figure are we doing 
args = commandArgs(trailingOnly=TRUE)
dir = args[1]
if (dir == 'pos'){
  figure = 'figure_s9'
  dir_label = 'Task Activated Regions'
} else if (dir == 'neg'){
  figure = 'figure_s10'
  dir_label = 'Task Deactivated Regions'
}

#Options for ggplot2
theme_set(theme_light())
gg_options = theme(axis.text.x=element_text(size=16),
                   axis.text.y=element_text(size=16),
                   axis.title.x=element_text(size=16, face="bold",
                                             margin = margin(t = 10, r = 0, 
                                                             b = 0, l = 0)),
                   axis.title.y=element_text(size=16, angle=90, 
                                             margin = margin(t = 0, r = 20, 
                                                             b = 0, l = 0), 
                                             face="bold"),
                   legend.text=element_text(size=20, face="bold"),
                   plot.title=element_text(size=24, face="bold", hjust=0.5),
                   legend.title=element_text(size=20, face="bold"),
                   legend.title.align=0.5)

#Define coordinates
pos_coords = c(c("z=0"), c("z=0"), c("z=27"), c("x=-6"))
neg_coords = c(c("x=-1.54"), c("z=4.3"), c("z=22"), c("z=7.8"))
coords = list(pos=pos_coords, neg=neg_coords)
roi_names = list(pos=c('Occipital',
                       'Ventral.Prefrontal',
                       'Lateral.Prefrontal',
                       'Frontal.Medial'),
                 neg=c('Parietal.Medial',
                       'Medial.Prefrontal',
                       'Parietal.Occipital',
                       'Auditory'))
roi_names = roi_names[[dir]]

#Main data directory
brain_path = file.path(
    '/opt/fsl', '/data/standard/MNI152_T1_2mm_brain.nii.gz'
)

#Load in data and get average runs
roi_data = read_excel('../common/data_values_file.xlsx', 'Task.Rois')
motion_sum = read_excel('../common/data_values_file.xlsx', 'Task.Motion') %>%
    filter(Metric=='MSE.FILT') %>% group_by(ID, Condition, Visit.Order, Task) %>%
    summarise(Value=mean(Value))

#Make "colormaps" for rois
colors = c("#b6006b", '#006bb6', '#6bb600', '#b6600b')
red_grad = scale_fill_gradient(low="#b6006b", high="#b6006b")
blue_grad = scale_fill_gradient(low='#006bb6', high='#006bb6')
green_grad = scale_fill_gradient(low='#6bb600', high='#6bb600')
brown_grad = scale_fill_gradient(low='#b6600b', high='#b6600b')
grads = c(red_grad, blue_grad, green_grad, brown_grad)

#Loop thorugh rois
plots = list()
idx = 1
for (roi in 1:4){
  
  #Loop through positive/negative response rois
  roi_path = str_glue('../common/rois/{dir}_roi_{roi}.nii.gz')
  
  #Define ggbrain stuff
  roi_obj = ggbrain(bg_color = "white", text_color = "black") + 
            geom_brain(definition="underlay") +
            geom_brain(definition="roi", limits=c(0, 1), show_legend=FALSE,
                       fill_scale=grads[[roi]]) + 
            images(c(underlay=brain_path, roi=roi_path)) + 
            slices(coords[[dir]][roi])
  roi_geom = roi_obj$render()
  plots[[idx]] = roi_geom + theme(plot.margin = margin(t=0, r=0, b=0, l=2, unit="cm"))
  plots[[idx + 1]] = NULL
  idx = idx + 2
  
  #Loop through tasks
  for (task in c('Encoding', 'Retrieval')){
    
    #Make design, with means for each condition
    motion_design = motion_sum %>%
      filter(Task == task) %>%
      group_by(ID, Condition) %>%
      summarise(Value=mean(Value))
    roi_design = roi_data %>%
      filter(Task == task & Region == roi_names[roi]) %>%
      group_by(ID, Condition) %>%
      summarise(Value=mean(Value))

    #Run model fit
    fit = lme(Value~Condition, random=~1|ID, data=roi_design)
    fit_p = summary(fit)$tTable[2,5]
    
    eq_val = abs(fixef(fit)[1] * 0.05)
    eq_less = summary(
      glht(fit, linfct = matrix(c(0, 1), nrow=1), rhs=eq_val, alternative='less')
    )
    eq_greater = summary(
      glht(fit, linfct = matrix(c(0, 1), nrow=1), rhs=-eq_val, alternative='greater')
    )
    
    print(eq_less)
    print(eq_greater)
    print(max(eq_less$test$pvalues[1], eq_greater$test$pvalues[1]))
    print(eq_less$test$pvalues[1] < 0.05 & eq_greater$test$pvalues[1] < 0.05)
    
    #Check for motion differences
    if(roi == 1){
      fit_abs = lme(Value~Condition, random=~1|ID, data=motion_design)
      print(summary(fit_abs))
      print(table(motion_design$Condition))
    }
    
    #Make roi plot
    roi_plot = ggplot(roi_design, aes(x=Condition, y=Value, group=ID))
    roi_point = geom_point(color=colors[[roi]])
    roi_line = geom_line(color=colors[[roi]])
    roi_x = scale_x_discrete('', labels=c('Eugly.', 'Hygly.'))
    roi_min = min(roi_design$Value)
    roi_max = max(roi_design$Value)
    if (roi_min > 0){
      roi_low = roi_min * 0.70
    } else{
      roi_low = roi_min * 1.20
    }
    if (dir == 'pos'){
      roi_high = roi_max * 1.30
      annot_y = roi_max * 1.2
    } else{
      roi_high = roi_max * 2.25
      annot_y = roi_max * 2
    }
    roi_limits = c(roi_low, roi_high)
    if (roi == 1){
      roi_y = scale_y_continuous("Reg. Coef.",
                                 limits=roi_limits)
    }
    else{
      roi_y = scale_y_continuous("", limits=roi_limits)
    }
    roi_annot = annotate("text",
                         x = 1.5,
                         y = annot_y,
                         label = deparse(bquote(italic(p)~"="~.(round(fit_p, 2)))), 
                         size=5, parse=TRUE)
    roi_fig = roi_plot + roi_point + roi_line + roi_x + roi_y +
              roi_annot + gg_options
    plots[[idx]] = roi_fig
    idx = idx + 1
    
  }
  
}

#Make base grid
grid = plot_grid(plotlist=plots, ncol=4, byrow=FALSE, 
                 rel_heights=c(0.8, 0.1, 1, 1),
                 labels=c("A)", "B)", "C)", "D)"),
                 label_size=22, label_x=0.05)

#Add task labels
grid = plot_grid(grid, NULL, rel_widths=c(1, 0.05))
enc_label = draw_label("Encoding", size = 20, angle = 90, x = 0.975, y = 0.56,
                       fontface = 'bold')
ret_label = draw_label("Retrieval", size = 20, angle = 90, x = 0.975, y = 0.22,
                       fontface = 'bold')
grid = grid + enc_label + ret_label

#Add title
grid = plot_grid(NULL, grid, rel_heights=c(0.125, 1), nrow=2)
fig_label = draw_label(dir_label, size=28, x=0.5, y=0.975,
                       fontface='bold')
grid = grid + fig_label

#Save figure
ggsave(str_glue("{figure}.tiff"), units="cm", height=20, width=35, plot=grid, bg='white')

