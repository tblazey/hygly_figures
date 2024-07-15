library(tidyverse)
library(readxl)
library(scales)
library(grid)
library(gridExtra)
rm(list=ls())

#Options for ggplot2
theme_set(theme_light())
gg_options = theme(axis.text.x=element_text(size=18),
                   axis.text.y=element_text(size=18),
                   axis.title.x=element_text(size=18, face="bold",
                                             margin = margin(t = 10, r = 0, b = 0, l = 0)),
                   axis.title.y=element_text(size=18, angle=90, 
                                             margin = margin(t = 0, r = 10, b = 0, l = 0), 
                                             face="bold"),
                   legend.text=element_text(size=16, face="bold"),
                   plot.title=element_text(size=24, face="bold", hjust=0.5),
                   legend.title=element_text(size=16, face="bold",
                                             margin = margin(t = 0, r = 10, b = 10, l = 0)),
                   legend.title.align=0.5,
                   plot.tag=element_text(size=24, face="bold"),
                   plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"))

# Constants and such
eu_val = 5.5
hy_val = 16.5

#Load in matrix data
mat = read_excel('../common/data_values_file.xlsx', "Sim.Matrix")

#Make plot
mat_plot = ggplot(mat, aes(x=Si_scale, y=Km_scale, fill=Per.Change))
mat_tile = geom_tile(color='black', linewidth=0.75)
mat_x = scale_x_continuous(trans = 'log10',
                           name = "Fraction of GM Intracellular Glucose",
                           breaks = trans_breaks('log10', function(x) 10^x),
                           labels = trans_format('log10', math_format(10^.x))) 
mat_y = scale_y_continuous(trans = 'log10',
                           name = bquote(bold(Fraction~of~HK1~K[m])),
                           breaks = trans_breaks('log10', function(x) 10^x),
                           labels = trans_format('log10', math_format(10^.x))) 
mat_fill = scale_fill_viridis_c(name='% Change', option='C', 
                                breaks=c(0, 20, 40, 60), limits=c(0, 60))
mat_annot = annotate("rect", xmin=0.077, xmax=4.13, ymin=4.35, ymax=7.25, 
                     colour="white", fill="transparent", linewidth=2)
mat_tag = labs(tag='A)') 
mat_fig = mat_plot + mat_tile + mat_x + mat_y + mat_fill + mat_annot + mat_tag + gg_options

#Load in curve data
curve_data = read_excel('../common/data_values_file.xlsx', "Sim.Curves")
wm_data = curve_data %>% filter(Tissue == 'White.Matter')
gm_data = curve_data %>% filter(Tissue == 'Gray.Matter')

#Normalize Vhex and change blood glucose units to mM
wm_data$Vhex = (wm_data$Vhex - 3.67) / 3.67 * 100
wm_data$So = wm_data$So * 1E3
gm_data$Vhex = (gm_data$Vhex - 6.5) / 6.5 * 100
gm_data$So = gm_data$So * 1E3

#Make plot at WM 5.6
wm_plot = ggplot(wm_data, aes(x=So, y=Vhex))
wm_line = geom_line(aes(group=Si, color=round(log10(Si), 2)), linewidth=1.25)
gm_line = geom_line(data=gm_data,  aes(y=Vhex), 
                    color='gray60', linewidth=2.25)
wm_color = scale_color_viridis_b(name=bquote(bold(10^{GM~S[i]~Fraction})), option='C',
                                 direction=-1,
                                 breaks=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5),
                                 limits=c(-1, 0.5),
                                 right=FALSE)
wm_x = scale_x_continuous(name = "Blood Glucose (mM)")
wm_y_label = expression(bold(paste(Delta, ' HK Flux (%)')))
wm_y = scale_y_continuous(name = wm_y_label)
wm_eu = geom_segment(aes(x=eu_val, xend=eu_val, y=-100, 
                             yend=70), 
                         linetype="dashed", color = "black", linewidth=.9)
wm_hy = geom_segment(aes(x=hy_val, xend=hy_val, y=-100, 
                         yend=70), 
                     linetype="dashed", color = "black", linewidth=.9)
wm_tag = labs(tag='B)')                    
wm_fig = wm_plot + wm_line + wm_color + wm_x + wm_y + gm_line + wm_eu + wm_hy + wm_tag +
         guides(color = guide_coloursteps()) + gg_options +
         theme(legend.text=element_text(size=12, face="bold", hjust=0))

# Combine plot and add arrows
tiff("figure_7.tiff", width = 18, height = 6, units='in', res=300)
grid_plot = grid.arrange(mat_fig, wm_fig, ncol=2)
grid.lines(x=unit(c(0.425, 0.575), "npc"), y=unit(c(0.725, 0.725), "npc"),
           gp=gpar(fill="black", lwd=4), 
           arrow=arrow(type="closed", length=unit(4,"mm")))
dev.off()

