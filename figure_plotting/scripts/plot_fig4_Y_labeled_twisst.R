#plot_fig4_Y_labeled_twisst.R
#plot twisst results from gene trees with Y labeled haplotypes
#input files are generated using y_assignment_trees_walkthough.txt


library(ggplot2)
library(cowplot)
library(scales)
source("~/gitreps/pungitius_sex_chromosome/gene_trees/scripts/twisst_plotting_functions.R")
setwd('~/gitreps/pungitius_sex_chromosome/figure_plotting/')
ll=load("Y_topo_weighting.Rdata")
plot(lp)


#plot simple line plot
dim(weights)
line.cols = hue_pal()(3) #switch to cowplot colors
LWD=2
weights_smooth$mb = window_data$mid / 1e6
norm.breaks = c(3.5, 18.9)
lp = ggplot(data= weights_smooth) + 
	geom_vline(xintercept= norm.breaks[1], lty=2, lwd=0.8, color='grey') + 
	geom_vline(xintercept= norm.breaks[2], lty=2, lwd=0.8,color='grey') +
	geom_line(aes(x=mb,y=topo3), col=line.cols[2], lwd=LWD) +
	geom_line(aes(x=mb,y=topo2), col=line.cols[3], lwd=LWD) +
	geom_line(aes(x=mb,y=topo1), col=line.cols[1], lwd=LWD) + 
	labs(y="Weight", x="Position (Mb)") + 
	scale_y_continuous(breaks=c(0,0.5, 1)) +
	scale_x_continuous(breaks=c(0,10,20))
plot(lp)


