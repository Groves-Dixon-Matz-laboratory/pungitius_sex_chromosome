#plot_fig2.R
#plot the Y-consistent gene trees figure

library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/figure_plotting")
source('~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R')


#------------ plot gene tree figure ------------#

#PLOT SDR-LIKE BOOLEAN DATA
ll=load("pun_Ylike.Rdata");SEX.CHROM = 'chrXII' #set this up with chrom_window_figuresV4_wrapped.R
dat$mb = dat$BIN_START / 1e6
bm = sdr_normalized(df=dat, xcol='norm.x', ycol='sexCon', sex.chrom= SEX.CHROM, YLAB="Y-consistent", sex.col ='black', other.col ='grey',LWD=0.5,PT.SIZE = 2, SHAPE=3, SPAN=.2, draw.legend=T)
tot.length = max(dat[dat$CHROM=='chrXII','mb'])




#add breaktpoints
norm.breaks = c(3.5, 18.9) / tot.length
bm=bm + geom_vline(xintercept= norm.breaks[1], lty=2, lwd=0.5, color='green') + 
	geom_vline(xintercept= norm.breaks[2], lty=2, lwd=0.5,color='green')


#set up coordinates for example trees
mx=max(dat$mb[dat$CHROM=='chrXII'])
Y1=c(0,0,1,1)
Y2=c(0.05,0.05,0.95,0.95)
labY = c(0.15,0.15,0.85,0.85)
labX = c(0.8, 2, 6.6, 14.7)/mx
lefts = c(1.1, 2, 6.6, 14.7)/mx
rights = c(1.2, 2.1, 6.7, 14.8)/mx
tlabs = c("B", "C", "D", "E")
tree.segs = data.frame(x1=lefts, x2=rights, y1=Y2, y2=Y2, tlab=tlabs, labY=labY, labX=labX)


plot(bm)
















