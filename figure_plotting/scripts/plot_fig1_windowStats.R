#plot_fig1_windowStats.R
#plot the 100 Kb window stats for each species

library(plotrix)
library(ggplot2)
library(gridExtra)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/figure_plotting")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")

#set some global variables
sex.col = 'black'
other.col = 'grey'
XCOL = 'norm.x'
SPAN=0.1
CUT= 0.99
XLAB=bquote(Pos.~.(sex.chrom.name)~(Mb))

#optionally don't set y limits
fst.ylim = "none"
r2.ylim = "none"
ss.ylim = "none"
dp.ylim = "none"

#assign y limits
fst.ylim = c(-.02, 0.3)
r2.ylim = c(0, 0.8)
ss.ylim = c(0, 85)
pi.ylim = c(-0.5, 2.5)
dp.ylim = c(-1,0.6)


#set up 3-spine panels
ll=load("po_assembledWindowStats.Rdata");SEX.CHROM = 'chrXIX'; sex.chrom.name="chr19";sex.specific = 'pctMaleSpecific'; specific.lab = '% Sex-Specific'
# img1=ggdraw() + draw_image("datasets/pacific_ocean.pdf", scale = 0.9)
fst1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, nox=T) #YLAB=bquote(Mean~F[ST])
r21 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, nox=T)
ss1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB= '', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, nox=T)
H = get_two_tailed(dat[,'bedtoolDepth'], CUT)
pi1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, nox=F)#;plot(pi1)
dp1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='bedtoolDepth', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= dp.ylim, verticals=F);#plot(dp1)

#set up pun panels
ll=load("pun_assembledWindowStats.Rdata");SEX.CHROM = 'chrXII';sex.chrom.name='chr12'
fst2 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, verticals=TRUE, nox=TRUE, noy=TRUE)#;plot(fst2)
r22 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, verticals=TRUE, nox=TRUE, noy=TRUE)
ss2 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, verticals=TRUE, nox=TRUE, noy=TRUE)
pi2 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, verticals=TRUE, nox=F, noy=TRUE)#;plot(pi2)
H = get_two_tailed(dat[,'bedtoolDepth'], CUT)
dp2 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='bedtoolDepth', sex.chrom= SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, draw.line=F, tailed = 2, horiz=H, YLIM=dp.ylim, verticals=TRUE, noy=F)#;plot(dp2)

#set up sin panels
ll=load("sin_assembledWindowStats.Rdata")#;sex.specific = 'pctFemaleSpecific';specific.lab = '% Sex-Specific'
fst3 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, nox=TRUE, noy=TRUE)
r23 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, nox=TRUE, noy=TRUE)
ss3 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, nox=TRUE, noy=TRUE);#plot(ss3)
pi3 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, nox=F, noy=TRUE)#;plot(pi3)
dp3 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='bedtoolDepth', sex.chrom= SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=dp.ylim, noy=T)

#set up tym panels
ll=load("tym_assembledWindowStats.Rdata")
fst4 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, nox=TRUE, noy=TRUE)
r24 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, nox=TRUE, noy=TRUE)
ss4 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, nox=TRUE, noy=TRUE)
pi4 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, nox=F, noy=TRUE)#;plot(pi4)
dp4 = multichrom_len_normalized(df=dat, xcol=XCOL, ycol='bedtoolDepth', sex.chrom=SEX.CHROM, YLAB='', XLAB=XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=dp.ylim, noy=TRUE)


#Fst, sex-specific, and pi
plot_grid(fst1,fst2,fst3,fst4,ss1,ss2,ss3,ss4,pi1,pi2,pi3,pi4, ncol=4, labels = '')


#plot depth for 3-spine and pungitius
plot_grid(dp1, dp2)
save(dp1, dp2, file="depthPlots.Rdata")



#other combinations:

#Fst, r2, and depth
plot_grid(fst1,fst2,fst3,fst4,r21,r22,r23,r24,dp1,dp2,dp3,dp4, ncol=4, labels = LETTERS[1:12])