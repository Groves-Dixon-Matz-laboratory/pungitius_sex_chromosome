#plot_figure6_degeneration.R

library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/figure_plotting")

#set up global plotting variables
male.col = 'dodgerblue'
female.col = 'firebrick'
sin.col = 'forestgreen'
tym.col = 'mediumseagreen'
xlabAngle=20

#load depth comparisons
ll=load('depthPlots.Rdata')
ll
dp1 = dp1 + labs(x="Position (Mb)", y='M:F coverage', subtitle='G. aculeatus') + theme(plot.subtitle=element_text(face="italic")) 
dp2 = dp2 + labs(x="Position (Mb)", y='M:F coverage', subtitle='P. pungitius') + theme(plot.subtitle=element_text(face="italic")) 
plot_grid(dp1, dp2, ncol=2)

#-------- REPETITIVE ELEMENT ENRICHMENT --------#
#LOAD DATA
ll=load("repetitiveElements.Rdata")
#rep.res = the DESeq2 results for repetitive elements
#s = the set of repetitive elements with sinificantly different fold coverage between males and females
#mt = the estimated proportion of reads mapped to repetitive elements
ll

#BUILD FIGURES

#density plot
den=density(s$log2FoldChange)
den.df = data.frame(x=den$x, y=den$y)
male.df = den.df[den.df$x >=max(den.df$x[den.df$x<0]),]
female.df = den.df[den.df$x <= 0,]
denRep=ggplot(data=female.df, aes(x=x, y=y)) + 
	geom_vline(aes(xintercept=0), linetype='dashed') +
	geom_line(col=female.col) + 
	geom_area(fill=female.col) +
	geom_line(data=male.df, aes(x=x,y=y), col=male.col) + 
	geom_area(data=male.df, aes(x=x,y=y), fill=male.col) +
	lims(x=c(-3.3, 3.3)) +
	# scale_x_continuous(limits = c(-0.03, 0.22)) +
	# geom_segment(data=d12, aes(x=d, y=y1, xend=d, yend=y2), lineend='round', lwd=1.2) +
	labs(x="", y="")
plot(denRep)


#stats
rep.res$sig = rep.res$padj < 0.1
m = rep.res[rep.res$log2FoldChange>0,]
osig = sum(m$sig, na.rm=T)
esig = round(sum(rep.res$sig, na.rm=T)/2, digits=0)
ons = sum(!m$sig, na.rm=T)
ens = round(sum(!rep.res$sig, na.rm=T)/2, digits=0)
sig = c(osig, esig)
notsig = c(ons, ens)
stbl = rbind(sig, notsig)
colnames(stbl) = c('observed', 'expected')
fisher.test(stbl, alternative='greater')

#BOXPLOT
mt$pct = mt$x / 1e6 * 100
repBox=ggplot(data=mt) + 
	geom_boxplot(aes(x=sex, y=pct), col=c(female.col, male.col), lwd=1) +
	labs(y="% reads", x="Sex") +
	labs(subtitle='Repetitive elements') 
	# lims(y=c(6.1, 6.65))
plot(repBox)

#stats
tapply(mt$x, INDEX=mt$sex, mean)
male = mt$x[mt$sex=="M"]
female = mt$x[mt$sex=="F"]
t.test(x=male, y=female, alternative='greater')


# PLOT VOLCANO PLOT
res.df = data.frame(rep.res)
res.df$logp = -log(res.df$pvalue, 10)
s.res.df = res.df[res.df$padj < 0.1,]
male = rep('male (10)', nrow(s.res.df))
male[s.res.df$log2FoldChange < 0]<-'female (3)'
s.res.df$significant = male
volcano = ggplot(data=res.df, aes(x=log2FoldChange, y=logp)) + 
	geom_point(col='black', alpha=0.3) +
	geom_point(data=s.res.df, aes(x=log2FoldChange, y=logp, col=significant), size=2.5) +
	scale_colour_manual(values=c('firebrick','dodgerblue')) +
	xlab(bquote(log[2]*"(M:F)")) +
	ylab(bquote(-log[10]*"(p)" )) +
	# lims(x=c(-3.2, 3.2), y=c(0, 20)) +
	labs(subtitle='Repetitive elements') +
	theme(legend.position="none")
plot(volcano)
volcanoLenged = ggplot(data=res.df, aes(x=log2FoldChange, y=logp)) + 
	geom_point(col='black', alpha=0.3) +
	geom_point(data=s.res.df, aes(x=log2FoldChange, y=logp, col=significant), size=2.5) +
	scale_colour_manual(values=c('firebrick','dodgerblue')) +
	xlab(bquote(log[2]*"(M:F)")) +
	ylab(bquote(-log[10]*"(p)" )) +
	labs(subtitle='Repetitive\nelements') +
	lims(x=c(-3.2, 3.2), y=c(0,17)) 
plot(volcanoLenged)

plot_grid(denRep, repBox, ncol=2)
plot_grid(volcano, repBox, ncol=2)

#====================================================#

#------------------ PAIRWISE dN/dS ------------------#
#The ggplot boxplots are strange, so went with violin plots, with the median values drawn manually

#LOAD
ll=load('pairwisedNdS.Rdata')
#ad = 'all data' these are the dNdS values for all variants
#pd = 'private data' these are the dNdS values based only on private variants
#sind = 'sin dnds' these are the pairwise comparisons with sinensis for all varants
pd$species=sub('private', '', pd$species)
ll

do.wilcox = function(sp1, sp2, df, stat){
	stat1 = df[df$species == sp1, stat]
	stat2 = df[df$species == sp2, stat]
	sub = df[df$species %in% c(sp1, sp2),]
	w=wilcox.test(x=stat1, y=stat2)
	boxplot(sub[,stat]~sub[,'species'], outline=F, main=paste("p =",w$p.value))
	return(w)
}


#ALL VARINATS VIOLIN PLOT
#set up data to plot in order
spp=c('tym', 'sin', 'punY', 'punX')
num = 1:4
names(num) = spp
num
#!!should just make these ordered factors
ad$sppNum = ad$species
ad$sppNum[ad$species=='tym']<-1
ad$sppNum[ad$species=='sin']<-2
ad$sppNum[ad$species=='Y']<-3
ad$sppNum[ad$species=='X']<-4

#set up data for line segments at medians
meds = data.frame(tapply(ad$dNdS, INDEX=ad$species, function(x) median(x, na.rm=T)))
colnames(meds) = c('y')
meds$x1 = 1:4 - 0.2
meds$x2 = 1:4 + 0.2
all.dnds.violin = ggplot(data=ad) + 
	geom_violin(aes(x=sppNum, y=dNdS, fill=sppNum, col=sppNum), na.rm=T) +
	# geom_point(data=pmeds, aes(x=x, y=y), pch="-", size=10) +
	geom_segment(data=meds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1.5) +
	scale_color_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	scale_fill_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	scale_x_discrete(breaks=1:4, labels=names(num)) +
	lims(y=c(-.1, 0.5)) +
	labs(y='dN/dS', x='', subtitle='All variants') +
	theme(legend.position="none",axis.text.x = element_text(angle= xlabAngle))
print(all.dnds.violin)

#standard boxplot to doublecheck
quartz()
boxplot(ad$dNdS~ad$species, outline=F)


#PRIVATE VARIANTS VIOLIN PLOT
pd$sppNum = pd$species
pd$sppNum[pd$species=='tym']<-1
pd$sppNum[pd$species=='sin']<-2
pd$sppNum[pd$species=='Y']<-3
pd$sppNum[pd$species=='X']<-4

pmeds = data.frame(tapply(pd$dNdS, INDEX=pd$sppNum, function(x) median(x, na.rm=T)))
colnames(pmeds) = c('y')
pmeds$sppNum = rownames(pmeds)
pmeds$spp = spp
pmeds$x1 = 1:4 - 0.2
pmeds$x2 = 1:4 + 0.2
pmeds

private.dnds.violin = ggplot(data=pd) + 
	geom_violin(aes(x=sppNum, y=dNdS, fill=sppNum, col=sppNum), na.rm=T) +
	geom_segment(data=pmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_x_discrete(breaks=1:4, labels=names(num)) +
	lims(y=c(-.1, 0.5)) +
	labs(y='dN/dS', x='', subtitle ='Private variants') +
	theme(legend.position="none",axis.text.x = element_text(angle= xlabAngle))
print(private.dnds.violin)

#boxplot to doublecheck
quartz()
boxplot(pd$dNdS~pd$species, outline=F, ylim=c(-.1, 0.5))

#stats
do.wilcox('X', 'tym', pd, 'dNdS')
do.wilcox('X', 'sin', pd, 'dNdS')
do.wilcox('Y', 'sin', pd, 'dNdS')
do.wilcox('X', 'Y', pd, 'dNdS')


#-------- PROVEAN RESULTS --------#
ll=load("provean.Rdata")
#d = the provean scores
#psub = the provean scores for mutations private to each species/chromosomes
#ratios = the ratio of bad to ok mutations (bad <= -2.5 score)
#private.ratios = same but for private mutations
#tbl = the table of bad and ok mutations to use for chi square test
#private.tbl = same but for private mutations
colnames(d)[7]='species'
colnames(psub)[7]='species'
d$score = d$score*-1
psub$score = psub$score*-1
head(d)
head(psub)

#ALL VARIANTS
#set up ordered species names
num
d$sppNum = d$species
d$sppNum[d$sppNum =='tym']<-1
d$sppNum[d$sppNum =='sin']<-2
d$sppNum[d$sppNum =='y']<-3
d$sppNum[d$sppNum =='x']<-4

#set up medians
dmeds = data.frame(tapply(d$score, INDEX=d$sppNum, function(x) median(x, na.rm=T)))
colnames(dmeds) = c('y')
dmeds$sppNum = rownames(dmeds)
dmeds$spp = spp
dmeds$x1 = 1:4 - 0.2
dmeds$x2 = 1:4 + 0.2
dmeds


#plot violin plot
all.provean.violin = ggplot(data=d) +
	geom_violin(aes(x=sppNum, y=score, col=sppNum, fill=sppNum), na.rm=T) +
	geom_segment(data=dmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_x_discrete(breaks=1:4, labels=names(num)) +
	labs(y='deleteriousness', x='', subtitle='All variants') + 
	theme(legend.position='none',axis.text.x = element_text(angle= xlabAngle))
plot(all.provean.violin)


#set up frequencies of bad (< -2.5) and ok (> -2.5)




#PRIVATE VARIANTS
#set up ordered species names
num
psub$sppNum = psub $species
psub$sppNum[psub$sppNum =='tym']<-1
psub$sppNum[psub$sppNum =='sin']<-2
psub$sppNum[psub$sppNum =='y']<-3
psub$sppNum[psub$sppNum =='x']<-4

#set up medians
psmeds = data.frame(tapply(psub$score, INDEX=psub$sppNum, function(x) median(x, na.rm=T)))
colnames(psmeds) = c('y')
psmeds$sppNum = rownames(psmeds)
psmeds$spp = spp
psmeds$x1 = 1:4 - 0.2
psmeds$x2 = 1:4 + 0.2
psmeds

#plot violin plot
private.provean.violin = ggplot(data=psub) +
	geom_violin(aes(x=sppNum, y=score, col=sppNum, fill=sppNum), na.rm=T) +
	geom_segment(data=psmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_x_discrete(breaks=1:4, labels=names(num)) +
	labs(y='deleteriousness', x='', subtitle="Private variants") + 
	theme(legend.position='none',
		axis.text.x = element_text(angle= xlabAngle))
plot(private.provean.violin)



#standard boxplot to doublecheck
boxplot(d$score~d$species, outline=F)


#DO STATS

siny = psub[psub$sample %in% c('tym', 'y'),]
siny$sample=as.factor(siny$sample)
t.test(siny$score~ siny$sample)

head(psub)
do.wilcox('y', 'sin', psub, 'score')
do.wilcox('y', 'tym', psub, 'score')
do.wilcox('y', 'x', psub, 'score')
do.wilcox('x', 'sin', psub, 'score')
do.wilcox('x', 'tym', psub, 'score')

#fisher's exact
head(psub)

ysubs = psub[psub$species=='y',]
sinsubs = psub[psub$species=='sin',]
bady = sum(ysubs$bad, na.rm=T)
badsin = sum(ysubs$bad, na.rm=T)
oky = sum(ysubs$ok, na.rm=T)
oksin = sum(sinsubs$ok, na.rm=T)
tab=rbind(c(bady, badsin), c(oky, oksin))
colnames(tab) = c('y', 'sin')
rownames(tab) = c('bad', 'ok')
fisher.test(tab)

by=ysubs[ysubs$bad==1,]
bs=sinsubs[sinsubs$bad==1,]
sum(bs$mut %in% by$mut)
bs[bs$mut %in% by$mut,]

#---- PLOT ALL TOGETHER ----#
#DEPTH
plot(dp1)
plot(dp2)

#REPETITIVE ELEMENTS
# plot(volcanoLenged)
plot(volcano)
plot(repBox)


#dNdS 
#all
plot(all.dnds.violin)
plot(private.dnds.violin)


#PROVEAN
plot(all.provean.violin)
plot(private.provean.violin)
LETTERS[1:8]


#build final plot
plot_grid(dp1, dp2, all.dnds.violin, private.dnds.violin, volcano, repBox, all.provean.violin, private.provean.violin, ncol=4, labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), label_size=18)


#alternative formats
plot_grid(dp1, dp2, all.dnds.violin, private.dnds.violin, volcanoLenged, repBox, all.provean.violin, private.provean.violin, ncol=4, labels=c('A', 'B', 'E', 'F', 'C', 'D', 'G', 'H'))
plot_grid(dp1, dp2, volcano, repBox, all.dnds.violin, private.dnds.violin, all.provean.violin, private.provean.violin, ncol=2)

