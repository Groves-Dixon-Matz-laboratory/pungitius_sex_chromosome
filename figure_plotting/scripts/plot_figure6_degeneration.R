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

#------------------ PAIRWISE dN/dS against 3-spine ------------------#
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
#set up data to plot species in order
spp=ad$species
spp[spp=="Y"]<-"punY"
spp[spp=="X"]<-"punX"
ad$species=spp
ad$species = factor(ad$species, levels=c("tym", "sin", "punY", "punX"), ordered=T)
levels(ad$species)

#set up data for line segments at medians
meds = data.frame(tapply(ad$dNdS, INDEX=ad$species, function(x) median(x, na.rm=T)))
colnames(meds) = c('y')
meds$x1 = 1:4 - 0.2
meds$x2 = 1:4 + 0.2
all.dnds.violin = ggplot(data=ad) + 
	geom_violin(aes(x=species, y=dNdS, fill=species, col=species), na.rm=T) +
	# geom_point(data=pmeds, aes(x=x, y=y), pch="-", size=10) +
	geom_segment(data=meds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1.5) +
	scale_color_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	scale_fill_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	# scale_x_discrete(breaks=1:4, labels=names(num)) +
	lims(y=c(-.1, 0.5)) +
	labs(y='dN/dS', x='', subtitle='All variants') +
	theme(legend.position="none",axis.text.x = element_text(angle= xlabAngle))
print(all.dnds.violin)

#standard boxplot to doublecheck the figure
quartz()
boxplot(ad$dNdS~ad$species, outline=F)


#make a similar figure for pairwise dS against 3-spine
meds = data.frame(tapply(ad$dS, INDEX=ad$species, function(x) median(x, na.rm=T)))
colnames(meds) = c('y')
meds$x1 = 1:4 - 0.2
meds$x2 = 1:4 + 0.2
all.dS.violin = ggplot(data=ad) + 
	geom_violin(aes(x=species, y=dS, fill=species, col=species), na.rm=T) +
	# geom_point(data=pmeds, aes(x=x, y=y), pch="-", size=10) +
	geom_segment(data=meds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1.5) +
	scale_color_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	scale_fill_manual(values=c(sin.col, tym.col, female.col, male.col)) +
	lims(y=c(-.1, 0.5)) +
	labs(y='pairwise dS vs G. aculeatus', x='', subtitle='All variants') +
	theme(legend.position="none",axis.text.x = element_text(angle= xlabAngle))
print(all.dS.violin)

#standard boxplot to doublecheck the figure
quartz()
boxplot(ad$dS~ad$species, outline=F)

#PRIVATE VARIANTS VIOLIN PLOT

#set up data to plot species in order
spp=pd$species
spp[spp=="Y"]<-"punY"
spp[spp=="X"]<-"punX"
pd$species=spp
pd$species = factor(pd$species, levels=c("tym", "sin", "punY", "punX"), ordered=T)
levels(pd$species)

#set up medians
pmeds = data.frame(tapply(pd$dNdS, INDEX=pd$species, function(x) median(x, na.rm=T)))
colnames(pmeds) = c('y')
pmeds$spp = rownames(pmeds)
pmeds$x1 = 1:4 - 0.2
pmeds$x2 = 1:4 + 0.2
pmeds

private.dnds.violin = ggplot(data=pd) + 
	geom_violin(aes(x=species, y=dNdS, fill=species, col=species), na.rm=T) +
	geom_segment(data=pmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	# scale_x_discrete(breaks=1:4, labels=names(num)) +
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


#REPEAT THE PRIVATE VIOLIN PLOT FOR dS against 3-spine

#set up medians
pmeds = data.frame(tapply(pd$dS, INDEX=pd$species, function(x) median(x, na.rm=T)))
colnames(pmeds) = c('y')
pmeds$spp = rownames(pmeds)
pmeds$x1 = 1:4 - 0.2
pmeds$x2 = 1:4 + 0.2
pmeds

private.dS.violin = ggplot(data=pd) + 
	geom_violin(aes(x=species, y=dS, fill=species, col=species), na.rm=T) +
	geom_segment(data=pmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	# scale_x_discrete(breaks=1:4, labels=names(num)) +
	lims(y=c(-.01, 0.1)) +
	labs(y='pairwise dS vs G. aculeatus', x='', subtitle ='Private variants') +
	theme(legend.position="none",axis.text.x = element_text(angle= xlabAngle))
print(private.dS.violin)

#boxplot to doublecheck
quartz()
boxplot(pd$dNdS~pd$species, outline=F, ylim=c(-.1, 0.5))



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
# d$score = d$score*-1
# psub$score = psub$score*-1
head(d)
head(psub)

#ALL VARIANTS

#set up ordered species names
spp=d$species
spp[spp=="y"]<-"punY"
spp[spp=="x"]<-"punX"
d$sppNum = factor(spp, levels=c("tym", "sin", "punY", "punX"), ordered=T)
levels(d$sppNum)

#set up medians
dmeds = data.frame(tapply(d$score, INDEX=d$sppNum, function(x) median(x, na.rm=T)))
colnames(dmeds) = c('y')
dmeds$sppNum = rownames(dmeds)
dmeds$x1 = 1:4 - 0.2
dmeds$x2 = 1:4 + 0.2
dmeds


#plot violin plot
all.provean.violin = ggplot(data=d) +
	geom_violin(aes(x=sppNum, y=score, col=sppNum, fill=sppNum), na.rm=T) +
	geom_segment(data=dmeds, aes(x=x1, y=y, xend=x2, yend=y), lwd=1) +
	scale_color_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	# scale_x_discrete(breaks=1:4, labels=names(num)) +
	labs(y='Provean score', x='', subtitle='All variants') + 
	theme(legend.position='none',axis.text.x = element_text(angle= xlabAngle))
plot(all.provean.violin)


#set up frequencies of bad (< -2.5) and ok (> -2.5)




#PRIVATE VARIANTS
#set up ordered species names
spp=psub$species
spp[spp=="y"]<-"punY"
spp[spp=="x"]<-"punX"
psub$sppNum = factor(spp, levels=c("tym", "sin", "punY", "punX"), ordered=T)
levels(psub$sppNum)

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
	labs(y='Provean score', x='', subtitle="Private variants") + 
	theme(legend.position='none',
		axis.text.x = element_text(angle= xlabAngle))
plot(private.provean.violin)



#standard boxplots to doublecheck
boxplot(d$score~d$sppNum, outline=F)
boxplot(psub$score~psub$sppNum, outline=F)


#PLOT THRESHOLD BASED WAY
prov.dat = psub
CUT=-2.5

get_ratios = function(prov.dat, cutoff){
	prov.dat$bad = as.numeric(prov.dat$score <= cutoff)
	prov.dat$ok = as.numeric(prov.dat$score > cutoff)
	sums = tapply(prov.dat$bad, INDEX= prov.dat$sppNum, function(x) sum(x, na.rm=T))
	oksums = tapply(prov.dat$ok, INDEX= prov.dat$sppNum, function(x) sum(x, na.rm=T))
	private.ratios = data.frame( (sums / (sums + oksums)))
	colnames(private.ratios) = c('pscore')
	private.ratios$sppNum = factor(rownames(private.ratios), ordered=T, levels = c('tym', 'sin', 'punY', 'punX'))
	private.tbl = rbind(sums, oksums)
	res = list(private.ratios, private.tbl)
	return(res)
}



all.res = get_ratios(d, CUT)
private.res = get_ratios(psub, CUT)
a.ratios =all.res[[1]]
p.ratios = private.res[[1]]
a.counts = all.res[[2]]
p.counts = private.res[[2]]

all.prov.barplot = ggplot(a.ratios) + 
	geom_bar(aes(y=pscore, x=sppNum, fill=sppNum), stat='identity') +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	labs(subtitle="All variants", y = "provean < -2.5", x='') +
	theme(legend.position='none',
		axis.text.x = element_text(angle= xlabAngle))
plot(all.prov.barplot)

p.prov.barplot = ggplot(p.ratios) + 
	geom_bar(aes(y=pscore, x=sppNum, fill=sppNum), stat='identity') +
	scale_fill_manual(values=c(tym.col, sin.col, male.col, female.col)) +
	labs(subtitle = "Private variants", y = "provean < -2.5", x='') +
	theme(legend.position='none',
		axis.text.x = element_text(angle= xlabAngle))
plot(p.prov.barplot)

#DO STATS

siny = psub[psub$sppNum %in% c('sin', 'punY'),]
siny$sample=as.factor(siny$sppNum)
t.test(siny$score~ siny$sample)

head(psub)
do.wilcox('y', 'sin', psub, 'score')
do.wilcox('y', 'tym', psub, 'score')
do.wilcox('y', 'x', psub, 'score')
do.wilcox('x', 'sin', psub, 'score')
do.wilcox('x', 'tym', psub, 'score')

#fisher's exact
chisq.test(p.counts[,c('punY','sin')])
chisq.test(p.counts[,c('punY','punX')])
chisq.test(p.counts[,c('sin','punX')])
chisq.test(private.tbl[,c('sin','y')])


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

plot_grid(all.prov.barplot, p.prov.barplot, labels=c("A", "B"))


#build final plot
plot_grid(dp1, dp2, all.dnds.violin, private.dnds.violin, volcano, repBox, all.provean.violin, private.provean.violin, ncol=4, labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), label_size=18)


#alternative formats
plot_grid(dp1, dp2, all.dnds.violin, private.dnds.violin, volcanoLenged, repBox, all.provean.violin, private.provean.violin, ncol=4, labels=c('A', 'B', 'E', 'F', 'C', 'D', 'G', 'H'))
plot_grid(dp1, dp2, volcano, repBox, all.dnds.violin, private.dnds.violin, all.provean.violin, private.provean.violin, ncol=2)


#
