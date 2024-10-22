#plot_figure7_exploratory.R
#This is a messy version of the plotting script that explores different ways of looking at the results
#For final clean copy see ~/gitreps/pungitius_sex_chromosome/figure_plotting/plot_figure7_expression.R

library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/")


#-------- MORE DIFFERENTIAL EXPRESSION ON 12 --------#

#load the ancestral comparison (white)
ll=load('~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/pungitiusAncCompare.Rdata')
#ddat = assembly of male:female differences with chromosomal region labeled 
#sdr.mr = the pungitius male vs 3spine male differences within the SDR


#format the data
ddat$absDiff = abs(ddat$diff)  #absolute M:F variation
ddat$type2=ddat$type            #set up region with rest of chr12 (not PAR)
ddat$type2[ddat$type==4]<-3    #
head(ddat)
ddat$sex.chrom = ddat$chr
ddat$sex.chrom[ddat$sex.chrom!='chrXII']<-'autosomes'
ddat$sex.chrom<-factor(ddat$sex.chrom, levels=c('chrXII', 'autosomes'), ordered=TRUE)

#M:F DIFF IN SEX VS AUTOSOME

#boxplot
mfChrom = ggplot(data=ddat) +
	geom_boxplot(aes(x= sex.chrom, y=absDiff), outlier.shape=26) + 
	lims(y=c(0,0.75)) +
	labs(x='', y="abs(M:F)")
plot(mfChrom)

#violin
mfViolin = ggplot(ddat) + 
	geom_violin(aes(x= sex.chrom, y=diff))
plot(mfViolin)


#UPLOAD THE MALE:FEMALE COMPARISON

ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/maleVfemalePungitius.Rdata")
ll
#gsdf = the male vs female in pungitius differences with SDR labeled


#format
sdr = gsdf[gsdf$type==3,]
den=density(sdr$log2FoldChange)
den.df = data.frame(x=den$x, y=den$y)
male.df = den.df[den.df$x >=max(den.df$x[den.df$x<0]),]
female.df = den.df[den.df$x <= 0,]
male.df$df = 'higher in males'
female.df$df = 'higher in females'
male.col = 'dodgerblue'
female.col = 'firebrick'



denRepFM=ggplot(data=female.df) + 
	# geom_vline(aes(xintercept=0), linetype='dashed') +
	geom_line(aes(x=x, y=y, color=df), na.rm=T) + 
	geom_area(aes(x=x, y=y, fill=df), na.rm=T) +
	geom_line(data=male.df, aes(x=x,y=y, color=df), na.rm=T) + 
	geom_area(data=male.df, aes(x=x,y=y, fill=df), na.rm=T) +
	scale_color_manual(values=c(female.col, male.col)) +
	scale_fill_manual(values=c(female.col, male.col)) +
	theme(legend.position=c(0.6, 0.75),
		legend.title=element_blank()) +
	labs(x=bquote(log[2]~"M:F pungitius"), y="density")
plot(denRepFM)



#look for degeneration of transcription from Y as preferential loss in male pun expression
head(sdr.mr)
den=density(sdr.mr$log2FoldChange)
den.df = data.frame(x=den$x, y=den$y)
pun.df = den.df[den.df$x >=max(den.df$x[den.df$x<0]),]
gacu.df = den.df[den.df$x <= 0,]
pun.df$df = 'higher in pungitius males'
gacu.df$df = ' higher in G.aculeatus males'


denRepMG=ggplot(data= pun.df) + 
	geom_line(aes(x=x, y=y, color=df), na.rm=T) + 
	geom_area(aes(x=x, y=y, fill=df), na.rm=T) +
	geom_line(data= gacu.df, aes(x=x,y=y, color=df), na.rm=T) + 
	geom_area(data= gacu.df, aes(x=x,y=y, fill=df), na.rm=T) +
	theme(legend.position=c(0.6, 0.75),
		legend.title=element_blank()) +
	labs(x=bquote(log[2]~" pungitius male : G.aculeatus male"), y="density")
plot(denRepMG)

#plot both
plot_grid(denRepFM, denRepMG, labels=c('A', 'B'), ncol=1, label_size =16)


#-------- DISTRIBUTIONS IN M:F FOLD COVERAGE --------#

#load m:f fold differences DNA
ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/maleVfemale_DNA.Rdata")
ll
dsr = data.frame(dna.sr$log2FoldChange)
colnames(dsr) = c('fc')
dsr$dset = 'dna'
dsr$geneId = rownames(dna.sr)
head(dsr)


#load differential expression by sex (von Hippel)
ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/maleVfemale_vonHippel.Rdata")
vsr = data.frame(s.r$log2FoldChange)
colnames(vsr) = c('fc')
vsr$dset = 'pelvic'
vsr$geneId = rownames(s.r)
head(vsr)

#load differential expression by sex (white)
ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/maleVfemalePungitius.Rdata")
wsr = data.frame(s.r$log2FoldChange)
colnames(wsr) = c('fc')
wsr$dset = 'brain'
wsr$geneId = rownames(s.r)
head(wsr)

#load differential expression by sex in 3-spine (white)
ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression/maleVfemale3spine.Rdata")
tsr = data.frame(s.r$log2FoldChange)
colnames(tsr) = c('fc')
tsr$dset = 'acu'
tsr$geneId = rownames(s.r)
head(tsr)


#CONTROL FOR MAPPING EFFICIENCY PEVLIC DATA
#to control for changes in mapping efficiency, substract
#the DNA fold differences from the RNA fold differences
dd = merge(vsr, dsr, by = 'geneId')
dd$fc = dd$fc.x - dd$fc.y
dd$dset = 'pelvic-dna'
pdsr=dd[,c('fc', 'dset', 'geneId')]
head(pdsr)

# see if the dna fold differences predict the RNA differences
ddg=merge(dd, gdat, by = 'geneId')

#for sdr
ddsdr = ddg[ddg$type==3,]
plot(ddsdr$fc.x~ddsdr$fc.y)
lm1=lm(ddsdr$fc.x~ddsdr$fc.y)
abline(lm1, col='red')
summary(lm1)

#for autosomes
ddauto = ddg[ddg$type==1,]
plot(ddauto$fc.x~ ddauto$fc.y)
lm2=lm(ddauto$fc.x~ ddauto$fc.y)
abline(lm2, col='red')
summary(lm2)


#CONTROL FOR MAPPING EFFICIENCY brain DATA
dd = merge(wsr, dsr, by = 'geneId')
dd$fc = dd$fc.x - dd$fc.y
dd$dset = 'brain-dna'
ddsr=dd[,c('fc', 'dset', 'geneId')]
head(ddsr)

# see if the dna fold differences predict the RNA differences
ddg=merge(dd, gdat, by = 'geneId')

#for sdr
#This shows that there is indeed an effect of mapping efficiency
ddsdr = ddg[ddg$type==3,]
plot(ddsdr$fc.x~ddsdr$fc.y)
lm1=lm(ddsdr$fc.x~ddsdr$fc.y)
abline(lm1, col='red')
summary(lm1)

#for autosomes
ddauto = ddg[ddg$type==1,]
plot(ddauto$fc.x~ ddauto$fc.y)
lm2=lm(ddauto$fc.x~ ddauto$fc.y)
abline(lm2, col='red')
summary(lm2)



#ASSEMBLE ALL THE MALE:FEMALE DATA TOGETHER
sr = rbind(dsr, vsr, wsr, tsr, ddsr, pdsr)
head(sr)


#add genomic locations
head(gdat)
msr = merge(sr, gdat, by='geneId', all.x=T)
msr$sex = 'auto'
# # msr$sex[msr$chr=='chrXII']<-'chrXII' #to use full chrXII
msr$sex[msr$type==3]<-'SDR' #to use SDR
msr$sex = factor(msr$sex, levels=c('SDR', 'auto'), ordered=T)
msr$sexdf = paste(msr$sex, msr$dset, sep="_")
msr$absDiff = abs(msr$fc)
msr$sexdf = factor(msr$sexdf, levels=c('SDR_brain', 'auto_brain', 'SDR_acu', 'auto_acu', 'SDR_pelvic', 'auto_pelvic', 'SDR_dna', 'auto_dna', 'SDR_brain-dna', 'auto_brain-dna', 'SDR_pelvic-dna', 'auto_pelvic-dna'), ordered=T)
boxplot(msr$absDiff~msr$sexdf, outline=F)


#subset
keep = c('SDR_pelvic-dna', 'auto_pelvic-dna', 'SDR_brain-dna', 'auto_brain-dna', 'SDR_acu', 'auto_acu')
ssr = msr[msr$sexdf %in% keep,]
ssr$sexdf = factor(ssr$sexdf, levels=keep, ordered=T)
smeds = tapply(ssr$absDiff, INDEX=ssr$sexdf, function(x) median(x, na.rm=T))
meddf = data.frame(y=smeds)
meddf$sexdf = names(smeds)
lineLen = 0.3
meddf$x1 = 1:6 - lineLen
meddf$x2 = 1:6 + lineLen

#plot main boxplot
absBox = ggplot(data=ssr) +
	geom_boxplot(aes(x=sexdf, y=absDiff, fill=sex), outlier.shape=26, lwd=0.75) +
	lims(y=c(-.1, 2)) + 
	labs(x='', y=bquote("|"*log[2]~"M:F|")) +
	theme(axis.text.x = element_blank(),
		axis.line.x = element_line(size=0),
		axis.ticks.x = element_blank(), 
		legend.title=element_blank(),
		legend.position=c(0.7, 0.8))
plot(absBox)


do.wilcox = function(sp1, sp2, df, col, stat){
	stat1 = df[df[,col] == sp1, stat]
	stat2 = df[df[,col] == sp2, stat]
	sub = df[df[,col] %in% c(sp1, sp2),]
	w=wilcox.test(x=stat1, y=stat2)
	boxplot(sub[,stat]~sub[,col], outline=F, main=paste("p =",w$p.value))
	return(w)
}
do.wilcox('SDR_brain-dna', 'auto_brain-dna', ssr, 'sexdf', 'absDiff')
do.wilcox('SDR_pelvic-dna', 'auto_pelvic-dna', ssr, 'sexdf', 'absDiff')
do.wilcox('SDR_acu', 'auto_acu', ssr, 'sexdf', 'absDiff')

brainp = do.wilcox('SDR_brain-dna', 'auto_brain-dna', ssr, 'sexdf', 'absDiff')$p.value
pelvicp = do.wilcox('SDR_pelvic-dna', 'auto_pelvic-dna', ssr, 'sexdf', 'absDiff')$p.value
acup = do.wilcox('SDR_acu', 'auto_acu', ssr, 'sexdf', 'absDiff')$p.value

pvals = c(brainp, pelvicp, acup)


#build the same with violin
absViolin = ggplot(data=ssr) +
	geom_violin(aes(x=sexdf, y=absDiff, fill=sex, col=sex), na.rm=T) +
	geom_segment(data=meddf, aes(x=x1, y=y, xend=x2, yend=y), lwd=1.5) +
	lims(y=c(-.1, 2)) + 
	labs(x='', y='abs(M:F)') + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_line(size=0))
plot(absViolin)

#---------- ALLELE-SPECIFIC EXPRESSION DATA ----------#

#load xy expression white
ll=load("~/gitreps/pungitius_sex_chromosome/gene_expression/allele_specific_expression/xy_white.Rdata")
#m.w = merged dataset with sex differences and XY within male differences for white dataset
ll


#scatterplot for Y:X vs Male:Female
xyscatter.w = ggplot(data=m.w) +  
	geom_smooth(data= m.w, aes(log2FoldChange.x,log2FoldChange.y), method='lm', formula=y~x, se=T) +
	geom_point(aes(x=log2FoldChange.x, y= log2FoldChange.y), alpha=1) +
	lims(y=c(-2, 1.5)) +
	labs(x=bquote(log[2]~"Y:X in males"), y=bquote(log[2]~"M:F") )
plot(xyscatter.w)
lmw=lm(m.w$log2FoldChange.y~m.w$log2FoldChange.x)
summary(lmw)


#repeat but with the DNA controlled values
dcont = ssr[ssr$dset=='brain-dna',]
xy = m.w[,c('Row.names', 'log2FoldChange.x')]
colnames(xy) = c('geneId', 'fcYX')
head(xy)
m=merge(dcont, xy, by = 'geneId')
head(m)
lmx = lm(m$fc~m$fcYX)
plot(m$fc~m$fcYX)
summary(lmx)


#keep the dna controlled one
xyscatter.w2 = ggplot(m) +  
	geom_smooth(aes(x=fcYX, y=fc), method='lm', formula=y~x, se=T) +
	geom_point(aes(x=fcYX, y= fc), alpha=1) +
	lims(y=c(-2, 1.5)) +
	labs(x=bquote(log[2]~"Y:X in males"), y=bquote(log[2]~"M:F") )
plot(xyscatter.w2)
lmw=lm(m$fc~m$fcYX)
summary(lmw)



#load xy expression von Hippel
ll=load("~/gitreps/stickle_back_sex_chromosomes/allele_specific_expression/xy_vonHippel.Rdata")
#m.vh = merged dataset with sex differences and XY within male differences for von Hippel dataset
ll

xyscatter.vh = ggplot(data=m.vh) +  
	geom_smooth(data= m.vh, aes(log2FoldChange.x,log2FoldChange.y), method='lm', formula=y~x, se=T) +
	geom_point(aes(x=log2FoldChange.x, y= log2FoldChange.y), alpha=1) +
	lims(y=c(-2, 1.5)) +
	labs(x=bquote(log[2]~"Y:X in males"), y=bquote(log[2]~"M:F") )
plot(xyscatter.vh)
head(m.vh)
lmvh=lm(m.vh$log2FoldChange.y~m.vh$log2FoldChange.x)
summary(lmvh)


#save data for plotting
save(absBox, ssr, m, m.vh, xyscatter.w2, xyscatter.vh, denRepMG, denRepFM, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/fig7_files.Rdata")





