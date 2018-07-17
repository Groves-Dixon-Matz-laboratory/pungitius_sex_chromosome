#deseqBaseline_white.R
#test for sex-biased expression in the White dataset

library(DESeq2)
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression")

#LOAD DATA (output from initialize_basicAlignment_counts.R)
ll=load("deseqBaselineInput.Rdata")
head(counts)

#ADD COLUMNS TO COLDATA
coldata$sppSex = paste(coldata$species, coldata$sex, sep="_")
malePuns = coldata$sample[coldata$sppSex=="pungitius_M"]
femalePuns = coldata$sample[coldata$sppSex=="pungitius_F"]
maleAcus = coldata$sample[coldata$sppSex=="aculeatus_M"]
femaleAcus = coldata$sample[coldata$sppSex=="aculeatus_F"]
coldata


#---- LOOK FOR EVIDENCE OF DEGENERATION ----#

#CHOOSE SPECIES YOU WANT TO LOOK AT
select.species = 'pungitius'
select.males = malePuns
select.females = femalePuns


#SET UP SEPARATE COUNTS DATA FOR MALES AND FEMALES 
m.data = coldata[coldata$sex=='M',]
f.data = coldata[coldata$sex=='F',]
m.counts = counts[,colnames(counts) %in% maleAcus | colnames(counts) %in% malePuns]
f.counts = counts[,colnames(counts) %in% femaleAcus | colnames(counts) %in% femalePuns]


#RUN DESEQ
#here we are testing for species differences within each sex
#if Y has degenerated trancription without compensation
#it will be lower than ancestral state seen in 3spine for males


#set up data
m=DESeqDataSetFromMatrix(m.counts,
	colData = m.data, 
	design = formula(~ species))

f=DESeqDataSetFromMatrix(f.counts,
	colData = f.data, 
	design = formula(~ species))

#run
m<-DESeq(m)
f<-DESeq(f)

#get results
m.r = results(m, contrast=c('species', 'pungitius', 'aculeatus'))
f.r = results(f, contrast=c('species', 'pungitius', 'aculeatus'))
head(m.r)
head(f.r)

#GET GTF DATA AND TAG GENES WITH THEIR REGION TYPE
gdat = read.table("gene_locations.txt", header = T, stringsAsFactors=F)

#set autosome, par, stratum2, and stratum1
subChr='chrXII'
gdat$type = '1' #autosomes
gdat$type[gdat$chr==subChr]<-'3' #inversion
gdat$type[gdat$chr==subChr & gdat$stop < 2.5e6]<-'2'   #PAR
gdat$type[gdat$chr==subChr & gdat$start > 18.9e6]<-'4' #low.recomb


#ASSEMBLE FOLD DIFFERENCES
#simply take the difference between male and female log2s

sum(rownames(m.r)==rownames(f.r))==nrow(m.r) #double-check genes are same
ldat = data.frame(m.r$log2FoldChange, f.r$log2FoldChange) #these are the species level differences divided by sex
rownames(ldat) = rownames(m.r)
colnames(ldat) = c('male', 'female')
plot(ldat$male~ldat$female)
ldat$diff = ldat$male - ldat$female
plot(density(m.r$log2FoldChange, na.rm=T)) #overall difference with ancestral (3spine) in males

#merge with expression data
ldat$geneId = rownames(ldat)
ddat = merge(gdat, ldat, by = 'geneId')
sdr.ddat = ddat[ddat$type==3,]
plot(sdr.ddat$male~sdr.ddat$female, ylab='Male Change', xlab='Female Change')

#merge up the species comparisons for males and females
#males 
mr.df = data.frame(m.r)
mr.df$geneId = rownames(mr.df)
g.mr = merge(mr.df, gdat, by = 'geneId')
sdr.mr = g.mr[g.mr$type==3,]
plot(density(sdr.mr$log2FoldChange, na.rm=T));abline(v=0)

#females
fr.df = data.frame(f.r)
fr.df$geneId = rownames(fr.df)
g.fr = merge(fr.df, gdat, by = 'geneId')
sdr.fr = g.fr[g.fr$type==3,]
plot(density(sdr.fr$log2FoldChange, na.rm=T));abline(v=0)




#---- SAVE/LOAD ----#
save(m.r, f.r, counts, coldata, gdat, ddat, sdr.mr, file='pungitiusAncCompare.Rdata')
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression")
ll=load('pungitiusAncCompare.Rdata')


#LOOK AT MALE-FEMALE EXPRESSION RELATIVE TO OUTGROUP

#boxplots for directional variation
par(mfrow=c(1,1))
boxplot(ddat$diff~ddat$type, main = 'log2(Male) - log2(Female)', outline=F, axes=F, ylab='log2(Male) - log2(Female)')
axis(1, at=1:4, labels=c('autosomes', 'PAR', 'inversion', 'low.recomb'))
axis(2)
abline(h=0, lty=2)
i=ddat$diff[ddat$type=='3']
a=ddat$diff[ddat$type=='1']
t.test(x=i, y=a, alternative='less')
i.lab = rep("inversion", length(i))
a.lab = rep("autosome", length(a))
d = c(i,a)
d.lab = c(i.lab, a.lab)
x = data.frame(d,d.lab)
wilcox.test(d ~ d.lab, data=x) 

#boxplots for absolute variation
par(mfrow=c(1,1))
head(ddat)
ddat$type2 = ddat$type
ddat$type2[ddat$type2==4]<-3
ddat$absDiff = abs(ddat$diff)
boxplot(ddat$absDiff~ddat$type2, main = 'log2(Male) - log2(Female)', outline=F, axes=F, ylab='log2(Male) - log2(Female)')
axis(1, at=1:3, labels=c('autosomes', 'PAR', 'SDR'))
axis(2)
abline(h=0, lty=2)

x=ddat[ddat$type2 %in% c(1, 3),]
t.test(x$absDiff~x$type2)
x=ddat[ddat$type2 %in% c(2, 3),]
t.test(x$absDiff~x$type2)
x=ddat[ddat$type2 %in% c(1, 2),]
t.test(x$absDiff~x$type2)

#not in correct direction. Mean inflated by a few really high values. Upregulated genes on Y?



#plot scatterplots
sectionScatter = function(num, chr.type){
	cx = -10:10
	cy = -10:10
	lmCnt = lm(cy~cx)
	sub=ddat[ddat$type== num,]
	plot(male~female, data=sub, main= chr.type)
	abline(lmCnt, lwd=0.5)
}

par(mfrow=c(2,2))
sectionScatter('1', 'autosomes')
sectionScatter('2', 'PAR')
sectionScatter('3', 'inversion')
sectionScatter('4', 'low.recomb')


#---- ENRICHMENT FOR DIFFERENTIAL EXPRESSION BY SEX ----#

#SUBSET THE DATA FOR JUST THE SPECIES OF INTEREST
head(counts)
head(coldata)
scounts = counts[,grep('DRR', colnames(counts))]
sdata = coldata[grep('DRR', coldata$sample),]

#RUN DESEQ
s=DESeqDataSetFromMatrix(scounts,
	colData = sdata, 
	design = formula(~ sex))

#run
s<-DESeq(s)

#get results
s.r = results(s, contrast=c('sex', 'M', 'F'), independentFiltering=F)
rld = rlog(s)


#VOLCANO PLOT
library(ggplot2)
library(cowplot)
sdf = na.omit(data.frame(s.r))
sdf$sig = sdf$padj < 0.1
g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
	geom_point()
plot(g)

#GET GENE LOCATIONS

head(gdat)
sdf$geneId = rownames(sdf)
gsdf = merge(gdat, sdf, by = 'geneId')


#LOOK AT DEGREE OF DIFF EXPRESSION ACCROSS REGIONS

#boxplot by type
par(mfrow=c(1,1))
boxplot(gsdf$log2FoldChange~gsdf$type, main="Male:Female Log2 Difference")
abline(h=0, lty=2)

#boxplot for selected chromosome
select = 'chrXII'
chrs = gsdf$chr
chrs[chrs!=select]<-'genome'
gsdf$abslog2 = abs(gsdf$log2FoldChange)
gsdf$chrType = chrs
boxplot(gsdf$log2FoldChange~chrs)
boxplot(gsdf$abslog2~chrs, outline=F, ylab="Abs(M:F log2 difference)", axes=F)
axis(1, at = c(1,2), labels = c("chrXII", "Autosomes"))
axis(2)
t.test(gsdf$abslog2~chrs, alternative='greater')


#look at frequency of significant differential expression
sig=gsdf[gsdf$sig,]
t=table(gsdf$chr)
barplot(t)


#SAVE FOR MAIN FIGURE PLOTTING
save(s.r, rld, gdat, gsdf, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/maleVfemalePungitius.Rdata")


