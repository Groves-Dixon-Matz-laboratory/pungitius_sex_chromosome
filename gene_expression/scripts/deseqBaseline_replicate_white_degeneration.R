#deseqBaseline_replicate_white.R
#purpose of this scrip is to replicate results shown in white 2014
#to make sure that the methods to be applied in reverse for pungitius
#are working correctly.


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
select.species = 'aculeatus'
select.males = maleAcus
select.females = femaleAcus

#SET UP SEPARATE COUNTS DATA FOR MALES AND FEMALES 
m.data = coldata[coldata$sppSex!=paste(select.species, "F", sep = "_"),]
f.data = coldata[coldata$sppSex!=paste(select.species, "M", sep = "_"),]
m.counts = counts[,!colnames(counts) %in% select.females]
f.counts = counts[,!colnames(counts) %in% select.males]


#RUN DESEQ
library(DESeq2)
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
m.r = results(m, contrast=c('species', 'aculeatus', 'pungitius'))
f.r = results(f, contrast=c('species', 'aculeatus', 'pungitius'))
head(m.r)
head(f.r)

#---- SAVE/LOAD ----#
save(m.r, f.r, counts, coldata, file='replicateWhite.Rdata')
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression")
ll=load('replicateWhite.Rdata')
#-------------------#

#ASSEMBLE FOLD DIFFERENCES
#simply take the difference between male and female log2s

sum(rownames(m.r)==rownames(f.r))==nrow(m.r) #double-check genes are same
ldat = data.frame(m.r$log2FoldChange, f.r$log2FoldChange)
rownames(ldat) = rownames(m.r)
colnames(ldat) = c('male', 'female')
plot(ldat$male~ldat$female)
ldat$diff = ldat$male - ldat$female
plot(density(ldat$diff, na.rm=T))


#GET GTF DATA AND TAG GENES WITH THEIR REGION TYPE
gdat = read.table("gene_locations.txt", header = T, stringsAsFactors=F)

#set autosome, par, stratum2, and stratum1
subChr='chrXIX'
gdat$type = 'autosomes'
gdat$type[gdat$chr==subChr]<-'stratum2'
gdat$type[gdat$chr==subChr & gdat$stop < 2.5e6]<-'PAR'
gdat$type[gdat$chr==subChr & gdat$start > 12e6]<-'stratum1'

#merge with expression data
ldat$geneId = rownames(ldat)
ddat = merge(gdat, ldat, by = 'geneId')


#LOOK AT MALE-FEMALE EXPRESSION RELATIVE TO OUTGROUP

#boxplots
par(mfrow=c(1,1))
boxplot(ddat$diff~ddat$type, main = 'log2(Male) - log2(Female)')
abline(h=0, lty=2)

i=ddat$diff[ddat$type=='stratum2']
a=ddat$diff[ddat$type=='autosomes']
i.lab = rep("st1", length(i))
a.lab = rep("autosome", length(a))
d = c(i,a)
d.lab = c(i.lab, a.lab)
x = data.frame(d,d.lab)
wilcox.test(d ~ d.lab, data=x) 




#plot scatterplots
sectionScatter = function(chr.type){
	cx = -10:10
	cy = -10:10
	lmCnt = lm(cy~cx)
	sub=ddat[ddat$type== chr.type,]
	plot(male~female, data=sub, main= chr.type)
	abline(lmCnt, lwd=0.5)
}

par(mfrow=c(2,2))
sectionScatter('autosomes')
sectionScatter('PAR')
sectionScatter('stratum2')
sectionScatter('stratum1')


#---- ENRICHMENT FOR DIFFERENTIAL EXPRESSION BY SEX ----#

#SUBSET THE DATA FOR JUST THE SPECIES OF INTEREST
head(counts)
head(coldata)
scounts = counts[,grep('SRR', colnames(counts))]
sdata = coldata[grep('SRR', coldata$sample),]

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
boxplot(gsdf$log2FoldChange~gsdf$type, ylab="Male:Female Log2 Difference", outline=F, axes=F)
axis(1, at = 1:4, labels=c('autosomes', 'PAR', 'stratum1', 'stratum2'))
axis(2)
abline(h=0, lty=2)

#boxplot for selected chromosome
select = 'chrXIX'
chrs = gsdf$chr
chrs[chrs!=select]<-'genome'
gsdf$abslog2 = abs(gsdf$log2FoldChange)
boxplot(gsdf$log2FoldChange~chrs)
boxplot(gsdf$abslog2~chrs, outline=F, ylab="Abs(M:F log2 difference)", axes=F, ylim=c(0,2.2))
axis(1, at = c(1,2), labels = c("chrXIX", "Autosomes"))
axis(2)
t.test(gsdf$abslog2~chrs, alternative='greater')


#look at frequency of significant differential expression
sig=gsdf[gsdf$sig,]
t=table(gsdf$chr)
barplot(t)

#SAVE FOR MAIN PLOTTING
save(s.r, rld, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/maleVfemale3spine.Rdata")


