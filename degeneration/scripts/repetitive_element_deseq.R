#repetitive_element_deseq.R

library('DESeq2')
library('ggplot2')
library('cowplot')
setwd("~/gitreps/pungitius_sex_chromosome/degeneration/repetitive_elements")

#LOAD DATA

raw=read.table("combinedFamilyCounts.tsv", sep="\t", comment.char="", strip.white=TRUE, header = T )


#FORMAT

#clean names
colnames(raw) = sub("_sorted.bam", "", colnames(raw))
head(raw)

#separate out counts and rep element info
counts = raw[,4:ncol(raw)]
families = as.character(raw$chr)

#some names are duplicates, tag to make rownames
dupCount = sum(duplicated(families))
families[duplicated(families)] = paste(families[duplicated(families)], "DUP2", sep="_")
rownames(counts) = families
counts = counts[,grep("pun", colnames(counts))]
counts = counts[!grepl("Unknown", rownames(counts)),] #remove unknown repetitive element types
head(counts)


#--------- get FPM ---------#
#upload overall counts
rdat = read.table("~/gitreps/pungitius_sex_chromosome/metadata/dna_raw_read_counts.tsv", header = T)
rdat = rdat[grep("pun", rdat$sample),]
m = rdat$readCount / 1e6
names(m) = rdat$sample
head(rdat)
print("Names match?")
sum(rdat$sample==colnames(counts))==ncol(counts)

#get fpm
fpm = sweep(counts, MARGIN=2, FUN="/", STATS=m)

#double-check
head(counts)
head(fpm)
m
counts[2,'pun11'] / m['pun11'] == fpm[2,'pun11']
counts[5,3] / m[3] == fpm[5,3]

#add project to coldata
sdat = read.table("~/gitreps/pungitius_sex_chromosome/metadata/genus_Pun_sex.txt")
colnames(sdat) = c('sample', 'sex')
rownames(sdat) = sdat$sample
samples = colnames(counts)
coldata=data.frame(sdat[samples,], row.names=samples)
head(coldata)
nrow(coldata)
ncol(counts)
sum(coldata$sample==colnames(fpm))==ncol(fpm)



#get ratio for mean FPKM between males and females
tdat = t(fpm)
head(tdat)
sum(coldata$sample==rownames(tdat))==ncol(fpm)
mdat = data.frame(t(apply(fpm, 1, function(x) tapply(x, INDEX=coldata$sex, mean))))
head(mdat)
hi = apply(mdat, 1, max)
m2=mdat
low = apply(m2, 1, min)
m2[m2==0] <- min(low[low>0])
m2$ratio = m2$M / m2$F
m2$log = log(m2$ratio, 2)


#LOOK AT CPM RESULTS
boxplot(m2$log)
boxplot(m2$log, outline=F);abline(h=0, lty=2, lwd=0.5)
plot(density(m2$log));abline(v=0, lty=2, lwd=0.5)


#RUN DEESEQ

#get size factors from PCR duplicate removed/mapped to whole genome dataset
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/genomeSizeFactors.Rdata") #this comes from bedtools_depth.R
print(genomeSizeFactors)
names(genomeSizeFactors) = tolower(names(genomeSizeFactors))
print("Names match?")
sum(names(genomeSizeFactors) == coldata$sample)==nrow(coldata)

#run DESEq
dds<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~sex))
sizeFactors(dds) <-genomeSizeFactors
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds)


#GET DESEQ RESULTS
res = results(dds, contrast = c('sex', 'M', 'F'), independentFiltering=F)
summary(res)
head(res[order(res$pvalue),], n=20)
rld = rlog(dds)
rld.df = assay(rld)



#VOLCANO PLOT

#set up dataframe
res.df = data.frame(res)
res.df$sig = res.df$padj < 0.1
res.df$sig[is.na(res.df$sig)]<-FALSE
res.df$logp = -log(res.df$pvalue, 10)
res.df$logBaseMean = log(res.df$baseMean)
head(res.df)

#plot
volcano = ggplot(data=res.df, aes(x=log2FoldChange, y=logp)) + 
	scale_colour_manual(values=c('black','red')) +
	geom_point(data=res.df, aes(color = sig), alpha=0.4) +
	xlab("log2 fold difference") +
	ylab("-log10 p-value") +
	lims(x=c(-3.2, 3.2)) 
plot(volcano)


#plot horizontally
g = ggplot(data=res.df, aes(y=log2FoldChange, x= logBaseMean)) + 
	scale_colour_manual(values=c('black','red')) +
	geom_point(data=res.df, aes(color = sig), alpha=0.4) +
	ylab("log2 fold difference") +
	xlab("Mean Coverage") +
	theme_bw() #+ ggtitle('Male vs Famale Fold Enrichment')
plot(g)


#high female repeats
w = res[res$log2FoldChange < -1,]
w

#significant repeats
s = res[!is.na(res$padj) & res$padj < 0.1,]
s

#look at density of significant
plot(density(s$log2FoldChange), xlim=c(-3.5, 3.5))


#log at overall log fold differences
boxplot(res.df$log2FoldChange)
boxplot(res.df$log2FoldChange, outline=F);abline(h=0, lty=2, lwd=0.5)
plot(density(res.df$log2FoldChange, na.rm=T))


#check agreement between DESEq and differences in mean FPKM
mgdat = merge(res.df, m2, by = 0)
plot(mgdat$log~mgdat$log2FoldChange, xlab="Deseq log2 diff", ylab='FPKM log2 diff')



#ASSESS THE TOTAL PROPORTION OF MALES AND FEMALES THAT ARE REPETITIVE ELEMENTS

#assemble total counter per million reads
tfpm = apply(fpm, 2, sum)
mt = merge(tfpm, coldata, by = 0)

#### PLOT RESULTS FIGURES ####
res.df = data.frame(res)
res.df$logp = -log(res.df$pvalue, 10)
s.res.df = res.df[res.df$padj < 0.1,]
male = rep('male', nrow(s.res.df))
male[s.res.df$log2FoldChange < 0]<-'female'
s.res.df$significant = male
volcano = ggplot(data=res.df, aes(x=log2FoldChange, y=logp)) + 
	geom_point(col='black', alpha=0.5) +
	geom_point(data=s.res.df, aes(x=log2FoldChange, y=logp, col=significant), size=2.25) +
	scale_colour_manual(values=c('firebrick','dodgerblue')) +
	xlab("log2 M:F") +
	ylab("-log10 p-value") +
	lims(x=c(-3.2, 3.2), y=c(0,17)) 
plot(volcano)

#DENSITY PLOT FOR SIGNIFICANT
den=density(s$log2FoldChange)
den.df = data.frame(x=den$x, y=den$y)
g1=ggplot(data=den.df, aes(x=x, y=y)) + 
	geom_vline(aes(xintercept=0), linetype='dashed') +
	geom_line(col='grey') + 
	geom_area(fill='grey') +
	lims(x=c(-3.3, 3.3)) +
	# scale_x_continuous(limits = c(-0.03, 0.22)) +
	# geom_segment(data=d12, aes(x=d, y=y1, xend=d, yend=y2), lineend='round', lwd=1.2) +
	labs(x="", y="")
plot(g1)
#stats
pos = s$log2FoldChange > 0
male = sum(pos)
female = sum(!pos)



#BOXPLOT
mt$pct = mt$x / 1e6 * 100
g=ggplot(data=mt) + 
	geom_boxplot(aes(x=sex, y=pct)) +
	labs(y="% reads", x="Sex") +
	lims(y=c(6.1, 6.65))
plot(g)
tapply(mt$x, INDEX=mt$sex, mean)
male = mt$x[mt$sex=="M"]
female = mt$x[mt$sex=="F"]
t.test(x=male, y=female, alternative='greater')



#save data for plotting with all the other degeration results
rep.res = res
save(rep.res, mt, s, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/repetitiveElements.Rdata")

