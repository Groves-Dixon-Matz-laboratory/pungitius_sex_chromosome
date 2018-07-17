#deseqBaseline_vonHippel.R


library('DESeq2')
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression")

#upload read counts
counts = read.table("pungitius_rnaseq_baseline_counts.tsv", header = T, row.names='geneID')
colnames(counts) = sub(".counts.txt", "", colnames(counts))
counts=counts[, grepl("SRR297", colnames(counts))]
head(counts)


#remove non-gene objects from counts data
to.remove=rownames(counts)[grep("__", rownames(counts))]
counts[to.remove,]
dim(counts)
counts= counts[!rownames(counts) %in% to.remove,]
dim(counts)

#get count totals
tots = apply(counts, 2, sum)
hist(tots)
write.table(data.frame(tots), file="counted_on_genes.tsv", sep="\t", quote=F) 
mean(tots)/1e6


#remove genes with low coverage
cut=2
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]



#set up coldata
sample = colnames(counts)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#set up coldata
sex = substrRight(sample, 1)
species=sample
species[grepl("SRR297", species)]<-'pungitius'
species[grepl("DRR02313", species)]<-'pungitius'
species[grepl("LS_2158", species)]<-'pungitius'
species[grepl("SRR19", species)]<-'aculeatus'
dataset=sample
dataset[grepl("SRR297", dataset)]<-'endocrine'
dataset[grepl("DRR02313", dataset)]<-'white.pun'
dataset[grepl("LS_2158", dataset)]<-'jun'
dataset[grepl("SRR19", dataset)]<-'white.acu'
tissue=sample
tissue[grepl("SRR297", tissue)]<-'pelvic'
tissue[grepl("DRR02313", tissue)]<-'brain'
tissue[grepl("LS_2158", tissue)]<-'liver'
tissue[grepl("SRR19", tissue)]<-'brain'
coldata = data.frame(sample, species, sex, dataset, tissue)
coldata


#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))

#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds)

#get variance stabilized counts and save them
rld = rlog(dds)
rld.df=assay(rld)
colnames(rld.df) = colnames(counts)

#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have geneswith too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
rld=rld[rownames(rld.df),]
dim(rld.df)
dim(rld)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#build sample heatmaps 
library(pheatmap)
quartz()
pheatmap(cor(rld.df))

#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
# 
#=====================================================================================

#Remove outliers by setting a branch cut threshold
# Plot a line to show the cut
cut.height = 150
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut.height, col = "red", lty = 2);
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 4)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
keepSampleNames = rownames(datExpr0)[keepSamples]
outlierNames = rownames(datExpr0)[clust==0]
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #number of samples left after outlier removal
print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
outlierNames
print(paste(nSamples, "samples were kept"))


#replot heatmap without outlier
rld.df = rld.df[, !colnames(rld.df) %in% outlierNames]
pheatmap(cor(rld.df))




#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')
counts=counts[,!colnames(counts) %in% outlierNames]
coldata=coldata[!coldata$sample %in% outlierNames,]


#---- DIFFERENTIAL EXPRESSION BY SEX ----#

#SUBSET THE DATA FOR JUST THE SPECIES OF INTEREST
head(counts)
head(coldata)


#RUN DESEQ
s=DESeqDataSetFromMatrix(counts,
	colData = coldata, 
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
gdat = read.table("gene_locations.txt", header = T, stringsAsFactors=F)
head(gdat)
sdf$geneId = rownames(sdf)
gsdf = merge(gdat, sdf, by = 'geneId')

#boxplot for selected chromosome
select = 'chrXII'
chrs = gsdf$chr
chrs[chrs!=select]<-'genome'
gsdf$abslog2 = abs(gsdf$log2FoldChange)
gsdf$type=chrs
boxplot(gsdf$log2FoldChange~chrs)
boxplot(gsdf$abslog2~chrs, outline=F, ylab="Abs(M:F log2 difference)", axes=F, ylim=c(0,2.2))
axis(1, at = c(1,2), labels = c("chrXII", "Autosomes"))
axis(2)
t.test(gsdf$abslog2~chrs, alternative='greater')

#look at frequency of significant differential expression
sig=gsdf[gsdf$sig,]
t=table(gsdf$type)
barplot(t)



#SAVE FOR MAIN PLOTTING
save(s.r, rld, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/maleVfemale_vonHippel.Rdata")












