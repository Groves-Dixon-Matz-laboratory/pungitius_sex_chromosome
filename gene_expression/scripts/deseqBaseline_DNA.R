#deseqBaseline_vonHippel.R


library(DESeq2)
library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/baseline_gene_expression")

#upload read counts
counts = read.table("dna_gene_counts.tsv", header = T, row.names='geneID')
colnames(counts) = sub(".counts.txt", "", colnames(counts))
counts=counts[, grepl("pun", colnames(counts))]
dirtyNames = colnames(counts)
cnames = sapply(dirtyNames, function(x) strsplit(x, "_")[[1]][3])
colnames(counts) = cnames
head(counts)
dim(counts)


#remove non-gene objects from counts data
to.remove=rownames(counts)[grep("__", rownames(counts))]
counts[to.remove,]
dim(counts)
counts= counts[!rownames(counts) %in% to.remove,]
dim(counts)


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


#set up sex data
sdat = read.table("~/gitreps/pungitius_sex_chromosome/metadata/genus_Pun_sex.txt")
colnames(sdat) = c('sample', 'sex')
sdat=sdat[grepl('pun', sdat$sample),]
sum(sdat$sample %in% colnames(counts))
rownames(sdat) = sdat$sample
coldata=sdat[colnames(counts),]
coldata
dim(coldata)
sum(coldata$sample==colnames(counts))==ncol(counts)


#RUN DESEQ
s=DESeqDataSetFromMatrix(counts,
	colData = coldata, 
	design = formula(~ sex))

#run
s<-DESeq(s)

#get results
dna.sr = results(s, contrast=c('sex', 'M', 'F'), independentFiltering=F)
head(dna.sr)
boxplot(dna.sr$log2FoldChange, outline=F)


#VOLCANO PLOT
sdf = na.omit(data.frame(dna.sr))
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
boxplot(gsdf$abslog2~chrs, outline=F, ylab="Abs(M:F log2 difference)", axes=F)
axis(1, at = c(1,2), labels = c("chrXII", "Autosomes"))
axis(2)
t.test(gsdf$abslog2~chrs, alternative='greater')

#look at frequency of significant differential fold coverage
sig=gsdf[gsdf$sig,]
t=table(gsdf$type)
barplot(t)

#SAVE FOR MAIN FIGURE PLOTTING
save(dna.sr, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/maleVfemale_DNA.Rdata")


