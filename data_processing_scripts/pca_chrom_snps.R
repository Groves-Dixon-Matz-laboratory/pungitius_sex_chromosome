#!/usr/bin/env Rscript
#pca_chrom_snps.R
#Groves Dixon
#9-11-17
#Use Adegenet to build PCAs from a set of vcf files

#requires file with list of file prefixes
#library(plotrix)
library(vcfR)
library(adegenet)


#parse arguments
args = commandArgs(trailingOnly=TRUE)
vcfInput = args[1]
sexIds = args[2]
N.CORES = as.numeric(args[3])


print("Running R script pca_chrom_snps.R...")
print(paste("Number of cores to use =", N.CORES))

#upload phenotype file
print(paste("Loading phenotype file", sexIds))
sdat = read.table(sexIds)



#upload the vcf
print(paste("Loading VCF file", vcfInput))
gll=vcfR2genlight(read.vcfR(vcfInput))
print("Done loading VCF.")

#look at genlight object
print("Converting data to matrix...")
x=as.matrix(gll)
gi=as.genind(x)


#assign sex
print("Matching up phenotypes")
names = gll@ind.names
print("Names from VCF:")
print(names)
sdat = read.table(sexIds)
colnames(sdat) = c('sample', 'SAMPLE', 'sex')
rownames(sdat) = sdat$sample
s=sdat[names,]

print("Sample names match between files?")
print(sum(rownames(s) == names) == nrow(s))
sex=s$sex
sexNum = s$sex
sex[sex==1]<-"Male"
sex[sex==0]<-"Female"
pop(gll)=sex


#set up colors based on ggplot standard
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col.set=gg_color_hue(2)
cols=col.set[sexNum]


#run PCA
print("Running PCA...")
pca=glPca(gll,nf=2, n.cores=N.CORES)
p=pca$scores
print("Done.")

#PLOT RESULTS

#point labels
pdfOut = sub('.vcf', '_labels.pdf', vcfInput)
pdf(pdfOut)
plot(p[,'PC2']~p[,'PC1'], bg=cols, pch=26, col='black', cex=1.5, xlab="PC1", ylab="PC2", main=vcfInput)
text(x=p[,'PC1'], y=p[,'PC2'], labels=names, col=cols)
dev.off()

#phenotype elipses
pdfOut = sub('.vcf', '_elipse.pdf', vcfInput)
pdf(pdfOut)
s.class(pca$scores,pop(gll),col=col.set,axesell=F,cstar=0,grid=F, label=unique(pop(gll)))
title(main=vcfInput)
dev.off()



#save results
outName=paste(sub(".vcf", "", vcfInput), 'pca.Rdata', sep='_')
save(p, file=outName)

