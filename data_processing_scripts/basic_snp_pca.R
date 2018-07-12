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
vcfInput = args[1] #table with sample labels
N.CORES = as.numeric(args[2])


print("Running R script pca_chrom_snps.R...")
print(paste("Number of cores to use =", N.CORES))


#upload the vcf
print(paste("Loading VCF file", vcfInput))
gll=vcfR2genlight(read.vcfR(vcfInput))
print("Done loading VCF.")

#look at genlight object
print("Converting data to matrix...")
x=as.matrix(gll)
gi=as.genind(x)


#run PCA
print("Running PCA...")
pca=glPca(gll,nf=6, n.cores=N.CORES)
p=pca$scores
print("Done.")

outFile = paste(vcfInput, "pca.Rdata", sep="_")
print(paste("Saving results as", outFile))
save(pca, file=outFile)