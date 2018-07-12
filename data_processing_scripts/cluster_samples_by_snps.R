#!/usr/bin/env Rscript
#pca_chrom_snps.R
#Groves Dixon
#9-11-17
#Use Adegenet to build PCAs from a set of vcf files

#requires file with list of file prefixes
# library(plotrix)
library(vcfR)
library(adegenet)


#parse arguments
args = commandArgs(trailingOnly=TRUE)
v = args[1]
sexIds = args[2]
outPrefix = args[3]



#upload the assigned sex
sdat = read.table(sexIds)




allRes=data.frame()
statRes=data.frame()


print(paste("Loading VCF file", v))
gll=vcfR2genlight(read.vcfR(v))
print("Done loading VCF.")
#look at genlight object
print("Converting data to matrix...")
x=as.matrix(gll)
gi=as.genind(x)


print("Matching up sexes")
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
sex[sex==2]<-"Female"
pop(gll)=sex

clus=find.clusters(gll, max.n.clus=3, n.pca=4, stat=c("BIC"), choose.n.clust=FALSE, criterion="diffNgroup") #check what stats are for different k
clusRes=find.clusters(gll, max.n.clus=3, n.pca=4, stat=c("BIC"), choose.n.clust=TRUE, n.clust=2) #coerse clustering into two groups

#assemble clustering results
res=data.frame(clusRes$grp)
res$sex.name=names
print("Do names in results match with original vcf?")
sum(rownames(res) == names) == nrow(res)
res$sex=sex
res$sexNum=sexNum
print("------")
print("results:")
print(res)
c = cor(x=as.numeric(res$clusRes.grp), y=res$sexNum)
dbic=clus$Kstat[2]-clus$Kstat[1]

print("Statsitcs assigned.")


#write out the clustering statsitcs
print("writing out stat results...")
statRes=data.frame()
statRes=rbind(statRes, c(clus$Kstat[1], clus$Kstat[2], clus$Kstat[3], dbic, c, v))
colnames(statRes) = c('bic1grp', 'bic2grp', 'bic3grp', 'deltaBIC_2v1_grps', 'corr', 'file')
statOut = paste(outPrefix, "stat.tsv", sep="_")
write.table(statRes, file=statOut, sep="\t", quote=F, row.names=F)

#write out the clustering results
print("writing out cluster results")
resOut=paste(outPrefix, "clusterRes.tsv", sep="_")
write.table(res, file=resOut, sep="\t", quote=F)



