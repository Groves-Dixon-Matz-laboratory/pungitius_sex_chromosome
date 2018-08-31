#plot_pcas.R
#This script plots PCA results generated with basic_snp_pca.R

library(ggplot2)
library(cowplot)
library(scales)
setwd("~/gitreps/pungitius_sex_chromosome/introgression/pca_results/weird_sin")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")


##########################
########## PCAs ##########
##########################
TEXTLABS=F
#PCA FOR MALES full chr12
ll=load("sample.vcf_pca.Rdata")
pcad = data.frame(pca$scores)
pcad$sample = rownames(pcad)
pcad$species = substr(rownames(pcad), start=1, stop=3)
head(pcad)


pcad$species[pcad$sample=="Sin8"]<-'sin8'
pcad$species[pcad$sample=="Sin20"]<-'sin20'
g = ggplot(data=pcad) + geom_point(aes(x=PC1, y=PC2, color=species, shape=species))
plot(g)

df.male = data.frame(pca$scores)
df.male = df.male[order(rownames(df.male)),]
sex = read.table("~/gitreps/pungitius_sex_chromosome/metadata/pun_pheno.txt", row.names = 1)
m = merge(df.male, sex, by = 0, sort=T, all.x=T)
sum(m$Row.names == rownames(df.male)) == nrow(df.male)
s = m$V3
s[s=='1']<-'F'
s[s=='2']<-'M'
s[is.na(s)]<-''
df.male$group = paste(substr(tolower(rownames(df.male)), start=1,stop=3), s, sep='')
shift.table = data.frame(y=c(0, 10, 0, 10), x=c(17, 5, 10, 5))
malepc = plot_snp_pca(df.male, pca, cx=1, cy=2, group='group', shift=T, shift.table=shift.table)
plot(malepc)

#PCA FOR MALES AUTOSOME (chr1)
ll=load("chrI_Pungitius_filt4.recode.vcf_pca.Rdata")
df.male = data.frame(pca$scores)
df.male = df.male[order(rownames(df.male)),]
m = merge(df.male, sex, by = 0, sort=T, all.x=T)
sum(m$Row.names == rownames(df.male)) == nrow(df.male)
s = m$V3
s[s=='1']<-'F'
s[s=='2']<-'M'
s[is.na(s)]<-''
df.male$group = paste(substr(tolower(rownames(df.male)), start=1,stop=3), s, sep='')
shift.table = data.frame(y=c(0, 10, 0, 10), x=c(17, 5, 10, 5))
male.auto = plot_snp_pca(df.male, pca, cx=1, cy=2, group='group', shift=T, shift.table=shift.table)
plot(male.auto)

#PCA FOR MALE PAR
ll=load("chrXII_Pungitius_filt4_PAR.recode.vcf_pca.Rdata")
df.male = data.frame(pca$scores)
df.male = df.male[order(rownames(df.male)),]
m = merge(df.male, sex, by = 0, sort=T, all.x=T)
sum(m$Row.names == rownames(df.male)) == nrow(df.male)
s = m$V3
s[s=='1']<-'F'
s[s=='2']<-'M'
s[is.na(s)]<-''
df.male$group = paste(substr(tolower(rownames(df.male)), start=1,stop=3), s, sep='')
shift.table = data.frame(y=c(-3, 2, -3, 1), x=c(2, 2, -1, 3))
male.par = plot_snp_pca(df.male, pca, cx=1, cy=2, group='group', shift=T, shift.table=shift.table)
plot(male.par)

#PCA FOR MALE SDR
ll=load("chrXII_Pungitius_filt4_SDR.recode.vcf_pca.Rdata")
df.male = data.frame(pca$scores)
df.male = df.male[order(rownames(df.male)),]
m = merge(df.male, sex, by = 0, sort=T, all.x=T)
sum(m$Row.names == rownames(df.male)) == nrow(df.male)
s = m$V3
s[s=='1']<-'F'
s[s=='2']<-'M'
s[is.na(s)]<-''
df.male$group = paste(substr(tolower(rownames(df.male)), start=1,stop=3), s, sep='')
shift.table = data.frame(y=c(0, 7, 0, 8), x=c(17, 0, 10, 0))
male.sdr = plot_snp_pca(df.male, pca, cx=1, cy=2, group='group', shift=T, shift.table=shift.table)
plot(male.sdr)

