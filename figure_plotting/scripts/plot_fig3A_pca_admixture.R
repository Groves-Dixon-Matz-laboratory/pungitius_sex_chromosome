#plot_figure3.R

library(ggplot2)
library(cowplot)
library(scales)
library(reshape)
setwd("~/gitreps/pungitius_sex_chromosome/figure_plotting")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")



show_col(hue_pal()(3))

##########################
########## PCAs ##########
##########################
TEXTLABS=F

#PCA FOR MALES AUTOSOME
ll=load("chrI_Pungitius_filt4.recode.vcf_pca.Rdata")
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
male.auto = plot_snp_pca(df.male, pca, cx=1, cy=2, group='group', shift=T, shift.table=shift.table)


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




###############################
########## ADMIXTURE ##########
###############################

tbl.auto=read.table("~/gitreps/pungitius_sex_chromosome/introgression/admixture/chrI_FULL.3.Q")
tbl.par=read.table("~/gitreps/pungitius_sex_chromosome/introgression/admixture/chrXII_PAR.3.Q")
tbl.sdr=read.table("~/gitreps/pungitius_sex_chromosome/introgression/admixture/chrXII_SDR.3.Q")
colnames(tbl.auto) = c('tym', 'pun', 'sin') #(set these based on the labeled plot below)
colnames(tbl.par) = c('tym', 'pun', 'sin')
colnames(tbl.sdr) = c('tym', 'pun', 'sin')
tbl.auto=tbl.auto[,c('pun', 'sin', 'tym')]
tbl.par=tbl.par[,c('pun', 'sin', 'tym')]
tbl.sdr=tbl.sdr[,c('pun', 'sin', 'tym')]
names = read.table("~/gitreps/pungitius_sex_chromosome/introgression/admixture/names.txt")
sex = read.table("~/gitreps/pungitius_sex_chromosome/metadata/pun_pheno.txt")
auto.in = format_admixture(tbl.auto, names, sex)
par.in = format_admixture(tbl.par, names, sex)
sdr.in = format_admixture(tbl.sdr, names, sex)
plot_admixture(auto.in)
plot_admixture(par.in)
plot_admixture(sdr.in)
legendPosition="top"
legendPosition="none"
auto = plot_admixture_bars_only(auto.in, legendPosition)
par = plot_admixture_bars_only(par.in, legendPosition)
sdr = plot_admixture_bars_only(sdr.in, legendPosition)

#plot clean admixture
plot_grid(auto, par,sdr,  ncol=3)


#plot all together
labs=LETTERS[1:6]
plot_grid(male.auto, male.par, male.sdr, auto, par,sdr,  ncol=3, labels = '', label_size=18)
