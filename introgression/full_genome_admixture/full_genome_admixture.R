#full_genome_admixture.R
#In addition to the single autosome and sex chromosome admixture
# data shown analysis, we ran admixture on all chromosomes.
#This script plots those results


library(ggplot2)
library(cowplot)
library(scales)
library(reshape)


#----------------- K3 -----------------#
#admixture results using K=3 ancestral populations

setwd("~/gitreps/pungitius_sex_chromosome/full_genome_admixture/k3_10Kb_thinned")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")


#upload data
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
names = read.table("names.txt")$V1
head(chromlist)
c=chromlist[1]
infile = paste(c, "10K_thinned.3.Q", sep="_")
tbl.auto = read.table(infile)
tbl.auto$chr = c
tbl.auto$name = names
for (c in chromlist[2:length(chromlist)]){
	infile = paste(c, "10K_thinned.3.Q", sep="_")
	to.add = read.table(infile)
	to.add$chr = c
	to.add$name=names
	tbl.auto = rbind(tbl.auto, to.add)
}


#format
colnames(tbl.auto) = c('tym', 'pun', 'sin', 'chr', 'name') #(set these based on the labeled plot below)
head(tbl.auto)
tbl.auto$spp = substr(tbl.auto$name, start=1, stop=3)
head(tbl.auto)

#look at sharing
nopun = tbl.auto[tbl.auto$spp !='Pun',]
boxplot(nopun$pun~nopun$spp)

#look at the outliers
some = nopun[nopun$pun>1e-5,]
table(some$spp)

nosin = tbl.auto[tbl.auto$spp !='Sin',]
boxplot(nosin$sin~nosin$spp)
head(nosin)
outliers = nosin[nosin$sin>0.1,]
nrow(outliers)

head(tbl.auto)

#plot barplot for each chromosome
sex = read.table("~/gitreps/pungitius_sex_chromosome/metadata/pun_pheno.txt")
names = read.table("names.txt")
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
chr='chrXII'
legendPosition="none"
tsub = tbl.auto[tbl.auto$chr==chr, c('pun', 'sin', 'tym')]
subin = format_admixture(tsub, names, sex)
plot_admixture(subin)
sp = plot_admixture_bars_only(subin, legendPosition)
sp = sp + ggtitle(chr)
plot(sp)

plist = list()

plot_grid(plotlist=plist)

#assemble list of plots from each chromosome
plist = list()
lindex = 0
for (chr in chromlist[2:length(chromlist)]){
	lindex = lindex + 1
	print(paste(chr, "..", sep="."))
	tsub = tbl.auto[tbl.auto$chr==chr, c('pun', 'sin', 'tym')]
	subin = format_admixture(tsub, names, sex)
	sp = plot_admixture_bars_only(subin, legendPosition)
	sp = sp + labs(subtitle=chr)
	plist[[lindex]] = sp
}

#plot them
plot_grid(plotlist=plist)

#----------------- K2 -----------------#

setwd("~/gitreps/pungitius_sex_chromosome/full_genome_admixture/k2")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")


#upload data
tbl.auto=read.table("chrXIII.2.Q")
tbl.auto$name = names
tbl.auto$spp = substr(tbl.auto$name, start=1, stop=3)
head(tbl.auto)

boxplot(tbl.auto$V1~tbl.auto$spp)
boxplot(tbl.auto$V2~tbl.auto$spp)
