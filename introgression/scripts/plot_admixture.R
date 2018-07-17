#plot_admixture.R


library(ggplot2)
library(cowplot)
library(scales)
library(reshape)
setwd("~/gitreps/pungitius_sex_chromosome/introgression/admixture")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")

###############################
########## ADMIXTURE ##########
###############################
#upload data
tbl.auto=read.table("chrI_FULL.3.Q")
tbl.par=read.table("chrXII_PAR.3.Q")
tbl.sdr=read.table("chrXII_SDR.3.Q")

#format
colnames(tbl.auto) = c('tym', 'pun', 'sin') #(set these based on the labeled plot below)
colnames(tbl.par) = c('tym', 'pun', 'sin')
colnames(tbl.sdr) = c('tym', 'pun', 'sin')
tbl.auto=tbl.auto[,c('pun', 'sin', 'tym')]
tbl.par=tbl.par[,c('pun', 'sin', 'tym')]
tbl.sdr=tbl.sdr[,c('pun', 'sin', 'tym')]
names = read.table("names.txt")
sex = read.table("~/gitreps/pungitius_sex_chromosome/metadata/pun_pheno.txt")
auto.in = format_admixture(tbl.auto, names, sex)
par.in = format_admixture(tbl.par, names, sex)
sdr.in = format_admixture(tbl.sdr, names, sex)

#build plots
plot_admixture(auto.in)
plot_admixture(par.in)
plot_admixture(sdr.in)
legendPosition="top"
legendPosition="none"
auto = plot_admixture_bars_only(auto.in, legendPosition)
par = plot_admixture_bars_only(par.in, legendPosition)
sdr = plot_admixture_bars_only(sdr.in, legendPosition)

#plot panel
plot_grid(auto, par,sdr,  ncol=3, label_size=20)


