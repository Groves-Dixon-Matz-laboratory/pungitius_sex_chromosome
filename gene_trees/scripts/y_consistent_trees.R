#y_consistent_trees.R
#examine phylogenetic results looking gene trees consistent with sex chromosomes

setwd("~/gitreps/pungitius_sex_chromosome/gene_trees/pre_Y_inference/")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")

#MAKE SELECTIONS FOR WHICH SPECIES/SEX TO ANALYZE

#pick the sex (male means looking for Ys, female for Ws)
sex='female'
sex='male'

#pick the species you want to look at
spp = 'po'
spp = 'tym'
spp = 'sin'
spp = 'pun'

#set up name
infile = paste(paste(spp, sex, sep="_"), 'allBiggestMonoRes.tsv', sep="_")
print(infile)

#read in data
dat = read.table(infile, header = T, stringsAsFactors=F)
dat = dat[!dat$CHROM %in% c("chrM", "chrUn"),]
head(dat)

#add normalized X values

ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
head(dat)
dat$norm.x = dat$BIN_START
for (chr in chromlist){
	len = max(dat$BIN_START[dat$CHROM==chr])
	dat$norm.x[dat$CHROM==chr]<-dat$BIN_START[dat$CHROM==chr] / len
}
head(dat)

#add numerics for sex chromosome consistent
sexCon = dat$SDRbool
sexCon[sexCon=='False']<-0
sexCon[sexCon=='True']<-1
dat$sexCon = as.numeric(sexCon)

#--------------- plot multipanel scatterplots ---------------#
x=dat$BIN_START
dat$mb = x/1e6
CUT = 0.99
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='mb', ycol='biggestMonoPhyletic', float= CUT, xlab="Position", ylab='Density Male-specific')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='mb', ycol='sexCon', float= CUT, xlab="Position", ylab='Density Male-specific')


#SAVE THE RESULTS FOR PLOTTING MAIN FIGURE
outName = '~/gitreps/pungitius_sex_chromosome/figure_plotting/pun_Ylike.Rdata'
if (spp=='pun'){
	save(dat, file=outName)
}




