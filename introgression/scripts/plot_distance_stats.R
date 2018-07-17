#plot_distances_stats.R

library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/introgression/distance_stats")
source("~/gitreps/pungitius_sex_chromosome//pun_sex_chrom_functions.R")


#load chromosome names
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
chroms = as.character(cdat$chroms)

#set up variables to examine chromosome regions from 4Mb to 17Mb to match SDR
to.remove = c("chrUn", "chrM", "PAR", "SDR"); selection = c("chrXII") #just subset from 12 instead of reading sdr and par
leftBound = 4
rightBound = 17

#----- PLOT COMPARISONS BETWEEN Y and SIN -----#
#these are intended to show that it's introgression and not ILS

ys = readinDistance("all_y_sin_distances.tsv", 'ys')
fst = plot_distance_density(ys, 'chrXII', 'fst', c(0, 0.15, 0.3), c(0, 0.3))
dxy = plot_distance_density(ys, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(0, 0.016))
dtilde = plot_distance_density(ys, 'chrXII', 'dtilde', c(0, 0.006, 0.012), c(0, 0.0125))
plot_grid(fst, dxy, dtilde, nrow=1, labels = LETTERS[1:3], label_size=16)


#----- PLOT COMPARISONS BETWEEN SIN, TYM, AND FEMALES PUNS -----#

#purpose here is to see what the distance looks like
#greater distance for the sex chromosome would support 12 as ancestral
fs = readinDistance("all_femalePun_sin_distances.tsv", 'fs')
ft = readinDistance("all_femalePun_tym_distances.tsv", 'ft')
fst1 = plot_distance_density(fs, 'chrXII', 'fst', c(0, 0.15, 0.3), c(0, 0.3))
fst2 = plot_distance_density(ft, 'chrXII', 'fst',  c(0, 0.15, 0.3), c(0, 0.3))
dxy1 = plot_distance_density(fs, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(-0.003, 0.018))
dxy2 = plot_distance_density(ft, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(-0.003, 0.018))
dtilde1 = plot_distance_density(fs, 'chrXII', 'dtilde', c(0, 0.006, 0.012), c(-0.003, 0.015))
dtilde2 = plot_distance_density(ft, 'chrXII', 'dtilde',  c(0, 0.006, 0.012), c(-0.003, 0.015))
plot_grid(fst1, dxy1, dtilde1, fst2, dxy2, dtilde2)


#Above are the main two figures
#below is stuff to plot other datasets


#--------------- 100 KB WINDOWS ---------------#

setwd("/Users/grovesdixon/gitreps/stickle_back_sex_chromosomes/results/my_dxy/100Kb_windows")
to.remove = c("chrUn", "chrM"); selection = c("PAR", "SDR", "chrXII")             #line for PAR and SDR
to.remove = c("chrUn", "chrM", "PAR", "chrXII"); selection = c("SDR")   #line for SDR only
to.remove = c("chrUn", "chrM", "PAR", "SDR"); selection = c("chrXII")  #line for chrXII only



#change to normalize to SDR-like region
leftBound = 4;rightBound = 17
#default for whole chromosomes
# leftBound = 0;rightBound = 1e3

#read in data
readinDistance = function(fileName, pair){
	x=read.table(fileName, sep="\t", header = T, stringsAsFactors=F)
	x=x[!x$chr %in% to.remove,]
	x$pair = pair
	x$mb = (x$lefts + x$rights) / 2 / 1e6
	x=x[x$mb > leftBound & x$mb < rightBound,]
	return(x)
}

#upload comparisons
yt = readinDistance("all_y_tym_distances.tsv", 'yt')
ys = readinDistance("all_y_sin_distances.tsv", 'ys')
yf = readinDistance("all_y_female_distances.tsv", 'yf')
ms = readinDistance("all_malePun_sin_distances.tsv", 'ms')
mt = readinDistance("all_malePun_tym_distances.tsv", 'mt')
fs = readinDistance("all_femalePun_sin_distances.tsv", 'fs')
ft = readinDistance("all_femalePun_tym_distances.tsv", 'ft')
st = readinDistance("all_sin_tym_distances.tsv", 'st')


ll=load("/Users/grovesdixon/gitreps/stickle_back_sex_chromosomes/metadata/chromList.Rdata")
ll


#select the dataset
df = ms
df = ys
df = yt
df = fs
df = ft
df = st
df=yf


#select the statistic to plot
stat = 'dtilde'
stat = 'hap.dtilde'
stat = 'fst'
stat = 'dxy'
stat = 'dminusx'


#plot density as before
df.auto = df[!df$chr %in% selection,]
den = density(df.auto[,stat]) #get densities for chromosomes only
den.df = data.frame(x=den$x, y=den$y)

sex = na.omit(df[df$chr %in% selection, stat])
den.sex = density(sex)
sex.df = data.frame(x= den.sex$x, y= den.sex$y)

mn.sex = data.frame(x=mean(sex, na.rm=T), y1=0, y2=max(den.sex$y))
mn.auto = data.frame(x=mean(df.auto[,stat]), y1=0, y2=max(den.sex$y))

g1=ggplot(data= sex.df, aes(x=x, y=y)) + 
	geom_area(data=sex.df, fill='black', alpha=0.75) +
	# geom_segment(data=mn.sex, aes(x=x, y=y1, xend=x, yend=y2), col='black', lty=2) +
	# geom_segment(data=mn.auto, aes(x=x, y=y1, xend=x, yend=y2), col='grey', lty=2) +
	geom_area(data=den.df, fill='grey', alpha=0.75) +
	labs(x=stat, y="")
plot(g1)




#----- PLOT COMPARISONS BETWEEN Y and SIN -----#
#these are intended to show that it's introgression and not ILS
fst = plot_distance_density(ys, 'chrXII', 'fst', c(0, 0.15, 0.3), c(0, 0.3))
dxy = plot_distance_density(ys, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(0, 0.016))
dtilde = plot_distance_density(ys, 'chrXII', 'dtilde', c(0, 0.006, 0.012), c(0, 0.0125))
plot_grid(fst, dxy, dtilde, nrow=1, labels = LETTERS[1:3], label_size=16)


#----- PLOT COMPARISONS BETWEEN SIN, TYM, AND FEMALES PUNS -----#

#purpose here is to see what the distance looks like
#greater distance for the sex chromosome would support 12 as ancestral


fst1 = plot_distance_density(fs, 'chrXII', 'fst', c(0, 0.15, 0.3), c(0, 0.3))
fst2 = plot_distance_density(ft, 'chrXII', 'fst',  c(0, 0.15, 0.3), c(0, 0.3))
dxy1 = plot_distance_density(fs, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(-0.003, 0.018))
dxy2 = plot_distance_density(ft, 'chrXII', 'dxy', c(0, 0.008, 0.016), c(-0.003, 0.018))
dtilde1 = plot_distance_density(fs, 'chrXII', 'dtilde', c(0, 0.006, 0.012), c(-0.003, 0.015))
dtilde2 = plot_distance_density(ft, 'chrXII', 'dtilde',  c(0, 0.006, 0.012), c(-0.003, 0.015))
plot_grid(fst1, dxy1, dtilde1, fst2, dxy2, dtilde2)




#get a legend
a=1:10
b=1:10
c=c(rep('SDR', 5), rep('autosomes', 5))
d=data.frame(a,b)
g=ggplot(data=d) + 
	geom_point(aes(x=a,y=b, color=c), pch=15, size=3, alpha=0.75) + 
	scale_color_manual(breaks=c("SDR","autosomes"), values=c('grey', 'black'))
	theme(legend.title=element_blank())
plot(g)







