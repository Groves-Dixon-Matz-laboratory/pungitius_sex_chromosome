#full_dS.R
#plot dS between species from autosomes to get idea of how closely related they are
library(ggplot2)
library(cowplot)
library(scales)
setwd("~/gitreps/pungitius_sex_chromosome/degeneration/full_genome_dS")


#upload gene location data
gdat = read.table("gene_locations.txt", header = T)
c12 = as.character(gdat[gdat$chr=='chrXII', 'geneId'])
head(gdat)
head(c12)

#upload pairwise dNdS vs 3spine
ddat = read.table("ref_pairwise_dNdS.tsv", header=T)
ddat$geneId = sapply(as.character(ddat$gene), function(x) return(strsplit(x, "_")[[1]][1]))
dim(ddat)
ddat=ddat[!ddat$geneId %in% c12,]
head(ddat)
dim(ddat)

boxplot(ddat$dS~ddat$species, outline=F, ylab = "pairwise dS vs G. aculeatus")
mns = tapply(ddat$dS, INDEX=ddat$species, function(x) mean(x, na.rm=T))
mns
length(unique(ddat$geneId))


#upload pairwise dNdS vs 3spine
ddats = read.table("sin_pairwise_dNdS.tsv", header=T, stringsAsFactors=F)
ddats=ddats[ddats$species!='ref',]
ddats$geneId = sapply(as.character(ddats$gene), function(x) return(strsplit(x, "_")[[1]][1]))
ddats=ddats[!ddats$geneId %in% c12,]
head(ddats)
dim(ddats)

boxplot(ddats$dS~ddats$species, outline=F, ylab = "pairwise dS vs P. sinensis")
mns = tapply(ddats$dS, INDEX=ddats$species, function(x) mean(x, na.rm=T))
table(na.omit(ddats)$species)
mns


#show the sinensis dS
ll=load("~/gitreps/pungitius_sex_chromosome/figure_plotting/pairwisedNdS.Rdata")
boxplot(sind[,'dS']~sind$species, outline=F, ylab="pairwise dS vs sinensis", main='vs sinensis')



#plot all three together
par(mfrow=c(1,1))
boxplot(ddat$dS~ddat$species, outline=F, ylab = "pairwise dS vs G. aculeatus")
boxplot(ddats$dS~ddats$species, outline=F, ylab = "pairwise dS vs P. sinensis")
boxplot(sind[,'dS']~sind$species, outline=F, ylab="pairwise dS vs sinensis", main='vs sinensis')


bp1 = ggplot(ddat) + 
	geom_boxplot(outlier.shape = 26, aes(x=species, y=dS)) +
	scale_y_continuous(limits=c(0,0.25)) +
	labs(x='')
plot(bp1)

bp2 = ggplot(ddats) + 
	geom_boxplot(outlier.shape = 26, aes(x=species, y=dS)) +
	scale_y_continuous(limits=c(0,0.08)) +
	labs(x='')
plot(bp2)

bp3 = ggplot(sind) + 
	geom_boxplot(outlier.shape = 26, aes(x=species, y=dS)) +
	scale_y_continuous(limits=c(0,0.3)) +
	labs(x='')
plot(bp3)


plot_grid(bp1, bp2, bp3, nrow=1, labels=c("A","B","C"))

#stats
table(na.omit(ddat)$species)
table(na.omit(ddats)$species)
table(na.omit(sind)$species)
tapply(sind$dS, INDEX=sind$species, function(x) median(x, na.rm=T))
tapply(sind$dS, INDEX=sind$species, function(x) mean(x, na.rm=T))
tapply(ddats$dS, INDEX= ddats$species, function(x) mean(x, na.rm=T))
x=sind$dS[sind$species=='Y']
y=sind$dS[sind$species=='tym']
t.test(x,y, alternative = 'less')
wilcox.test(x=x, y=y)
wilcox.test(x=ddat$dS[ddat$species=='sin'], y=ddat$dS[ddat$species=='pun'])
