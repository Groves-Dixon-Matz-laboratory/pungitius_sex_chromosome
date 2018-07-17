#compare_pairwise_dnds.R
#plot and do stats on pairwise dN/dS values from PAML


setwd("~/gitreps/pungitius_sex_chromosome/degeneration/dnds")


#function to do signed rank test between two chromosomes
do.wilcox = function(sp1, sp2, df, stat){
	stat1 = df[df$species == sp1, stat]
	stat2 = df[df$species == sp2, stat]
	sub = df[df$species %in% c(sp1, sp2),]
	w=wilcox.test(x=stat1, y=stat2)
	boxplot(sub[,stat]~sub[,'species'], outline=F, main=paste("p =",w$p.value))
	return(w)
}


#compare for pairwise against threespine
d=read.table('ref_chrXII_pairwise_dNdS.tsv', header = T, stringsAsFactors=F)

#plot for all mutations
ad = d[!grepl('private', d$species),]
boxplot(ad$dN~ad$species, outline=F, ylab='dN')
boxplot(ad$dS~ad$species, outline=F, ylab='dS')
boxplot(ad$dNdS~ad$species, outline=F, ylab="dN/dS")
do.wilcox('X', 'tym', ad, 'dNdS')


#plot for "private" mutations
pd = d[grepl('private', d$species),]
boxplot(pd$dN~pd$species, outline=F, ylab='dN')
boxplot(pd$dS~pd$species, outline=F, ylab='dS')
boxplot(pd$dNdS~pd$species, outline=F, ylab='dN/dS')
do.wilcox('Xprivate', 'tymprivate', pd, 'dNdS')
do.wilcox('Xprivate', 'sinprivate', pd, 'dNdS')
do.wilcox('Yprivate', 'sinprivate', pd, 'dNdS')
do.wilcox('Xprivate', 'Yprivate', pd, 'dNdS')



#look at X
xd=read.table('X_pairwise_dNdS.tsv', header = T, stringsAsFactors=F)
xd=sind[!grepl('private', xd$species),]
boxplot(xd[,'dN']~xd$species, outline=F, ylab="pairwise dN vs sinensis", main='vs sinensis')
boxplot(xd[,'dS']~xd$species, outline=F, ylab="pairwise dS vs sinensis", main='vs sinensis')
boxplot(xd[,'dNdS']~xd$species, outline=F, ylab="pairwise dN/dS vs sinensis", main='vs sinensis')
do.wilcox('X', 'Y', d, 'dS')
do.wilcox('tym', 'Y', d, 'dS')


#look at sinensis
sind=read.table('sin_pairwise_dNdS.tsv', header = T, stringsAsFactors=F)
sind=sind[!grepl('private', sind$species),]
boxplot(sind[,'dN']~sind$species, outline=F, ylab="pairwise dN vs sinensis", main='vs sinensis')
boxplot(sind[,'dS']~sind$species, outline=F, ylab="pairwise dS vs sinensis", main='vs sinensis')
boxplot(sind[,'dNdS']~sind$species, outline=F, ylab="pairwise dN/dS vs sinensis", main='vs sinensis')
do.wilcox('X', 'Y', d, 'dS')
do.wilcox('tym', 'Y', d, 'dS')



#save objects for plotting with the rest of the degeneration figures
save(ad, pd, sind, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/pairwisedNdS.Rdata")


