#plot_fig3B_chrom_df_ABBA_BABA_stats.R

library(ggplot2)
library(cowplot)
setwd("~/gitreps/pungitius_sex_chromosome/introgression/abbababba")

#load chromosome names
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
cdat
chroms = as.character(cdat$chroms)


#UPLOAD DATA

readin = function(fileName){
	df = read.table(fileName, header = T)
	df$mb = (df$lefts + df$rights) / 2 / 1e6
	df = df[df$mb > left & df$mb < right,]
	return(df)
}

left = 4
right = 17

fms = readin("all_FMS.tsv"); xlim=c(-0.48, 0.48)
fmt = readin("all_FMT.tsv")
tsp = readin("all_TSP.tsv"); xlim=c(-0.4, 0.4)



#SELECT WHICH DATASET TO PLOT
df = fms
df = fmt
df = tsp

#format
selection = c("chrXII")
df.auto = df[!df$chr %in% selection,]
sex = df[df$chr %in% selection,]

#setup density df
den = density(df.auto[,'fd'], na.rm=T) #get densities for chromosomes only
den.df = data.frame(x=den$x, y=den$y)
den.sex = density(sex[,'fd'], na.rm=T)
sex.df = data.frame(x= den.sex$x, y= den.sex$y)

#plot
g1=ggplot(data= sex.df, aes(x=x, y=y)) + 
	geom_area(data=sex.df, fill='black', alpha=0.75) +
	# geom_segment(data=mn.sex, aes(x=x, y=y1, xend=x, yend=y2), col='black', lty=2) +
	# geom_segment(data=mn.auto, aes(x=x, y=y1, xend=x, yend=y2), col='grey', lty=2) +
	geom_area(data=den.df, fill='grey', alpha=0.75) +
	lims(x=xlim) +
	labs(x="D", y="")
plot(g1)



