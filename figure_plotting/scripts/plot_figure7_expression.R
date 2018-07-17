#plot_figure7_expression.R
#This is a cleaned version of plot_figure7_exploratory.R
#That one has a lot more in case these break or aren't clear
setwd("~/gitreps/pungitius_sex_chromosome/figure_plotting/")
ll=load("fig7_files.Rdata")


#---------- PLOT ALL TOGETHER ----------#
#Figure 7:
quartz()
plot_grid(absBox, labels='')
quartz()
plot_grid(xyscatter.w2, xyscatter.vh, ncol=2, label_size =16, axis='b')

#supplemental showing its not degeneration
plot(denRepMG)
plot(denRepFM)
plot_grid(denRepFM, denRepMG, labels=c('A', 'B'), ncol=1, label_size =16)



#---------- STATS ----------#

#STATS FOR BOXPLOTS
do.wilcox = function(sp1, sp2, df, col, stat){
	stat1 = df[df[,col] == sp1, stat]
	stat2 = df[df[,col] == sp2, stat]
	sub = df[df[,col] %in% c(sp1, sp2),]
	w=wilcox.test(x=stat1, y=stat2)
	boxplot(sub[,stat]~sub[,col], outline=F, main=paste("p =",w$p.value))
	return(w)
}
do.wilcox('SDR_brain-dna', 'auto_brain-dna', ssr, 'sexdf', 'absDiff')
do.wilcox('SDR_pelvic-dna', 'auto_pelvic-dna', ssr, 'sexdf', 'absDiff')
do.wilcox('SDR_acu', 'auto_acu', ssr, 'sexdf', 'absDiff')


#STATS FOR SCATTERPLOTS

#for dna controlled brain scatterplot
plot(m$fc~m$fcYX)
lmw=lm(m$fc~m$fcYX)
abline(lmw, col='red')
summary(lmw)


#for von Hippel scatterplot
plot(m.vh$log2FoldChange.y~m.vh$log2FoldChange.x)
lmvh=lm(m.vh$log2FoldChange.y~m.vh$log2FoldChange.x)
abline(lmvh, col='red')
summary(lmvh)

