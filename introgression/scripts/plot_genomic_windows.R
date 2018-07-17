#plot_genomic_windows.R
library(ggplot2)
library(cowplot)
source("~/gitreps/pungitius_sex_chromosome//pun_sex_chrom_functions.R")

#---- ABBA BABA ----#
setwd("~/gitreps/pungitius_sex_chromosome/introgression/abbababba")


#load chromosome names
ll=load("~/gitreps/pungitius_sex_chromosome/metadata/chromList.Rdata")
chroms = as.character(cdat$chroms)

#pick boundaries for SDR
left = 0
right = 1e10
span = 0.2

#upload the fd data
#funciton to read in and format the fd data


#
tsp = readin_fd("all_TSP.tsv")
fms = readin_fd("all_FMS.tsv")



#select the dataset
df=tsp
df=fms


#build plots
YLIM=c(0,0.25)
clist = list()
counter = 0
for (chr in chroms){
	counter = counter + 1
	main = paste("Chr", counter, sep='')
	main=counter
	csub = na.omit(df[df$chr==chr,])
	lastTick = round(max(csub$mb))
	cplot = ggplot(data=csub) + 
		geom_smooth(aes(x=mb, y=fd), span = span, se=F) +
		lims(y=c(0,0.25)) +
		labs(x='', title=main) +
		scale_x_continuous(breaks = c(0, floor(lastTick/2), lastTick)) +
		theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
    if (chr=='chrXII'){
    	cplot = cplot + 
    	geom_vline(xintercept= 3.8, lty=2, lwd=1, color='green') +
    	geom_vline(xintercept= 18.9, lty=2, lwd=1, color='green') +
    	geom_smooth(aes(x=mb, y=fd), span = span, se=F)
    }
	clist[[chr]]<-cplot
}
length(clist)

#plot
plot_grid(plotlist = clist, ncol=length(clist))


#---- FST AND DTILDE ----#
setwd("~/gitreps/pungitius_sex_chromosome/introgression/distance_stats")
to.remove = c("chrUn", "chrM", "PAR", "SDR"); selection = c("chrXII") #just subset from 12 instead of reading sdr and par
leftBound = 0
rightBound = 1e10
ys = readinDistance("all_y_sin_distances.tsv", 'ys')


#build plots
flist = list()
dlist=list()
counter = 0
for (chr in chroms){
	counter = counter + 1
	main = paste("Chr", counter)
	main = ""
	csub = na.omit(ys[ys$chr==chr,])
	lastTick = round(max(csub$mb))
	
	
	#plot fst
	fstplot = ggplot(data=csub) + 
		geom_smooth(aes(x=mb, y=fst), span = span, se=F) +
		lims(y=c(0,0.32)) +
		labs(x='', title=main) +
		scale_x_continuous(breaks = c(0, floor(lastTick/2), lastTick)) +
		theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        # axis.ticks.x=element_blank()
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
          
    #plot dtilde
    dplot = ggplot(data=csub) + 
		geom_smooth(aes(x=mb, y=dtilde), span = span, se=F) +
		lims(y=c(0,0.012)) +
		labs(x='', title=main) +
		scale_x_continuous(breaks = c(0, floor(lastTick/2), lastTick), labels=c('0', '', lastTick)) +
		theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        # axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
    
    if (chr=='chrXII'){
    	fstplot = fstplot + 
    	geom_vline(xintercept= 3.8, lty=2, lwd=1, color='green') +
    	geom_vline(xintercept= 18.9, lty=2, lwd=1, color='green') +
    	geom_smooth(aes(x=mb, y=fst), span = span, se=F)
    	
    	dplot = dplot + 
    	geom_vline(xintercept= 3.8, lty=2, lwd=1, color='green') +
    	geom_vline(xintercept= 18.9, lty=2, lwd=1, color='green') +
    	geom_smooth(aes(x=mb, y=dtilde), span = span, se=F)
    }

    flist[[chr]]<-fstplot
	dlist[[chr]]<-dplot
}
length(dlist)
length(flist)

#plot individually
# plot_grid(plotlist = flist, nrow=1)
# plot_grid(plotlist = dlist, nrow=1)


#plot all together
full = append(append(clist, flist), dlist)
plot_grid(plotlist=full, nrow=3)


