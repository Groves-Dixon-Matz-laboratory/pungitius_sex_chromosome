#plot_provean.R
#format and plot provean score variation
#save formatted data for final plotting

library(plotrix)
setwd("~/gitreps/pungitius_sex_chromosome/degeneration/provean")


#read in the data
readin = function(fileName, sampleName){
	d=read.table(fileName, header = T)
	d$sample = sampleName
	return(d)
}

x=readin('x_provean_res.tsv', 'x')
y= readin('y_provean_res.tsv', 'y')
s= readin('sin_provean_res.tsv', 'sin')
t= readin('tym_provean_res.tsv', 'tym')
d=rbind(x,y,s,t)


#plot overall scores
boxplot(d$score~d$sample, outline=F)
mns = tapply(d$score, INDEX=d$sample, function(x) mean(x, na.rm=T))
ses = tapply(d$score, INDEX=d$sample, function(x) std.error(x, na.rm=T))
plotCI(x=1:length(mns), y=mns, uiw=ses)


#frequency of bad mutations
CUT=-2.5
d$bad = as.numeric(d$score <= CUT)
d$ok = as.numeric(d$score > CUT)
sums = tapply(d$bad, INDEX=d$sample, function(x) sum(x, na.rm=T))
oksums = tapply(d$ok, INDEX=d$sample, function(x) sum(x, na.rm=T))
ratios = sums / (sums + oksums)
tbl = rbind(sums, oksums)
barplot(ratios)
chisq.test(tbl)


#SUBSET FOR PRIVATE MUTATIONS
#first accross full dataset
#these will work for sin, tym, and Y, but not X, 
#since Y is likely to have stray X variants
p = d
p$mut = paste(p$gene, p$var, sep="_")
counts = table(p$mut)
plot(density(counts)) #most of them are shared
private = counts[counts==1]
pmuts= names(private)

#to get private alleles for X, remove the Y mutations and repeat
no.y = p[p$sample !='y',]
counts.noy = table(no.y$mut)
private.x = counts.noy[counts.noy==1]
pmuts.x = names(private.x)
psub.nonx = p[p$mut %in% pmuts & p$sample != 'x',]
px = p[p$mut %in% pmuts.x & p$sample == 'x',]

#now combine the two together
psub = rbind(psub.nonx, px)


#plot overall scores
boxplot(psub$score~psub$sample, outline=F)
mns = tapply(psub$score, INDEX=psub$sample, function(x) mean(x, na.rm=T))
ses = tapply(psub$score, INDEX=psub$sample, function(x) std.error(x, na.rm=T))
plotCI(x=1:length(mns), y=mns, uiw=ses, axes=F, xlim=c(0.5,4.5))
axis(1, at = c(1:4), labels=names(mns))
axis(2)


#frequency of bad mutations
CUT=-2.5
psub$bad = as.numeric(psub$score <= CUT)
psub$ok = as.numeric(psub$score > CUT)
sums = tapply(psub$bad, INDEX=psub$sample, function(x) sum(x, na.rm=T))
oksums = tapply(psub$ok, INDEX=psub$sample, function(x) sum(x, na.rm=T))
private.ratios = sums / (sums + oksums)
private.tbl = rbind(sums, oksums)
barplot(sums)
barplot(private.ratios)
chisq.test(private.tbl)
chisq.test(private.tbl[,c('x','y')])
chisq.test(private.tbl[,c('sin','y')])

#save data for full figure plotting
save(d, psub, ratios, private.ratios, tbl, private.tbl, file="~/gitreps/pungitius_sex_chromosome/figure_plotting/provean.Rdata")


