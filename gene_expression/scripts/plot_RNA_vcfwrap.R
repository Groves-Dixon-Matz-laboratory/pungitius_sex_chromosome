#plot_RNA_vcfwrap.R
#same as script for plotting figure 1 but for the RNA data

library(plotrix)
library(ggplot2)
library(gridExtra)
library(cowplot)

#set global variables for script
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/rna_vcf_wrapper")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")
chromsToIgnore = c("chrUn", "chrM")
male.color = 'dodgerblue'
female.color = 'firebrick'


#choose focal species
spp="DRR"
spp="SRR"


#upload the VCF wrapper results
wrapFile = paste(spp, "allWindows.tsv", sep="_")
dat = read.table(wrapFile, header=T, sep="\t", na.strings = c("NA", "no_snps_in_this_window"))
dat=dat[!dat$CHROM %in% chromsToIgnore,]
dat$mb = dat$BIN_START / 1e6
head(dat)
dim(dat)
dat$chr_pos = paste(dat$CHROM, dat$BIN_START, sep="_")
dat$mfPi = log( (dat$SNP_COUNT_male/dat$SNP_COUNT_female), 2)
dat$fmPi = log( (dat$PI_female/dat$PI_male), 2)
dat$mfDens = log( (dat$SNP_COUNT_male/dat$SNP_COUNT_female), 2)
dat$fmDens = log( (dat$SNP_COUNT_female/dat$SNP_COUNT_male), 2)
dat$pctFemaleSpecific = dat$pFemaleSpecific*100
dat$pctMaleSpecific = dat$pMaleSpecific*100



#build chromosome set
chroms = unique(as.character(dat$CHROM))
chroms = chroms[!chroms %in% c("chrM", "chrUn")]
chroms = chroms[order(chroms)]
nums = c(1, 2, 3, 4, 9, 5, 6, 7, 8, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 20, 21)
cdat = data.frame(chroms, nums)
cdat=cdat[order(cdat$nums),]
print(cdat)
chromlist = as.character(cdat$chroms)



#ADD NORMALIZED X VALUES
head(dat)
dat$norm.x = dat$BIN_START
for (chr in chromlist){
	len = max(dat$BIN_START[dat$CHROM==chr])
	dat$norm.x[dat$CHROM==chr]<-dat$BIN_START[dat$CHROM==chr] / len
}
head(dat)


#--------- SAVE/LOAD
datFileName = paste(spp, 'windowData.Rdata', sep="_")
datFileName
save(dat, file= datFileName)



#-------- plot the supplemental figure for the RNA datasets --------#

setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/rna_vcf_wrapper")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")

#set up plotting vars
sex.col = 'black'
other.col = 'grey'
XCOL = 'norm.x'
XLAB=""
SPAN=0.1
SEX.CHROM = 'chrXII'
sex.specific = 'pSexSpecific'; specific.lab = '% Sex-Specific'
fst.ylim = c(-.2, 0.75)
r2.ylim = c(0, 0.8)
ss.ylim = c(0, 1)
pi.ylim = c(-2, 4)


#White Dataset
ll=load(paste('SRR', 'windowData.Rdata', sep="_"))
fst1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, nox=T)#;plot(fst1)
r21 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, nox=T)
ss1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB= '', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, nox=T)#;plot(ss1)
pi1 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, nox=F)#;plot(pi1)



#vonHippel Dataset
ll=load(paste("DRR", 'windowData.Rdata', sep="_"))
fst2 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='MEAN_FST', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=fst.ylim, nox=T, noy=F)#;plot(fst2)
r22 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='male_r2', sex.chrom=SEX.CHROM, YLAB='', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=r2.ylim, nox=T)
ss2 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol= sex.specific, sex.chrom=SEX.CHROM, YLAB= '', XLAB= XLAB,  sex.col = sex.col, other.col = other.col,SPAN=SPAN, YLIM=ss.ylim, nox=T)#;plot(ss2)
pi2 = multichrom_len_normalized(df=dat, xcol= XCOL, ycol='mfPi', sex.chrom= SEX.CHROM, YLAB='', XLAB= XLAB, sex.col = sex.col, other.col = other.col,SPAN=SPAN, , draw.line=F, tailed = 2, horiz=H, YLIM= pi.ylim, nox=F)#;plot(pi2)


plot_grid(fst1, fst2, ss1, ss2, pi1, pi2, ncol=2)



#--------------- plot multipanel scatterplots ---------------#
x=dat$BIN_START
dat$BIN_START = x/1e6
CUT = 0.99
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='pMaleSpecific', float= CUT, xlab="Position", ylab='Prop. Male-specific')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='pctMaleSpecific', float= CUT, xlab="Position", ylab='Density Male-specific')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='pFemaleSpecific', float= CUT, xlab="Position", ylab='Prop. Female-specific')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='pctFemaleSpecific', float= CUT, xlab="Position", ylab='Density Male-specific')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='mfDens', float= CUT, xlab="Position", ylab='M:F SNP Density')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='fmDens', float= CUT, xlab="Position", ylab='F:M SNP Density')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='mfPi', float= CUT, xlab="Position", ylab='M:F Pi')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='fmPi', float= CUT, xlab="Position", ylab='F:M Pi')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='MEAN_FST', float= CUT, xlab="Position", ylab='Mean Fst')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='Dxy', float= CUT, xlab="Position", ylab='Mean Dxy')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='male_r2', float= CUT, xlab="Position", ylab='Male r2')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='female_r2', float= CUT, xlab="Position", ylab='Female r2')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='m.bigMono', float= CUT, xlab="Position", ylab='Biggest Monophyletic Male')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='m.sdrBool', float= CUT, xlab="Position", ylab='Y-like toplogy')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='f.bigMono', float= CUT, xlab="Position", ylab='Biggest monophyletic female')
plot_multipanel_chroms(df0=dat, chromlist=chromlist, xcol='BIN_START', ycol='f.sdrBool', float= CUT, xlab="Position", ylab='W-like toplogy')

#------------ plot length standardized line plots ------------#
sex.col = 'black'
other.col = 'grey'
XCOL = 'norm.x'
XLAB=""
SPAN=0.1
SEX.CHROM = 'chrXII'
sex.specific = 'pSexSpecific'; specific.lab = '% Sex-Specific'



#optionally don't set y limits
fst.ylim = "none"
r2.ylim = "none"
ss.ylim = "none"
dp.ylim = "none"
pi.ylim = "none"

#assign y limits
fst.ylim = c(-.2, 0.75)
r2.ylim = c(0, 0.8)
ss.ylim = c(0, 1)
pi.ylim = c(-2, 4)








plot_grid(fst, ss, ncol=1)

#------------ plot gene tree figure ------------#
bm = multichrom_len_normalized(df=dat, xcol='BIN_START', ycol='m.sdrBool', sex.chrom= SEX.CHROM, YLAB="Y-like Topology", sex.col = sex.col, other.col = other.col)
plot(bm)

#--------------- plot single scatterplots ---------------#
#choose chr
CHR="chrXII"
CHR="chrVII"
CHR="CM002717.1"
CHR="CM002720.1"


CUT=0.99
XLAB = 'Position (Mb)'
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='pMaleSpecific', float= CUT, xlab= XLAB, ylab='Prop. Male-specific')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='pFemaleSpecific', float= CUT, xlab= XLAB, ylab='Prop. Female-specific')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='mfDens', float= CUT, xlab= XLAB, ylab='M:F SNP Density')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='fmDens', float= CUT, xlab= XLAB, ylab='F:M SNP Density')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='fmPi', float= CUT, xlab= XLAB, ylab='F:M Pi')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='MEAN_FST', float= CUT, xlab= XLAB, ylab='Mean Fst')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='male_r2', float= CUT, xlab= XLAB, ylab='Male r2')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='female_r2', float= CUT, xlab= XLAB, ylab='Female r2')
plot_sexChrom_plus_dist(df0=dat, chromlist=chromlist, chr=CHR, xcol='BIN_START', ycol='piRatio', float= CUT, xlab= XLAB, ylab='M:F Pi')


#LOOK AT PEAKS
ycol='pFemaleSpecific'
ycol='pMaleSpecific'
ycol='female_r2'
ycol='MEAN_FST'

pk = dat[dat[,ycol] > 0.15 & dat$CHROM==CHR,]
plot(dat[,ycol]~dat[,'mb'], data=dat[dat$CHROM==CHR,])
abline(v=14.95, lwd=2)


#------------- HMM for finding strata -------------#
library(depmixS4)
chr='chrXII'
xcol='BIN_START'
cdat = dat[dat$CHROM==chr,]
x=cdat[,xcol]

#set up state parameters
Nstates = 3
defaultTrans = 0
stayTrans = 5
fTrans = 1

instart = c(1, rep(0, Nstates-1))
trs0 = rep(defaultTrans, Nstates^2)
states = paste("S", 1:Nstates, sep="")
trstart = matrix(trs0, nrow=Nstates, ncol=Nstates)

for (i in 1:Nstates){
	trstart[i,i] = stayTrans}
for (i in 1:(Nstates-1)){
	trstart[i,(i+1)] = fTrans
	}
tt = t(trstart)




#Multiresponse
response = list(MEAN_FST~1, piRatio~1, pMaleSpecific~1, snpDensRat~1)
y = cdat[,c('MEAN_FST', 'piRatio', 'pMaleSpecific', 'snpDensRat')]
family=list()
for (i in 1:length(response)){
	family[[i]]<-gaussian()
}
hmm = depmix(response=response, family=family, nstates=Nstates, instart=instart, trstart=tt, data=y)
summary(hmm)
f = fit(hmm)
r=posterior(f)
res=cbind(y, x, r)
par(mfrow=c(2,2))
XLAB
plot(res$MEAN_FST ~res$x, col=res$state, pch=19, xlab = XLAB, ylab="Fst(Male vs Female)")
plot(res$piRatio~res$x, col=res$state, pch=19, xlab = XLAB, ylab = "M:F Pi")
plot(res$pMaleSpecific*100 ~res$x, col=res$state, pch=19, xlab = XLAB, ylab = "% SNPs Male-specific")
plot(res$snpDensRat ~res$x, col=res$state, pch=19, xlab = XLAB, ylab="M:F SNP Density")




#Single response
response = list(pMaleSpecific~1)
y=data.frame(pMaleSpecific=cdat[,c('pMaleSpecific')])
family=list()
for (i in 1:length(response)){
	family[[i]]<-gaussian()
}
hmm = depmix(response=response, family=family, nstates=3, instart=c(1, 0, 0), trstart=tt, data=y)




summary(hmm)
f = fit(hmm)
r=posterior(f)
res=cbind(y, x, r)
plot(res$pMaleSpecific~res$x, col=res$state, pch=19)


#Need your variables in a data frame
tXY = data.frame(var=x)
#The main depmix call, simplest version for 3 states.
response = list()
mXY =depmix(var~1,family=gaussian(),nstates=3,data=tXY)
summary(mXY)
fXY = fit(mXY)
r=posterior(fXY)
res=data.frame(r$state, x)
plot()


#A more advanced call with constraints. instart gives proabilities of starting in one of the states, trstart gives the starting transition matrix (n x n, n is your number of states), and repstart is the starting values for your gaussian parameters (mean, stdev for each state). Any 0s in the trstart are fixed rather than optimized, for some reason, but that's how you make certain states be non-transitionable from or to.
mXY = depmix(var~1,family=gaussian(),nstates=3,instart=c(1,0,0),trstart = c(0.5,0.25,0.25,0.25,0.5,0.25,0,0,1),repstart=c(4,2,0.5,1,1,1),data =tXY)

#Actually fit the model, add verbose=FALSE if you're doing lots of these.
fXY = fit(mXY)

#The results are given in the posterior(fXY). The first column is the state call, the other n are probabilities of each of the states for each observation.

#And to get the overall parameters of the optimized model
summary(fXY)

Lemme know if you have questions,







