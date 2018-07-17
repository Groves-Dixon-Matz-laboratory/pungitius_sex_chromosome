#bedtools_depth.R
#This script compares fold coverages from bedtools


#set global variables for script
setwd("~/gitreps/pungitius_sex_chromosome/male_female_comparisons/dna_100Kb_window_results")
source("~/gitreps/pungitius_sex_chromosome/pun_sex_chrom_functions.R")

###############################
####### WINDOW DEPTH ##########
###############################
library('DESeq2')

#pick one of the three Pungitius species
SPP="Pun"; spp='pun'
SPP="Tym"; spp='tym'
SPP="Sin"; spp='sin'

#or pick Pacific Ocean (3-spine)
SPP="po"; spp="po"




#upload and format data
ddat = read.table("all_depth.tsv", header = T, stringsAsFactors=F)
windows = paste(paste(ddat$chr,ddat$start,sep="_"), ddat$end, sep="_")
wdat = ddat[,1:3]
chroms=wdat$chr
mids=apply(wdat[,2:3], 1, mean)
chrList = unique(chroms)
chrList = chrList[!chrList %in% c('chrM', 'chrUn')]
chrNum = c(1, 2, 3, 4, 9, 5, 6, 7, 8, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 20, 21)
chr.df = data.frame(chrList, chrNum)
chr.df=chr.df[order(chr.df$chrNum),]
chrList=as.character(chr.df$chrList)



#now subset the counts 
counts=ddat[, grep(SPP, colnames(ddat))]
colnames(counts) = sub("_chrI.bam", "", colnames(counts))
head(counts)
dim(counts)




#UNCOMMENT CHUNK BELOW FOR PACIFIC OCEAN 3-SPINE
# #now subset the counts 
# t = read.table("~/gitreps/pungitius_sex_chromosome/metadata/po_sex.txt")
# colnames(t) = c('sample', 'sex')
# t=t[grep('DRS', t$sample),]
# counts=ddat[,colnames(ddat) %in% t$sample]
# dim(counts)
# coldata=t



#set up trait data (double-check you have all samples still)
t = read.table("~/gitreps/pungitius_sex_chromosome/metadata/multispecies_fish_info.tsv", header = T, sep="\t", row.names='sample')
coldata = t[colnames(counts),]
head(coldata)
dim(coldata)

#double-check names align
sum(rownames(coldata) == colnames(counts)) == ncol(counts)


#-------------- FPM WAY ----------------#
m = apply(counts, 2, sum) / 1e6
fpm = sweep(counts, MARGIN=2, FUN='/', STATS=m)
head(fpm)
male.samples = rownames(coldata)[coldata$sex=='M']
female.samples = rownames(coldata)[coldata$sex=='F']
m.fpm = fpm[,male.samples]
f.fpm = fpm[,female.samples]
m.mn = apply(m.fpm, 1, mean)
f.mn = apply(f.fpm, 1, mean)
res=data.frame(chr=ddat$chr, start=ddat$start, end=ddat$end, m.fpm=m.mn, f.fpm=f.mn)
res$ratio=res$m.fpm/res$f.fpm
res$logmf = log(res$ratio, 2)

select = res[res$chr=='chrXII',]
# # select = res[res$chr=='chrXIX',] #use for 3-spine
plot(select$logmf~ select$start, ylim=c(-1,1))


#-------------- DESEQ WAY ----------------#
#run DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~sex))
	
dds = DESeq(ddsHTSeq, fitType='mean')
resultsNames(dds)
res=results(dds, contrast = c('sex', 'M', 'F'))
head(res)
rld=rlog(ddsHTSeq, fitType='mean')
colnames(rld) = colnames(counts)
rld.df = assay(rld)


#organize the results
res.df=cbind(wdat, data.frame(res))


#plot log2 fold diff female vs male
par(mfrow=c(5,5))
for (chr in chrList){
	cres = res.df[res.df$chr==chr,]
	plot(cres$log2FoldChange~cres$start, cex=0.5, main=chr, ylab = "M:F Log2 Ratio", xlab="Position (kb)")
	loess_fit <- loess(cres$log2FoldChange~cres$start, span = 0.2)
	lines(cres$start, predict(loess_fit), col='red', lwd=2)
	abline(h=0,lty=2, col='grey')
}


#write out the depth results
head(res.df)
outName=paste(spp, 'MvsF_foldDiff.Rdata', sep="_")
save(res.df, file=outName)

### optionally save the size factors for another analysis
# #Save for pungitius read counts against the whole genome
# #This will be used later for running DESEq on reads mapped to repetitive elements
# genomeSizeFactors = sizeFactors(dds)
# save(genomeSizeFactors, file="~/gitreps/pungitius_sex_chromosome/metadata/genomeSizeFactors.Rdata") ## saved for other analyses


