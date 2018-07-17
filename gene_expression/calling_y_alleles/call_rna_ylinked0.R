#get_putative_fixed0.R
#Use this script to estimate false positive rates for Y-assignment criteria
#Works for either RNAseq datasets
#delete the # in front of main VCF header line before importing


setwd("/Users/grovesdixon/lab_files/projects/rna_sex/call_y_linked")


dat = read.table("all_whitePun.vcf", header = T )
dat = read.table("all_endocrine.vcf", header = T)
head(dat)

#double-check the weird samples
should.be.male = colnames(dat)[grep("DRR023134", colnames(dat))]
should.be.male
print(paste("Male switched correctly =", grepl("_M", should.be.male)))
dat=dat[,colnames(dat) != 'DRR023137_M']



gdat = dat[,10:ncol(dat)]
head(gdat)

#grab genotypes and separate by sex
g = apply(gdat, c(1,2), function(x) substr(x, start=1, stop=3))
g[g=="./."]<-NA
m = g[,grep("_M", colnames(g))]
f = g[,grep("_F", colnames(g))]

is.full.het = function(x, threshold){
	l = length(na.omit(x))
	t.het = table(x)["0/1"]
	if(is.na(t.het)){
		return(FALSE)
	}
	else{
		return(t.het>=l*threshold)
	}
}

is.fixed = function(x){
	l = length(na.omit(x))
	c = table(x)
	fix.check = c(c["0/0"]==l, c["1/1"]==l)
	fix.check[is.na(fix.check)]<-FALSE
	fixed = sum(fix.check)==TRUE
	return(fixed)
}


is.ref = function(x){
	xstr = na.omit(unique(as.character(x)))[1]
	return(xstr=='0/0')
}


tset = seq(0.0, 1, 0.2)
fps = c()
totals = c()
for (threshold in tset){
	#select putative fixed
	print(paste(threshold, "...", sep=""))
	m.het = apply(m, 1, function(x) is.full.het(x, threshold))
	f.fixed = apply(f, 1, function(x) is.fixed(x))
	put.y = m.het & f.fixed
	f.fixed.ref = apply(f, 1, function(x) is.ref(x))
	sum(put.y)
	length(m.het)
	length(f.fixed)
	
	
	#get proportions on chromosomes
	ccount = table(dat$CHROM)
	pydat = dat[put.y,]
	scount = table(pydat$CHROM)
	scount
	false = 1-scount['chrXII']/sum(scount)
	fps = append(fps, false)
	totals = append(totals, scount['chrXII'])
}
x=tset*100
plot(fps~x, ylab="False Positive", xlab="Heterozygosity  Cutoff")
# abline(v=80, lty=2, col='red')
plot(totals~x, ylab="N loci", xlab="Heterozygosity Cutoff")
# abline(v=80, lty=2, col='red')
real = totals - fps
ratio = totals / fps
plot(fps~totals, ylab="False Positive", xlab="total loci")



## CHOSE 1 AS CUTOFF FOR WHITE DATASET (same as 0.8)
## CHOSE .8 AS CUTOFF FOR ENDOCRINE