#!/usr/bin/env Rscript
#rewrite_ylab_vcf.R
#Groves Dixon
#2-19-18


#parse arguments
args = commandArgs(trailingOnly=TRUE)
vcfIn = args[1]
yIn = args[2]
yOut = args[3]



#READ IN DATASETS
#vcf
print(paste("Reading in VCF file", vcfIn))
vdat = read.table(vcfIn, header =T)
print("VCF header:")
head(vdat)

#Ycalls 
print("Reading in Y-assignment data")
ydat = read.table(yIn, header = T)
print("Y-assignment header:")
print(head(ydat))

#subset the VCF for the male and female chromosomes
print("Separating VCF into males and females based on samples present in Y-assignment file")
males.s = unique(as.character(ydat$s))
print("Male sample Names:")
print(males.s)
males.a = paste(males.s, "A", sep = "_")
males.b = paste(males.s, "B", sep = "_")
males.c = append(males.a, males.b)
print("Male chromosome labeles as they should appear in VCF:")
print(males.c)
mdat = vdat[,colnames(vdat) %in% males.c]
fdat = vdat[,!colnames(vdat) %in% males.c]
pos = vdat$POS
print("Female VCF header:")
print(head(fdat))
print("Dimentions:")
print(dim(fdat))
print("Male VCF header:")
print(head(mdat))
print("Dimentions:")
print(dim(mdat))


#assign y genotypes
yvdat = fdat[,1:9]
head(yvdat)
count=0
for (ms in males.s){
	count=count+1
	print(paste("Calling Ys for male sample ", count, sep = "#"))
	print(ms)
	ysub = ydat[ydat$s==ms,]
	print("Positions match expectation:")
	print(sum(ysub$p==pos) == length(pos))
	ycall = ysub$y
	print("Summary of Y-assignment:")
	print(table(ycall))
	score=table(ycall)[c('A', 'B')]
	if (score[1] != score[2]){
		topY = names(score[score==max(score)])[1]
		lowY = names(score[score==min(score)])[1]
		}
	else{
		topY='A'
		lowY='B'
		}
	print(paste("topY =", topY))
	print(paste("bottomY =", lowY))
	alab = paste(ms, "A", sep="_")
	blab = paste(ms, "B", sep="_")
	ylab = paste(ms, "Y", sep="_")
	xlab = paste(ms, "X", sep="_")
	mostYlab = paste(ms, topY, sep="_")
	mostXlab = paste(ms, lowY, sep="_")
	a=vdat[,alab]
	b=vdat[,blab]
	mostY = vdat[, mostYlab]
	mostX = vdat[, mostXlab]
	
	#assemble the Y genotypes
	y=a                            #assign the Y as the 'A' column. This is correct only for rows when 'A' was called as Y
	y[ycall=='B']<-b[ycall=='B']   #when B is called as Y, replace the Y genotype with the call from column 'B'
	y[ycall=='ambiguous']<-mostY[ycall=='ambiguous'] #add the best guess for the Y for the ambigiuous calls
	
	#assemble the Y genotypes
	x=b                            #assign the X as the 'B' column. This is correct only for rows when 'A' was called Y
	x[ycall=='B']<-a[ycall=='B'] #when B is called as Y, replace the X genotype with call from column 'A'
	x[ycall=='ambiguous']<-mostX[ycall=='ambiguous'] #add the best guess for the Y for the ambigiuous calls
	#Here the abmiguous calls are left as the arbitrary calls from BEAGLE
	yvdat[,ylab] = y
	yvdat[,xlab] = x
	print("Done.")
	# print("Current Header:")
	# print(head(yvdat))
}


#now add the females back to the vcf
head(fdat)
print("Appending the male and female data...")
print("Positions match:")
print(sum(fdat$POS==yvdat$POS)==nrow(fdat))
rdat = cbind(fdat, yvdat[,10:ncol(yvdat)])
colnames(rdat)[1]='#CHROM'
print("Final dataframe:")
print(head(rdat))
print("Writing out results to following file:")
print(yOut)
write.table(rdat, file=yOut, row.names=F, quote=F, sep = "\t")
