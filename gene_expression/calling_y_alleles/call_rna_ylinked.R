#get_putative_fixed.R
#used varaints called from RNAseq data to build a 
#revised VCF for N-masking a refence with SNPsplit

#putative X/Y linked follow these criteria:
#1: biallelic
#2: heterozygous in all genotyped males
#3: homozygous in all genotyped females

setwd("~/gitreps/stickle_back_sex_chromosomes/call_y_linked_rna")

#pick which dataset to run for
dat = read.table("all_whitePun.vcf", header = T, ); dataset="whitePun"   #White 2014
dat = read.table("all_endocrine.vcf", header = T, ); dataset="endocrine" #von Hippel 2018
head(dat)

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

## CHOSE 4 AS CUTOFF (CAN HAVE 1 SAMPLE THAT MISSED ITS HETEROZYGOSITY)


#select putative fixed
threshold=0.79   #this is based off of figures built in get_putative_fixed0.R
m.het = apply(m, 1, function(x) is.full.het(x, threshold))
print(sum(m.het))
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
barplot(scount)
barplot(scount[names(scount)!='chrXII'])
false = 1-scount['chrXII']/sum(scount)


#build the SNPs_file that SNPsplit needs
#recording the Y-linked allele as the alternative
dat$f.fixed.ref = f.fixed.ref
arbitrary.male = colnames(dat)[grep("_M", colnames(dat))][1]
sdat = dat[put.y, c('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', arbitrary.male, 'f.fixed.ref')]
switch = !sdat$f.fixed.ref


#write out the sex-linked snps with the reference as X and alternative as Y
#*note this is different from how the VCF is output below
snpId = 1:nrow(sdat)
chr = sdat$CHROM
pos = sdat$POS
strand = 1
refSnp = paste(sdat$REF, sdat$ALT, sep='/')
refSnp[switch]<-paste(sdat$ALT, sdat$REF, sep='/')[switch]
fres = data.frame(snpId, chr, pos, strand, refSnp)


#double-check
chrPos = paste(dat$CHROM, dat$POS, sep="_")
sdat$chrPos = paste(fres$chr, fres$pos, sep="_")
refDat = sdat[sdat$f.fixed.ref,]
altDat = sdat[sdat$f.fixed.ref==FALSE,]
c=cbind(f, m)
femaleRef = c[chrPos %in% refDat$chrPos,]
femaleAlt = c[chrPos %in% altDat$chrPos,]
x=dat[chrPos %in% refDat$chrPos,]
y=dat[chrPos %in% altDat$chrPos,]
z=dat[chrPos %in% sdat$chrPos,]
sum(c(nrow(x), nrow(y)))
#can look at these manually to insure that the altDat has females fixed for the alternative allele, 
#as they should be, and that 
write.table(y, file="~/junk/altDat.tsv", quote=F, row.names=F)
write.table(z, file="~/junk/subDat.tsv", quote=F, row.names=F)



#write out sex-linked SNPs
colnames(fres) = c('SNP-ID', 'Chromosome', 'Position', 'Strand', 'Ref/SNP')
head(fres)
sexLinkedSnpFile=paste(dataset, "newSex_linked_snps.tsv",  sep="_")
write.table(fres, file=sexLinkedSnpFile, row.names = F, quote=F, sep = "\t")



##Write out a fake VCF to prepare genome with
header = c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'N_MASKING_DUMMY_SAMPLE')
id = rep(".", length(chr))
ref = sdat$REF
alt = sdat$ALT
qual = sdat$QUAL
filter = sdat$FILTER
info = sdat$INFO
format = paste(sdat$FORMAT, 'FI', sep=":")

#add the vcf header text back
format=rep('GT:FI', length(alt))
format=rep("GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI", length(alt))
gt = rep("1/1", length(alt)) ##set them all to the alternative allele, so that all are N-masked (don't care whether on X or Y)
fi = rep("1", length(gt))
fake.dat = rep("22:6:0.166667:152,22,0:137,18,0:2:36:6:0,0,6,0:0:-0.616816:.:1", length(gt))
geno = paste(gt, fake.dat,sep=":")



vdat = data.frame(chr, pos, id, ref, alt, qual, filter, info, format, geno)
colnames(vdat) = header
head(vdat)

vcf.header = '##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chrI,length=29641652>
##contig=<ID=chrII,length=23708229>
##contig=<ID=chrIII,length=17807731>
##contig=<ID=chrIV,length=34153411>
##contig=<ID=chrIX,length=20593793>
##contig=<ID=chrM,length=15742>
##contig=<ID=chrUn,length=26725976>
##contig=<ID=chrV,length=15563594>
##contig=<ID=chrVI,length=18854982>
##contig=<ID=chrVII,length=30850397>
##contig=<ID=chrVIII,length=20538625>
##contig=<ID=chrX,length=18035287>
##contig=<ID=chrXI,length=17646579>
##contig=<ID=chrXII,length=20772122>
##contig=<ID=chrXIII,length=20752762>
##contig=<ID=chrXIV,length=16171761>
##contig=<ID=chrXIX,length=20612724>
##contig=<ID=chrXV,length=17323106>
##contig=<ID=chrXVI,length=19522594>
##contig=<ID=chrXVII,length=20254136>
##contig=<ID=chrXVIII,length=15990693>
##contig=<ID=chrXX,length=20460780>
##contig=<ID=chrXXI,length=17357772>
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="Whether a sample was a Pass(1) or fail (0) based on FILTER values">
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

outPath = paste(dataset, "vcf_for_genome_prep.vcf", sep="_")
write(vcf.header, file= outPath, ncolumns=1)
write(colnames(vdat), file= outPath, ncolumns=ncol(vdat), append=T, sep = "\t")
write(t(as.matrix(vdat)), file= outPath, ncolumns=ncol(vdat), append=T, sep = "\t")


#------- done outputting putative Y-linked SNPs -------#

#check overlap between the two RNAseq datasets
edat = read.table("endocrine_newSex_linked_snps.tsv", header = T)
wdat = read.table("whitePun_newSex_linked_snps.tsv", header = T)
edat$chr_pos = paste(edat$Chromosome, edat$Position)
wdat$chr_pos = paste(wdat$Chromosome, wdat$Position)
rownames(edat) = edat$chr_pos
rownames(wdat) = wdat$chr_pos
head(edat)
overlap = edat$chr_pos[edat$chr_pos %in% wdat$chr_pos]
e=edat[overlap,]
w=wdat[overlap,]
nrow(w)/nrow(wdat)
nrow(e)/nrow(edat)
sum(e$Ref.SNP==w$Ref.SNP) / nrow(w)


#check the putative fixed vars against what we see in dna
ff = read.table("datasets/femalePunChr12_Rin.freq", sep = "\t", header = T)
mf = read.table("datasets/malePunChr12_Rin.freq", sep = "\t", header = T)


#remove NA positions
naf = is.na(ff$f1)
nam = is.na(mf$f1)
remove = naf | nam

ff = ff[!remove,]
mf = mf[!remove,]
dim(ff)
dim(mf)

#check they always sum to 1
t=ff$f1 + ff$f2
table(t)

#get fixed
CUT=1
ffr = ff$f1>=CUT
ffr[is.na(ffr)]<-TRUE
ffa = ff$f2>=CUT
ffa[is.na(ffa)]<-TRUE

ff$dna.f.fixed.ref = TRUE
ff$dna.f.fixed.ref[ffa]<-FALSE


#get male freq
yMin = 0.3   #minimum frequency an allele should be at in males to potentially be Y-linked
yMax = 0.7   #max frequency an allele should be at in males to potentially be Y-linked

yr = mf$f1 > yMin & mf$f1 < yMax
yr[is.na(yr)]<-FALSE
ya = mf$f2 > yMin & mf$f2 < yMax
ya[is.na(ya)]<-FALSE



#putative Y-reference alleles
pyR = ffa & ya
length(pyR)
sum(pyR)

#putative Y-alternative alleles
pyA = ffr & yr
length(pyA)
sum(pyA)

#set up full dataframe to look at both together
pymdat = mf[pyR | pyA,]
pyfdat = ff[pyR | pyA,]
pydat = pyfdat
pydat$male.f1 = pymdat$f1
pydat$male.f2 = pymdat$f2
head(pydat, n=100)

#check overlap
sub = sdat[sdat$CHROM=='chrXII',]
pos = sub$POS
overlap = pos[pos %in% pydat$POS]
length(overlap)
length(overlap) / length(pos)


#merge up 
m=merge(sub, pydat, by = 'POS')
fagree = m$f.fixed.ref == m$dna.f.fixed.ref



#subset for the genic dna variants
gtf = read.table("datasets/chrom12.gtf", sep = "\t")
g = gtf[gtf$V3=='CDS',c('V4', 'V5')]
colnames(g) = c('left', 'right')


check.in.exon = function(x){
	return(sum(x>= g$left & x <= g$right) > 0)
}
pydat$in.exon = sapply(pydat$POS, function(x) check.in.exon(x))
sum(pydat$in.exon)
sum(pydat$in.exon) / nrow(pydat)
genic.pos = pydat$POS[pydat$in.exon]
length(genic.pos)



#summarize
nrow(sub) #total X-linked loci called with RNAseq
nrow(m)   #total that had matched calls in DNA data
nrow(m) / nrow(sub)  #proportion of RNA Y-linked calls found in DNA data
nrow(m) / length(genic.pos)
sum(fagree) / nrow(m)  #do they always agree on which allele is X-linked ? -- yes


#output high-confidence Y-alleles
hi.conf = data.frame(m$CHROM.x, m$POS)
write.table(hi.conf, file="datasets/high_confidence_Y_loci.tsv", quote=F, row.names = F, sep="\t", col.names = F)

#------------------------------------------
#at this point, go back to TACC and output the DNA data for these hig-confidence Y-linked loci



#-----------------------------------------
setwd("/Users/grovesdixon/lab_files/projects/rna_sex/")
depth = read.table("results/Y_calling/high_confidenceY_male.idepth", header = T)
het = read.table("results/Y_calling/high_confidenceY_male.het", header = T)
c = cbind(depth, het)
c$het = 1-c$O.HOM. / c$N_SITES
plot(c$het~c$MEAN_DEPTH, xlab = 'Mean Depth', ylab = "Male Ho", main = "High-Confidence Y loci")


