#xy_white_deseq.R
#use read counts from snpSplit to look at X:Y differences in transcription


library('DESeq2')
library('ggplot2')
library('cowplot')
setwd("~/gitreps/pungitius_sex_chromosome/gene_expression/allele_specific_expression")

#load data
raw.counts = read.table('readCounts_white_SNPsplit.tsv', header = T, row.names='geneID', sep="\t")

#revise mismatched and remove outlier
dim(raw.counts)
to.remove = colnames(raw.counts)[grep("DRR023137", colnames(raw.counts))]
raw.counts = raw.counts[,!colnames(raw.counts) %in% to.remove]
dim(raw.counts)

cn0 = colnames(raw.counts)
cn = sub("DRR023134_F", "DRR023134_M", cn0)
colnames(raw.counts) = cn



#isolate the X and Y reads for males
drr = colnames(raw.counts)[grep("DRR023", colnames(raw.counts))]
males = drr[grep("_M.genome", drr)]
counts0=raw.counts[,colnames(raw.counts) %in% males]
colnames(counts0) = sub(".counts.txt", "", colnames(counts0))
st = counts0[,grep("_st", colnames(counts0))]
non = counts0[,-grep("_st", colnames(counts0))]
head(st)
head(non)
sum(sub("_st", "", colnames(st)) == colnames(non)) == ncol(non)
colnames(st) = sub("_st", "", colnames(st))
counts = st + non
head(counts)
dim(counts)
#double-check that worked
sum(counts$DRR023131_M.genome2 == counts0$DRR023131_M.genome2 + counts0$DRR023131_M.genome2_st) == nrow(counts0)


#remove non-gene objects from counts data
to.remove=rownames(counts)[grep("__", rownames(counts))]
counts[to.remove,]
dim(counts)
counts= counts[!rownames(counts) %in% to.remove,]
dim(counts)

#remove genes with low coverage
cut=1
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]
dim(counts)



#add project to coldata
sample = colnames(counts)
chr = sample
chr[grep('genome1', chr)]<-'X'
chr[grep('genome2', chr)]<-'Y'
coldata = data.frame(row.names=sample, chrom=chr)
rownames(coldata)==colnames(counts)


#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~chrom))

#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds, contrast = c('chrom', 'Y', 'X'), independentFiltering=T)
summary(res)
head(res)



#check X-linked
gdat = read.table("gene_locations.txt", header = T)
rownames(gdat) = gdat$geneId
chr12genes = gdat[gdat$chr=="chrXII",'geneId']
r12 = rownames(res)[rownames(res) %in% chr12genes]
length(r12) / nrow(res)



#volcano plot
res.df = data.frame(res)
res.df$sig = res.df$padj < 0.1
res.df$sig[is.na(res.df$sig)]<-FALSE
res.df$logp = -log(res.df$pvalue, 10)
res.df$chr12 = rownames(res.df) %in% chr12genes
res.df=res.df[rownames(res.df) %in% chr12genes,]
head(res.df)
dim(res.df)

g = ggplot(data=res.df, aes(x=log2FoldChange, y=logp)) + 
	scale_colour_manual(values=c('black','red')) +
	geom_point(data=res.df, aes(color = sig), alpha=0.4) +
	xlab("log2 fold difference") +
	ylab("-log10 p-value") +
	theme_bw()
plot(g)

plot(density(res.df$log2FoldChange), main='', xlab='log2 Y:X')


#compare with basic alignment results
ll=load("~/gitreps/pungitius_sex_chromosome/figure_plotting/maleVfemalePungitius.Rdata")
head(s.r)
s.df = data.frame(s.r)
plot(density(s.df$log2FoldChange, na.rm=T), main='', xlab='log2 Male:Female')


m=merge(res.df, s.df, by = 0)
plot(m$log2FoldChange.y~m$log2FoldChange.x, xlab="log2 Y:X in males", ylab="log2 Male:Female")
abline(h=0)
abline(v=0)
lm1=lm(m$log2FoldChange.y~m$log2FoldChange.x)
abline(lm1, col='red')
summary(lm1)


#plot relationship
g = ggplot(data=m) +  
	labs(x="Y:X within males", y="Male:Female") +
	geom_smooth(data=m, aes(log2FoldChange.x,log2FoldChange.y), method='lm', formula=y~x, se=T) +
	geom_point(aes(x=log2FoldChange.x, y= log2FoldChange.y), alpha=1)
plot(g)

#save objects
m.w = m
sr.w = s.r
xy.w = res.df
save(m.w, sr.w, xy.w, file='~/gitreps/pungitius_sex_chromosome/figure_plotting/xy_white.Rdata')


