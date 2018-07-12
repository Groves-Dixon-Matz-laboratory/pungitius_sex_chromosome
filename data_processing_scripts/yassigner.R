#!/usr/bin/env Rscript
#yassinger.R

#parse arguments
args = commandArgs(trailingOnly=TRUE)
infiles=args[1:length(args)]


print(infiles)

print("Reading in first file:")
print(infiles[1])
d = read.table(infiles[1], header = TRUE)
print(head(d))


adat = data.frame(d$AisY)
bdat = data.frame(d$BisY)
colnames(adat) = c('1')
colnames(bdat) = c('1')

print(head(adat))
print(head(bdat))


for (i in seq(2, length(infiles))) {
	infile=infiles[i]
	print("Adding data for next file:")
	print(infile)
	dnew = read.table(infile, header = T)
	adat[as.character(i)] = dnew$AisY
	bdat[as.character(i)] = dnew$BisY
}
print("Done adding data:")
print(head(adat))
print(head(bdat))


print("Counting up Y calls...")
acount = apply(adat, 1, sum)
bcount = apply(bdat, 1, sum)
print(head(acount))
print(head(bcount))

print("Making assignments based on counts...")
a=acount > bcount
b=bcount > acount

print(head(a))
print(head(b))


print("Formatting calls")
res=rep('ambiguous', length(acount))
res[a==TRUE]<-'A'
res[b==TRUE]<-'B'
namb = sum(res=='ambiguous')
pamb = namb/length(res)
pct=round((1-pamb), digits=4)*100
pctAmb = 100-pct


print("proportion ambiguous:")
print(paste(pctAmb, '%', sep=''))
print("proportion called:")
print(paste(pct, '%', sep=''))



print("Outputting results...")
cdat = data.frame(d$p, d$s, res)
colnames(cdat) = c('p', 's', 'y')
print(head(cdat))

write.table(cdat, file="concensusY.tsv", row.names=F, quote=F)
