#!/usr/bin/env Rscript
#concencus_vcf.R
#script takes a vcf and returns a new vcf with a single 
#individual with genotypes called as the modal genotype from the original

#parse arguments
args = commandArgs(trailingOnly=TRUE)
infileName = args[1]
concensusName = args[2]
outName = args[3]

print("Reading in infile...")
vdat = read.table(infileName)
print("Header:")
print(head(vdat))
ndat = vdat[,1:9]
colnames(ndat) = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
gdat = vdat[,10:ncol(vdat)]


#funciton for Mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

c = apply(gdat, 1, function(x) Mode(x))
ndat[, concensusName] = c


print("Writing out results")
print("Outfile:")
print(outName)
write.table(ndat, file=outName, row.names=F, quote=F, sep = "\t")