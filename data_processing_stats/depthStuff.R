

setwd("~/gitreps/pungitius_sex_chromosome/data_processing_stats")
d=read.csv("depth_data.csv")
d=d[,!grepl("X", colnames(d))]
head(d)

#get estimate for genome length
cl = read.table("chromLengths.txt")
colnames(cl) = c('chr', 'len')
totl = sum(cl$len)



#get fold coverage
readLen = 150
pe = 2
d$cov = (d$raw * readLen * pe) / totl

#get overall mean coverage
mean(d$cov)

d$spp = substr(d$sample.1, start=1, stop=1)
sppMeans = tapply(d$cov, INDEX=d$spp, mean)