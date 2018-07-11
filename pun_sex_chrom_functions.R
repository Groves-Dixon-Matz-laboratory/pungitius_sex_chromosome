#pun_sex_chrom_functions.R


#Functions:


#reads in a 'biggest monophyletic' results file
#see gene_treesV2.txt for how to get these
read_big_mono = function(filePath){
	gdat = read.table(filePath, header = T)
	gdat$BIN_START = gdat$BIN_START + 1
	gdat$BIN_END = gdat$BIN_END + 1
	gdat$SDRbool = as.numeric(gdat$SDRbool) - 1
	gdat$chr_pos = paste(gdat$CHROM, gdat$BIN_START, sep="_")
	return(gdat)
}


add.column = function(d1, d2, mby, col2add){
	m=merge(d1, d2, by = mby, all.x=T, sort=T)
	print("Everything Match?")
	print(sum(m[,mby] == d1[,mby]) == nrow(d1))
	if (sum(m[,mby] == d1[,mby]) == nrow(d1)){
		return(m[,col2add])
		}
}




#get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#### pull_species
#funciton to split strings in a vector
split_string_vector = function(fileNameVector, delim, pos){
	subset = c()
	for (f in fileNameVector){
		s=strsplit(f, delim)[[1]][pos]
		subset =append(subset, s)
	}
	return(subset)
}


#### get_conf_bounds
get_conf_bounds = function(v, float){
	q=quantile(v, probs=seq(0, 1, (1-float)), na.rm=T)
	low = paste((1-float)*100, "%", sep="")
	high = paste((float)*100, "%", sep="")
	l=q[low]
	h=q[high]
	return(h)
}


get_upper_tail = function(v, float){
	q=quantile(v, probs=seq(0, 1, (1-float)), na.rm=T)
	high = paste((float)*100, "%", sep="")
	h=q[high]
	return(list(h))
}

get_lower_tail = function(v, float){
	q=quantile(v, probs=seq(0, 1, (1-float)), na.rm=T)
	low = paste((1-float)*100, "%", sep="")
	l=q[low]
	return(list(l))
}

get_two_tailed = function(v, float){
	q=quantile(v, probs=seq(0, 1, ((1-float)/2)), na.rm=T)
	h = q[(length(q) - 1)]
	l = q[2]
	return(list(l, h))
}


plot_conf_bounds = function(df, column, float){
	v=df[,column]
	h=get_conf_bounds(v, float)
	print("Confidence Boundary:")
	print(h)
	
	
g=ggplot(df, aes_string(column)) +
	geom_density() +
	theme_minimal()
	dpb = ggplot_build(g)
	x2 <- min(which(dpb$data[[1]]$x >= h))
	x1 <-which(dpb$data[[1]]$x == max(dpb$data[[1]]$x))
	x3 <-which(dpb$data[[1]]$x == min(dpb$data[[1]]$x))
	g = g + geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
	        y=dpb$data[[1]]$y[x1:x2]),
	        aes(x=x, y=y), fill=gg_color_hue(2)[2])
	g = g + geom_area(data=data.frame(x=dpb$data[[1]]$x[x3:x2],
	        y=dpb$data[[1]]$y[x3:x2]),
	        aes(x=x, y=y), fill=gg_color_hue(2)[1])
	return(g)
}


#### single_chrom_scatterplot
single_chrom_scatterplot = function(df0, chr, xcol, ycol, float=0.95, span=0.2, main=chr, xlab=xcol, ylab=ycol){
	v=df0[,ycol]
	h=get_conf_bounds(v, float)
	df=df0[df0$CHROM==chr,]
	df$sig = df[,ycol] > h
	g = ggplot(df, aes_string(x = xcol, y = ycol)) + 
       geom_point(aes(color=sig)) +
       stat_smooth(span = span, se=FALSE, color='black', lwd=.5, method="loess") +
       geom_hline(yintercept=h, lwd=0.5, lty=2) +
       theme_minimal() + ggtitle(main) + labs(y=ylab, x = xlab)
    return(g)
}

plot_multipanel_chroms = function(df0, chromlist, xcol, ycol, float=0.95, span=0.2, xlab=xcol, ylab=ycol){
	splist = list()
	count=0
	for (i in seq(1, length(chromlist))){
		count=count + 1
		chr=chromlist[i]
		print(chr)
		g=single_chrom_scatterplot(df0, chr, xcol, ycol, float)
		g = g + ggtitle(chr) + labs(y=ylab, x = xlab)
		splist[[i]]<-g
	}
	print(length(splist))
	gd = plot_conf_bounds(df0, ycol, float)
	splist[[count+1]]<-gd
	print(length(splist))
	grid.arrange(grobs=splist, ncol=5)
}


plot_sexChrom_plus_dist = function(df0, chromlist, chr, xcol, ycol, float=0.95, span=0.2, xlab=xcol, ylab=ycol){
	splist = list()
	count=0
	chromlist=as.character(chromlist[chromlist==chr])
	for (i in seq(1, length(chromlist))){
		count=count + 1
		chr=chromlist[i]
		print(chr)
		g=single_chrom_scatterplot(df0, chr, xcol, ycol, float)
		g = g + ggtitle(chr) + labs(y=ylab, x = xlab)
		splist[[i]]<-g
	}
	print(length(splist))
	gd = plot_conf_bounds(df0, ycol, float)
	splist[[count+1]]<-gd
	print(length(splist))
	grid.arrange(grobs=splist, ncol=2)
}

# get_tails = function(v, tail.type, alpha){
	# if (tail.type=="upper"){
		# q=quantile(v, probs=seq(0,1,alpha))
		# pct=paste( (1-alpha)*100, "%", sep='')
		# cut = q[pct]
		# }
	# if (tail.type=="two.tailed"){
		# q=quantile(v, probs=seq(0,1,alpha/2))
		# pct1=paste( (1-alpha/2)*100, "%", sep='')
		# pct1=paste( (alpha/2)*100, "%", sep='')
		# cut = c(q[pct1], q[pc2]
	# }
	# return(list(pct, cut))
# }


#similar to multichrom_len_normalized but for boolean gene tree data
multichrom_len_normalized = function(df, xcol, ycol, sex.chrom, YLAB=ycol, XLAB='', sex.col = 'black', other.col = 'grey', ln.col='blue', SPAN=.1, SE=F, draw.legend=F, SHAPE=19,LWD=0.75,PT.SIZE = 0.25, draw.line=F, tailed = 2, horiz=FALSE, YLIM='none', verticals=F, breakpoints=c(3.5, 18.9), nox=F, noy=F){
	# XLAB = bquote(Pos.~.(sex.chrom)~(Mb))
	df$chrom.type = as.character(df$CHROM)
	df$chrom.type[df$CHROM != sex.chrom]<-'autosomes'
	df$chrom.type[df$CHROM == sex.chrom]<-sex.chrom
	g = ggplot(data=df, aes_string(x =xcol, y = ycol, group="CHROM", colour="chrom.type"))
	if (verticals){
		print("Outlining breakpoints.")
		tot.length = max(df[df$CHROM==sex.chrom,'mb'])
		norm.breaks = breakpoints / tot.length
		g=g+ geom_vline(xintercept= norm.breaks[1], lty=2, lwd=0.8, color='green') + geom_vline(xintercept= norm.breaks[2], lty=2, lwd=0.8, 			color='green')
	}
	g=g + scale_colour_manual(values=c(other.col, sex.col)) +
		geom_smooth(lwd=LWD, span = SPAN, se=SE) +
		geom_point(size= PT.SIZE, shape=SHAPE) +
		labs(y=YLAB, x = XLAB) +
		theme(legend.title=element_blank())
	df.sub = df[df$CHROM==sex.chrom,]
	g = g + geom_smooth(data=df.sub, col=ln.col, lwd=LWD, span=SPAN, se=SE) + geom_point(data=df.sub, size=PT.SIZE, col=sex.col, shape=SHAPE)
	if (draw.line){
		if (tailed==2){
			g = g + geom_hline(aes(yintercept=horiz[[1]]), colour="black", linetype="dashed")
			g = g + geom_hline(aes(yintercept=horiz[[2]]), colour="black", linetype="dashed")
		}
	}
	if (!draw.legend){
		g = g + theme(legend.position="none")
	}
	if (nox){
		g=g+theme(axis.text.x=element_blank())
	}
	if (noy){
		g=g+theme(axis.text.y=element_blank())
	}
	#set the x-tickmarks for the putative sex chromosome
	df.sub = df[df$CHROM==sex.chrom,]
	breaks = c(seq(0, 1, by = 0.25))
	labs = round(breaks*max(df.sub$mb), digits=0)
	g = g + scale_x_continuous(breaks=c(seq(0, 1, by = 0.25)), labels= labs)
	if (YLIM[1] != "none"){
		print("Scaling Y axis")
		g = g + coord_cartesian(ylim = YLIM) 
	}
	return(g)
}

#Funtion to plot a statistic for all chromosomes on single plot with sex chrom overlaid
sdr_normalized = function(df, xcol, ycol, sex.chrom, YLAB=ycol, XLAB='Position (Mb)', sex.col = 'black', other.col = 'grey', SPAN=.1, SE=F, draw.legend=T, SHAPE=19,LWD=0.5,PT.SIZE = 0.25, draw.line=F, tailed = 2, horiz=FALSE){
	df$chrom.type = as.character(df$CHROM)
	df$chrom.type[df$CHROM != sex.chrom]<-'autosomes'
	df$chrom.type[df$CHROM == sex.chrom]<-sex.chrom
	df.sub = df[df$CHROM==sex.chrom,]
	absTicks = c(0, 5, 10, 15, 20)
	xticks = c(0, 5, 10, 15, 20) / max(df.sub$mb)
	yticks = c(0,1)
	ylabs = c("False", "True")
	tot.length = max(df[df$CHROM=='chrXII','mb'])
	norm.breaks = c(3.5, 18.9) / tot.length
	bpDat = data.frame(x=norm.breaks, y1=c(0,0), y2=c(1,1))
	g = ggplot() + 
		geom_segment(data=bpDat, aes(x=x, y=y1, xend=x, yend=y2), lineend='round', lwd=0.5, lty=2, color='green') +
		scale_colour_manual(values=c(other.col, sex.col)) +
		geom_point(data=df, aes_string(x =xcol, y = ycol, group="CHROM", colour="chrom.type"), size= PT.SIZE, shape=SHAPE) +
		labs(y=YLAB, x = XLAB) +
		theme(legend.title=element_blank()) + 
		geom_smooth(data=df.sub, aes_string(x =xcol, y = ycol), col='blue', lwd=LWD, span=SPAN, se=SE) +
		geom_point(data=df.sub, aes_string(x =xcol, y = ycol), size=PT.SIZE, col=sex.col, shape=SHAPE) +
		scale_x_continuous(breaks=c(seq(0, 1, by = 0.25)), labels= absTicks) +
		scale_y_continuous(breaks=yticks, labels= ylabs) + 
	if (!draw.legend){
		g = g + theme(legend.position="none")
	}
	return(g)
}

#biplots
dopairvar = function(df, c1, c2, ci){
	par(mfrow=c(2,1))
	c1mns =tapply(df[,c1], df[,ci], mean)
	c1stders = tapply(df[,c1], df[,ci], std.error)
	c2mns = tapply(df[,c2], df[,ci], mean)
	c2stders = tapply(df[,c2], df[,ci], std.error)
	minx = min(c1mns)-max(c1stders)
	maxx = max(c1mns)+max(c1stders)
	miny = min(c2mns)-max(c2stders)
	maxy = max(c2mns)+max(c2stders)
	xlim=c(minx, maxx)
	ylim=c(miny, maxy)
	col.set = gg_color_hue(length(unique(as.character(df[,ci]))))
	plotCI(x=c1mns, y=c2mns, uiw=c1stders, col=col.set[1:length(df[,ci])], pch=19, err="x", xlim=xlim, ylim=ylim)
	plotCI(x=c1mns, y=c2mns, uiw=c2stders, err="y", add=T, col=col.set[1:length(df[,ci])], pch=19)
}


#plot snp pca
format_y_snp_pca_df = function(pca.object){
	df = data.frame(pca.object$scores)
	df_out = df
	df_out$group <- substr(tolower(rownames(df)), start = 1, stop=3)
	df_out$group = paste(df_out$spp, df_out$chrom, sep="")
	return(df_out)
}

#format y-labeled snps for pca
format_y_snp_pca_df = function(pca.object){
	df = data.frame(pca.object$scores)
	df_out = df
	df_out$chrom <- sapply( strsplit(as.character(row.names(df)), "_"), "[[", 2 )
	df_out$spp <- substr(sapply( strsplit(tolower(as.character(row.names(df))), "_"), "[[", 1 ), start = 1, stop=3)
	df_out[df_out =='A']<-''
	df_out[df_out =='B']<-''
	df_out$group = paste(df_out$spp, df_out$chrom, sep="")
	return(df_out)
}





plot_snp_pca = function(pca.df, pca.object, cx=1, cy=2, group='group', draw.legend=F, shift=F, shift.table=NULL, invertX=1, invertY=1, text.labs=TEXTLABS){
	# pct1 = paste(round(pca.object$eig[cx] / sum(pca$eig), digits=3)*100, "%", sep='')
	# pct2 = paste(round(pca.object$eig[cy] / sum(pca$eig), digits=3)*100, "%", sep='')
	color.set0 = c(hue_pal()(3), hue_pal()(4)[4])
	color.set = color.set0[c(1,4,2,3)]
	pcx = paste("PC", cx, sep='')
	pcy = paste("PC", cy, sep='')
	pca.df[,pcx]=pca.df[,pcx]*invertX
	pca.df[,pcy]=pca.df[,pcy]*invertY
	p<-ggplot(pca.df,aes_string(x=pcx,y=pcy,color=group))
	p<-p+geom_point() + scale_color_manual(values= color.set)
	mn1 = tapply(pca.df[,pcx], INDEX=pca.df[,group], mean)
	mn2 = tapply(pca.df[,pcy], INDEX=pca.df[,group], mean)
	mn.df = data.frame(mn1, mn2)
	mn.df$sexmn = rownames(mn.df)
	print("Means:")
	print(mn.df)
	if (text.labs){
		if (shift){
			mn.df$mn1 = mn.df$mn1 + shift.table$x
			mn.df$mn2 = mn.df$mn2 + shift.table$y
			print("Shifted means:")
			print(mn.df)
		}
		p = p + geom_text(data=mn.df, size=5, aes(x=mn1, y=mn2, label= mn.df$sexmn)) + labs(x=paste(pcx), y=paste(pcy))
	}
	if (!draw.legend){
			p = p + theme(legend.position="none")
		}
	p=p+theme(axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.y = element_blank(),
	axis.text.x = element_blank(),
	axis.ticks.y = element_blank(),
	axis.ticks.x = element_blank() )
	return(p)
}


format_admixture = function(tbl, names, sex){
	m = merge(names, sex, by = "V1", all.x=T)
	m$V3[is.na(m$V3)]<-''
	m$V3[m$V3==1]<-'F'
	m$V3[m$V3==2]<-'M'
	rownames(m) = m$V1
	om = m[names$V1,]
	sum(om$V1==names$V1) == length(names$V1)
	fnames = paste(om$V3, tolower(om$V1), sep='')
	rownames(tbl) = fnames
	otbl = tbl[order(rownames(tbl)),]
	otbl$names = rownames(otbl)
	tbl2=melt(otbl, id.vars="names")
	return(tbl2)
}

plot_admixture = function(tbl){
	p<-ggplot(data=tbl, aes(x=names, y=value, fill=variable)) +
	geom_bar(stat="identity") +
	labs(y="Ancestry")
	return(p)
}


plot_admixture_bars_only = function(tbl, legend.pos="none"){
	p<-ggplot(data=tbl, aes(x=names, y=value, fill=variable, width=1)) +
	geom_bar(stat="identity") +
	theme( legend.position = legend.pos,
	legend.title=element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.y = element_blank(),
	axis.text.x = element_blank(),
	axis.ticks.x=element_blank() ) +
	labs(y='') 
	return(p)
}
		
plot_distance_density = function(ddf, chr.selection, dstat, xticks, xlim){
	fdf.auto = ddf[!ddf$chr %in% chr.selection,]
	fden = density(fdf.auto[,dstat]) #get densities for chromosomes only
	fden.df = data.frame(x=fden$x, y=fden$y)
	
	fsex = na.omit(ddf[ddf$chr %in% chr.selection, dstat])
	fden.sex = density(fsex)
	fsex.df = data.frame(x= fden.sex$x, y= fden.sex$y)
	

	gg=ggplot(data= fsex.df, aes(x=x, y=y)) + 
		geom_area(data=fsex.df, fill='black', alpha=0.75) +
		geom_area(data=fden.df, fill='grey', alpha=0.75) +
		scale_x_continuous(breaks = xticks, limits=xlim) +
		labs(x=dstat, y="") 
	return(gg)
}	
		
		
		