#twisst_plotting_functions.R

simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights)
    y.predict <- predict(y.loess,x)
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth_df <- function(x, df, span, col.names=NULL, weights=NULL, min=NULL, max=NULL){
    smoothed <- df
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, max = max, min = min, weights = weights)
        }
    smoothed
    }


stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }


plot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (is.matrix(x)==TRUE) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) p=plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty, axes=F)
    
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA)
            lines(x,y, type = "l", col = line_cols[n])
            }
        }
    }


ggplot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE, draw.legend=F){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (is.matrix(x)==TRUE) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) p=plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty, axes=F)
    
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            # polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], border=NA)
            }
        #build dataframe for ggplot
        upper = data.frame(y_stacked[["upper"]])
        lower = data.frame(y_stacked[["lower"]])
        df_stacked = data.frame()
		for (n in 1:ncol(upper)){
			add.df = data.frame(x = c(x, rev(x)), y = c(upper[,n], rev(lower[,n])), group = n)
			df_stacked = rbind(df_stacked, add.df)
		}
		df_stacked$group = as.factor(df_stacked$group)
		head(df_stacked)
		g = ggplot() + geom_polygon(data= df_stacked, mapping=aes(x=x, y=y, group=group, fill=group)) +
		scale_fill_manual(values=fill_cols) +
		labs(x = "Position (Mb)", y = "Topology weight")
		if (!draw.legend==TRUE){
			g = g + theme(legend.position="none")
		}
		plot(g)
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA)
            lines(x,y, type = "l", col = line_cols[n])
            }
        }
    return(g)
    }



    
cols = c(
# "#F0A3FF", #Amethyst 
"#0075DC", #Blue
"#993F00", #Caramel
"#4C005C", #Damson
"#8F7C00", #Khaki
"#2BCE48", #Green
"#191919", #Ebony
"#005C31", #Forest
"#FFCC99", #Honeydew
"#808080", #Iron
"#94FFB5", #Jade
"#8F7C00", #Khaki
"#9DCC00", #Lime
"#C20088", #Mallow
"#003380", #Navy
"#FFA405", #Orpiment
"#FFA8BB", #Pink
"#426600", #Quagmire
"#FF0010", #Red
"#5EF1F2", #Sky
"#00998F", #Turquoise
"#E0FF66", #Uranium
"#740AFF", #Violet
"#990000", #Wine
"#FFFF80", #Xanthin
"#FFFF00", #Yellow
"#FF5005" #Zinnia
)

#semi-transparent version of each colour
trans_cols = paste0(cols, "25")

#get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
