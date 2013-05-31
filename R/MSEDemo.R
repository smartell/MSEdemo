# MSEDemo.R
require(reshape2)
require(ggplot2)
require(PBSmodelling)
# setwd('/Users/stevenmartell1/Documents/MSEdemo/demo/')
# load('demo.Rdata')

setwd('~/Documents/IPHC/MSEdemo/R/')

ifiles   <- read.table("MSEView.txt",header=TRUE,sep=",")
createWin("MSEDemoWin.txt")
hdr      <- ifiles[ ifiles$Select, ]
fn       <- as.character(hdr$File.Path)
getObj   <- function(fn){load(fn);return(sims)}
M        <- lapply(fn,getObj)
names(M) <- hdr$Label
IDX       <- sample(1:length(M[[1]]),3)
.FONTSIZE <- 18
#load('fixedCCC.Rdata')

.viewPlot <- function()
{
	guiInfo <- getWinVal(scope="L")
	mdlname <- guiInfo$ifiles$Label[guiInfo$ifiles$Select]
	names(M) <- mdlname
	if( plotType=="biomass" )
	{
		.plotTrueBiomass( M )
	}
	if( plotType=="catch" )
	{
		.plotCatch( M )
	}
}

.plotTrueBiomass <- function( S )
{
	n  <- length(S)
	mdf<- NULL
	qdf<- NULL
	for(i in 1:n)
	{
		fn  <- function(X){return(X$t_bt)}
		df  <- sapply(S[[i]],fn)
		qf  <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.95)))
		# if(i==1)
		# IDX  <- sample(1:dim(df)[2],5)
		df  <- data.frame(Model=names(S)[i],Year=(1964+1:dim(qf)[1]),df[,IDX])
		mdf <- rbind(mdf,melt(df,id.vars=c("Model","Year")))
		qf  <- data.frame(Model=names(S)[i],Year=(1964+1:dim(qf)[1]),qf)
		colnames(qf) = c("Model","Year","lci","median","uci")
		qdf <- rbind(qdf,qf)
	}



	eb  <- aes(ymin=lci,ymax=uci)
	ggp <- ggplot(qdf,aes(x=Year,y=median))
	ggp <- ggp + geom_ribbon(eb,alpha=0.40)
	ggp <- ggp + geom_line(data=mdf,aes(x=Year,y=value,col=variable),alpha=1)
	ggp <- ggp + geom_line(size=2,alpha=0.5) 
	ggp <- ggp + labs(x="Year",y="Biomass (t)")
	ggp <- ggp + theme_bw(.FONTSIZE) + theme(legend.position="none")
	print(ggp + facet_wrap(~Model) )
}

.plotCatch <- function( S )
{
	n  <- length(S)
	mdf <- NULL
	qdf <- NULL
	for(i in 1:n)
	{
		fn <- function(S){return(as.vector(S$ct))}
		df <- sapply(S[[i]],fn)
		qr <- t(apply(df,1,quantile,probs=c(0.00,1.0)))
		qf <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.9)))
		# if(i==1)
		# IDX  <- sample(1:dim(df)[2],5)
		df  <- data.frame(Model=names(S)[i],Year=(1964+1:dim(qf)[1]),df[,IDX])
		mdf <- rbind(mdf,melt(df,id.vars=c("Model","Year")))
		qf <- data.frame(Model=names(S)[i],Year=(1964+1:dim(qf)[1]),qf,qr)
		colnames(qf) = c("Model","Year","lci","median","uci","lrange","urange")
		qdf <- rbind(qdf,qf)
	}

	eb  <- aes(ymin=lci,ymax=uci)
	rb  <- aes(ymin=lrange,ymax=urange)
	ggp <- ggplot(qdf,aes(x=Year,y=median))
	ggp <- ggp + geom_ribbon(rb,alpha=0.15)
	ggp <- ggp + geom_ribbon(eb,alpha=0.25)
	ggp <- ggp + geom_line(data=mdf,aes(x=Year,y=value,col=variable),alpha=1)
	ggp <- ggp + geom_line(size=2,alpha=0.5) 
	ggp <- ggp + labs(x="Year",y="Catch (t)") 
	ggp <- ggp + theme_bw(.FONTSIZE) + theme(legend.position="none")
	print(ggp+facet_wrap(~Model))
}

# .plotTrueBiomass(sims)
# .plotCatch(sims)