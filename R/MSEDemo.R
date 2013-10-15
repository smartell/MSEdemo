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
.NSAMP    <- 5
IDX       <- sample(1:length(M[[1]]),.NSAMP)
.FONTSIZE <- 18
#load('fixedCCC.Rdata')

.saveDataFrame <- function(M)
{
	n <- length(M)
	S_HCR <- strsplit(names(M),"_")
	year  <- 1964:(1964+length(M[[1]][[1]]$t_bt)-1)
	mdf   <- qdf  <- NULL
	for (i in 1:n)
	{
		cat("Scenario ",S_HCR[[i]][1],"\t")
		cat("HCR ",S_HCR[[i]][2],"\n")
		fn  <- function(X){return(X$t_bt)}
		df  <- sapply(M[[i]],fn)
		qf  <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.95)))

		fnc <- function(X){return(X$ct)}
		cdf <- sapply(M[[i]],fnc)
		cqf <- t(apply(cdf,1,quantile,probs=c(0.05,0.5,0.95)))
		#cmu <- apply(cdf,1,mean)
		#csd <- apply(cdf,1,sd)

		# depletion
		fnd <- function(X)(return(X$t_bt/X$t_bo))
		ddf <- sapply(M[[i]],fnd)
		dmed <- t(apply(ddf,1,quantile,probs=c(0.05,0.5,0.95)))
		#cdf <- (apply(ddf,1,mean))
		#sdf <- (apply(ddf,1,sd))

		# perfect information
		fnp <- function(X)(return(X$p_bt))
		tmp <- sapply(M[[i]],fnp)
		mdf <- apply(tmp,1,median)

		# AAV
		fnaav <- function(X)(return(X$t_aav))
		tmp <- sapply(M[[i]],fnaav)
		aav <- apply(tmp,1,median)

		# No of closures
		fnclose <- function(X)(return(X$ct != 0))
		tmp <- sapply(M[[i]],fnclose)
		closed <- apply(tmp,1,FUN=function(tmp){sum(tmp==FALSE)})
		closed <- closed/prod(dim(tmp)) * 100
		# closed <- apply(tmp,1,sum)

		# HarvestRate
		fnhr <- function(X)(return(X$ct/X$t_bt))
		dfhr <- sapply(M[[i]],fnhr)
		dhr  <- t(apply(dfhr,1,quantile,probs=c(0.05,0.5,0.95)))

		# Biomass traces
		fntb <- function(X)(return(X$t_bt))
		tmp  <- sapply(M[[i]],fntb)
		btrc <- tmp[,IDX]

		# Catch Traces
		fntc <- function(X)(return(X$ct))
		tmp  <- sapply(M[[i]],fntc)
		ctrc <- tmp[,IDX]
		# print(i)
		qf  <- data.frame(Scenario=S_HCR[[i]][1],MP=S_HCR[[i]][2],
		                  Year=year,qf,cqf,dmed,mdf,aav,dhr,closed,
		                  btrc,ctrc)
		colnames(qf) <- c("Scenario","MP","Year","Bt.lci","Biomass","Bt.uci",
		                  "Ct.lci","Catch","Ct.uci",
		                  "Dt.lci","Depletion","Dt.uci",
		                  "PerfectInfo.Bt","AAV",
		                  "Hr.lci","HarvestRate","Hr.uci",
		                  "Closures",
		                  paste0("btrc",1:length(IDX)),
		                  paste0("ctrc",1:length(IDX)))	
		qdf <- rbind(qdf,qf)
	}
	scen <- qdf$Scenario; sn<-paste0("S",1:length(unique(scen))); levels(scen)=sn
	manp <- qdf$MP;       pn<-paste0("P",1:length(unique(manp))); #levels(manp)=pn
	qdf <- cbind(qdf,SP=paste(scen,manp,sep=""))

	save(qdf,file="QDF.Rdata")
}

.viewPlot <- function()
{
	guiInfo <- getWinVal(scope="L")
	mdlname <- guiInfo$ifiles$Label[guiInfo$ifiles$Select]
	ii      <- names(M) %in% mdlname
	M       <- M[ii]
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
	cat(n)
	if(n > 1) print(ggp + facet_wrap(~Model) ) else print(ggp)
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
		qf <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.95)))
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