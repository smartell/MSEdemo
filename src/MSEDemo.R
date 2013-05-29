# MSEDemo.R
require(reshape2)
require(ggplot2)
# setwd('/Users/stevenmartell1/Documents/MSEdemo/demo/')
# load('demo.Rdata')

setwd('~/Documents/IPHC/MSEdemo/fortyten/')
load('fortyten.Rdata')

.plotTrueBiomass <- function( S )
{
	fn  <- function(S){return(S$t_bt)}
	df  <- sapply(S,fn)
	qf  <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.95)))
	ii  <- sample(1:dim(df)[2],5)
	df  <- data.frame(Year=(1964+1:dim(qf)[1]),df[,ii])
	mdf <- melt(df,id.vars=c("Year"))
	qf  <- data.frame(Year=(1964+1:dim(qf)[1]),qf)
	colnames(qf) = c("Year","lci","median","uci")

	eb  <- aes(ymin=lci,ymax=uci)
	ggp <- ggplot(qf,aes(x=Year,y=median))
	ggp <- ggp + geom_ribbon(eb,alpha=0.40)
	ggp <- ggp + geom_line(data=mdf,aes(x=Year,y=value,col=variable),alpha=1)
	ggp <- ggp + geom_line(size=2,alpha=0.5) 
	ggp <- ggp + labs(x="Year",y="Biomass (t)") + theme(legend.position="none")
	print(ggp)
}

.plotCatch <- function( S )
{
	fn <- function(S){return(S$ct)}
	df <- sapply(S,fn)
	qr <- t(apply(df,1,quantile,probs=c(0.00,1.0)))
	qf <- t(apply(df,1,quantile,probs=c(0.05,0.5,0.9)))
	ii  <- sample(1:dim(df)[2],5)
	df  <- data.frame(Year=(1964+1:dim(qf)[1]),df[,ii])
	mdf <- melt(df,id.vars=c("Year"))
	qf <- data.frame(Year=(1964+1:dim(qf)[1]),qf,qr)
	colnames(qf) = c("Year","lci","median","uci","lrange","urange")

	eb  <- aes(ymin=lci,ymax=uci)
	rb  <- aes(ymin=lrange,ymax=urange)
	ggp <- ggplot(qf,aes(x=Year,y=median))
	ggp <- ggp + geom_ribbon(rb,alpha=0.15)
	ggp <- ggp + geom_ribbon(eb,alpha=0.25)
	ggp <- ggp + geom_line(data=mdf,aes(x=Year,y=value,col=variable),alpha=1)
	ggp <- ggp + geom_line(size=2,alpha=0.5) 
	ggp <- ggp + labs(x="Year",y="Catch (t)") + theme(legend.position="none")
	print(ggp)
}

.plotTrueBiomass(sims)
.plotCatch(sims)