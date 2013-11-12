

.saveDataFrame <- function(M)
{
	n <- length(M)
	.NSAMP    <- 5
	IDX       <- sample(1:length(M[[1]]),.NSAMP)
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
		closed <- closed/dim(tmp)[2]
		
		# P(Bt < 0.20)
		fnp <- function(X)(return(X$t_bt/X$t_bo <= 0.2))
		tmp <- sapply(M[[i]],fnp)
		pb20 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb20 <- pb20/dim(tmp)[2]

		# P(Bt < 0.30)
		fnp <- function(X)(return(X$t_bt/X$t_bo <= 0.3))
		tmp <- sapply(M[[i]],fnp)
		pb30 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb30 <- pb30/dim(tmp)[2]

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
		                  btrc,ctrc,pb20,pb30)
		colnames(qf) <- c("Scenario","MP","Year","Bt.lci","Biomass","Bt.uci",
		                  "Ct.lci","Catch","Ct.uci",
		                  "Dt.lci","Depletion","Dt.uci",
		                  "PerfectInfo.Bt","AAV",
		                  "Hr.lci","HarvestRate","Hr.uci",
		                  "Closures",
		                  paste0("btrc",1:length(IDX)),
		                  paste0("ctrc",1:length(IDX)),
		                  "P(SSB<0.2)","P(SSB<0.3)")	
		qdf <- rbind(qdf,qf)
	}
	scen <- qdf$Scenario; sn<-paste0("S",1:length(unique(scen))); levels(scen)=sn
	manp <- qdf$MP;       pn<-paste0("P",1:length(unique(manp))); #levels(manp)=pn
	qdf <- cbind(qdf,SP=paste(scen,manp,sep=""))

	save(qdf,file="QDF.Rdata")
}