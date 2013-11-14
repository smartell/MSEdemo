# MSEDemo.R
require(reshape2)
require(ggplot2)
require(PBSmodelling)
# setwd('/Users/stevenmartell1/Documents/MSEdemo/demo/')
# load('demo.Rdata')

setwd('~/Documents/IPHC/MSEdemo/R/')
source("global.R")

.RFILES     <- list.files("./R-scripts/",pattern="\\.[Rr]$")
for(nm in .RFILES) source(file.path("./R-scripts", nm), echo=FALSE)




.FONTSIZE <- 12
#load('fixedCCC.Rdata')


.guiSetup <- function()
{
	ifiles   <- read.table("MSEView.txt",header=TRUE,sep=",")
	createWin("MSEDemoWin.txt")
	hdr      <- ifiles[ ifiles$Select, ]
	fn       <- as.character(hdr$File.Path)
	getObj   <- function(fn){load(fn);return(A)}
	M        <<- lapply(fn,getObj)
	names(M) <<- hdr$Label
}

guiView <- function()
{
	.guiSetup()
}


.viewPlot <- function()
{
	guiInfo <- getWinVal(scope="L")
	mdlname <- guiInfo$ifiles$Label[guiInfo$ifiles$Select]
	ii      <- names(M) %in% mdlname
	M       <- M[ii]

	# ASS Plots
	if( plotType=="ass.biomass" )
	{
		.plotASSbiomass( M )
	}
	if( plotType=="ass.ut" )
	{
		.plotASSexploitation( M )
	}
	if( plotType=="ass.catch" )
	{
		.plotASScatch( M )
	}
	if( plotType=="ass.cpue" )
	{
		.plotASScpue( M )
	}
	if( plotType=="ass.cpue.res" )
	{
		.plotCPUEresiduals( M )
	}
	if( plotType=="ass.rec.dev" )
	{
		.plotRecDevs( M )
	}


	# MSE Plots
	if( plotType=="mse.biomass" )
	{
		.plotMSETrueBiomass( M )
	}
	if( plotType=="mse.catch" )
	{
		.plotMSECatch( M )
	}
}


.plotCPUEresiduals <- function( M )
{
	cat(".plotCPUEresiduals\n")
	n <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		df  <- M[[i]]$it_data
		ep  <-na.omit(as.vector(t(M[[i]]$epsilon)))
		df  <- data.frame(Model=names(M)[i],df,Residual=ep)
		names(df) = c("Model","Year","Epoch","CPUE","CV","Residual")
		mdf <- rbind(mdf,df)
	}

	p <- ggplot(mdf,aes(x=Year,y=Residual,fill=factor(Epoch))) 
	p <- p + geom_bar(stat='identity',position='dodge')
	p <- p + facet_wrap(~Model)
	print(p)

}

.plotRecDevs <- function( M )
{
	cat("plotRecDevs.\n")
	n   <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		yr  <- M[[i]]$yr
		wt  <- M[[i]]$wt
		df  <- data.frame(Model=names(M)[i],Year=yr,wt=wt)
		mdf <- rbind(mdf,df)
	}

	p <- ggplot(mdf,aes(x=Year,y=wt)) + geom_bar(stat='identity')
	p <- p + facet_wrap(~Model)
	print(p)
}





# .plotTrueBiomass(sims)
# .plotCatch(sims)