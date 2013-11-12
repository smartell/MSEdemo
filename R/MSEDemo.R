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








# .plotTrueBiomass(sims)
# .plotCatch(sims)