.plotASSbiomass <- function( M ) 
{
	n <- length(M)
	cat(".plotASSBiomass\n")
	mdf <- NULL
	for( i in 1:n )
	{
		bt  <- M[[i]]$bt
		yrs <- M[[i]]$yrs

		df  <- data.frame(Model=names(M)[i],Year=yrs,Biomass=bt)
		mdf <- rbind(mdf,melt(df,id.vars=c("Model","Year")))

	}
	p  <- ggplot(mdf,aes(Year,value,linetype=Model)) + geom_line()
	p  <- p + labs(x="Year",y="Biomass (million net pounds)")
	print(p)
}