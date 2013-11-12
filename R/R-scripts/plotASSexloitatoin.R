.plotASSexploitation <- function( M ) 
{
	n <- length(M)
	cat(".plotASSExploitation\n")
	mdf <- NULL
	for( i in 1:n )
	{
		ft  <- M[[i]]$ft
		yr <- M[[i]]$yr

		df  <- data.frame(Model=names(M)[i],Year=yr,Exploitation=ft)
		mdf <- rbind(mdf,melt(df,id.vars=c("Model","Year")))

	}
	p  <- ggplot(mdf,aes(Year,value,linetype=Model)) + geom_line()
	p  <- p + labs(x="Year",y="Exploitation (Catch/Biomass)")
	print(p)
}
