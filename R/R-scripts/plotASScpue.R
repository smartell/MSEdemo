.plotASScpue <- function( M )
{
	cat(".plotASScpue\n")
	n   <- length(M)
	mdf <- NULL

	for(i in 1:n)
	{
		df  <- M[[i]]$it_data
		df  <- data.frame(Model=names(M)[i],df)
		names(df) = c("Model","Year","Epoch","CPUE","CV")
		mdf <- rbind(mdf,df)
	}

	p <- ggplot(mdf,aes(x=Year,y=CPUE,col=factor(Epoch))) + geom_point()
	p <- p + geom_linerange(aes(Year,ymin=exp(log(CPUE)-1.96*CV),
	                                 ymax=exp(log(CPUE)+1.96*CV)))
	p <- p + facet_wrap(~Model)
	print(p)
}
