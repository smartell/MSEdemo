.plotASScatch <- function( M )
{
	cat("plotASScatch\n")
	n <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		yr <- M[[i]]$yr
		ct <- M[[i]]$ct
		df <- data.frame(Model=names(M)[i],Year=yr,Catch=ct)
		mdf <- rbind(mdf,melt(df,id.vars=c("Model","Year")))
	}

	p <- ggplot(mdf,aes(x=factor(Year),y=value,fill=variable))
	p <- p + geom_bar(stat='identity') + facet_wrap(~Model)
	p <- p + scale_x_discrete(breaks=pretty(mdf$Year))
	print(p)
}