#server.R
library(shiny)
library(ggplot2)
library(plyr)
load("QDF.Rdata")
raw.data <- qdf
mi <- c(-4,-6,-7,-9,-11:-13)
argv <- names(raw.data)[mi]


shinyServer(function(input,output){

	data <- reactive({
		a <- subset(raw.data, MP %in% input$iclr & Year <= input$x)
		# a <- cbind(a,SMP=paste(a$MP))
		a <- droplevels(a)
		return(a)
	})


	captionText <- reactive({
		paste(input$y)
	})

	# Return caption
	output$caption <- renderText({
		captionText()
	})


	# TAB 1
	mse_plotinput = function(...){
		tagList(
				selectInput('yy','Variable',argv[-1:-3],select='Biomass'),
				selectInput('facet_irow','Facet Row',c(None='.',argv[c(1,2)])),
				selectInput('facet_jcol','Facet Column',c(None='.',argv[c(1,2)]),select='Scenario')
			)
	}

	output$gVisoutput <- renderUI({
		
		# selectInput("performance.measure","Peformance Measure",
		            # c(Depletion="Depletion"),select='Depletion')

		wellPanel({
			mse_plotinput()
		})

	})

	output$msePlot <- renderPlot({
		p <- ggplot(data(),aes_string(x="Year",y=input$yy))

		if(input$yy == 'Biomass')
		{
			eb <- aes(ymin=Bt.lci,ymax=Bt.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		if(input$yy == 'Landings')
		{
			eb <- aes(ymin=Ct.lci,ymax=Ct.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		if(input$yy == 'Depletion.Median')
		{
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		# if(input$iclr != 'none')
		{
			 p <- p + aes_string(color="MP",fill="MP")
		}

		

		facets <- paste(input$facet_irow, '~', input$facet_jcol)
		if(facets != '. ~ .')
		{
			p <- p + facet_grid(facets)
		}

		print(p + theme_bw())
	})

	# TAB 2

	output$motionchart <- renderGvis({
		# Combine Scenarios
		mdf <- melt(data(),id=c("Scenario","MP","Year"))
		tmp <- cast(mdf,MP+Year~variable,mean)
		tmp <- cbind(tmp,Procedure=tmp$MP)
		
		M1 <- gvisMotionChart(tmp, idvar="MP", timevar="Year",
		                      xvar="Depletion.Median",yvar="Landings",
		                      sizevar="Catch.CV",colorvar="Procedure")
		return(M1)
	})


})