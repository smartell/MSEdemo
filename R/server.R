#server.R
library(shiny)
library(ggplot2)
library(plyr)
load("QDF.Rdata")
raw.data <- qdf



shinyServer(function(input,output){

	data <- reactive({
		a <- subset(raw.data, MP %in% input$iclr & Year <= input$x)
		a <- cbind(SMP=paste(a$Scenario,a$MP,sep="."),a)
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
	output$msePlot <- renderPlot({
		p <- ggplot(data(),aes_string(x="Year",y=input$y))

		if(input$y == 'Biomass')
		{
			eb <- aes(ymin=Bt.lci,ymax=Bt.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$y),alpha=1)
		}
		if(input$y == 'Landings')
		{
			eb <- aes(ymin=Ct.lci,ymax=Ct.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$y),alpha=1)
		}

		# if(input$iclr != 'none')
		{
			 p <- p + aes_string(color="MP",fill="MP")
		}

		

		facets <- paste(input$facet_row, '~', input$facet_col)
		if(facets != '. ~ .')
		{
			p <- p + facet_grid(facets)
		}

		print(p + theme_bw())
	})

	# TAB 2
	output$gVisoutput <- renderUI({
		print("Hello Bob")
		selectInput("performance.measure","Peformance Measure",
		            c(Depletion="Depletion"),select='Depletion')
	})

	output$motionchart <- renderGvis({
		M1 <- gvisMotionChart(data(), idvar="SMP", timevar="Year",
		                      xvar="Depletion.Mean",yvar="Landings",
		                      colorvar="MP",sizevar="Catch.CV")
		return(M1)
	})


})