#server.R
library(shiny)
library(ggplot2)
load("QDF.Rdata")
dataset <- qdf



shinyServer(function(input,output){

	captionText <- reactive({
		paste(input$y)
	})

	# Return caption
	output$caption <- renderText({
		captionText()
	})


	output$msePlot <- reactivePlot(function(){
		p <- ggplot(qdf,aes_string(x=input$x,y=input$y)) + geom_line()

		if(input$y == 'Biomass')
		{
			eb <- aes(ymin=Bt.lci,ymax=Bt.uci)
			p <- p + geom_ribbon(eb,alpha=0.05)
		}

		if(input$clr != 'none')
		{
			 p <- p + aes_string(color=input$clr)
		}

		facets <- paste(input$facet_row, '~', input$facet_col)
		if(facets != '. ~ .')
		{
			p <- p + facet_grid(facets)
		}

		print(p + theme_bw())
	})

})