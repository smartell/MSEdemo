#server.R
# library(shiny)
# library(ggplot2)
# library(plyr)
# load("QDF.Rdata")
# raw.data <- qdf
# mi <- c(-4,-6,-7,-9,-11:-13)
# argv <- names(raw.data)[mi]


shinyServer(function(input,output){

	data <- reactive({
		a <- subset(raw.data, 
		            MP %in% input$iclr           &
		            Scenario %in% input$scenario &
		            Year >= input$syr            &
		            Year <= input$nyr)
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
		    withTags(div(class='row-fluid',    
			  div(class='span3',	selectInput('yy','Variable',argv[-1:-3],select='Biomass')),
			  div(class='span3',	selectInput('facet_irow','Facet Row',c(None='.',argv[c(1,2)]))),
			  div(class='span3',	selectInput('facet_jcol','Facet Column',c(None='.',argv[c(1,2)]),select='Scenario'))
			))
			)
	}

	output$gVisoutput <- renderUI({
		
		# selectInput("performance.measure","Peformance Measure",
		            # c(Depletion="Depletion"),select='Depletion')

		wellPanel({
			mse_plotinput()
		})

	})

	# This is the ggplot on the opening page.
	output$msePlot <- renderPlot({
		p <- ggplot(data(),aes_string(x="Year",y=input$yy))

		if(input$yy == 'Biomass')
		{
			eb <- aes(ymin=Bt.lci,ymax=Bt.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		if(input$yy == 'Catch')
		{
			eb <- aes(ymin=Ct.lci,ymax=Ct.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		if(input$yy == 'Depletion')
		{
			p <- p + ylim(0,1)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
		}
		# print(length(input$iclr))
		if(length(input$iclr) > 1)
		{
			 p <- p + aes_string(color="MP",fill="MP")
		}

		facets <- paste(input$facet_irow, '~', input$facet_jcol)
		if( length(input$iclr)==1 && facets =='. ~ .')
		{
			p <- p + aes_string(color="Scenario",fill="Scenario")	
		}
		

		if(facets != '. ~ .')
		{
			p <- p + facet_grid(facets)
		}

		print(p + theme_bw())
	})

	# TAB 2
	# Function for the motion chart on the second tab
	output$motionchart <- renderGvis({

		# Combine Scenarios
		print(input$integrate)
		if( input$integrate )
		{
			mdf <- melt(data(),id=c("Scenario","MP","Year","SP"))
			tmp <- cast(mdf,MP+Year~variable,mean)
			tmp <- cbind(tmp,Procedure=tmp$MP)			
		
		
			M1 <- gvisMotionChart(tmp, idvar="MP", timevar="Year",
		                      xvar="Depletion",yvar="Catch",
		                      sizevar="AAV",colorvar="Procedure",
		                      options=list(height=700,width=900))
		}
		if( !input$integrate )
		{  #NOT WORKING YET
			tmp <- data()
			tmp <- cbind(tmp,Procedure=tmp$MP)

			M1 <- gvisMotionChart(tmp, idvar="SP", timevar="Year",
		                      xvar="Depletion",yvar="Catch",
		                      sizevar="AAV",colorvar="Scenario",
		                      options=list(height=700,width=900))
		}


		return(M1)
	})


})