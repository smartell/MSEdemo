#ui.R
library(googleVis)
library(reshape)
library(reldist)
library(shiny)
library(ggplot2)
load("QDF.Rdata")
raw.data <- qdf
mi <- c(-4,-6,-7,-9)
argv <- names(raw.data)[mi]

shinyUI(pageWithSidebar(

	# Application title
	headerPanel("Management Strategy Evaluation Explorer"),

	sidebarPanel(
		
	    numericInput('x','Maximum year',value=2010,
	                 min=min(qdf$Year),max=max(qdf$Year),step=5),

		# selectInput('y','Variable',argv[-1:-3],select='Biomass'),

		h4('Management Procedure:'),
		checkboxGroupInput('iclr','Harvest Control Rule',
		                   levels(raw.data$MP),selected="FortyTen"),

		selectInput('facet_row','Facet Row',c(None='.',argv)),
		selectInput('facet_col','Facet Column',c(None='.',argv),select='Scenario')

		
		),

	mainPanel(
		tabsetPanel(
			tabPanel("Plot Variables",uiOutput("gVisoutput"),plotOutput("msePlot")),
			tabPanel("Performance Measures",htmlOutput("motionchart")),
			tabPanel("Table",   tableOutput("table"))
			)	          
	    )
))