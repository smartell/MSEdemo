#ui.R
# library(googleVis)
# library(reshape)
# library(reldist)
# library(shiny)
# library(ggplot2)
# load("QDF.Rdata")
# raw.data <- qdf
# mi <- c(-4,-6,-7,-9)
# argv <- names(raw.data)[mi]

shinyUI(pageWithSidebar(

	# Application title
	headerPanel("Management Strategy Evaluation Explorer"),

	sidebarPanel(
		
		
		numericInput('syr','Start Year',value=1980,
                 min=min(qdf$Year),max=max(qdf$Year),step=1),
    	numericInput('nyr','Terminal year',value=2015,
                 min=min(qdf$Year),max=max(qdf$Year),step=1),
	    
		# selectInput('y','Variable',argv[-1:-3],select='Biomass'),
	    h4('Scenarios:'),
	    selectInput('scenario','Select multiple scenarios with shift key'
	                ,levels(raw.data$Scenario),
	                multiple=TRUE,selected='Deterministic.N.Theft'),

		checkboxInput('integrate','Integrate Scenarios',value=TRUE),

		h4('Management Procedure:'),
		checkboxGroupInput('iclr','Harvest Control Rule',
		                   levels(raw.data$MP),selected="FortyTen")


		# selectInput('facet_row','Facet Row',c(None='.',argv)),
		# selectInput('facet_col','Facet Column',c(None='.',argv),select='Scenario')

		
		),

	mainPanel(
		tabsetPanel(
			tabPanel("Worm Plots",uiOutput("gVisoutput"),plotOutput("msePlot")),
			tabPanel("MSE Player",htmlOutput("motionchart")),
			tabPanel("Table",   tableOutput("table"))
			)	          
	    )
))