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
	    #
	    selectInput('scenario','Scenarios'
	                ,levels(raw.data$Scenario),
	                multiple=TRUE,selected=c('DET.N.Theft','PDO.N.Theft')),
	    helpText('Note: Select multiple scenarios using the shift or control key.'),

		checkboxInput('integrate','Integrate Scenarios',value=TRUE),

		wellPanel( 
			h5('Harvest Control Rule'),
			sliderInput('hcr1', 'SSB Limit'    , min=0, max=1.0, value=0.2, step = 0.05),
			sliderInput('hcr2', 'SSB Threshold', min=0, max=1.0, value=0.3, step = 0.05),
			sliderInput('hcr3', 'SSB Target   ', min=0, max=1.0, value=0.4, step = 0.05),
			numericInput('hcr4', 'F target',value=0.2,min=0.05,max=0.50,step=0.05)
		),

		h4('Management Procedure:'),
		checkboxGroupInput('iclr','Harvest Control Rule',
		                   levels(raw.data$MP),selected="ThirtyTwenty")


		# selectInput('facet_row','Facet Row',c(None='.',argv)),
		# selectInput('facet_col','Facet Column',c(None='.',argv),select='Scenario')

		
		),

	mainPanel(
		tabsetPanel(
			tabPanel("Worm Plots",uiOutput("gVisoutput"),plotOutput("msePlot")),
			tabPanel("MSE Player",htmlOutput("motionchart")),
			tabPanel("Table",
			         h4("Median Depletion"), tableOutput("viewDepletionTable"),
			         h4("Median Catch"),tableOutput("viewCatchTable"),
			         h4("5-year Average Annual Catch Variation"),tableOutput("viewAAVTable"),
			         h4("Probability of closing the fishery"),tableOutput("viewClosedTable"),
			         h4("Probability of SSB < 0.20 (Objective P < 0.05)"),tableOutput("viewSSBlimit"),
			         h4("Probability of SSB < 0.30 (Objective P < 0.25)"),tableOutput("viewSSBthreshold") ),
		    tabPanel("Harvest control Rule",plotOutput("gHCR"),plotOutput("gCt"))
			)	          
	    )
))