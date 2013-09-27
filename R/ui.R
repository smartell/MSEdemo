#ui.R
library(shiny)
library(ggplot2)
load("QDF.Rdata")
dataset <- qdf
mi <- c(-4,-6,-7,-9)
argv <- names(dataset)[mi]

shinyUI(pageWithSidebar(

	# Application title
	headerPanel("Management Strategy Evaluation Explorer"),

	sidebarPanel(
		selectInput('x','X',argv,select='Year'),
		selectInput('y','Y',argv,select='Biomass'),
		selectInput('clr','Color',c('none',argv),select='MP'),

		selectInput('facet_row','Facet Row',c(None='.',argv)),
		selectInput('facet_col','Facet Column',c(None='.',argv),select='Scenario')

		# selectInput("scn","Recruitment Scenario:",
		#             list("Stationary distribution"="S1",
		#                  "Cyclic (PDO-like)"="S2"))

		# selectInput("hcr","Harvest Control Rule:",
		#             list("Fixed Escapement"="fixE",
		#                  "Fixed Escapement Cap"="fixcap",
		#                  "Fixed Exploitation"="fixU",
		#                  "Fourty Ten Rule"="fourten",
		#                  "Thirty Twenty Rule"="thirtytwenty",
		#                  "Conditional Constant Catch"="ccc"),
		#             selected=TRUE,
		#             multiple=TRUE),

		# radioButtons("plotType","Output:",
		#              list("Biomass"="sbio",
		#                   "Total catch"="catch")),

		# numericInput("txtNSAM","Number of traces",3,min=0,max=10)
		),

	mainPanel(
		tabsetPanel(
			tabPanel("Plot",  h3(textOutput("caption")),  plotOutput("msePlot")),
			tabPanel("Summary", verbatimTextOutput("caption")),
			tabPanel("Table",   tableOutput("table"))
			)	          
	    )
))