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

	output$gHCR <- renderPlot({
		x <- seq(0,1,b=0.01)
		y <- rep(0,length(x))
		sblim = input$hcr1
		sbthr = input$hcr2
		sbtar = input$hcr3
		frate = input$hcr4
		y[x<=sblim] <- 0
		y[x>=sbthr] <- frate
		ii          <- x>sblim&x<sbthr
		y[ii]       <- frate * (x[ii]-min(x[ii]))/(diff(range(x[ii])))

		plot(x,y,type="l")
	})

	output$gCt <- renderPlot({
		x <- seq(0,1,b=0.01)
		y <- rep(0,length(x))
		sblim = input$hcr1
		sbthr = input$hcr2
		sbtar = input$hcr3
		frate = input$hcr4
		y[x<=sblim] <- 0
		y[x>=sbthr] <- frate
		ii          <- x>sblim&x<sbthr
		y[ii]       <- frate * (x[ii]-min(x[ii]))/(diff(range(x[ii])))
		fe = y
		be = 1
		rk = 55
		s  = 0.85
		a  = rk*(1-s)
		b  = (rk-1)/be
		ce = rep(0,length=length(fe))
		bmsy = be*(-1+sqrt(k))/(k-1)
		msy  = (bmsy*(-1+s)*(-1+k)*(-be-bmsy)) / (be+(k-1)*bmsy)
		for(i in 1:length(fe))
		{
			be = ((a/(1-s+fe[i]))-1)/b
			ce[i] = be * fe[i]
			cat(fe[i]," ",be," ",ce[i],"\n")
			# be    = s*be +(a*be)/(1+b*be)-ce[i]
		}
		plot(x,ce,type="l")
	})

	# TAB 1
	mse_plotinput = function(...){
		tagList(
		    withTags(div(class='row-fluid',    
			  div(class='span3', selectInput('yy','Variable',argv[-1:-3],select='Biomass')),
			  div(class='span3', selectInput('facet_irow','Facet Row',c(None='.',argv[c(1,2)]))),
			  div(class='span3', selectInput('facet_jcol','Facet Column',c(None='.',argv[c(1,2)]),select='Scenario'))
			)),
			withTags(div(class='row-fluid',
			  div(class='span5', checkboxInput('savePNG','Save (CurrentImage.PNG)',value=FALSE)),
			  div(class='span3', selectInput('.FONTSIZE','Font size',c(10,12,16,18),select=18))
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
		if( is.null(input$yy)) return(NULL)
		if(input$yy == 'Biomass')
		{
			eb <- aes(ymin=Bt.lci,ymax=Bt.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
			if(length(input$iclr)==1)
			{
				d <- data()[,c(1:3,grep("btrc",names(qdf)))]
				d <- melt(d,id=c("Scenario","MP","Year"))
				# print(head(d))

				p <- p + geom_line(data=d,aes_string(x="Year",y='value',col='variable'))
			}
		}
		if(input$yy == 'Catch')
		{
			eb <- aes(ymin=Ct.lci,ymax=Ct.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
			if(length(input$iclr)==1)
			{
				d <- data()[,c(1:3,grep("ctrc",names(qdf)))]
				d <- melt(d,id=c("Scenario","MP","Year"))
				# print(head(d))

				p <- p + geom_line(data=d,aes_string(x="Year",y='value',col='variable'))
			}
		}
		if(input$yy == 'HarvestRate')
		{
			eb <- aes(ymin=Hr.lci,ymax=Hr.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)	
		}

		if(input$yy == 'Depletion')
		{
			eb <- aes(ymin=Dt.lci,ymax=Dt.uci)
			p <- p + geom_ribbon(eb,alpha=0.15)
			p <- p + geom_line(data=data(),aes_string(x="Year",y=input$yy),alpha=1)
			p <- p + ylim(0,1)
		}


		# print(length(input$iclr))
		if(length(input$iclr) > 1)
		{
			.LEGPOS <- 'top'
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
		print(input$.FONTSIZE)
		.FONTSIZE = as.integer(input$.FONTSIZE)
		print(p + theme_bw(.FONTSIZE) + theme( legend.position = .LEGPOS ))
		# print(p)
		if(input$savePNG)
		{
			png("CurrentImage.png",width=9.08,height=6.81,res=600,units="in")
			print(p + theme_bw(.FONTSIZE) + theme( legend.position = .LEGPOS ))
			dev.off()			
		}
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
		                      options=list(sizeAxis=FALSE))
		}
		if( !input$integrate )
		{  
			tmp <- data()
			tmp <- cbind(tmp,Procedure=tmp$MP)

			M1 <- gvisMotionChart(tmp, idvar="SP", timevar="Year",
		                      xvar="Depletion",yvar="Catch",
		                      sizevar="AAV",colorvar="Scenario"
		                      )
		}


		return(M1)
	})


	# Table panel
	output$viewDepletionTable <- renderTable({
		dt<-data()
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		# tmp <- cast(subset(mdf,variable=="Catch"),MP~Scenario,mean,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="Depletion"))
		return(tmp)
	})

	output$viewCatchTable <- renderTable({
		dt<-data()
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		# tmp <- cast(subset(mdf,variable=="Catch"),MP~Scenario,mean,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="Catch"))
		return(tmp)
	})

	output$viewAAVTable <- renderTable({
		dt<-data()
		print(head(dt))
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		#tmp <- cast(subset(mdf,variable=="AAV"),MP~Scenario,mean,na.rm=TRUE,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="AAV"))
		return(tmp)
	})

	output$viewClosedTable <- renderTable({
		dt<-data()
		print(head(dt))
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		#tmp <- cast(subset(mdf,variable=="AAV"),MP~Scenario,mean,na.rm=TRUE,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="Closures"))
		return(tmp)
	})

	output$viewSSBlimit <- renderTable({
		dt<-data()
		print(head(dt))
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		#tmp <- cast(subset(mdf,variable=="AAV"),MP~Scenario,mean,na.rm=TRUE,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="P(SSB<0.2)"))
		return(tmp)
	})

	output$viewSSBthreshold <- renderTable({
		dt<-data()
		print(head(dt))
		mdf <- melt(dt,id=c("Scenario","MP","SP","Year"))
		#tmp <- cast(subset(mdf,variable=="AAV"),MP~Scenario,mean,na.rm=TRUE,margins=TRUE)
		if(  input$integrate ) mar = "Scenario"
		if( !input$integrate ) mar = FALSE
		
		tmp <- dcast(mdf,MP~Scenario,mean,na.rm=TRUE,margins=mar,subset=.(variable=="P(SSB<0.3)"))
		return(tmp)
	})
})