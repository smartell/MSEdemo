require(googleVis)
require(reshape)
require(reshape2)
require(reldist)
require(shiny)
require(ggplot2)
require(plyr)
load("QDF.Rdata")
raw.data <- qdf
mi <- c(-4,-6,-7,-9,-10,-12,-13,-14,-15,-17,-18)
mi <- c(mi,-grep("btrc",names(raw.data)),
           -grep("SP",names(raw.data)),
           -grep("ctrc",names(raw.data)),
           -grep("SSB<0.2",names(raw.data)),
           -grep("SSB<0.3",names(raw.data)) )
argv      <- names(raw.data)[mi]
.LEGPOS   <- 'none'
.FONTSIZE <- 12