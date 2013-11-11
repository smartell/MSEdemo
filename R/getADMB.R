#getADMB.R
source(file.path("../..","R","read.admb.r"));

A <- read.admb("om")
nf <- paste(basename(getwd()),".Rdata",sep="");
save(A,file=nf)