source(file.path("..","src","r","read.admb.r"));

readOutput <- function(d){
  return(read.admb(file.path(d,"om")));
}

dn <- dir(pattern="^[[:digit:]]");
nf <- paste(basename(getwd()),".Rdata",sep="");
sims <- lapply(dn,readOutput);
save(sims,file=nf)
