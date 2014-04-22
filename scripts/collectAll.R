source(file.path("../../../..dist/R/lib","read.admb.r"));

readOutput <- function(d){
	cat(d)
  return(read.admb(file.path(d,"om")));
}

dn <- dir(pattern="^[[:digit:]]");
nf <- paste(basename(getwd()),".Rdata",sep="");
sims <- lapply(dn,readOutput);
save(sims,file=nf)
