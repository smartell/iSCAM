source(file.path("../../../../dist/R/lib","read.admb.r"));
readOutput <- function(d){
  return(read.admb(file.path(d,"iscam")));
}
wd <- getwd()
nf <- paste(basename(getwd()),".rda",sep="");
setwd("./reps")
dn <- dir(pattern="^[[:digit:]]");
sims <- lapply(dn,readOutput);
setwd(wd)
save(sims,file=nf)