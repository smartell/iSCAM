## Libraries 
## This script is used to construct the *.Rdata files 
## for use in Management Strategy Evaluation.
##
## It is called from the makefile using:
## 	make collect
source(file.path("../../../dist/R/lib","read.admb.r"));
wd <- getwd()


readOutput <- function(d)
{
	A <- read.admb(file.path(d,"iscam"))
	B <- read.rep(file.path(d,"milka.rep"))
	C <- c(A,B)
	return( C );
}

mse.dirs <- dir("../DATA",recursive=FALSE,pattern="^mse_",full.names=TRUE)

for(dir in mse.dirs)
{
	setwd(dir)
	file.name <- paste(basename(dir),".Rdata",sep="")
	dn <- dir(pattern="^[[:digit:]]");
	sims <- lapply(dn,readOutput);
	save(sims,file=file.name)
	setwd(wd)
}

