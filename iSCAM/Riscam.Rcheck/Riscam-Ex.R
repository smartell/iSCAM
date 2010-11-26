pkgname <- "Riscam"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Riscam')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("A")
### * A

flush(stderr()); flush(stdout())

### Name: A
### Title: Pacific hake results
### Aliases: A
### Keywords: datasets

### ** Examples

data(A)
## maybe str(A) ; plot(A) ...



cleanEx()
nameEx("read.admb")
### * read.admb

flush(stderr()); flush(stdout())

### Name: read.admb
### Title: Data input
### Aliases: read.admb read.fit read.rep read.psv
### Keywords: ADMB iscam

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(ifile)
{	
	ret=read.fit(ifile)
	
	fn=paste(ifile,'.rep', sep='')
	A=read.rep(fn)
	A$fit=ret
	
	pfn=paste(ifile,'.psv',sep='')
	if(file.exists(pfn))
		A$mc=read.psv(pfn)
	
	return(A)
}



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
