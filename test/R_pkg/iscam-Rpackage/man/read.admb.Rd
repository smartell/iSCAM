\name{read.admb}
\alias{read.admb}
\alias{read.fit}
\alias{read.rep}
\alias{read.psv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Data input
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reads the parameter files, correlation file, and report file outputs from an
ADMB program.
}
\usage{
read.admb(ifile)
read.fit(ifile)
read.rep(fn)
read.psv(fn,nsamples=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ifile}{the file name (without the extension) of the ADMB program}
  \item{fn}{the file name with the extension, use for read.rep & read.psv}
  \item{nsamples}{the number of posterior samples to read in from the binary *.psv file}
%%     ~~Describe \code{ifile} here~~

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function read.admb reads the contents of the ADMB output parameter file (*.par)
and the correlation file (*.cor) as well as the report file (*.rep).  The function
read.rep, just reads the report file only.

The read.psv function reads the binary output from the -mcsave command line option 
in ADMB.  The default number of samples (rows) is 10,000. Note that this function is
also called from read.admb if, and only if, the *.psv file exists.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
S. Martell, 
Anders Nielsen, and the original contribution came from George Watters
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ADMB }
\keyword{ iscam }% __ONLY ONE__ keyword per line
