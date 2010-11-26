\name{fig1}
\alias{fig1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fig1(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function(A)
{	#Age-schedule information Length, Growth/maturity
	op=par(no.readonly=T)
	par(mfcol=c(2, 1))
	plot(A$age, A$la, xlab="Age", ylab="Length (cm)",type="o",  
	ylim=c(0,max(A$la)))
	plot(A$age, A$wa, xlab="Age", ylab="Weight (kg)/fecundity")
	lines(A$age, A$fa)
	legend("bottom",c("Weight-at-age", "Fecundity-at-age"), 
	lty=c(-1, 1), pch=c(1, -1), bty="n")
	
	fn=paste(prefix,"Fig1.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=4)
	par(op)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
