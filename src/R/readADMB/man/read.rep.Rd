\name{read.rep}
\alias{read.rep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.rep(fn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
%%     ~~Describe \code{fn} here~~
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
function (fn) 
{
    options(warn = -1)
    ifile = scan(fn, what = "character", flush = TRUE, blank.lines.skip = FALSE, 
        quiet = TRUE)
    idx = sapply(as.double(ifile), is.na)
    vnam = ifile[idx]
    nv = length(vnam)
    A = list()
    ir = 0
    for (i in 1:nv) {
        ir = match(vnam[i], ifile)
        if (i != nv) 
            irr = match(vnam[i + 1], ifile)
        else irr = length(ifile) + 1
        dum = NA
        if (irr - ir == 2) 
            dum = as.double(scan(fn, skip = ir, nlines = 1, quiet = TRUE, 
                what = ""))
        if (irr - ir > 2) {
            dum = as.matrix(read.table(fn, skip = ir, nrow = irr - 
                ir - 1, fill = TRUE, row.names = NULL))
        }
        if (is.numeric(dum)) {
            A[[vnam[i]]] = dum
        }
    }
    options(warn = 0)
    return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
