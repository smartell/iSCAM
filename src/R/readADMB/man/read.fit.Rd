\name{read.fit}
\alias{read.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
read parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.fit(ifile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ifile}{
%%     ~~Describe \code{ifile} here~~
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
function (ifile) 
{
    ret <- list()
    parfile <- as.numeric(scan(paste(ifile, ".par", sep = ""), 
        what = "", n = 16, quiet = TRUE)[c(6, 11, 16)])
    ret$nopar <- as.integer(parfile[1])
    ret$nlogl <- parfile[2]
    ret$maxgrad <- parfile[3]
    file <- paste(ifile, ".cor", sep = "")
    if (file.exists(file)) {
        lin <- readLines(file)
        ret$npar <- length(lin) - 2
        ret$logDetHess <- as.numeric(strsplit(lin[1], "=")[[1]][2])
        sublin <- lapply(strsplit(lin[1:ret$npar + 2], " "), 
            function(x) x[x != ""])
        ret$names <- unlist(lapply(sublin, function(x) x[2]))
        ret$est <- as.numeric(unlist(lapply(sublin, function(x) x[3])))
        ret$std <- as.numeric(unlist(lapply(sublin, function(x) x[4])))
        ret$cor <- matrix(NA, ret$npar, ret$npar)
        corvec <- unlist(sapply(1:length(sublin), function(i) sublin[[i]][5:(4 + 
            i)]))
        ret$cor[upper.tri(ret$cor, diag = TRUE)] <- as.numeric(corvec)
        ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)]
        ret$cov <- ret$cor * (ret$std \%o\% ret$std)
    }
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
