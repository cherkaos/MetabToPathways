\name{FindPathwayToCsv}
\alias{FindPathwayToCsv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FindPathwayToCsv(puchem, name, kegg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{puchem}{
%%     ~~Describe \code{puchem} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{kegg}{
%%     ~~Describe \code{kegg} here~~
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
function (puchem, name, kegg) 
{
    result <- FindPathway(puchem, name, kegg)
    for (i in 1:dim(result)[2]) {
        result[, i]$KeggPathway <- paste(result[, i]$KeggPathway, 
            collapse = "__")
        result[, i]$SmpdbPathway <- paste(result[, i]$SmpdbPathway, 
            collapse = "__")
        result[, i]$HmdbPathway <- paste(result[, i]$HmdbPathway, 
            collapse = "__")
        result[, i]$BiocycPathway <- paste(result[, i]$BiocycPathway, 
            collapse = "__")
    }
    write.csv(t(as.matrix(result)), "result.csv", row.names = TRUE, 
        col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
