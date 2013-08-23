\name{TopPathwayKegg}
\alias{TopPathwayKegg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TopPathwayKegg(name, pubchem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{pubchem}{
%%     ~~Describe \code{pubchem} here~~
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
\author{Sarah Cherkaoui
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
function (name, pubchem) 
{
    Pathways <- matrix(ncol = 3)
    for (i in 1:length(puchem)) {
        s1 <- SearchKEGG(puchem[i], "pubchem")
        if (!is.list(s1)) {
            s1 <- SearchKEGG(name[i], "name")
        }
        if (is.list(s1)) {
            KeggPathway <- s1$KeggPathway
            for (j in 1:length(KeggPathway)) {
                if (is.character(KeggPathway[j]) && length(KeggPathway[j]) > 
                  0) {
                  s2 <- Mmatch(KeggPathway[j], Pathways[, 1])
                  if (length(s2) > 0) {
                    position <- s2[1]
                    Pathways[position, 2] <- as.numeric(Pathways[position, 
                      2]) + 1
                  }
                  else {
                    Pathways <- rbind(Pathways, c(KeggPathway[j], 
                      "1", s1$Name))
                  }
                }
            }
        }
    }
    Pathways <- Pathways[!is.na(Pathways[, 1]), ]
    Pathways <- Pathways[order(Pathways[, 2], decreasing = TRUE), 
        ]
    colnames(Pathways) <- c("Pathway", "Number", "Example")
    return(Pathways)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
