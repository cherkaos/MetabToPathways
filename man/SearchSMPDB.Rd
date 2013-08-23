\name{SearchSMPDB}
\alias{SearchSMPDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SearchSMPDB(id, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (id, type) 
{
    position <- vector()
    if (type == "name") {
        position <- Mmatch(id, DBSMPDB[, 4])
    }
    if (type == "hmdb") {
        position <- Mmatch(id, DBSMPDB[, 3])
    }
    if (type == "kegg") {
        position <- Mmatch(id, DBSMPDB[, 5])
    }
    if (!length(position) == 0) {
        Name <- fixlc(DBSMPDB[position[1], 4])
        HMDB <- fixlc(DBSMPDB[position[1], 3])
        KeggID <- fixlc(DBSMPDB[position[1], 5])
        SmpdbPathway <- fixlc(DBSMPDB[position, 2])
        return(list(Name = Name, HMDB = HMDB, KeggID = KeggID, 
            SmpdbPathway = SmpdbPathway))
    }
    else (return("error"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
