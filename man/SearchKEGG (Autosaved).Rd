\name{SearchKEGG}
\alias{SearchKEGG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to find a compound information using his id by looking in database which has kegg pathways
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SearchKEGG(id, type)
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
1. 4Creek DJ, Jankevics A, Burgess KE, Breitling R, Barrett MP. IDEOM: an Excel
interface for analysis of LC-MS-based metabolomics data. Bioinformatics. 2012 Apr
1;28(7):1048-9. 

}
\author{Sarah Cherkaoui}
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
        position <- Mmatch(id, DBKEGG[, 3])
    }
    if (type == "pubchem") {
        position <- Mmatch(id, DBKEGG[, 9])
    }
    if (type == "hmdb") {
        position <- Mmatch(id, DBKEGG[, 11])
    }
    if (type == "kegg") {
        position <- Mmatch(id, DB1DBKEGG[, 8])
    }
    if (type == "biocyc") {
        position <- Mmatch(id, DBKEGG[, 12])
    }
    if (type == "inchikey") {
        position <- Mmatch(id, DBKEGG[, 15])
    }
    if (length(position) > 0) {
        Name <- fixlc(DBKEGG[position, 3])
        Pubchem <- fixlc(DBKEGG[position, 9])
        HMDB <- fixlc(DBKEGG[position, 11])
        KeggID <- fixlc(DBKEGG[position, 8])
        BioCyc <- fixlc(DBKEGG[position, 12])
        InchiKey <- fixlc(DBKEGG[position, 15])
        Map <- fixlc(DBKEGG[position, 5])
        KeggPathway <- fixlc(strsplit(DBKEGG[position, 6], "__"))
        KeggPathID <- fixlc(DBSMPDB[position, 7])
        return(list(Name = Name, Pubchem = Pubchem, HMDB = HMDB, 
            KeggID = KeggID, BioCyc = BioCyc, InchiKey = InchiKey, 
            Map = Map, KeggPathway = KeggPathway, KeggPathID = KeggPathID))
    }
    else (return("error"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
