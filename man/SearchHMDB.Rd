\name{SearchHMDB}
\alias{SearchHMDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SearchHMDB(hmdb, info = "Pathways")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmdb}{
%%     ~~Describe \code{hmdb} here~~
}
  \item{info}{
%%     ~~Describe \code{info} here~~
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
function (hmdb, info = "Pathways") 
{
    switch(info, Pathways = "pathway/name", Biocyc = "biocyc_id", 
        Synonyms = "synonym", CellularLocation = "cellular_location", 
        BiofluidLocation = "biofluid", TissueLocation = "tissue")
    if (!hmdb == "error" || !is.null(hmdb)) {
        url <- paste0("http://www.hmdb.ca/metabolites/", hmdb, 
            ".xml")
        xmlhttp <- getURL(url)
        if (grepl("metabolite", xmlhttp, ignore.case = TRUE) == 
            TRUE) {
            if (is.null(xmlhttp)) {
                return("error")
            }
            else {
                doc <- (xmlParse(xmlhttp))
                src <- tryCatch(xpathSApply(doc, paste0("//", 
                  info), xmlValue), error = function(e) {
                  NULL
                })
                if (length(src) > 0) {
                  src <- unlist(src)
                  return(src)
                }
                else (return("error"))
            }
        }
        else {
            return("error")
        }
    }
    else {
        return("error")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
