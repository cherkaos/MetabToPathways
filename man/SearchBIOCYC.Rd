\name{SearchBIOCYC}
\alias{SearchBIOCYC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SearchBIOCYC(biocyc, organism = "Multiorganism")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biocyc}{
%%     ~~Describe \code{biocyc} here~~
}
  \item{organism}{
%%     ~~Describe \code{organism} here~~
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
function (biocyc, organism = "Multiorganism") 
{
    switch(organism, Multiorganism = "META", Human = "HUMAN", 
        Ecoli = "ECOLI", Basillus = "BSUB", )
    if (!OBJECTID == "error") {
        url <- paste0("http://websvc.biocyc.org/apixml?fn=pathways-of-compound&id=", 
            organism, ":", biocyc)
        xmlhttp <- getURL(url)
        if (grepl("ptools-xml", xmlhttp, ignore.case = TRUE) == 
            TRUE) {
            if (is.null(xmlhttp)) {
                return("error")
            }
            else {
                doc <- xmlParse(xmlhttp)
                src <- tryCatch(xpathApply(doc, "//common-name", 
                  xmlValue), error = function(e) {
                  NULL
                })
                html <- c("&", "<i>", "</i>")
                for (i in 1:length(html)) {
                  src <- sub(html[i], "", src, fixed = FALSE)
                }
                src <- as.vector(fixlc(src))
                return(src)
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
