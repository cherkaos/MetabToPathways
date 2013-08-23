\name{TopPathwayKeggToLatex}
\alias{TopPathwayKeggToLatex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TopPathwayKeggToLatex(puchem, name, top = 20, index = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{puchem}{
%%     ~~Describe \code{puchem} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
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
function (puchem, name, top = 20, index = 1) 
{
    Pathways <- TopPathKegg(puchem, name)
    Pathways <- Pathways[1:top, ]
    output <- capture.output(stargazer(as.data.frame(Pathways), 
        title = "Pathway Analysis", font.size = "scriptsize", 
        type = "latex", summary = FALSE, label = paste0("table:tab", 
            index + 1, seq = ""), table.placement = "h"))
    cat(paste(output, collapse = "\n"), "\n", file = "pathways.tex", 
        append = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
