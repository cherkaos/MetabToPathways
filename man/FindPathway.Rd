\name{FindPathway}
\alias{FindPathway}
\title{Find, for every compound, his pathways using multiple databases}
\description{Function to find, for every compound, his pathways using the 4 databases, using soecific id. Returns different id and pathway of every compound to get the best coverage.}
\usage{
FindPathway(name, puchem, kegg)
}
\arguments{
  \item{name}{List of name of every compounds to query.Should be the same length as puchem and kegg. Accept names that are blank.
}
  \item{puchem}{List of puchem id of every compounds to query.Should be the same length as name and kegg. Accept puchem id that are blank.
}
  \item{kegg}{List of kegg id of the compounds to query.Should be the same length as puchem and name. Accept kegg id that are blank.
}
}
\details{
	The same position in the 3 lists should repesent the same compound but with different id. A compound can be find with his puchem, name or kegg. Can accept id's that are blanks.
}
\value{
Returns a list of information for every compound.
\item{Name}{Name of compound as character}
\item{Pubchem}{Pubchem id of compound as character}
\item{HMDB}{HMDB id of compound as character}
\item{KeggID}{Kegg id of compound as character}
\item{BioCyc}{Biocy id of compound as character}
\item{InchiKey}{Inchikey of compound as character}
\item{Map}{Map of the 1st pathway of the list}
\item{KeggPathway}{List of pathway from Kegg where the compound is found}
\item{KeggPathID}{Kegg id of the 1st pathway of the list}
\item{SmpdbPathway}{List of pathway from SMPDB where the compound is found}
\item{HmdbPathway}{List of pathway from HMDB where the compound is found}
\item{BiocycPathway}{List of pathway from Biocyc where the compound is found}
}
\references{
Kegg: Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic
Acids Res. 2000 Jan 1;28(1):27-30.

SMPDB: Wishart DS, Frolkis A, Knox C, et al., SMPDB: The Small Molecule Pathway Database. Nucleic Acids Res. 2010 Jan;38(Database issue):D480-7

HMDB: Wishart DS, Tzur D, Knox C, et al., HMDB: the Human Metabolome Database. Nucleic Acids Res. 2007 Jan;35(Database issue):D521-6

Biocyc: Caspi R, Altman T, Dale JM, Dreher K, Fulcher CA, Gilham F, Kaipa P, Karthikeyan AS, Kothari A,Krummenacker M, Latendresse M, Mueller LA, Paley S, Popescu L, Pujar A, Shearer AG, Zhang P, Karp PD. The MetaCyc database of metabolic pathways and enzymes and the BioCyc collection of pathway/genome databases. Nucleic Acids Res. 2010 Jan;38(Database issue):D473-9.

}
\author{Sarah Cherkaoui}
\note{
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


## The function is currently defined as
FindPathway (list$name, list$puchem, list$kegg) 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pathways }
\keyword{ ~kegg }
\keyword{ ~smpdb }
\keyword{ ~hmdb }
\keyword{ ~biocyc }

