\encoding{UTF-8}
\name{df2genind}
\alias{df2genind}
\alias{genind2df}
\title{ Convert a data.frame of genotypes to a genind object, and conversely.}
\description{
  The function \code{df2genind} converts a data.frame (or a matrix) into
  a \linkS4class{genind} object. The data.frame must meet the following
  requirements:\cr
  - genotypes are in row (one row per genotype)\cr
  - markers are in columns\cr
  - each element is a string of characters coding alleles with or without
  separator. If no separator is used, the function tries to find how
  many characters code each genotypes at a locus, but it is safer to state it
  (\code{ncode} argument). Uncomplete strings are filled with "0" at the
  begining. \cr

  The function \code{genind2df} converts a \linkS4class{genind} back to
 such a data.frame. Alleles of a given locus can be coded as a single
 character string (with specified separators), or provided on different
 columns (see \code{oneColPerAll} argument).
 }
\usage{
df2genind(X, sep=NULL, ncode=NULL, ind.names=NULL, loc.names=NULL,
 pop=NULL, missing=NA, ploidy=2, type=c("codom","PA"))
genind2df(x,pop=NULL, sep="", usepop=TRUE, oneColPerAll=FALSE)
}
\arguments{
  \item{X}{a matrix or a data.frame (see decription)}
  \item{sep}{a character string separating alleles. See details.}
  \item{ncode}{an optional integer giving the number of characters used for coding one
    genotype at one locus. If not provided, this is determined from data.}
  \item{ind.names}{an optional character vector giving the individuals
    names; if NULL, taken from rownames of X.}
  \item{loc.names}{an optional character vector giving the markers
    names; if NULL, taken from colnames of X.}
  \item{pop}{an optional factor giving the population of each individual.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{ploidy}{an integer indicating the degree of ploidy of the genotypes.}
  \item{type}{ a character string indicating the type of
      marker: 'codom' stands for 'codominant' (e.g. microstallites,
      allozymes); 'PA' stands for 'presence/absence' markers (e.g. AFLP, RAPD).}
  \item{x}{a \linkS4class{genind} object}
  \item{usepop}{a logical stating whether the population (argument \code{pop}
    or \code{x@pop} should be used (TRUE, default) or not (FALSE)).}
  \item{oneColPerAll}{a logical stating whether alleles of one locus should be
    provided on separate columns (TRUE) rather than as a single character
    string (FALSE, default).}
}
   
\details{=== There are 3 treatments for missing values === \cr
  - NA: kept as NA.\cr
  
  - 0: allelic frequencies are set to 0  on all alleles of the concerned
  locus. Recommended for a PCA on compositionnal data.\cr
  
  - "mean": missing values are replaced by the mean frequency of the
  corresponding allele, computed on the whole set of
  individuals. Recommended for a centred PCA.\cr\cr

  === Details for the \code{sep} argument ===\cr
  this character is directly used in reguar expressions like
  \code{gsub}, and thus require some characters to be preceeded by
  double backslashes. For instance, "/" works but "|" must be coded as "\\\\|".
}
\value{an object of the class \linkS4class{genind} for \code{df2genind};
  a matrix of biallelic genotypes for \code{genind2df}}
\seealso{\code{\link{import2genind}}, \code{\link{read.genetix}},
  \code{\link{read.fstat}}, \code{\link{read.structure}}}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## simple example
df <- data.frame(locusA=c("11","11","12","32"),
locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
row.names(df) <- .genlab("genotype",4)
df

obj <- df2genind(df, ploidy=2)
obj
truenames(obj)

## converting a genind as data.frame 
genind2df(obj)
genind2df(obj, sep="/")
genind2df(obj, oneColPerAll=TRUE)
}
\keyword{manip}
