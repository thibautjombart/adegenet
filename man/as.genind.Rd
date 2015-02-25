\encoding{UTF-8}
\name{genind constructor}
\alias{genind-methods}
\alias{genind}
\alias{as.genind}
\alias{is.genind}
\title{genind constructor}
\description{Constructor for \linkS4class{genind} objects.\cr
  The function \code{genind} creates a \linkS4class{genind} object from
  a matrix of allelic frequency where genotypes are in rows and alleles
  in columns. This table must have correct names for rows and
  columns.\cr

  The function \code{as.genind} is an alias for \code{genind} function.\cr
 
  \code{is.genind} tests if an object is a valid genind object.\cr

  Note: to get the manpage about \linkS4class{genind}, please type
  'class ? genind'.
}

\usage{
genind(tab,pop=NULL,prevcall=NULL, ploidy=2, type=c("codom","PA"))
as.genind(tab,pop=NULL,prevcall=NULL, ploidy=2, type=c("codom","PA"))
is.genind(x)
}

\arguments{
  \item{tab}{A table corresponding to the @tab slot of a genind object,
  with individuals in rows and alleles in columns.
    Its content depends on \code{type} (type of marker).\cr
    - 'codom': table contains allele frequencies (numeric values summing
  to 1).\cr
    - 'PA': table contains binary values, which indicate
  presence(1)/absence(0) of alleles.\cr
  }
  \item{pop}{a factor giving the population of each genotype in 'x'}
  \item{prevcall}{call of an object}
  \item{ploidy}{an integer indicating the degree of ploidy of
    the genotypes. Beware: 2 is not an integer, but as.integer(2) is.}
  \item{type}{ a character string indicating the type of
    marker: 'codom' stands for 'codominant' (e.g. microstallites,
    allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
  \item{x}{an object}
}
\value{
  For \code{genind} and \code{as.genind}, a genind object. For \code{is.genind}, a logical.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{
  \code{\linkS4class{genind} class}, and \code{\link{import2genind}} for
  importing from various types of file.\cr

  Related classes:\cr
  - \linkS4class{genpop} for storing data per populations\cr

  - \linkS4class{genlight} for an efficient storage of binary SNPs genotypes\cr
}
\examples{
data(nancycats)
nancycats@loc.names

# isolate one marker, fca23
obj <- seploc(nancycats)$"fca23"
obj
}
\keyword{manip}
