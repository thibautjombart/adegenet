\encoding{UTF-8}
\name{genpop constructor}
\alias{genpop-methods}
\alias{genpop}
\alias{as.genpop}
\alias{is.genpop}
\title{genpop constructor}
\description{Constructor for \linkS4class{genpop} objects.\cr
  The function \code{genpop} creates a \linkS4class{genpop} object from
  a matrix of alleles counts where genotypes are in rows and alleles
  in columns. This table must have correct names for rows and
  columns.\cr

  The function \code{as.genpop} is an alias for \code{genpop} function.\cr
 
  \code{is.genpop} tests if an object is a valid genpop object.\cr

   Note: to get the manpage about \linkS4class{genpop}, please type
  'class ? genpop'.
}
\usage{
genpop(tab,prevcall=NULL, ploidy=as.integer(2), type=c("codom","PA"))
as.genpop(tab, prevcall=NULL, ploidy=as.integer(2), type=c("codom","PA"))
is.genpop(x)
}
\arguments{
  \item{tab}{a pop x alleles matrix which terms are numbers of alleles,
    i.e. like in a genpop object}  
  \item{prevcall}{call of an object}
  \item{ploidy}{an integer indicating the degree of ploidy of
    the genotypes. Beware: 2 is not an integer, but as.integer(2) is.}
  \item{type}{a character string indicating the type of
    marker: 'codom' stands for 'codominant' (e.g. microstallites,
    allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP, RAPD).}
  \item{x}{an object}
}
\value{
  For \code{genpop} and \code{as.genpop}, a genpop object. For \code{is.genpop}, a logical.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{
  \code{\linkS4class{genpop} class}, and \code{\link{genind2genpop}} for
  conversion from a genind to a genpop object.
}
\examples{
data(nancycats)
obj <- genind2genpop(nancycats)

# isolate one locus, fca77
obj <- seploc(obj)$"fca77"
obj
}
\keyword{manip}
