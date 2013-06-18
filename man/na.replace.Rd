\name{na.replace-methods}
\docType{methods}
\alias{na.replace}
\alias{na.replace-methods}
\alias{na.replace,genind-method}
\alias{na.replace,genpop-method}
\title{ Replace missing values (NA) from an object}
\description{
  The generic function \code{na.replace} replaces NA in an object by
  appropriate values as defined by the argument \code{method}.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop} objects.
}
\usage{
\S4method{na.replace}{genind}(x,method, quiet=FALSE)
\S4method{na.replace}{genpop}(x,method, quiet=FALSE)
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{method}{a character string: can be "0" or "mean" for
    \linkS4class{genind} objects, and "0" or "chi2" for
    \linkS4class{genpop} objects.}
  \item{quiet}{logical stating whether a message should be
    printed (TRUE,default) or not (FALSE).} 
 }
 \value{
   A \linkS4class{genind} and \linkS4class{genpop} object without
   missing values.
}
\details{The argument "method" have the following effects:\cr
  - "0": missing values are set to "0". An entity (individual or
  population) that is not type on a locus has zeros for all alleles of
  that locus.\cr

    - "mean": missing values are set to the mean of the concerned
    allele, computed on all available observations (without distinction
    of population).\cr

    - "chi2": if a population is not typed for a marker, the
    corresponding count is set to that of a theoretical count in of a
    Chi-squared test. This is obtained by the product of the sums of
    both margins divided by the total number of alleles.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)

obj1 <- genind2genpop(nancycats)
# note missing data in this summary
summary(obj1)

# NA are all in pop 17 and marker fca45
which(is.na(obj1$tab),TRUE)
truenames(obj1)[17,]

# replace missing values
obj2 <- na.replace(obj1,"chi2")
obj2$loc.names

# missing values where replaced
truenames(obj2)[,obj2$loc.fac=="L4"]
}
}
\keyword{methods}
\keyword{manip}
