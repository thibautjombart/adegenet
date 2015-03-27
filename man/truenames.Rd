\encoding{UTF-8}
\docType{methods}
\name{truenames}
\alias{truenames}
\alias{truenames-methods}
\alias{truenames,ANY-method}
\alias{truenames,genind-method}
\alias{truenames,genpop-method}
\title{ Restore true labels of an object}
\description{
  The function \code{truenames} returns some elements of an object
  (\linkS4class{genind} or \linkS4class{genpop}) using true names
  (as opposed to generic labels) for individuals, markers, alleles, and
  population.\cr

  Important: as of adegenet_1.5-0, these functions are deprecated as
  true labels are used whenever possible. Please use the function
  \code{\link{tab}} instead.
}
\usage{
\S4method{truenames}{genind}(x)
\S4method{truenames}{genpop}(x)
}
\arguments{
  \item{x}{a \linkS4class{genind} or a \linkS4class{genpop} object}
}
\value{
  If x\$pop is empty (NULL), a matrix similar to the x\$tab slot but
  with true labels.
  
  If x\$pop exists, a list with this matrix (\$tab) and a population
  vector with true names (\$pop).\cr
  
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{\code{\link{tab}}}