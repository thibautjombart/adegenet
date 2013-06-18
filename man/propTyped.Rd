\name{propTyped-methods}
\docType{methods}
\alias{propTyped}
\alias{propTyped-methods}
\alias{propTyped,genind-method}
\alias{propTyped,genpop-method}
\title{ Compute the proportion of typed elements }
\description{
  The generic function \code{propTyped} is devoted to investigating the
  structure of missing data in adegenet objects.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
  objects. They can return the proportion of available
  (i.e. non-missing) data per individual/population, locus, or the
  combination of both in with case the matrix indicates which entity
  (individual or population) was typed on which locus.
}
\usage{
\S4method{propTyped}{genind}(x,  by=c("ind","loc","both"))
\S4method{propTyped}{genpop}(x,  by=c("pop","loc","both"))
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{by}{a character being "ind","loc", or "both" for
  \linkS4class{genind} object and "pop","loc", or "both" for
  \linkS4class{genpop} object. It specifies whether proportion of typed
  data are provided by entity ("ind"/"pop"), by locus ("loc") or both
  ("both"). See details.}
 }
 \value{
   A vector of proportion (when \code{by} equals "ind", "pop", or
  "loc"), or a matrix of binary data (when \code{by} equals "both")
}
\details{
  When \code{by} is set to "both", the result is a matrix of binary data
  with entities in rows (individuals or populations) and markers in
  columns. The values of the matrix are 1 for typed data, and 0 for NA.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)
propTyped(nancycats,by="loc")
propTyped(genind2genpop(nancycats),by="both")
}
}
\keyword{methods}
\keyword{manip}
