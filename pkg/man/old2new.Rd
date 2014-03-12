\name{old2new}
\docType{methods}
\alias{old2new}
\alias{old2new,ANY-method}
\alias{old2new-methods}
\alias{old2new,genind-method}
\alias{old2new,genpop-method}
\title{Convert objects with obsolete classe into new objects}
\description{
  Adegenet classes changed from S3 to S4 types starting from version
  1.1-0. \code{old2new} has two methods for genind and genpop objects,
  so that old adegenet objects can be retrieved and used in recent versions.
}
\usage{
\S4method{old2new}{genind}(object)
\S4method{old2new}{genpop}(object)
}
\arguments{
  \item{object}{a genind or genpop object in S3 version, i.e. prior
    adegenet\_1.1-0}
}
\details{
  Optional content but \code{$pop} and \code{$pop.names} will not be
  converted. These are to be coerced into a list and set in the
  \code{@other} slot of the new object.  
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\keyword{methods}
\keyword{classes}
\keyword{manip}