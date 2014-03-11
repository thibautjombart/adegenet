\name{as methods in adegenet}
\docType{methods}
\alias{as-method}
\alias{as,genind,data.frame-method}
\alias{as,genpop,data.frame-method}
\alias{as,genind,matrix-method}
\alias{as,genpop,matrix-method}
\alias{as,genind,genpop-method}
\alias{ktab-class}
\alias{as,genind,ktab-method}
\alias{as,genpop,ktab-method}
\alias{coerce,genind,data.frame-method}
\alias{coerce,genpop,data.frame-method}
\alias{coerce,genind,matrix-method}
\alias{coerce,genpop,matrix-method}
\alias{coerce,genind,genpop-method}
\alias{coerce,genind,ktab-method}
\alias{coerce,genpop,ktab-method}
\alias{as.data.frame.genind}
\alias{as.data.frame.genpop}
\alias{as.matrix.genind}
\alias{as.matrix.genpop}
\alias{as.genpop.genind}
\alias{as.ktab.genind}
\alias{as.ktab.genpop}
\title{Converting genind/genpop objects to other classes}
\description{
 These S3 and S4 methods are used to coerce \linkS4class{genind} and
 \linkS4class{genpop} objects to matrix-like objects. In most cases,
 this is equivalent to calling the \code{@tab} slot. An exception to
 this is the convertion to \code{\link[ade4]{ktab}} objects used in the
 ade4 package as inputs for K-tables methods (e.g. Multiple Coinertia
 Analysis).\cr
}
\section{Usage}{
\code{as(object, Class)}
}
\section{Arguments}{
\describe{
  \item{\code{object}}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
  \item{\code{Class}}{the name of the class to which the object should
    be coerced, for instance \code{"data.frame"} or \code{"matrix"}.}
}
}
\section{Methods}{
\describe{
  \item{coerce}{from one object class to another using
    \code{as(object,"Class")}, where the \code{object} is of the old
    class and the returned object is of the new class \code{"Class"}.}
}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
\dontrun{
data(microbov)
x <- na.replace(microbov,method="0")
as(x[1:3],"data.frame")

## dudi functions attempt to convert their first argument
## to a data.frame; so they can be used on genind/genpop objects.
## perform a PCA
pca1 <- dudi.pca(x, scale=FALSE, scannf=FALSE)
pca1

x <- genind2genpop(microbov,miss="chi2")
x <- as(x,"ktab")
class(x)
## perform a STATIS analysis
statis1 <- statis(x, scannf=FALSE)
statis1
plot(statis1)

}
}
\keyword{methods}