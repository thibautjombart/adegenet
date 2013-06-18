\name{isPoly-methods}
\docType{methods}
\alias{isPoly}
\alias{isPoly-methods}
\alias{isPoly,genind-method}
\alias{isPoly,genpop-method}
\title{ Assess polymorphism in genind/genpop objects }
\description{
  The simple function \code{isPoly} can be used to check which loci are
  polymorphic, or alternatively to check which alleles give rise to polymorphism.
}
\usage{
\S4method{isPoly}{genind}(x, by=c("locus","allele"), thres=1/100)
\S4method{isPoly}{genpop}(x, by=c("locus","allele"), thres=1/100)
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{by}{a character being "locus" or "allele", indicating whether
    results should indicate polymorphic loci ("locus"), or alleles giving
    rise to polymorphism ("allele").}
  \item{thres}{a numeric value giving the minimum frequency of an allele
  giving rise to polymorphism (defaults to 0.01).}
 }
 \value{
   A vector of logicals.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)
isPoly(nancycats,by="loc", thres=0.1)
isPoly(nancycats[1:3],by="loc", thres=0.1)
genind2df(nancycats[1:3])
}
}
\keyword{methods}
\keyword{manip}
