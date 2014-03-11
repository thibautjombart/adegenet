\encoding{UTF-8}
\name{propShared}
\alias{propShared}
\title{Compute proportion of shared alleles}
\description{The function \code{propShared} computes the proportion of
  shared alleles in a set of genotypes (i.e. from a \linkS4class{genind}
  object). Current implementation works for any level of ploidy.
}
\usage{
propShared(obj)
}
\arguments{
  \item{obj}{a \linkS4class{genind} object.}
 }
 \details{
   Computations of the numbers of shared alleles are done in C.
   Proportions are computed from all available data, i.e. proportion can
   be computed as far as there is at least one typed locus in common
   between two genotypes.
}
\value{Returns a matrix of proportions}
\seealso{\code{\link{dist.genpop}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{

\dontrun{
## make a small object
data(microbov)
obj <- microbov[1:5,microbov@loc.fac \%in\% c("L01","L02")]

## verify results
propShared(obj)
genind2df(obj,sep="|")

## Use this similarity measure inside a PCoA
##  ! This is for illustration only !
## the distance should be rendered Euclidean before
## (e.g. using cailliez from package ade4).
matSimil <- propShared(microbov)
matDist <- exp(-matSimil)
D <- cailliez(as.dist(matDist))
pcoa1 <- dudi.pco(D,scannf=FALSE,nf=3)
s.class(pcoa1$li,microbov$pop,lab=microbov$pop.names)

}
}
\keyword{manip}
\keyword{multivariate}
