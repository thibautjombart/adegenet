\encoding{UTF-8}
\name{F statistics}
\alias{fstat}
\alias{FST}
\alias{fst}
\alias{pairwise.fst}
\title{F statistics for genind objects}
\description{
  \code{pairwise.fst} computes Nei's pairwise Fst between all pairs of
  populations using a \linkS4class{genind} object. Heretozygosities are
  weighted by group sizes (see details).

  The function \code{fstat} is a wrapper for \code{varcomp.glob} of the
  package \code{hierfstat}. For Fst, Fis and Fit, an alternative is
  offered by \code{Fst} from the \code{pagas} package (see example).
}
\usage{
fstat(x, pop=NULL, fstonly=FALSE)
pairwise.fst(x, pop=NULL, res.type=c("dist","matrix"), truenames=TRUE)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
    pop is seeked from \code{pop(x)}. Note that the term population refers in
    fact to any grouping of individuals'.}
  \item{res.type}{the type of result to be returned: a \code{dist} object, or a
    symmetric matrix}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
  \item{fstonly}{a logical stating whether only the Fst should be returned.}
}
\value{
  A vector, a matrix, or a dist object containing F statistics.
}
\seealso{\code{\link{Hs}}
}
\references{
  Nei, M. (1973) Analysis of gene diversity in subdivided
  populations. Proc Natl Acad Sci USA, 70: 3321-3323
}
\details{
  Let \eqn{A} and \eqn{B} be two populations of population sizes \eqn{n_A} and
  \eqn{n_B}, with expected heterozygosity (averaged over loci) \eqn{Hs(A)} and \eqn{Hs(B)},
  respectively. We denote \eqn{Ht} the expected heterozygosity of a population
  pooling \eqn{A} and \eqn{B}. Then, the pairwise \eqn{Fst} between \eqn{A} and \eqn{B} is computed
  as:\cr

 \eqn{  Fst(A,B) = \frac{(Ht - (n_A Hs(A) + n_B Hs(B))/(n_A + n_B) )}{Ht}} \cr
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
data(nancycats)

\dontrun{
## pairwise Fst
mat.fst <- pairwise.fst(nancycats, res.type="matrix")
mat.fst
}

## Fst, Fis, Fit
## using hierfstat
if(require(hierfstat)){
fstat(nancycats)
}

## using pegas
if(require(pegas)){
data(nancycats)

## conversion to pegas's format
as.loci(nancycats)

## use Fst from pegas
fsttab <- Fst(as.loci(nancycats))

## average over loci
apply(fsttab, 2, mean)
}
}
\keyword{multivariate}