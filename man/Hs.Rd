\encoding{UTF-8}
\name{Hs}
\alias{Hs}
\title{Expected heterozygosity}
\description{
  This function computes the expected heterozygosity (Hs) within
  populations of a \linkS4class{genpop} object. This function is
  available for codominant markers (\code{@type="codom"}) only. Hs is
  commonly used for measuring within population genetic diversity (and
  as such, it still has sense when computed from haploid data).
}
\usage{
Hs(x, truenames=TRUE)
}
\arguments{
  \item{x}{an object of class \linkS4class{genpop}.}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
}
\details{
  Let \emph{m(k)} be the number of alleles of locus \emph{k}, with a
  total of \emph{K} loci. We note \eqn{f_i} the allele frequency of
  allele \emph{i} in a given population. Then, \eqn{Hs} is given for a
  given population by:\cr

  \eqn{\frac{1}{K} \sum_{k=1}^K (1 - \sum_{i=1}^{m(k)} f_i^2)} \cr
}
\value{A vector of Hs values (one value per population).}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)
Hs(genind2genpop(nancycats))
}
}
\keyword{multivariate}