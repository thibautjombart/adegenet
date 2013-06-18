\encoding{UTF-8}
\name{HWE.test.genind}
\alias{HWE.test.genind}

\title{Hardy-Weinberg Equilibrium test for multilocus data}
\description{The function \code{HWE.test} is a generic function to
  perform Hardy-Weinberg Equilibrium tests defined by the
  \code{genetics} package. adegenet proposes a method for \code{genind}
  objects.\cr
   
  The output can be of two forms:\cr
  - a list of tests (class \code{htest}) for each locus-population
  combinaison \cr
  - a population x locus matrix containing p-values of the tests  
}
\usage{
\method{HWE.test}{genind}(x,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{pop}{a factor giving the population of each individual. If NULL,
  pop is seeked from x\$pop.}
  \item{permut}{a logical passed to \code{HWE.test} stating whether
    Monte Carlo version (TRUE) should be used or not (FALSE, default).}
  \item{nsim}{number of simulations if Monte Carlo is used (passed to \code{HWE.test}).}
  \item{hide.NA}{a logical stating whether non-tested loci (e.g., when
    an allele is fixed) should be hidden in the results (TRUE, default)
    or not (FALSE).}
  \item{res.type}{a character or a character vector whose only first
    argument is considered giving the type of result to display. If
    "full", then a list of complete tests is returned. If "matrix",
    then a matrix of p-values is returned.}
}
\details{Monte Carlo procedure is quiet computer-intensive when large
  datasets are involved. For more precision on the performed test, read
  \code{HWE.test} documentation (\code{genetics} package).
}
\value{Returns either a list of tests or a matrix of p-values. In the
  first case, each test is designated by locus first and then by
  population. For instance if \code{res} is the "full" output of the
  function, then the test for population "PopA" at locus "Myloc" is
  given by res$Myloc$PopA. If \code{res} is a matrix of p-values,
  populations are in rows and loci in columns. P-values are given for
  the upper-tail: they correspond to the probability that an oberved
  chi-square statistic as high as or higher than the one observed
  occured under H0 (HWE).\cr  
  
  In all cases, NA values are likely to appear in fixed loci, or
  entirely non-typed loci.
}
\seealso{\code{\link[genetics]{HWE.test}},\code{\link[stats]{chisq.test}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)
obj <- nancycats
if(require(genetics)){
obj.test <- HWE.test(obj)

# pvalues matrix to have a preview
HWE.test(obj,res.type="matrix")

#more precise view to...
obj.test$fca90$P10
}
}
}
\keyword{manip}
\keyword{multivariate}