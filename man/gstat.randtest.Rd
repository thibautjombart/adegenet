\encoding{UTF-8}
\name{gstat.randtest}
\alias{gstat.randtest}
\title{Goudet's G-statistic Monte Carlo test for genind object}
\description{
  The function \code{gstat.randtest} implements Goudet's
  G-statistic Monte Carlo test (\code{g.stats.glob}, package
  \code{hierfstat}) for \code{genind} object. \cr
  The output is an object of the class \code{randtest} (package
  \code{ade4}) from a \code{genind} object.\cr
    
  This procedure tests for genetic structuring of individuals using 3 different
  schemes (see details).
}
\usage{
gstat.randtest(x,pop=NULL, method=c("global","within","between"),
sup.pop=NULL, sub.pop=NULL, nsim=499)
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
  pop is seeked from x\@pop. Note that the term population refers in
  fact to any grouping of individuals'.}
  \item{method}{a character (if a vector, only first argument is kept)
    giving the method to be applied: 'global', 'within' or 'between'
    (see details).}
  \item{sup.pop}{a factor indicating any grouping of individuals at a
    larger scale than 'pop'. Used in 'within' method.}
  \item{sub.pop}{a factor indicating any grouping of individuals at a
    finer scale than 'pop'. Used in 'between' method.}
  \item{nsim}{number of simulations to be used for the randtest.}
 }
 \details{
   This G-statistic Monte Carlo procedure tests for population
   structuring at different levels. This is determined by the argument
   'method':\cr
   
   - "global": tests for genetic structuring given 'pop'.\cr
   - "within": tests for genetic structuring within 'pop' inside each
   'sup.pop' group (i.e., keeping sup.pop effect constant).\cr
   - "between": tests for genetic structuring between 'pop' keeping
   individuals in their 'sub.pop' groups (i.e., keeping sub.pop effect constant).

}
\value{Returns an object of the class \code{randtest} (package \code{ade4}).
}
\seealso{%\code{\link[hierfstat]{g.stats.glob}},
  \code{\link{fstat}},
  %\code{\link[hierfstat]{test.g}}%,\code{\link[hierfstat]{test.within}},
  %\code{\link[hierfstat]{test.between}},\code{\link[ade4]{as.randtest}},
  \code{\link{genind2hierfstat}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
if(require(hierfstat)){
# here the example of g.stats.glob is taken using gstat.randtest
data(gtrunchier)
x <- df2genind(X=gtrunchier[,-c(1,2)],pop=gtrunchier$Patch)

# test in hierfstat
gtr.test<- g.stats.glob(gtrunchier[,-1])
gtr.test

# randtest version
x.gtest <- gstat.randtest(x,nsim=99)
x.gtest
plot(x.gtest)

# pop within sup.pop test
gstat.randtest(x,nsim=99,method="within",sup.pop=gtrunchier$Locality)

# pop test with sub.pop kept constant
gstat.randtest(x,nsim=99,pop=gtrunchier$Locality,method="between",sub.pop=gtrunchier$Patch)
}
}
}
\keyword{manip}
\keyword{multivariate}