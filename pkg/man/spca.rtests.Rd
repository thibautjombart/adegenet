\name{global.rtest}
\alias{global.rtest}
\alias{local.rtest}
\title{Global and local tests}
\description{
  These two Monte Carlo tests are used to assess the existence of global
  and local spatial structures. They can be used as an aid to interprete
  global and local components of spatial Principal Component Analysis
  (sPCA).\cr
  
  They rely on the decomposition of a data matrix X into global and local
  components using multiple regression on Moran's Eigenvector Maps (MEMs).
  They require a data matrix (X) and a list of weights derived from a
  connection network. X is regressed onto global MEMs (U+) in the global
  test and on local ones (U-) in the local test. One mean \eqn{R^2}{R^2}
  is obtained for each MEM, the k highest being summed to form the test
  statistic.
  
  The reference distribution of these statistics are obtained by
  randomly permuting the rows of X.
}
\usage{
global.rtest(X, listw, k = 1, nperm = 499)
local.rtest(X, listw, k = 1, nperm = 499)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data matrix, with variables in columns}
  \item{listw}{a list of weights of class \code{listw}. Can be obtained
    easily using the function \code{chooseCN}.}
  \item{k}{integer: the number of highest \eqn{R^2}{R^2} summed to form the test statistics}
  \item{nperm}{integer: the number of randomisations to be performed.}
}
\details{
 This test is purely R code. A C or C++ version will be developped soon.
}
\value{
An object of class \code{randtest}.
}
\references{
  Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
  Revealing cryptic spatial patterns in genetic variability by a new
  multivariate method. \emph{Heredity}, \bold{101}, 92--103.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{ \code{\link{chooseCN}}, \code{\link{spca}}, \code{\link{monmonier}}}
\examples{
\dontrun{
 data(sim2pop)
if(require(spdep)){
cn <- chooseCN(sim2pop@other$xy,ask=FALSE,type=1,plot=FALSE,res="listw")

# global test
Gtest <- global.rtest(sim2pop@tab,cn)
Gtest

# local test
Ltest <- local.rtest(sim2pop@tab,cn)
Ltest
}
}
}
\keyword{multivariate}
\keyword{spatial}
