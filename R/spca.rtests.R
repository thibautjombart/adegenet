########################################
#
# Tests for global and local structures
#
# Thibaut Jombart 2007
# t.jombart@@imperial.ac.uk
#
########################################



###############
# global.rtest
###############


#' Global and local tests
#' 
#' These two Monte Carlo tests are used to assess the existence of global and
#' local spatial structures. They can be used as an aid to interpret global
#' and local components of spatial Principal Component Analysis (sPCA).\cr
#' 
#' They rely on the decomposition of a data matrix X into global and local
#' components using multiple regression on Moran's Eigenvector Maps (MEMs).
#' They require a data matrix (X) and a list of weights derived from a
#' connection network. X is regressed onto global MEMs (U+) in the global test
#' and on local ones (U-) in the local test. One mean \eqn{R^2}{R^2} is
#' obtained for each MEM, the k highest being summed to form the test
#' statistic.
#' 
#' The reference distribution of these statistics are obtained by randomly
#' permuting the rows of X.
#' 
#' This test is purely R code. A C or C++ version will be developed soon.
#' 
#' @aliases global.rtest local.rtest
#' @param X a data matrix, with variables in columns
#' @param listw a list of weights of class \code{listw}. Can be obtained easily
#' using the function \code{chooseCN}.
#' @param k integer: the number of highest \eqn{R^2}{R^2} summed to form the
#' test statistics
#' @param nperm integer: the number of randomisations to be performed.
#' @return An object of class \code{randtest}.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{chooseCN}}, \code{\link{spca}}, \code{\link{monmonier}}
#' @references Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
#' Revealing cryptic spatial patterns in genetic variability by a new
#' multivariate method. \emph{Heredity}, \bold{101}, 92--103.
#' @keywords multivariate spatial
#' @examples
#' 
#' \dontrun{
#'  data(sim2pop)
#' if(require(spdep)){
#' cn <- chooseCN(sim2pop@@other$xy,ask=FALSE,type=1,plot=FALSE,res="listw")
#' 
#' # global test
#' Gtest <- global.rtest(sim2pop@@tab,cn)
#' Gtest
#' 
#' # local test
#' Ltest <- local.rtest(sim2pop@@tab,cn)
#' Ltest
#' }
#' }
#' 
#' @export global.rtest
global.rtest <- function(X, listw, k=1, nperm=499){
  if (!require(spdep)) stop("Package spdep is required.")
  if (!inherits(listw, "listw")) stop("object of class 'listw' expected")
  if (listw$style != "W") stop("object of class 'listw' with style 'W' expected")
  if(any(is.na(X))) stop("NA entries in X")

  n <- nrow(X)
  X <- scalewt(X)

  # computation of U+
  temp <- orthobasis.listw(listw)
  val <- attr(temp,"values")
  U <- as.matrix(temp)
  Upos <-  U[,val > -1/(n-1)]

  # test statistic
  calcstat <- function(X,k){
    R <- ( t(X) %*% Upos ) / n
    R2 <- R*R
    temp <- sort(apply(R2,2,mean),decreasing=TRUE)
    stat <- sum(temp[1:k])
    return(stat)
  }

  ini <- calcstat(X,k)

  sim <- sapply(1:nperm, function(i) calcstat( X[sample(1:n),], k ) )

  res <- as.randtest(sim=sim, obs=ini, alter="greater")
  res$call <- match.call()

  return(res)
} #end global.rtest





###############
# local.rtest
###############
local.rtest <- function(X, listw, k=1, nperm=499){
  if (!require(spdep)) stop("Package spdep is required.")
  if (!inherits(listw, "listw")) stop("object of class 'listw' expected")
  if (listw$style != "W") stop("object of class 'listw' with style 'W' expected")
  if(any(is.na(X))) stop("NA entries in X")

  n <- nrow(X)
  X <- scalewt(X)

  # computation of U-
  temp <- orthobasis.listw(listw)
  val <- attr(temp,"values")
  U <- as.matrix(temp)
  Uneg <-  U[,val < -1/(n-1)]

  X <- scalewt(X)

  # test statistic
  calcstat <- function(X,k){
    R <- ( t(X) %*% Uneg ) / n
    R2 <- R*R
    temp <- sort(apply(R2,2,mean),decreasing=TRUE)
    stat <- sum(temp[1:k])
    return(stat)
  }

  ini <- calcstat(X,k)

  sim <- sapply(1:nperm, function(i) calcstat( X[sample(1:n),], k ) )

  res <- as.randtest(sim=sim, obs=ini, alter="greater")
  res$call <- match.call()

  return(res)
} #end local.rtest

