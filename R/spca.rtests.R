########################################
#
# Tests for global and local structures
#
# Thibaut Jombart 2007
# t.jombart@imperial.ac.uk
#
########################################



###############
# global.rtest
###############
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

