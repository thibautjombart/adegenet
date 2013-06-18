## pctest <- function(x, nperm=99, center=TRUE, scale=TRUE, method=c("sigma", "binom"), quiet=FALSE, plot=TRUE){
##     ## check x
##     if(!is.genind(x) & !is.genpop(x)){
##         stop("x is not a genind or a genpop object")
##     }

##     ## a few general variables
##     N <- nrow(x@tab)
##     P <- ncol(x@tab)

##     ## make tables of allele frequencies
##     X <- scaleGen(x, center=center, scale=scale, method=method, missing="mean")
##     fac.loc <- factor(sub("[.][^.]*$","",colnames(X)))
##     lX <- lapply(levels(fac.loc), function(id) X[,fac.loc==id,drop=FALSE])

##     ## auxil function to compute the first eigenvalue
##     if(N > P){ # N > P
##         f1 <- function(A){
##             Z <- t(A) %*% A / N
##             return(eigen(Z, symmetric=TRUE, only.values=TRUE)$values[1])
##         }
##     } else { #p <= n
##         f1 <- function(A){
##             Z <- A %*% t(A) / N
##             return(eigen(Z, symmetric=TRUE, only.values=TRUE)$values[1])
##         }
##     }


##     ## Monte Carlo procedure
##     makeOnePerm <- function(listX){
##         return(as.matrix(data.frame(lapply(listX, function(e) e[sample(N),,drop=FALSE]))))
##     }

##     if(quiet){
##         sim <- sapply(1:nperm, function(i) f1(makeOnePerm(lX)))
##     } else {
##         cat("\n Computing", nperm, "simulated eigenvalues ")
##         sim <- sapply(1:nperm, function(i) {cat(ifelse(i%%10==0, i, "."));return(f1(makeOnePerm(lX)))} )
##         cat(" done.\n")
##     }
##     ini <- f1(X)

##     ## return res
##     myCall <- match.call()
##     res <- as.randtest(sim=sim, obs=ini, alter="greater", call=myCall)
##     if(plot) {
##         plot(res, nclass=NULL, main="1st eigenvalue vs simulated eigenvalues (histogram)")
##     }
##     return(res)
## }
