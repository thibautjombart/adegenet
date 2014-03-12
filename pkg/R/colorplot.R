##
## COLOR PLOT
##
## used to plot up to 3 variables in space using RGB system
##
## all coded in S3 method (arguments vary largely)
##


##########
# generic
##########
colorplot <- function(...){
    UseMethod("colorplot")
}



#################
# default method
#################
colorplot.default <- function(xy, X, axes=NULL, add.plot=FALSE, defaultLevel=0, transp=FALSE, alpha=.5, ...){

    ## some checks
    if(any(is.na(xy))) stop("NAs exist in xy")
    xy <- as.matrix(xy)
    if(!is.numeric(xy)) stop("xy is not numeric")
    if(nrow(xy) != nrow(X)) stop("xy and X have different row numbers")
    if(is.null(axes)) {
        axes <- 1:min(ncol(X),3)
    }
    X <- as.matrix(X[,axes,drop=FALSE])
    if(any(is.na(X))) stop("NAs exist in X")
    if(!is.numeric(X)) stop("X is not numeric")
    if(defaultLevel < 0 | defaultLevel>1) stop("defaultLevel must be between 0 and 1")

    ## function mapping x to [0,+inf[
    f1 <- function(x){
        if(any(x<0)) {
            x <- x + abs(min(x))
        }
        return(x)
    }

    ## apply f1 to X
    X <- apply(X, 2, f1)

    v1 <- X[,1]
    if(ncol(X)>=2) {v2 <- X[,2]} else {v2 <- defaultLevel}
    if(ncol(X)>=3) {v3 <- X[,3]} else {v3 <- defaultLevel}

    ## make colors
      if(transp){
        col <- rgb(v1/max(X), v2/max(X), v3/max(X), alpha)
    } else {
        col <- rgb(v1, v2, v3, maxColorValue=max(X))
    }

    ## handle ...
    listArgs <- list(...)
    if(is.null(listArgs$pch)) {listArgs$pch <- 20}

    ## build list of arguments
    listArgs$x <- xy
    listArgs$col <- col

    ## plot data
    if(!add.plot) {
        do.call(plot,listArgs)
    } else {
        do.call(points,listArgs)
    }

    ##return(invisible(match.call()))
    return(invisible(col))
} # end colorplot.default
