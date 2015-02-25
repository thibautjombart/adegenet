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


#' Represents a cloud of points with colors
#' 
#' The \code{colorplot} function represents a cloud of points with colors
#' corresponding to a combination of 1,2 or 3 quantitative variables, assigned
#' to RGB (Red, Green, Blue) channels. For instance, this can be useful to
#' represent up to 3 principal components in space. Note that the property of
#' such representation to convey multidimensional information has not been
#' investigated.\cr
#' 
#' \code{colorplot} is a S3 generic function. Methods are defined for
#' particular objects, like \code{\link{spca}} objects.
#' 
#' 
#' @aliases colorplot colorplot.default
#' @param xy a numeric matrix with two columns (e.g. a matrix of spatial
#' coordinates.
#' @param X a matrix-like containing numeric values that are translated into
#' the RGB system. Variables are considered to be in columns.
#' @param axes the index of the columns of X to be represented. Up to three
#' axes can be chosen. If null, up to the first three columns of X are used.
#' @param add.plot a logical stating whether the colorplot should be added to
#' the existing plot (defaults to FALSE).
#' @param defaultLevel a numeric value between 0 and 1, giving the default
#' level in a color for which values are not specified. Used whenever less than
#' three axes are specified.
#' @param transp a logical stating whether the produced colors should be
#' transparent (TRUE) or not (FALSE, default).
#' @param alpha the alpha level for transparency, between 0 (fully transparent)
#' and 1 (not transparent); see \code{?rgb} for more details.
#' @param \dots further arguments to be passed to other methods. In
#' \code{colorplot.default}, these arguments are passed to plot/points
#' functions. See \code{?plot.default} and \code{?points}.
#' @return Invisibly returns a vector of colours used in the plot.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords multivariate hplot
#' @examples
#' 
#' # a toy example
#' xy <- expand.grid(1:10,1:10)
#' df <- data.frame(x=1:100, y=100:1, z=runif(100,0,100))
#' colorplot(xy,df,cex=10,main="colorplot: toy example")
#' 
#' \dontrun{
#' # a genetic example using a sPCA
#' if(require(spdep)){
#' data(spcaIllus)
#' dat3 <- spcaIllus$dat3
#' spca3 <- spca(dat3,xy=dat3$other$xy,ask=FALSE,type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)
#' colorplot(spca3, cex=4, main="colorplot: a sPCA example")
#' text(spca3$xy[,1], spca3$xy[,2], dat3$pop)
#' mtext("P1-P2 in cline\tP3 random \tP4 local repulsion")
#' }
#' }
#' 
#' @export colorplot
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
