##############
# loadingplot
##############


#' Represents a cloud of points with colors
#' 
#' The \code{loadingplot} function represents positive values of a vector and
#' identifies the values above a given threshold. It can also indicate groups
#' of observations provided as a factor. \cr
#' 
#' Such graphics can be used, for instance, to assess the weight of each
#' variable (loadings) in a given analysis.
#' 
#' 
#' @aliases loadingplot loadingplot.default
#' @param x either a vector with numeric values to be plotted, or a matrix-like
#' object containing numeric values. In such case, the \code{x[,axis]} is used
#' as vector of values to be plotted.
#' @param at an optional numeric vector giving the abscissa at which loadings
#' are plotted. Useful when variates are SNPs with a known position in an
#' alignement.
#' @param threshold a threshold value above which values of x are identified.
#' By default, this is the third quartile of x.
#' @param axis an integer indicating the column of x to be plotted; used only
#' if x is a matrix-like object.
#' @param fac a factor defining groups of observations.
#' @param byfac a logical stating whether loadings should be averaged by groups
#' of observations, as defined by \code{fac}.
#' @param lab a character vector giving the labels used to annotate values
#' above the threshold; if NULL, names are taken from the object.
#' @param cex.lab a numeric value indicating the size of annotations.
#' @param cex.fac a numeric value indicating the size of annotations for groups
#' of observations.
#' @param lab.jitter a numeric value indicating the factor of randomisation for
#' the position of annotations. Set to 0 (by default) implies no randomisation.
#' @param main the main title of the figure.
#' @param xlab the title of the x axis.
#' @param ylab the title of the y axis.
#' @param srt rotation of the labels; see ?text.
#' @param adj adjustment of the labels; see ?text.
#' @param \dots further arguments to be passed to the plot function.
#' @return Invisibly returns a list with the following components:\cr -
#' threshold: the threshold used\cr - var.names: the names of observations
#' above the threshold\cr - var.idx: the indices of observations above the
#' threshold\cr - var.values: the values above the threshold\cr
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords multivariate hplot
#' @examples
#' 
#' x <- runif(20)
#' names(x) <- letters[1:20]
#' grp <- factor(paste("group", rep(1:4,each=5)))
#' 
#' ## basic plot
#' loadingplot(x)
#' 
#' ## adding groups
#' loadingplot(x,fac=grp,main="My title",cex.lab=1)
#' 
#' @export loadingplot
loadingplot <- function (x, ...) UseMethod("loadingplot")


loadingplot.default <- function(x, at=NULL, threshold=quantile(x,0.75), axis=1, fac=NULL, byfac=FALSE,
                        lab=NULL, cex.lab=0.7, cex.fac=1, lab.jitter=0,
                        main="Loading plot", xlab="Variables", ylab="Loadings", srt=0, adj=NULL, ...){
    ## some checks
    if(is.data.frame(x) | is.matrix(x)){
        if(is.null(lab)) {lab <- rownames(x)}
        x <- x[,axis]
    } else {
        if(is.null(lab)) {lab <- names(x)}
    }

    names(x) <- lab <- rep(lab, length=length(x))

    if(!is.numeric(x)) stop("x is not numeric")
    if(any(is.na(x))) stop("NA entries in x")
    if(any(x<0)) {
        warning("Some values in x are less than 0\n Using abs(x) instead, but this might not be optimal.")
        x <- abs(x)
    }
    if(is.null(at)){
        at <- 1:length(x)
    } else {
        if(length(at) != length(x)) stop("x and at do not have the same length.")
    }

    ## preliminary computations
    y.min <- min(min(x),0)
    y.max <- max(max(x),0)
    y.offset <- (y.max-y.min)*0.02
    if(is.null(lab)) {lab <- 1:length(x)}

    if(!is.null(fac)){
        if(byfac){
            x <- tapply(x, fac, mean)
            if(length(lab) != length(x)) lab <- names(x)
        } else {
            fac <- factor(fac, levels=unique(fac))
            grp.idx <- cumsum(table(fac)) + 0.5
            grp.lab.idx <- tapply(1:length(x), fac, mean)
            grp.lab <- names(grp.idx)
            grp.idx <- grp.idx[-length(grp.idx)]
    }
    } # end fac handling


    ## start the plot
    dat <- cbind(at, x)
    plot(dat, type="h", xlab=xlab, ylab=ylab,
         main=main, xaxt="n", ylim=c(y.min,y.max*1.2), ...)

    ## add groups of variables (optional)
    if(!is.null(fac) & !byfac) {
        abline(v=grp.idx,lty=2) # split groups of variables
        text(x=grp.lab.idx,y=y.max*1.15, labels=grp.lab, cex=cex.fac) # annotate groups
    }

    ## annotate variables that are above the threshold
    if(sum(x > threshold)>0){
        x.ann <- at[x > threshold]
        x.ann <- jitter(x.ann,factor=lab.jitter)

        y.ann <- x[x > threshold] + y.offset
        y.ann <- jitter(y.ann,factor=lab.jitter)

        txt.ann <- lab[x > threshold]
        text(x=x.ann, y=y.ann, label=txt.ann, cex=cex.lab, srt=srt, adj=adj)
    
    ## indicate the threshold
    abline(h=threshold, col="grey")

    ## build the result
    
    res <- list(threshold=threshold,
                var.names=txt.ann,
                var.idx=which(x > threshold),
                var.values=x[x > threshold])
    return(invisible(res))
    }

    return(NULL) # if no point above threshold
} # end loadingplot
