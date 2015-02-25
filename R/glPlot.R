###########
## glPlot
############


#' Plotting genlight objects
#' 
#' \linkS4class{genlight} object can be plotted using the function
#' \code{glPlot}, which is also used as the dedicated \code{plot} method. These
#' functions relie on \code{\link{image}} to represent SNPs data. More
#' specifically, colors are used to represent the number of second allele for
#' each locus and individual.
#' 
#' 
#' @aliases glPlot plot.genlight plot,genlight-method
#' @param x a \linkS4class{genlight} object.
#' @param col an optional color vector; the first value corresponds to 0
#' alleles, the last value corresponds to the ploidy level of the data.
#' Therefore, the vector should have a length of (\code{ploidy(x)+1}).
#' @param legend a logical indicating whether a legend should be added to the
#' plot.
#' @param posi a character string indicating where the legend should be
#' positioned. Can be any concatenation of "bottom"/"top" and "left"/"right".
#' @param bg a color used as a background for the legend; by default,
#' transparent white is used; this may not be supported on some devices, and
#' therefore background should be specified (e.g. \code{bg="white"}).
#' @param \dots further arguments to be passed to \code{\link{image}}.
#' @param y ununsed argument, present for compatibility with the \code{plot}
#' generic.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{\link{genlight}}: class of object for storing massive
#' binary SNP data.
#' 
#' - \code{\link{glSim}}: a simple simulator for \linkS4class{genlight}
#' objects.
#' 
#' - \code{\link{glPca}}: PCA for \linkS4class{genlight} objects.
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' ## simulate data
#' x <- glSim(100, 1e3, n.snp.struc=100, ploid=2)
#' 
#' ## default plot
#' glPlot(x)
#' plot(x) # identical plot
#' 
#' ## disable legend
#' plot(x, leg=FALSE)
#' 
#' ## use other colors
#' plot(x, col=heat.colors(3), bg="white")
#' }
#' 
#' @export glPlot
glPlot <- function(x, col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5),...) {

    ## get plotted elements ##
    X <- t(as.matrix(x))
    X <- X[,ncol(X):1]
    ylabpos <- pretty(1:nInd(x),5)
    if(is.null(col)) {
        myCol <- colorRampPalette(c("royalblue3", "firebrick1"))(max(X,na.rm=TRUE)+1)
    } else {
        myCol <- col
    }

    ## draw the plot ##
    ## main plot
    image(x=1:nLoc(x), y=1:nInd(x), z=X, xlab="SNP index", ylab="Individual index", yaxt="n", col=myCol, ...)

    ## add y axis
    axis(side=2, at=nInd(x)-ylabpos+1, labels=ylabpos)

    ## add legend
    if(legend){
        legend(posi, fill=myCol, legend=0:max(X,na.rm=TRUE), horiz=TRUE, bg=bg, title="Number of 2nd allele")
    }

    return(invisible())
} # end plot for glPlot



## hack to remove the NOTE in R CMD check about:
## "plot,genlight: no visible binding for global variable ‘y’"
if(getRversion() >= "2.15.1")  utils::globalVariables("y")

## plot method
setMethod("plot", signature(x="genlight", y="ANY"), function(x, y=NULL, col=NULL, legend=TRUE,
                                          posi="bottomleft", bg=rgb(1,1,1,.5),...) {
    glPlot(x, col=col, legend=legend, posi=posi, bg=bg, ...)
})
