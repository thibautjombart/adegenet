
###########
## glPlot
############
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
