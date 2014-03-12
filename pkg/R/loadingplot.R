##############
# loadingplot
##############
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
