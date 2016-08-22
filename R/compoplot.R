#' Genotype composition plot
#'
#' The compoplot uses a barplot to represent the group assignment probability of individuals to several groups. It is a generic with methods for the following objects:
#'
#' \itemize{
#' \item \code{matrix}: a matrix with individuals in row and genetic clusters in column, each entry being an assignment probability of the corresponding individual to the corresponding group
#' \item \code{dapc}: the output of the \code{dapc} function; in this case, group assignments are based upon geometric criteria in the discriminant space
#' \item \code{genclust.em}: the output of the \code{genclust.em} function; in this case, group assignments are based upon the likelihood of genotypes belonging to their groups
#' }
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @rdname compoplot
#' @aliases compoplot
#'
#' @param x an object to be used for plotting (see description)
#' @param ... further arguments to be passed to \code{barplot}
compoplot <- function(x, ...){
    UseMethod("compoplot", x)
}



#' @rdname compoplot
#'
#' @aliases compoplot.matrix
#' @export
#'
#' @param col.pal a color palette to be used for the groups; defaults to \code{funky}
#' @param show.lab a logical indicating if individual labels should be displayed
#' @param lab a vector of individual labels; if NULL, row.names of the matrix are used
#' @param legend a logical indicating whether a legend should be provided for the colors
#' @param txt.leg a character vector to be used for the legend
#' @param ncol the number of columns to be used for the legend
#' @param posi the position of the legend
#' @param cleg a size factor for the legend
#' @param bg the background to be used for the legend
#'
compoplot.matrix <- function(x, col.pal = funky, show.lab = FALSE,
                             lab = rownames(x), legend = TRUE,
                             txt.leg = colnames(x), ncol = 4,
                             posi = NULL, cleg = .8, bg = transp("white"),
                             ...) {

    ## generate colors, process arguments
    col <- col.pal(ncol(x))

    ## individual labels
    if (!show.lab) {
        lab <- rep("", nrow(x))
    }

    ## position of the legend
    if (is.null(posi)) {
        posi <- list(x=0, y=-.01)
    }

    ## make the plot ##
    barplot(t(x), border = NA, col = col, ylab = "membership probability",
            names.arg = lab, las = 3, ...)

    if (legend) {
        oxpd <- par("xpd")
        par(xpd=TRUE)
        legend(posi, fill=col, legend = txt.leg, cex=cleg, ncol=ncol, bg=bg)
        on.exit(par(xpd=oxpd))
    }

    return(invisible())
}





#' @rdname compoplot
#' @aliases compoplot.dapc
#' @export
compoplot.dapc <- function(x, only.grp=NULL, subset=NULL, new.pred=NULL, col=NULL, lab=NULL,
                      legend=TRUE, txt.leg=NULL, ncol=4, posi=NULL, cleg=.8, bg=transp("white"), ...){

    ## HANDLE ARGUMENTS ##
    ngrp <- length(levels(x$grp))

    ## HANDLE DATA FROM PREDICT.DAPC ##
    if(!is.null(new.pred)){
        n.new <- length(new.pred$assign)
        x$grp <- c(as.character(x$grp), rep("unknown", n.new))
        x$assign <- c(as.character(x$assign), as.character(new.pred$assign))
        x$posterior <- rbind(x$posterior, new.pred$posterior)
        lab <- c(lab, rownames(new.pred$posterior))
    }


    ## TREAT OTHER ARGUMENTS ##
    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
        lab <- lab[only.grp==ori.grp]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
        lab <- lab[subset]
    }

    ## call matrix method
    compoplot(t(x$posterior))

    return(invisible(match.call()))
} # end compoplot









#' @rdname genclust.em
#' @export
compoplot.genclust.em <- function(x, col.pal = funky, show.lab=TRUE, ...) {
    if (!show.lab) {
        lab <- rep("", nrow(x$proba))
    } else {
        lab <- NULL
    }
    barplot(t(x$proba), col=col.pal(ncol(x$proba)), border=NA, las=3,
            names.arg = lab, ylab = "Group assignment probability", ...)
}
