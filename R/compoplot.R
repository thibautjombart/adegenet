#' Genotype composition plot
#'
#' The compoplot uses a barplot to represent the group assignment probability of
#' individuals to several groups. It is a generic with methods for the following
#' objects:
#'
#' \itemize{
#' 
#' \item \code{matrix}: a matrix with individuals in row and genetic clusters in
#' column, each entry being an assignment probability of the corresponding
#' individual to the corresponding group
#' 
#' \item \code{dapc}: the output of the \code{dapc} function; in this case,
#' group assignments are based upon geometric criteria in the discriminant space
#'
#' \item \code{snapclust}: the output of the \code{snapclust} function; in
#' this case, group assignments are based upon the likelihood of genotypes
#' belonging to their groups
#' 
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
#' 
#' @param ... further arguments to be passed to \code{barplot}
#' 
compoplot <- function(x, ...){
    UseMethod("compoplot", x)
}


#' Palette parser
#'
#' @param inPAL a palette function OR a character vector (named or unnamed)
#' @param npop number of populations/colors desired
#' @param pnames names of these populations if `inPal` is a function
#'
#' @md
#' @noRd
#' @return a named character vector specifying the colors for the palette.
#' @keywords internal
#' 
#' @note This was originally from the poppr package [commit a0818eed6](https://github.com/grunwaldlab/poppr/commit/a0818eed6a72d9145e46da73715dc22be0640b0c)
#'
#' @examples
#' .palette_parser(rainbow, 5, letters[1:5])
#' .palette_parser(colors()[1:5], 5, letters[1:5])
.palette_parser <- function(inPAL, npop, pnames){
  PAL <- try(match.fun(inPAL, descend = FALSE), silent = TRUE)
  if ("try-error" %in% class(PAL)){
    if (all(pnames %in% names(inPAL))){
      color <- inPAL[pnames]
    } else if (npop == length(inPAL)){
      color <- stats::setNames(inPAL, pnames)
    } else if (npop < length(inPAL)){
      warning("Number of populations fewer than number of colors supplied. Discarding extra colors.")
      color <- stats::setNames(inPAL[1:npop], pnames)
    } else {
      warning("insufficient color palette supplied. Using funky().")
      color <- stats::setNames(funky(npop), pnames)
    }
  } else {
    color   <- stats::setNames(PAL(npop), pnames)
  }
  return(color)
}


#' @rdname compoplot
#'
#' @aliases compoplot.matrix
#' @export
#'
#' @param col.pal a color palette to be used for the groups; defaults to \code{funky}
#'
#' @param border a color for the border of the barplot; use \code{NA} to
#' indicate no border. 
#' 
#' @param show.lab a logical indicating if individual labels should be displayed
#' 
#' @param lab a vector of individual labels; if NULL, row.names of the matrix are used
#' 
#' @param legend a logical indicating whether a legend should be provided for the colors
#' 
#' @param txt.leg a character vector to be used for the legend
#' 
#' @param n.col the number of columns to be used for the legend
#' 
#' @param posi the position of the legend
#' 
#' @param cleg a size factor for the legend
#' 
#' @param bg the background to be used for the legend
#' 
#' @param subset a subset of individuals to retain
#'
compoplot.matrix <- function(x, col.pal = funky, border = NA,
                             subset = NULL, show.lab = FALSE,
                             lab = rownames(x), legend = TRUE,
                             txt.leg = colnames(x), n.col = 4,
                             posi = NULL, cleg = .8,
                             bg = transp("white"), ...) {

    ## individual labels
    if (!show.lab) {
        lab <- rep("", nrow(x))
    }

    ## handle subset
    if (!is.null(subset)) {
        names(lab) <- rownames(x)
        x <- x[subset, , drop=FALSE]
        lab <- lab[rownames(x)]
    }

    ## group labels
    if (is.null(txt.leg)) {
        txt.leg <- colnames(x)
    }
    if (is.null(txt.leg)) {
        txt.leg <- seq_len(ncol(x))
    }
    ## generate colors, process arguments
    col <- .palette_parser(col.pal, ncol(x), txt.leg)
    
    ## position of the legend
    if (is.null(posi)) {
        posi <- list(x=0, y=-.01)
    }

    ## make the plot: we need to suppress warnings because '...' could contain
    ## arguments from other methods not meant to be used by 'barplot'

    suppressWarnings(
        out <- barplot(t(x), col = col,
                       ylab = "membership probability",
                       names.arg = lab, las = 3,
                       border = border, ...) )

    if (legend) {
        oxpd <- par("xpd")
        par(xpd=TRUE)
        legend(posi, fill=col, legend = txt.leg,
               cex = cleg, ncol = n.col, bg = bg)
        on.exit(par(xpd=oxpd))
    }

    return(invisible(out))
}





#' @rdname compoplot
#' @aliases compoplot.dapc
#' @export
#' @param only.grp a subset of groups to retain

## The compoplot for DAPC is basically a compoplot.matrix on the predicted group membership
## probabilities. Only extra features related to keeping a subset of groups or individuals.

compoplot.dapc <- function(x, only.grp=NULL, border = NA, ...){
    ## get predictions and subset if needed
    pred <- predict(x)$posterior

    ## handle group subsetting
    if (!is.null(only.grp)) {
        if(is.numeric(only.grp) || is.logical(only.grp)) {
            only.grp <- levels(x$grp)[only.grp]
        }
        to.keep <- as.character(x$grp) %in% only.grp
        pred <- pred[to.keep, , drop=FALSE]
        lab <- lab[to.keep]
    }

    ## call matrix method
    compoplot(pred, border = border, ...)

} # end compoplot









#' @rdname compoplot
#' @export
compoplot.snapclust <- function(x, border = NA, ...) {
    compoplot(x$proba, border = border, ...)
}
