#'
#' Expected heterozygosity (Hs)
#'
#' This function computes the expected heterozygosity (Hs) within
#' populations of a \linkS4class{genpop} object. This function is
#' available for codominant markers (\code{@@type="codom"}) only. Hs is
#' commonly used for measuring within population genetic diversity (and
#' as such, it still has sense when computed from haploid data).
#'
#' @aliases Hs
#' @aliases Hs.test
#'
#' @export
#'
#' @seealso \code{\link{spca}}
#' @details
#' Let \emph{m(k)} be the number of alleles of locus \emph{k}, with a
#' total of \emph{K} loci. We note \eqn{f_i} the allele frequency of
#' allele \emph{i} in a given population. Then, \eqn{Hs} is given for a
#' given population by:\cr
#'
#' \eqn{\frac{1}{K} \sum_{k=1}^K (1 - \sum_{i=1}^{m(k)} f_i^2)}
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a \linkS4class{genind} or \linkS4class{genpop} object.
#' @param pop only used if x is a \linkS4class{genind}; an optional factor to be used as population; if not provided, pop(x) is used.
#'
#' @return a vector of Hs values (one value per population)
#' @examples
#' \dontrun{
#' data(nancycats)
#' Hs(genind2genpop(nancycats))
#' }
#'
Hs <- function(x, pop=NULL) {

    ## CHECKS
    if(is.genind(x)){
        if(!is.null(pop)) pop(x) <- pop
        x <- genind2genpop(x, quiet=TRUE)
    }
    if(!is.genpop(x)) stop("x is not a valid genpop object")
    if(x@type=="PA") stop("not implemented for presence/absence markers")


    ## MAIN COMPUTATIONS
    ## get number of typed loci per pop
    if(any(is.na(tab(x)))){
        temp <- apply(tab(x),1,tapply, locFac(x), function(e) !any(is.na(e)))
        nLoc.typed <- apply(temp,2,sum)
    } else {
        nLoc.typed <- nLoc(x)
    }

    F <- tab(x, freq=TRUE, NA.method="asis")
    res <- 1-(apply(F^2,1,sum, na.rm=TRUE))/nLoc.typed

    return(res)
} # end Hs


