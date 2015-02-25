############################
# Hs (expected heterozygosity)
############################


#' Expected heterozygosity
#' 
#' This function computes the expected heterozygosity (Hs) within populations
#' of a \linkS4class{genpop} object. This function is available for codominant
#' markers (\code{@@type="codom"}) only. Hs is commonly used for measuring
#' within population genetic diversity (and as such, it still has sense when
#' computed from haploid data).
#' 
#' Let \emph{m(k)} be the number of alleles of locus \emph{k}, with a total of
#' \emph{K} loci. We note \eqn{f_i} the allele frequency of allele \emph{i} in
#' a given population. Then, \eqn{Hs} is given for a given population by:\cr
#' 
#' \eqn{\frac{1}{K} \sum_{k=1}^K (1 - \sum_{i=1}^{m(k)} f_i^2)} \cr
#' 
#' @param x an object of class \linkS4class{genpop}.
#' @param truenames a logical indicating whether true labels (as opposed to
#' generic labels) should be used to name the output.
#' @return A vector of Hs values (one value per population).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' data(nancycats)
#' Hs(genind2genpop(nancycats))
#' }
#' 
#' @export Hs
Hs <- function(x, truenames=TRUE) {

    ## CHECKS
    if(is.genind(x)){
        x <- genind2genpop(x, quiet=TRUE)
    }
    if(!is.genpop(x)) stop("x is not a valid genpop object")
    if(x@type=="PA") stop("not implemented for presence/absence markers")


    ## MAIN COMPUTATIONS
    x.byloc <- seploc(x, truenames=truenames)

    lX <- lapply(x.byloc, function(e) makefreq(e, quiet=TRUE, truenames=truenames)$tab)
    lres <- lapply(lX, function(X) 1- apply(X^2,1,sum))
    res <- apply(as.matrix(data.frame(lres)),1,mean)

    return(res)
} # end Hs
