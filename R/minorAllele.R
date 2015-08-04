#' Compute minor allele frequency
#'
#' This function computes the minor allele frequency for each locus in a \linkS4class{genind} object. To test if loci are polymorphic, see the function  \code{\link{isPoly}}.
#'
#' @param x  a \linkS4class{genind} object
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## LOAD DATA
#' data(nancycats)
#'
#' ## COMPUTE ALLELE FREQUENCIES
#' x <- nancycats
#' apply(tab(x, freq=TRUE),2,mean, na.rm=TRUE)
#'
#' ## GET MINOR ALLELE FREQUENCY
#' m.freq <- minorAllele(x)
#' m.freq
#' }
#'
#' @seealso \code{\link{isPoly}}
#'
minorAllele <- function(x){
    ## CHECK INPUT
    if(!is.genind(x)) stop("x is not a valid genind object")

    ## AUXIL FUNCTION
    f1 <- function(vec){
        vec <- sort(vec, decreasing=TRUE)
        res <- vec[min(2,length(vec))]
        return(res)
    }

    ## GET ALLELE FREQUENCIES
    freq <- apply(tab(x, freq=TRUE),2,mean, na.rm=TRUE)

    ## GET MINOR ALLELE
    out <- tapply(freq, locFac(x), f1)

    ## RETURN OUTPUT
    return(out)
}
