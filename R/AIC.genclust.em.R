#' Compute Akaike Information Criterion (AIC) for genclust.em
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @param object An object returned by the function \code{\link{genclust.em}}.
#'
#' @param ... Further arguments for compatibility with the \code{AIC} generic
#'     (currently not used).
#'
#' @seealso  \code{\link{genclust.em}} to generate clustering solutions.
#'
#'
AIC.genclust.em <- function(object, ...) {

    ## The number of parameters is defined as:
    ## (number of independent allele frequencies) x (nb clusters).

  -2 * object$ll + 2 * object$n.param
}
