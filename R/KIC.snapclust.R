#' Compute Akaike Information Criterion for small samples (AICc) for snapclust
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @param object An object returned by the function \code{\link{snapclust}}.
#'
#' @param ... Further arguments for compatibility with the \code{AIC} generic
#'     (currently not used).
#'
#' @seealso  \code{\link{snapclust}} to generate clustering solutions.
#'
#' @rdname KIC
#'
KIC <- function(object, ...) {
    UseMethod("KIC", object)
}





#' @export
#' @aliases KIC.snapclust
#' @rdname KIC
KIC.snapclust <- function(object, ...) {

    ## The number of parameters is defined as:
    ## (number of independent allele frequencies) x (nb clusters).
    k <- object$n.param
    -2 * object$ll + 3 * (k + 1)

}
