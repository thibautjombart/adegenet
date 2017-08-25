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
#' @rdname AICc
#'
AICc <- function(object, ...) {
    UseMethod("AICc", object)
}





#' @export
#' @aliases AICc.snapclust
#' @rdname AICc
AICc.snapclust <- function(object, ...) {

    ## The number of parameters is defined as:
    ## (number of independent allele frequencies) x (nb clusters).
    k <- object$n.param
    n <- length(x$group)
    AIC(object) + (2 * k * (k + 1)) / (n - k -1)
}
