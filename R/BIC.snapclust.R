#' Compute Bayesian Information Criterion (BIC) for snapclust
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @param object An object returned by the function \code{\link{snapclust}}.
#'
#' @param ... Further arguments for compatibility with the \code{BIC} generic
#'     (currently not used).
#'
#' @seealso  \code{\link{snapclust}} to generate clustering solutions.
#'
#'
BIC.snapclust <- function(object, ...) {

  ## The number of parameters is defined as:
  ## (number of independent allele frequencies) x (nb clusters).

  n <- length(object$group)

  -2 * object$ll + log(n) * object$n.param
}
