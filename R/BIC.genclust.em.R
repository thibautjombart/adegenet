#' Compute Bayesian Information Criterion (BIC) for genclust.em
#'
#' Do not use. We work on that stuff. Contact us if interested.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @param object An object returned by the function \code{\link{genclust.em}}.
#'
#' @param ... Further arguments for compatibility with the \code{BIC} generic
#'     (currently not used).
#'
#' @seealso  \code{\link{genclust.em}} to generate clustering solutions.
#'
#'
BIC.genclust.em <- function(object, ...) {

  ## The number of parameters is defined as:
  ## (number of independent allele frequencies) x (nb clusters).

  n <- length(object$group)

  -2 * object$ll + log(n) * object$n.param
}
