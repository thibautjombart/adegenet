#' Export analysis for webapp visualisation
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#'
#' @param x The analysis to be exported.
#'
#' @param pop A factor indicating the population of each individual, in the same
#'   order as in the analysis.
#'
#' @param info A data.frame with additional information containing at least the
#'   following columns: \code{key} (individual label), \code{latitude}, and
#'   \code{longitude}.
#'
#' @param ... Further arguments to pass to other methods.
#'
#' @export
#'
#' @rdname export_to_webapp

export_to_webapp <- function(x, ...) {
  UseMethod("export_to_webapp")
}





#' @export
#' @rdname export_to_webapp

export_to_webapp.default <- function(x, ...) {
  msg <- sprintf("No method available for the class %s",
                 paste(class(x), collapse = ", "))
  stop(msg)
}






## All method will consist in merging output from the analysis with extra info
## containing latitude and longitude, stored in 'info'.

#' @export
#' @rdname export_to_webapp

export_to_webapp.dapc <- function(x, info, ...) {

  ## Extract principal components, groups, assigned groups and the corresponding
  ## probability.

  pcs <- x$ind.coord
  colnames(pcs) <- paste0("PC", 1:ncol(pcs))
  key <- rownames(pcs)
  grp <- x$grp
  assigned_grp <- x$assign
  support <- apply(x$posterior, 1, max)


  analysis <- cbind.data.frame(key, pcs,
                               grp,
                               assigned_grp,
                               support)

  ## process 'info' (checks that required columns are there)
  info <- check_info(info, key)

  out <- merge(analysis, info, by = "key")
  return(out)
}






## This internal function will merely check the content of the extra 'info'
## being provided, making sure key, latitude and longitude are provided.

check_info <- function(info, ref_keys) {
  info <- as.data.frame(info)

  if (!"key" %in% names(info)) {
    stop("'info' is missing a 'key' column")
  }

  if (!"lat" %in% names(info)) {
    stop("'info' is missing a 'lat' column")
  }

  if (!"lon" %in% names(info)) {
    stop("'info' is missing a 'lon' column")
  }

  nb_missing <- sum(!ref_keys %in% info$key)
  if (nb_missing > 0L) {
    msg <- sprintf("%d individuals are not documented in 'info'",
                   nb_missing)
    warning(msg)
  }

  return(info)
}
