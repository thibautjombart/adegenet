#' Export analysis for webapp visualisation
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#'
#' @param x The analysis to be exported.
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
#'

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
#' @examples
#'
#' data(sim2pop)
#'
#' dapc1 <- dapc(sim2pop, n.pca = 10, n.da = 1)
#'
#' info <- data.frame(key = indNames(sim2pop),
#'                    lat = other(sim2pop)$xy[,2],
#'                    lon = other(sim2pop)$xy[,1],
#'                    Population = pop(sim2pop))
#'
#' out <- export_to_webapp(dapc1, info)
#' head(out)

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
  info <- .check_info(info, key)

  out <- merge(analysis, info, by = "key")
  return(out)
}






#' @export
#' @rdname export_to_webapp

export_to_webapp.dudi <- function(x, info, ...) {

  ## Extract principal components, groups, assigned groups and the corresponding
  ## probability.

  pcs <- x$li
  colnames(pcs) <- paste0("PC", 1:ncol(pcs))
  key <- rownames(pcs)

  analysis <- cbind.data.frame(key, pcs)

  ## process 'info' (checks that required columns are there)
  info <- .check_info(info, key)

  out <- merge(analysis, info, by = "key")
  return(out)
}






#' @export
#' @rdname export_to_webapp
#' @examples
#'
#' data(rupica)
#'
#' spca1 <- spca(rupica, type=5, d1 = 0, d2 = 2300,
#'               plot = FALSE, scannf = FALSE,
#'               nfposi = 2,nfnega = 0)
#'
#' info <- data.frame(key = indNames(rupica),
#'                    lat = rupica$other$xy[,2],
#'                    lon = rupica$other$xy[,1])
#'
#' out <- export_to_webapp(spca1, info)
#' head(out)
#'

export_to_webapp.spca <- function(x, info, ...) {

  ## Extract principal components, groups, assigned groups and the corresponding
  ## probability.

  pcs <- x$li
  colnames(pcs) <- paste0("PC", 1:ncol(pcs))
  lag_pcs <- x$ls
  colnames(lag_pcs) <- paste0("Lag_PC", 1:ncol(pcs))
  key <- rownames(pcs)

  analysis <- cbind.data.frame(key, pcs, lag_pcs)

  ## process 'info' (checks that required columns are there)
  info <- .check_info(info, key)

  out <- merge(analysis, info, by = "key")
  return(out)
}






## This internal function will merely check the content of the extra 'info'
## being provided, making sure key, latitude and longitude are provided.

.check_info <- function(info, ref_keys,
                       look_for = c("key", "lat", "lon")) {

  info <- as.data.frame(info)

  if (length(look_for) > 0L) {
    for (e in look_for) {
      if (!e %in% names(info)) {
        msg <- sprintf("'info' is missing a '%s' column", e)
        stop(msg)
      }
    }
  }

  nb_missing <- sum(!ref_keys %in% info$key)
  if (nb_missing > 0L) {
    msg <- sprintf("%d individuals are not documented in 'info'",
                   nb_missing)
    warning(msg)
  }

  return(info)
}
