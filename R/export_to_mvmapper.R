#' Export analysis for mvmapper visualisation
#'
#' \code{mvmapper} is an interactive tool for visualising outputs of a
#' multivariate analysis on a map from a web browser. The function
#' \code{export_to_mvmapper} is a generic with methods for several standard
#' classes of analyses in \code{adegenet} and \code{ade4}. Information on
#' individual locations, as well as any other relevant data, is passed through
#' the second argument \code{info}. By default, the function returns a formatted
#' \code{data.frame} and writes the output to a .csv file.\cr
#'
#' \code{mvmapper} can be found at:
#' \url{https://github.com/genomeannotation/mvMapper}
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#'
#' @param x The analysis to be exported. Can be a \code{dapc}, \code{spca}, or a
#'   \code{dudi} object.
#'
#' @param info A \code{data.frame} with additional information containing at
#' least the following columns: \code{key} (unique individual identifier),
#' \code{lat} (latitude), and \code{lon} (longitude). Other columns will be
#' exported as well, but are optional.
#'
#'
#' @param write_file A \code{logical} indicating if the output should be written
#'   out to a .csv file. Defaults to \code{TRUE}.
#'
#' @param out_file A character string indicating the file to which the output
#'   should be written. If NULL, the file used will be named
#'   \code{'mvmapper_data_[date and time].csv'}
#'
#' @param ... Further arguments to pass to other methods.
#'
#' @return A \code{data.frame} which can serve as input to \code{mvmapper},
#' containing at least the following columns:
#'
#' \itemize{
#'
#' \item \code{key}: unique individual identifiers
#'
#' \item \code{PC1}: first principal component; further principal components are
#' optional, but if provided will be numbered and follow \code{PC1}.
#'
#' \item \code{lat}: latitude for each individual
#'
#' \item \code{lon}: longitude for each individual
#'
#' }
#'
#' In addition, specific information is added for some analyses:
#'
#' \itemize{
#'
#' \item \code{spca}: \code{Lag_PC} columns contain the lag-vectors of the
#' principal components; the lag operator computes, for each individual, the
#' average score of neighbouring individuals; it is useful for clarifying
#' patches and clines.
#'
#' \item \code{dapc}: \code{grp} is the group used in the analysis;
#' \code{assigned_grp} is the group assignment based on the discriminant
#' functions; \code{support} is the statistical support (i.e. assignment
#' probability) for \code{assigned_grp}.
#'
#' }
#'
#'
#'
#' @export
#'
#' @rdname export_to_mvmapper
#'
#' @seealso
#'
#' \code{mvmapper} is available at:
#' \url{https://popphylotools.github.io/mvMapper/}
#'


export_to_mvmapper <- function(x, ...) {
  UseMethod("export_to_mvmapper")
}





#' @export
#' @rdname export_to_mvmapper

export_to_mvmapper.default <- function(x, ...) {
  msg <- sprintf("No method available for the class %s",
                 paste(class(x), collapse = ", "))
  stop(msg)
}






## All method will consist in merging output from the analysis with extra info
## containing latitude and longitude, stored in 'info'.

#' @export
#' @rdname export_to_mvmapper
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
#' out <- export_to_mvmapper(dapc1, info)
#' head(out)

export_to_mvmapper.dapc <- function(x, info, write_file = TRUE, out_file = NULL, ...) {

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
  .write_mvmapper_output(out, write_file, out_file)
  return(out)
}






#' @export
#' @rdname export_to_mvmapper

export_to_mvmapper.dudi <- function(x, info, write_file = TRUE, out_file = NULL, ...) {

  ## Extract principal components, groups, assigned groups and the corresponding
  ## probability.

  pcs <- x$li
  colnames(pcs) <- paste0("PC", 1:ncol(pcs))
  key <- rownames(pcs)

  analysis <- cbind.data.frame(key, pcs)

  ## process 'info' (checks that required columns are there)
  info <- .check_info(info, key)

  out <- merge(analysis, info, by = "key")
  .write_mvmapper_output(out, write_file, out_file)
  return(out)
}






#' @export
#' @rdname export_to_mvmapper
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
#' out <- export_to_mvmapper(spca1, info)
#' head(out)
#'

export_to_mvmapper.spca <- function(x, info, write_file = TRUE, out_file = NULL, ...) {

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
  .write_mvmapper_output(out, write_file, out_file)
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






## This internal function writes results to a csv file if needed, and does
## nothing otherwise.
##
## 'x' is the data.frame output from the export function
## other arguments as documented
##
.write_mvmapper_output <- function(x, write_file = TRUE, out_file = NULL) {
  if (write_file) {
    if (is.null(out_file)) {
      out_file <- paste0("mvmapper_data_",
                         gsub(" ", "_", Sys.time()),
                         ".csv")
    }
    message("Writing output to the file: ", out_file)
    write.csv(x, out_file, row.names = FALSE)
  }
}
