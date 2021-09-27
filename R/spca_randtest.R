#' Monte Carlo test for sPCA
#'
#' The function \code{spca_randtest} implements Monte-Carlo tests for the
#' presence of significant spatial structures in a sPCA object. Two tests are
#' run, for global (positive autocorrelation) and local (negative
#' autocorrelation) structures, respectively. The test statistics used are the
#' sum of the absolute values of the corresponding eigenvalues.
#'
#' @export
#'
#' @author Original code by Valeria Montano adapted by Thibaut Jombart.
#'
#' @param x A \code{\link{spca}} object.
#'
#' @param nperm The number of permutations to be used for the test.
#'
#' @return
#'
#' A list with two objects of the class 'randtest' (see
#' \code{\link[ade4]{as.randtest}}), the first one for 'global' structures
#' (positivie autocorrelation) and the second for 'local' structures (negative
#' autocorrelation).
#'
#' @examples
#'
#' \dontrun{
#' ## Load data
#' data(sim2pop)
#'
#' ## Make spca
#' spca1 <- spca(sim2pop, type = 1, scannf = FALSE, plot.nb = FALSE)
#'
#' spca1
#' plot(spca1)
#'
#' ## run tests (use more permutations in practice, e.g. 999)
#' tests <- spca_randtest(spca1, nperm = 49)
#'
#' ## check results
#' tests
#' plot(tests[[1]]) # global structures
#  plot(tests[[2]]) # local structures
#'
#' }
#'
spca_randtest <-function(x, nperm = 499){

  if (!requireNamespace("adespatial", quietly=TRUE)) {
    install <- paste0('install.packages(', shQuote("adespatial"), ")")
    msg <- c("The adespatial package is required. Please use `", install, "` to install it")
    stop(paste(msg, collapse = ""))
  }

  if(!inherits(x, "spca")){
    stop("x must be an spca object")
  }


  ## This function compute the test statistics for a given data object. Two test
  ## statistics are computed, from the eigenvalues of the sPCA, called 'lambda':

  ## sum(lambda >= 0)
  ## sum(lambda < 0)

  get_stats <- function(obj){
    obj_pca <- ade4::dudi.pca(obj, center = FALSE, scale = FALSE,
                              scannf = FALSE)
    obj_spca <- adespatial::multispati(dudi = obj_pca,
                                 listw = x$lw, scannf = FALSE,
                                 nfposi = 1, nfnega = 1)
    lambda <- obj_spca$eig
    lambda_pos <- lambda[lambda >= 0]
    lambda_neg <- lambda[lambda < 0]
    stats <- c(pos = sum(lambda_pos),
               neg = sum(abs(lambda_neg)))

    return(stats)
  }


  ## This function permutes individuals (rows) in the dataset.

  perm_data <- function(obj = x$tab){
    obj[sample(1:nrow(obj)), , drop = FALSE]
  }

  sims <- vapply(seq_len(nperm),
                 function(i) get_stats(perm_data()),
                 double(2))

  obs <- get_stats(x$tab)

  pos_test <- as.randtest(sim = sims[1,], obs = obs[1], alter = "greater")
  neg_test <- as.randtest(sim = sims[2,], obs = obs[2], alter = "greater")

  list(global = pos_test, local = neg_test)
}

