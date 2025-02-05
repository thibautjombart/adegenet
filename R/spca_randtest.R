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
spca_randtest <-function(x, nperm = 499, p=.05){

  if(!inherits(x, "spca")){
    stop("x must be an spca object")
  }
  if (!requireNamespace("adespatial", quietly = TRUE)) {
    stop("the adespatial package is required for this funciton to work")
  }

  ## This function compute the eigenvalues for a given data object.

  get_stats <- function(obj){
    obj_pca <- ade4::dudi.pca(obj, center = FALSE, scale = FALSE,
                              scannf = FALSE)
    obj_spca <- adespatial::multispati(dudi = obj_pca,
                                 listw = x$lw, scannf = FALSE,
                                 nfposi = 1, nfnega = 1)
    lambda <- obj_spca$eig
 
    return(lambda)
  }


  ## This function permutes individuals (rows) in the dataset.

  perm_data <- function(obj = x$tab){
    obj[sample(1:nrow(obj)), , drop = FALSE]
  }

  ## Retains all simulated eigenvalues

  sims <- vapply(seq_len(nperm),
                 function(i) get_stats(perm_data()),
                 double(length(x$eig)))

  obs <- get_stats(x$tab)

  ## sum(eigenvalues >= 0) for global structure
  ## sum(eigenvalues < 0) for local structure

  obs_pos <- obs[obs >= 0]
  obs_neg <- obs[obs < 0]
  stats <- c(pos = sum(obs_pos), 
               neg = sum(abs(obs_neg)))

  pos_test <- as.randtest(sim = sum(sims[sims >= 0]), obs = stats[1], alter = "greater")
  neg_test <- as.randtest(sim = sum(abs(sims[sims < 0])), obs = stats[2], alter = "greater")

  ## Single eigenvalue observed p-value and Bonferroni correction
 
  eigen_p<-sapply(1:length(obs), function(e) as.randtest(abs(sims[e,]), abs(obs[e]), alter="greater")[[5]])

  bon_corr<-sapply(1:length(obs), function(e) if(obs[e] >= 0){p/e
			}else{p/(length(obs)-e+1)})

  eigentest<-cbind(round(obs, 4), round(eigen_p, 4), round(bon_corr, 4))
 
  colnames(eigentest)<-c("eigen", "sim_p", "bonf_p")

  list(global = pos_test, local = neg_test, eigentest = eigentest)

}
	
