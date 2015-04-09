
####################
# Function makefreq
####################

#' Compute allelic frequencies
#'
#' The function \code{makefreq} is a generic to compute allele frequencies.
#' These can be derived for \linkS4class{genind} or \linkS4class{genpop} objects.
#' In the case of \linkS4class{genind} objects, data are kept at the individual level, but standardised so that allele frequencies sum up to 1.
#'
#' There are 3 treatments for missing values: \cr - NA: kept as NA.\cr - 0:
#' missing values are considered as zero. Recommended for a PCA on
#' compositionnal data.\cr - "mean": missing values are given the mean
#' frequency of the corresponding allele. Recommended for a centred PCA.
#'
#' Note that this function is now a simple wrapper for the accessor \code{\link{tab}}.
#'
#' @param x a \linkS4class{genind} or \linkS4class{genpop} object.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @param missing treatment for missing values. Can be NA, 0 or "mean" (see
#' details)
#' @param truenames deprecated; there for backward compatibility
#' @param ... further arguments (curently unused)
#'
#' @return Returns a list with the following components: \item{tab}{matrix of
#' allelic frequencies (rows: populations; columns: alleles).}
#' \item{nobs}{number of observations (i.e. alleles) for each population x
#' locus combinaison.} \item{call}{the matched call}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{genpop}}
#' @keywords manip multivariate
#'
#'#' @examples
#'
#' \dontrun{
#' data(microbov)
#' obj1 <- microbov
#' obj2 <- genind2genpop(obj1)
#'
#' # perform a correspondance analysis on counts data
#' Xcount <- tab(obj2, NA.method="zero")
#' ca1 <- dudi.coa(Xcount,scannf=FALSE)
#' s.label(ca1$li,sub="Correspondance Analysis",csub=1.2)
#' add.scatter.eig(ca1$eig,nf=2,xax=1,yax=2,posi="topleft")
#'
#' # perform a principal component analysis on frequency data
#' Xfreq <- makefreq(obj2, missing="mean")
#' Xfreq <- tab(obj2, NA.method="mean") # equivalent to line above
#' pca1 <- dudi.pca(Xfreq,scale=FALSE,scannf=FALSE)
#' s.label(pca1$li,sub="Principal Component Analysis",csub=1.2)
#' add.scatter.eig(pca1$eig,nf=2,xax=1,yax=2,posi="top")
#' }
#'
#' @export makefreq
#' @docType methods
#' @rdname makefreq
#'
setGeneric("makefreq", function(x, ...) standardGeneric("makefreq"))

#' @export
#' @rdname makefreq
#' @aliases makefreq,genind-methods
#' @aliases makefreq.genind
setMethod ("makefreq", signature(x="genind"), function(x, quiet=FALSE, missing=NA, truenames=TRUE, ...){
  if(!quiet) cat("\n Finding allelic frequencies from a genpop object... \n")

  ## treat 'missing' argument for compatibility with 'tab
  ## NA -> "asis"
  ## 0 -> "zero"
  ## "mean" -> "mean"
  if(is.na(missing)) {
      NA.method <- "asis"
  } else if(is.numeric(missing) && missing==0){
      NA.method <- "zero"
  } else {
      NA.method <- "mean"
  }

  out <- tab(x, freq=TRUE, NA.method=NA.method)

  if(!quiet) cat("\n...done.\n\n")

  return(out)
}) #end makefreq for genind




## genpop method
#' @export
#' @rdname makefreq
#' @aliases makefreq,genpop-methods
#' @aliases makefreq.genpop
setMethod ("makefreq", signature(x="genpop"), function(x, quiet=FALSE, missing=NA, truenames=TRUE, ...){

  if(!quiet) cat("\n Finding allelic frequencies from a genpop object... \n")

   ## treat 'missing' argument for compatibility with 'tab
  ## NA -> "asis"
  ## 0 -> "zero"
  ## "mean" -> "mean"
  if(is.na(missing)) {
      NA.method <- "asis"
  } else if(is.numeric(missing) && missing==0){
      NA.method <- "zero"
  } else {
      NA.method <- "mean"
  }

  out <- tab(x, freq=TRUE, NA.method=NA.method)

  if(!quiet) cat("\n...done.\n\n")

  return(out)
}) #end makefreq for genpop

