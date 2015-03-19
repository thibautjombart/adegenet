
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
#'
#' obj2 <- genind2genpop(obj1)
#'
#' Xfreq <- makefreq(obj2,missing="mean")
#'
#'
#' # perform a correspondance analysis on counts data
#'
#' Xcount <- genind2genpop(obj1,missing="chi2")
#' ca1 <- dudi.coa(as.data.frame(Xcount@@tab),scannf=FALSE)
#' s.label(ca1$li,sub="Correspondance Analysis",csub=1.2)
#' add.scatter.eig(ca1$eig,nf=2,xax=1,yax=2,posi="topleft")
#'
#' # perform a principal component analysis on frequency data
#' pca1 <- dudi.pca(Xfreq$tab,scale=FALSE,scannf=FALSE)
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

  res <- list()

  out <- x@tab / x@ploidy

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

  f1 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(rep(NA,length(v)))
    return(v/(sum(v,na.rm=TRUE)))
  }

  res <- list()

  tabcount <- x@tab

  eff.pop <- t(apply(tabcount,1,function(r) tapply(r,x@loc.fac,sum,na.rm=TRUE)))
  if(nLoc(x)==1){ # fix for nloc==1
      eff.pop <- t(eff.pop)
  }

  # out is a pop x loci table of allelic frequencies
  out <- t(apply(tabcount,1,function(r) unlist(tapply(r,x@loc.fac,f1))))
  if(length(x@loc.nall)==1 && x@loc.nall[1]==1) out <- t(out) # matrix is transposed by apply if there's a single allele
  colnames(out) <- colnames(x@tab)

  # NA treatment
  # NA can be kept as is, or replaced 0 or by the mean frequency of the allele.
  if(!is.na(missing)){
    if(missing==0) out[is.na(out)] <- 0
    if(toupper(missing)=="MEAN") {
      moy <- apply(out,2,function(c) mean(c,na.rm=TRUE))
      for(j in 1:ncol(out)) {out[,j][is.na(out[,j])] <- moy[j]}
    }
  }

  if(!quiet) cat("\n...done.\n\n")

  return(out)
}) #end makefreq for genpop

