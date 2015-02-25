####################
# Function makefreq
####################


#' Function to generate allelic frequencies
#' 
#' The function \code{makefreq} generates a table of allelic frequencies from
#' an object of class \code{genpop}.
#' 
#' There are 3 treatments for missing values: \cr - NA: kept as NA.\cr - 0:
#' missing values are considered as zero. Recommended for a PCA on
#' compositionnal data.\cr - "mean": missing values are given the mean
#' frequency of the corresponding allele. Recommended for a centred PCA.
#' 
#' @param x an object of class \code{genpop}.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @param missing treatment for missing values. Can be NA, 0 or "mean" (see
#' details)
#' @param truenames a logical indicating whether true labels (as opposed to
#' generic labels) should be used to name the output.
#' @return Returns a list with the following components: \item{tab}{matrix of
#' allelic frequencies (rows: populations; columns: alleles).}
#' \item{nobs}{number of observations (i.e. alleles) for each population x
#' locus combinaison.} \item{call}{the matched call}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{genpop}}
#' @keywords manip multivariate
#' @examples
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
#' ca1 <- dudi.coa(as.data.frame(Xcount@tab),scannf=FALSE)
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
makefreq <- function(x,quiet=FALSE,missing=NA,truenames=TRUE){

  if(!is.genpop(x)) stop("x is not a valid genpop object")
  ##if(x@type=="PA") stop("frequencies not computable for presence/asbence data")
  checkType(x)

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

  # tabfreq is a pop x loci table of allelic frequencies
  tabfreq <- t(apply(tabcount,1,function(r) unlist(tapply(r,x@loc.fac,f1))))
  if(length(x@loc.nall)==1 && x@loc.nall[1]==1) tabfreq <- t(tabfreq) # matrix is transposed by apply if there's a single allele
  colnames(tabfreq) <- colnames(x@tab)

  # NA treatment
  # NA can be kept as is, or replaced 0 or by the mean frequency of the allele.
  if(!is.na(missing)){
    if(missing==0) tabfreq[is.na(tabfreq)] <- 0
    if(toupper(missing)=="MEAN") {
      moy <- apply(tabfreq,2,function(c) mean(c,na.rm=TRUE))
      for(j in 1:ncol(tabfreq)) {tabfreq[,j][is.na(tabfreq[,j])] <- moy[j]}
    }
  }

  if(!quiet) cat("\n...done.\n\n")

  res$tab <- tabfreq
  res$nobs <- eff.pop
  res$call <- match.call()

  ## handle truenames
  if(truenames){
      temp <- rep(x@loc.names,x@loc.nall)
      colnames(res$tab) <- paste(temp,unlist(x@all.names),sep=".")
      rownames(res$tab) <- x@pop.names

      colnames(res$nobs) <- x@loc.names
      rownames(res$nobs) <- x@pop.names
  }

  return(res)
} #end makefreq

