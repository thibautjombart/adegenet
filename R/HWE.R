##################
# HWE.test.genind
##################

HWE.test.genind <- function(x,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix")){

  if(!is.genind(x)) stop("x is not a valid genind object")
  if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
  checkType(x)

  if(!require(genetics)) stop("genetics package is required. Please install it.")
  if(is.null(pop)) pop <- x@pop
  if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))
  res.type <- tolower(res.type[1])
  if(res.type != "full" && res.type != "matrix") stop("unknown res.type specified.")

  kGen <- genind2genotype(x,pop=pop,res.type="list")

  # ftest tests HWE for a locus and a population
  ftest <- function(vec,permut=permut,nperm=nsim){
    temp <- unique(vec)
    temp <- temp[!is.na(temp)]
    if(length(temp) < 2) return(NA)
    if(res.type=="full") {
      res <- HWE.chisq(vec, simulate.p.value=permut, B=nperm)
    } else {
      res <- HWE.chisq(genotype(vec), simulate.p.value=permut, B=nperm)$p.value
    }
    return(res)
  }

  res <- lapply(kGen,function(e) lapply(e,ftest,permut,nsim))

  # clean non-tested elements in the results list
  if(hide.NA && res.type=="full"){
    newres=list()
    tokeep <- which(unlist(lapply(res,function(e) !all(is.na(e)))))
    if(length(tokeep) > 0) for(i in 1:length(tokeep)) {newres[[i]] <- res[[tokeep[i]]]}
    newres <- lapply(newres,function(e) {e[!is.na(e)] })
    names(newres) <- names(res)[tokeep]
    res <- newres
  }

  if(res.type=="matrix"){
    res <- as.data.frame(lapply(res,unlist))
    rnam <- rownames(res)
    rownames(res) <- gsub(".X-squared","",rnam)
    res <- as.matrix(res)
  }

  return(res)
}
