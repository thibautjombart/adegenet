####################
# Function makefreq
####################
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

