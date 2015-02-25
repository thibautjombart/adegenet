#########################
# Function genind2genpop
#########################
genind2genpop <- function(x,pop=NULL,missing=c("NA","0","chi2"),quiet=FALSE,
                          process.other=FALSE, other.action=mean){

  if(!is.genind(x)) stop("x is not a valid genind object")
  checkType(x)

  if(is.null(x@pop) && is.null(pop)) {
      if(!quiet) warning("\npop is not provided either in x or in pop - assuming one single group")
      pop <- factor(rep(1, nrow(x@tab)))
  }

  missing <- match.arg(missing)

  if(!quiet) cat("\n Converting data from a genind to a genpop object... \n")

  ## choose pop argument over x@pop
   if(!is.null(pop)) {
    if(length(pop) != nrow(x@tab)) stop("inconsistent length for factor pop")
    # keep levels in order of appearance
    pop <- as.character(pop)
    pop <- factor(pop, levels=unique(pop))
  } else {
    pop <- pop(x)
    # keep levels in order of appearance
    ##pop <- as.character(pop)
    ##pop <- factor(pop, levels=unique(pop))
    ##if(!is.null(x@pop.names)) levels(pop) <- x@pop.names # restore real names
  }

  # make generic pop labels, store real pop names
  # pop.names <- levels(pop) ## no longer used

  # tabcount is a matrix pop x alleles, counting alleles per pop
  # *ploidy to have alleles counts
  f1 <- function(v){
    if(all(is.na(v))) return(NA) else return(sum(v,na.rm=TRUE))
  }

  f2 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(NA)
    return(v/(sum(v,na.rm=TRUE)))
  }

  tabcount <- x@ploidy * apply(x@tab,2,function(c) tapply(c,pop,f1))
  tabcount <- round(tabcount,digits=0)
  # restitute matrix class when only one pop
  if(is.null(dim(tabcount))) {
    lab.col <- names(tabcount)
    tabcount <- matrix(tabcount,nrow=1)
    colnames(tabcount) <- lab.col
    rownames(tabcount) <- levels(pop)[1]
  }

  ## make final object
  if(x@type=="codom"){
      temp <- paste(rep(x@loc.names,x@loc.nall),unlist(x@all.names),sep=".")
  } else{
      temp <- x@loc.names
  }
  colnames(tabcount) <- temp

  prevcall <- match.call()

  res <- genpop(tab=tabcount, prevcall=prevcall, ploidy=x@ploidy, type=x@type)

  ## handle @other here
  res@other <- x@other
  if(process.other){
      ## auxiliary function doing the job
      fOther <- function(e){
          N <- nrow(x@tab)
          if(is.vector(e) && is.numeric(e) && length(e)==N){ # numeric vector
              res <- tapply(e, pop, other.action)
              return(res)
          } else if(is.matrix(e) && is.numeric(e) && nrow(e)==N){ # numeric matrix
              res <- apply(e, 2, function(vec) tapply(vec, pop, other.action))
              colnames(res) <- colnames(e)
              return(res)
          } else if(is.data.frame(e) && nrow(e)==N && all(sapply(e,is.numeric)) ){ # df of numeric vectors
              res <- lapply(e, function(vec) tapply(vec, pop, other.action))
              res <- data.frame(res)
              names(res) <- names(e)
              return(res)
          } else return(e)
      } # end fOther

      res@other <- lapply(res@other, fOther)
  } # end if(process.other)

  if(missing != "NA"){
      res <- na.replace(res, method=missing, quiet=quiet)
  }

  if(!quiet) cat("\n...done.\n\n")

  return(res)

} # end genind2genpop
