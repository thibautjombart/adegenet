

setMethod("$","genpop",function(x,name) {
  return(slot(x,name))
})


setMethod("$<-","genpop",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})


.drop_alleles <- function(x, toKeep){
  all.vec <- unlist(x@all.names, use.names = FALSE)[toKeep]
  loc.fac <- factor(x@loc.fac[toKeep])
  
  x@all.names <- split(all.vec, loc.fac)
  x@loc.n.all  <- setNames(tabulate(loc.fac), levels(loc.fac))
  x@loc.fac   <- loc.fac
  return(x)
}


###############
# '[' operator
###############
## genind
setMethod("[", signature(x="genind", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., pop=NULL, loc=NULL, treatOther=TRUE, quiet=TRUE, drop=FALSE) {
  
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  
  ## HANDLE 'POP'
  if(!is.null(pop) && !is.null(pop(x))){
    if(is.factor(pop)) pop <- as.character(pop)
    if(!is.character(pop)) pop <- popNames(x)[pop]
    temp <- !pop %in% pop(x)
    if (any(temp)) { # if wrong population specified
      warning(paste("the following specified populations do not exist:", pop[temp]))
    }
    i <- pop(x) %in% pop
  }
  
  ## handle population factor
  if(!is.null(x@pop)) {
    pop <- factor(pop(x)[i])
  } else {
    pop <- NULL
  }
  
  tab       <- x@tab
  old.other <- other(x)
  hier      <- x@strata
  prevcall  <- match.call()
  
  if (x@type == "codom"){
    ## handle loc argument
    if(!is.null(loc)){
      if(is.factor(loc)) loc <- as.character(loc)
      if(!is.character(loc)) loc <- locNames(x)[loc]
      temp <- !loc %in% x@loc.fac
      if (any(temp)) { # if wrong loci specified
        warning(paste("the following specified loci do not exist:", loc[temp]))
      }
      j <- x$loc.fac %in% loc
    } # end loc argument
    if (drop){
      tab    <- tab[i, , ..., drop = FALSE]
      allNb  <- colSums(tab, na.rm=TRUE) # allele absolute frequencies
      toKeep <- (allNb > 1e-10)
      j      <- j & toKeep
      tab    <- tab[, j, ..., drop=FALSE]
    } else {
      tab <- tab[i, j, ..., drop=FALSE]
    }
  } else { # PA case
    tab <- tab[i, j, ..., drop = FALSE]
  }
  
  
  ## handle 'other' slot
  nOther <- length(x@other)
  namesOther <- names(x@other)
  counter <- 0
  if(treatOther){
    f1 <- function(obj,n=nrow(x@tab)){
      counter <<- counter+1
      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
        obj <- obj[i,,drop=FALSE]
      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
        obj <- obj[i]
        if(is.factor(obj)) {obj <- factor(obj)}
      } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}
      
      return(obj)
    } # end f1
    
    x@other <- lapply(x@other, f1) # treat all elements
    
  } else {
    other(x) <- old.other
  } # end treatOther
  
  x@tab    <- tab
  x@pop    <- pop
  x@call   <- prevcall
  x@type   <- x@type
  
  # Treat sample and strata
  x@ploidy    <- x@ploidy[i]
  x@hierarchy <- x@hierarchy
  x@strata    <- hier[i, , drop = FALSE]
  
  if (x@type == "codom"){
    # Treat locus items
    x <- .drop_alleles(x, j)
  }
  return(x)
})





## genpop
setMethod("[", "genpop", function(x, i, j, ..., loc=NULL, treatOther=TRUE, drop=FALSE) {
  
  if (missing(i)) i <- TRUE
  if (missing(j)) j <- TRUE
  
  tab <- x@tab
  old.other <- other(x)
  
  
  ## handle loc argument
  if(!is.null(loc)){
    if(is.factor(loc)) loc <- as.character(loc)
    if(!is.character(loc)) loc <- locNames(x)[loc]
    temp <- !loc %in% x@loc.fac
    if(any(temp)) { # si mauvais loci
      warning(paste("the following specified loci do not exist:", loc[temp]))
    }
    j <- x$loc.fac %in% loc
  } # end loc argument
  
  prevcall <- match.call()
  # tab <- tab[i, j, ...,drop=FALSE]
  
  # if(drop){
  #     allNb <- apply(tab, 2, sum, na.rm=TRUE) # allele absolute frequencies
  #     toKeep <- (allNb > 1e-10)
  #     tab <- tab[,toKeep, drop=FALSE]
  # }
  if (drop){
    tab    <- tab[i, , ..., drop = FALSE]
    allNb  <- colSums(tab, na.rm=TRUE) # allele absolute frequencies
    toKeep <- (allNb > 1e-10)
    j      <- j & toKeep
    tab    <- tab[, j, ..., drop=FALSE]
  } else {
    tab <- tab[i, j, ..., drop=FALSE]
  }
  
  # res <- genpop(tab,prevcall=prevcall,ploidy=x@ploidy)
  
  ## handle 'other' slot
  nOther <- length(x@other)
  namesOther <- names(x@other)
  counter <- 0
  if(treatOther){
    f1 <- function(obj,n=nrow(x@tab)){
      counter <<- counter+1
      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
        obj <- obj[i,,drop=FALSE]
      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
        obj <- obj[i]
        if(is.factor(obj)) {obj <- factor(obj)}
      } else {warning(paste("cannot treat the object",namesOther[counter]))}
      
      return(obj)
    } # end f1
    
    x@other <- lapply(x@other, f1) # treat all elements
    
  } else {
    other(x) <- old.other
  } # end treatOther
  
  x@tab    <- tab
  x@call   <- prevcall
  x@type   <- x@type
  
  # Treat populations
  x@ploidy    <- x@ploidy
  
  # Treat locus items
  x <- .drop_alleles(x, j)
  
  return(x)
})


##########################
## Method show for genind
##########################
setMethod ("show", "genind", function(object){
  ## GET USEFUL VARIABLES
  x <- object
  indTxt <- ifelse(nInd(x)>1, "individuals;", "individual;")
  locTxt <- ifelse(nLoc(x)>1, "loci;", "locus;")
  allTxt <- ifelse(ncol(tab((x)))>1, "alleles;", "allele;")
  
  ## HEADER
  cat("/// GENIND OBJECT /////////")
  
  cat("\n\n //", format(nInd(x), big.mark=","), indTxt,
      format(nLoc(x), big.mark=","), locTxt,
      format(ncol(tab(x)), big.mark=","), allTxt,
      "size:", format(object.size(x), units="auto"))
  
  ## BASIC CONTENT
  cat("\n\n // Basic content")
  p <- ncol(x@tab)
  len <- 7
  
  cat("\n   @tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of allele counts" )
  
  if(!is.null(x@loc.n.all)){
    alleletxt <- paste("(range: ", paste(range(x@loc.n.all), collapse="-"), ")", sep="")
    cat("\n   @loc.n.all: number of alleles per locus", alleletxt)
  }
  
  if(!is.null(x@loc.fac)){
    cat("\n   @loc.fac: locus factor for the", ncol(x@tab), "columns of @tab")
  }
  if(!is.null(x@all.names)){
    cat("\n   @all.names: list of allele names for each locus")
  }
  
  ploidytxt <- paste("(range: ", paste(range(x@ploidy), collapse="-"), ")", sep="")
  cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  cat("\n   @type: ",x@type)
  cat("\n   @call: ")
  print(x@call)
  
  ## OPTIONAL CONTENT
  cat("\n // Optional content")
  optional <- FALSE
  if(!is.null(x@pop)){
    optional <- TRUE
    poptxt <- paste("(group size range: ", paste(range(table(x@pop)), collapse="-"), ")", sep="")
    cat("\n   @pop:", paste("population of each individual", poptxt))
  }
  
  if (!is.null(x@strata)){
    optional <- TRUE
    cat("\n   @strata: ")
    levs <- names(x@strata)
    if (length(levs) > 6){
      levs <- paste(paste(head(levs), collapse = ", "), "...", sep = ", ")
    } else {
      levs <- paste(levs, collapse = ", ")
    }
    cat("a data frame with", length(x@strata), "columns (", levs, ")")
  }
  
  if (!is.null(x@hierarchy)){
    optional <- TRUE
    cat("\n   @hierarchy:", paste(x@hierarchy, collapse = ""))
  }
  
  if(!is.null(x@other)){
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  }
  
  if(!optional) cat("\n   - empty -")
  
  cat("\n")
}
) # end show method for genind




##########################
## Method show for genpop
##########################
setMethod ("show", "genpop", function(object){
  ## GET USEFUL VARIABLES
  x <- object
  popTxt <- ifelse(nPop(x)>1, "populations;", "population;")
  locTxt <- ifelse(nLoc(x)>1, "loci;", "locus;")
  allTxt <- ifelse(ncol(tab((x)))>1, "alleles;", "allele;")
  
  ## HEADER
  cat("/// GENPOP OBJECT /////////")
  
  cat("\n\n //", format(nPop(x), big.mark=","), popTxt,
      format(nLoc(x), big.mark=","), locTxt,
      format(ncol(tab(x)), big.mark=","), allTxt,
      "size:", format(object.size(x), units="auto"))
  
  ## BASIC CONTENT
  cat("\n\n // Basic content")
  p <- ncol(x@tab)
  len <- 7
  
  cat("\n   @tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of allele counts" )
  
  if(!is.null(x@loc.n.all)){
    alleletxt <- paste("(range: ", paste(range(x@loc.n.all), collapse="-"), ")", sep="")
    cat("\n   @loc.n.all: number of alleles per locus", alleletxt)
  }
  
  if(!is.null(x@loc.fac)){
    cat("\n   @loc.fac: locus factor for the", ncol(x@tab), "columns of @tab")
  }
  if(!is.null(x@all.names)){
    cat("\n   @all.names: list of allele names for each locus")
  }
  
  ploidytxt <- paste("(range: ", paste(range(x@ploidy), collapse="-"), ")", sep="")
  cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  cat("\n   @type: ",x@type)
  cat("\n   @call: ")
  print(x@call)
  
  ## OPTIONAL CONTENT
  cat("\n // Optional content")
  optional <- FALSE
  if(!is.null(x@other)){
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  }
  
  if(!optional) cat("\n   - empty -")
  
  cat("\n")
  
}
) # end show method for genpop





############################
# Method summary for genind
############################
if(!isGeneric("summary")){
  setGeneric("summary", function(object, ...) standardGeneric("summary"))
}
setMethod ("summary", signature(object="genind"), function(object, verbose = TRUE, ...){
  x <- object
  if(!is.genind(x)) stop("Provided object is not a valid genind.")
  
  
  if(is.null(x@pop)){
    x@pop <- factor(rep(1,nrow(x@tab)))
    x@pop.names <- ""
    names(x@pop.names) <- "P1"
  }
  
  ## BUILD THE OUTPUT ##
  ## type-independent stuff
  res <- list()
  
  res$N <- nrow(x@tab)
  
  res$pop.eff <- as.numeric(table(x@pop))
  names(res$pop.eff) <- x@pop.names
  
  ## PA case ##
  if(x@type=="PA"){
    ## % of missing data
    res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab))
    
    ## display and return
    listlab <- c("# Total number of genotypes: ",
                 "# Population sample sizes: ",
                 "# Percentage of missing data: ")
    
    if (verbose) {
      cat("\n",listlab[1],res[[1]],"\n")
      for(i in 2:3){
        cat("\n",listlab[i],"\n")
        print(res[[i]])
      }
    }
    
    return(invisible(res))
  }
  
  
  ## codom case ##
  res$loc.nall <- x@loc.nall
  
  temp <- genind2genpop(x,quiet=TRUE)@tab
  
  res$pop.nall <- apply(temp,1,function(r) sum(r!=0,na.rm=TRUE))
  
  ##  res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab)) <- wrong
  res$NA.perc <- 100*(1-mean(propTyped(x,by="both")))
  
  ## handle heterozygosity
  if(any(x@ploidy > 1)){
    ## auxiliary function to compute observed heterozygosity
    temp=lapply(seploc(x),tab, freq=TRUE)
    f1 <- function(tab){
      H <- apply(tab, 1, function(vec) any(vec > 0 & vec < 1))
      H <- mean(H,na.rm=TRUE)
      return(H)
    }
    
    res$Hobs <- unlist(lapply(temp,f1))
    
    ## auxiliary function to compute expected heterozygosity
    ## freq is a vector of frequencies
    f2 <- function(freq){
      H <- 1-sum(freq*freq,na.rm=TRUE)
      return(H)
    }
    
    temp <- genind2genpop(x,pop=rep(1,nInd(x)),quiet=TRUE)
    temp <- tab(temp, freq=TRUE, quiet=TRUE)
    res$Hexp <-tapply(temp^2, x@loc.fac, function(e) 1-sum(e, na.rm=TRUE))
  } else { # no possible heterozygosity for haploid genotypes
    res$Hobs <- 0
    res$Xexp <- 0
  }
  
  ## print to screen
  listlab <- c("# Total number of genotypes: ",
               "# Population sample sizes: ",
               "# Number of alleles per locus: ",
               "# Number of alleles per population: ",
               "# Percentage of missing data: ",
               "# Observed heterozygosity: ",
               "# Expected heterozygosity: ")
  
  if (verbose) {
    cat("\n",listlab[1],res[[1]],"\n")
    for(i in 2:7){
      cat("\n",listlab[i],"\n")
      print(res[[i]])
    }
  }
  
  return(invisible(res))
})  # end summary.genind





############################
# Method summary for genpop
############################
setMethod ("summary", signature(object="genpop"), function(object, verbose = TRUE, ...){
  x <- object
  if(!inherits(x,"genpop")) stop("To be used with a genpop object")
  
  ## BUILD THE OUTPUT ##
  ## type-independent stuff
  res <- list()
  
  res$npop <- nrow(x@tab)
  
  ## PA case ##
  if(x@type=="PA"){
    ## % of missing data
    res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab))
    
    ## display and return
    listlab <- c("# Total number of genotypes: ",
                 "# Percentage of missing data: ")
    
    if (verbose) {
      cat("\n",listlab[1],res[[1]],"\n")
      for(i in 2){
        cat("\n",listlab[i],"\n")
        print(res[[i]])
      }
    }
    
    return(invisible(res))
  }
  
  
  ## codom case ##
  res$loc.n.all <- x@loc.n.all
  
  res$pop.nall <- apply(x@tab,1,function(r) sum(r>0,na.rm=TRUE))
  
  ##  res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab)) <- old version
  mean.w <- function(x,w=rep(1/length(x),length(x))){
    x <- x[!is.na(x)]
    w <- w[!is.na(x)]
    w <- w/sum(w)
    return(sum(x*w))
  }
  
  w <- apply(x@tab,1,sum,na.rm=TRUE) # weights for populations
  res$NA.perc <- 100*(1-mean.w(propTyped(x), w=w))
  ## res$NA.perc <- 100*(1-mean(propTyped(x,by="both"))) <- old
  
  # print to screen
  listlab <- c("# Number of populations: ",
               "# Number of alleles per locus: ",
               "# Number of alleles per population: ",
               "# Percentage of missing data: ")
  
  if (verbose) {
    cat("\n",listlab[1],res[[1]],"\n")
    for(i in 2:4){
      cat("\n",listlab[i],"\n")
      print(res[[i]])
    }
  }
  
  return(invisible(res))
  
}
)# end summary.genpop



#} # end .initAdegenetClasses()






###############
# Methods "is"
###############
is.genind <- function(x){
  res <- ( is(x, "genind") & validObject(x))
  return(res)
}

is.genpop <- function(x){
  res <- ( is(x, "genpop") & validObject(x))
  return(res)
}



.hasUniquePloidy <- function(x){
  return(length(unique(x@ploidy))==1)
}
