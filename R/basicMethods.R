

setMethod("$","genpop",function(x,name) {
  return(slot(x,name))
})


setMethod("$<-","genpop",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})


.drop_alleles <- function(x, toKeep){
  all.vec <- unlist(alleles(x), use.names = FALSE)[toKeep]
  loc.fac <- factor(locFac(x)[toKeep])
  present_alleles <- colSums(tab(x), na.rm = TRUE) > 0L

  x@all.names <- split(all.vec, loc.fac)
  x@loc.n.all <- setNames(tabulate(loc.fac), levels(loc.fac)) #vapply(split(present_alleles, loc.fac), sum, integer(1))
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

  ## HANDLE 'i' as character
  if(is.factor(i)) i <- as.character(i)
  if(is.character(i)){
      old.i <- i
      i <- match(i, indNames(x))
      if(any(is.na(i))){
          warning(paste("the following specified individuals do not exist:", paste0(old.i[is.na(i)], collapse = ", ")), call. = FALSE)
          i <- i[!is.na(i)]
          if(length(i)==0) {
              warning("no individual selected - ignoring", call. = FALSE)
              i <- TRUE
          }
      }
  }

  ## HANDLE 'POP'
  if(!is.null(pop) && !is.null(pop(x))){
    if(is.factor(pop)) pop <- as.character(pop)
    if(!is.character(pop)) pop <- popNames(x)[pop]
    temp <- !pop %in% pop(x)
    if (any(temp)) { # if wrong population specified
      warning(paste("the following specified populations do not exist:", paste0(pop[temp], collapse = ", ")), call. = FALSE)
      if (all(temp)){
        warning("no populations selected - ignoring", call. = FALSE)
        pop <- pop(x)
      }
    }
    i <- pop(x) %in% pop
  }

  ## handle population factor
  if(!is.null(pop(x))) {
    pop <- factor(pop(x)[i])
  } else {
    pop <- NULL
  }

  tab       <- tab(x)
  old.other <- other(x)
  hier      <- x@strata
  prevcall  <- match.call()

  if (x@type == "codom"){
    ## handle loc argument
    if(!is.null(loc)){
      if(is.factor(loc)) loc <- as.character(loc)
      if(!is.character(loc)) loc <- locNames(x)[loc]
      temp <- !loc %in% locFac(x)
      if (any(temp)) { # if wrong loci specified
        warning(paste("the following specified loci do not exist:", paste0(loc[temp], collapse = ", ")), call. = FALSE)
        if (all(temp)){
          warning("no loci selected - ignoring", call. = FALSE)
          loc <- x@loc.fac
        }
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
  nOther <- length(other(x))
  namesOther <- names(other(x))
  counter <- 0
  if(treatOther){
    f1 <- function(obj,n=nrow(tab(x))){
      counter <<- counter+1
      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
        obj <- obj[i,,drop=FALSE]
      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
        obj <- obj[i]
        if(is.factor(obj)) {obj <- factor(obj)}
      } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}

      return(obj)
    } # end f1

    x@other <- lapply(other(x), f1) # treat all elements

  } else {
    other(x) <- old.other
  } # end treatOther

  x@tab    <- tab
  x@pop    <- pop
  x@call   <- prevcall
  x@type   <- x@type

  # Treat sample and strata
  x@ploidy    <- ploidy(x)[i]
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

  tab <- tab(x)
  old.other <- other(x)


  ## HANDLE 'i' as character
  if(is.factor(i)) i <- as.character(i)
  if(is.character(i)){
      old.i <- i
      i <- match(i, popNames(x))
      if(any(is.na(i))){
          warning(paste("the following specified populations do not exist:", paste0(old.i[is.na(i)], collapse = ", ")), call. = FALSE)
          i <- i[!is.na(i)]
          if(length(i)==0) {
              warning("no population selected - ignoring", call. = FALSE)
              i <- TRUE
          }
      }
  }

  ## handle loc argument
  if(!is.null(loc)){
    if(is.factor(loc)) loc <- as.character(loc)
    if(!is.character(loc)) loc <- locNames(x)[loc]
    temp <- !loc %in% locFac(x)
    if(any(temp)) { # si mauvais loci
      warning(paste("the following specified loci do not exist:", paste0(loc[temp], collapse = ", ")), call. = FALSE)
    }
    if (all(temp)){
      warning("no loci selected - ignoring", call. = FALSE)
      loc <- x@loc.fac
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
  nOther <- length(other(x))
  namesOther <- names(other(x))
  counter <- 0
  if(treatOther){
    f1 <- function(obj,n=nrow(tab(x))){
      counter <<- counter+1
      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
        obj <- obj[i,,drop=FALSE]
      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
        obj <- obj[i]
        if(is.factor(obj)) {obj <- factor(obj)}
      } else {warning(paste("cannot treat the object",namesOther[counter]), call. = FALSE)}

      return(obj)
    } # end f1

    x@other <- lapply(other(x), f1) # treat all elements

  } else {
    other(x) <- old.other
  } # end treatOther

  x@tab    <- tab
  x@call   <- prevcall
  x@type   <- x@type

  # Treat populations
  x@ploidy    <- ploidy(x)

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
  p <- ncol(tab(x))
  len <- 7

  cat("\n   @tab: ", nrow(tab(x)), "x", ncol(tab(x)), "matrix of allele counts" )

  if (!is.null(nAll(x))){
    alleletxt <- paste("(range: ", paste(range(nAll(x, onlyObserved = FALSE)), collapse="-"), ")", sep="")
    cat("\n   @loc.n.all: number of alleles per locus", alleletxt)
  }

  if(!is.null(locFac(x))){
    cat("\n   @loc.fac: locus factor for the", ncol(tab(x)), "columns of @tab")
  }
  if(!is.null(alleles(x))){
    cat("\n   @all.names: list of allele names for each locus")
  }

  ploidytxt <- paste("(range: ", paste(range(ploidy(x)), collapse="-"), ")", sep="")
  cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  cat("\n   @type: ",x@type)
  cat("\n   @call: ")
  print(x@call)

  ## OPTIONAL CONTENT
  cat("\n // Optional content")
  optional <- FALSE
  if(!is.null(pop(x))){
    optional <- TRUE
    poptxt <- paste("(group size range: ", paste(range(table(pop(x))), collapse="-"), ")", sep="")
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

  if(!is.null(other(x))){
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(is.null(names(other(x))), "elements without names", paste(names(other(x)), collapse= "  ")), "\n")
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
  p <- ncol(tab(x))
  len <- 7

  cat("\n   @tab: ", nrow(tab(x)), "x", ncol(tab(x)), "matrix of allele counts" )

  if(!is.null(nAll(x))){
    alleletxt <- paste("(range: ", paste(range(nAll(x, onlyObserved = FALSE)), collapse="-"), ")", sep="")
    cat("\n   @loc.n.all: number of alleles per locus", alleletxt)
  }

  if(!is.null(locFac(x))){
    cat("\n   @loc.fac: locus factor for the", ncol(tab(x)), "columns of @tab")
  }
  if(!is.null(alleles(x))){
    cat("\n   @all.names: list of allele names for each locus")
  }

  ploidytxt <- paste("(range: ", paste(range(ploidy(x)), collapse="-"), ")", sep="")
  cat("\n   @ploidy: ploidy of each individual ", ploidytxt)
  cat("\n   @type: ",x@type)
  cat("\n   @call: ")
  print(x@call)

  ## OPTIONAL CONTENT
  cat("\n // Optional content")
  optional <- FALSE
  if(!is.null(other(x))){
    optional <- TRUE
    cat("\n   @other: ")
    cat("a list containing: ")
    cat(ifelse(is.null(names(other(x))), "elements without names", paste(names(other(x)), collapse= "  ")), "\n")
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


  if(is.null(pop(x))){
    pop(x) <- rep("P1", nInd(x))
  }

  ## BUILD THE OUTPUT ##
  ## type-independent stuff
  res <- list()

  res$n <- nrow(tab(x))

  res$n.by.pop <- as.numeric(table(pop(x)))
  names(res$n.by.pop) <- popNames(x)

  ## PA case ##
  if(x@type=="PA"){
    ## % of missing data
    res$NA.perc <- 100*sum(is.na(tab(x)))/prod(dim(tab(x)))

    return(invisible(res))
  }


  ## codom case ##
  res$loc.n.all <- nAll(x, onlyObserved = TRUE)

  temp <- tab(genind2genpop(x,quiet=TRUE))

  res$pop.n.all <- apply(temp,1,function(r) sum(r!=0,na.rm=TRUE))

  res$NA.perc <- 100*(1-mean(propTyped(x,by="both")))

  ## handle heterozygosity
  if(any(ploidy(x) > 1)){
    ## auxiliary function to compute observed heterozygosity
    temp <- lapply(seploc(x),tab, freq=TRUE)
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
    res$Hexp <-tapply(temp^2, locFac(x), function(e) 1-sum(e, na.rm=TRUE))
  } else { # no possible heterozygosity for haploid genotypes
    res$Hobs <- 0
    res$Xexp <- 0
  }

  ## add class and return
  class(res) <- "genindSummary"
  return(res)
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

  res$n.pop <- nrow(tab(x))

  ## PA case ##
  if(x@type=="PA"){
    ## % of missing data
    res$NA.perc <- 100*sum(is.na(tab(x)))/prod(dim(tab(x)))

    ## add class and return
    class(res) <- "genpopSummary"
    return(res)
  }


  ## codom case ##
  res$loc.n.all <- nAll(x, onlyObserved = TRUE)

  res$pop.n.all <- apply(tab(x),1,function(r) sum(r>0,na.rm=TRUE))

  ##  res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab)) <- old version
  mean.w <- function(x,w=rep(1/length(x),length(x))){
    x <- x[!is.na(x)]
    w <- w[!is.na(x)]
    w <- w/sum(w)
    return(sum(x*w))
  }

  w <- apply(tab(x),1,sum,na.rm=TRUE) # weights for populations
  res$NA.perc <- 100*(1-mean.w(propTyped(x), w=w))

  ## add class and return
  class(res) <- "genpopSummary"
  return(res)
}
)# end summary.genpop




#######################
## print for summaries
#######################
print.genindSummary <- function(x, ...){
    if(!is.null(x$n)) cat("\n// Number of individuals:", x$n)
    if(!is.null(x$n.by.pop)) cat("\n// Group sizes:", x$n.by.pop)
    if(!is.null(x$loc.n.all)) cat("\n// Number of alleles per locus:", x$loc.n.all)
    if(!is.null(x$pop.n.all)) cat("\n// Number of alleles per group:", x$pop.n.all)
    if(!is.null(x$NA.perc)) cat("\n// Percentage of missing data:", round(x$NA.perc,2), "%")
    if(!is.null(x$Hobs)) cat("\n// Observed heterozygosity:", round(x$Hobs,2))
    if(!is.null(x$Hexp)) cat("\n// Expected heterozygosity:", round(x$Hexp,2))
    cat("\n")
} # end print.genindSummary


print.genpopSummary <- function(x, ...){
    if(!is.null(x$n.pop)) cat("\n// Number of populations:", x$n.pop)
    if(!is.null(x$loc.n.all)) cat("\n// Number of alleles per locus:", x$loc.n.all)
    if(!is.null(x$pop.n.all)) cat("\n// Number of alleles per group:", x$pop.n.all)
    if(!is.null(x$NA.perc)) cat("\n// Percentage of missing data:", round(x$NA.perc,2), "%")
    cat("\n")
} # end print.genpopSummary



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
  return(length(unique(ploidy(x)))==1)
}
