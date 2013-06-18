

setMethod("$","genpop",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","genpop",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})





###############
# '[' operator
###############
## genind
setMethod("[", signature(x="genind", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., loc=NULL, treatOther=TRUE, quiet=TRUE, drop=FALSE) {

    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    pop <- NULL
    if(is.null(x@pop)) { tab <- truenames(x) }
    if(!is.null(x@pop)) {
        temp <- truenames(x)
        tab <- temp$tab
        pop <- temp$pop
        pop <- factor(pop[i])
    }

    old.other <- other(x)

    ## handle loc argument
    if(!is.null(loc)){
        loc <- as.character(loc)
        temp <- !loc %in% x@loc.fac
        if(any(temp)) { # si mauvais loci
            warning(paste("the following specified loci do not exist:", loc[temp]))
        }
        j <- x$loc.fac %in% loc
    } # end loc argument

    prevcall <- match.call()

    tab <- tab[i, j, ...,drop=FALSE]

    if(drop){
        allNb <- apply(tab, 2, sum, na.rm=TRUE) # allele absolute frequencies
        toKeep <- (allNb > 1e-10)
        tab <- tab[,toKeep, drop=FALSE]
    }

    res <- genind(tab,pop=pop,prevcall=prevcall, ploidy=x@ploidy, type=x@type)

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

        res@other <- lapply(x@other, f1) # treat all elements

    } else {
        other(res) <- old.other
    } # end treatOther

    return(res)
})





## genpop
setMethod("[", "genpop", function(x, i, j, ..., loc=NULL, treatOther=TRUE, drop=FALSE) {

    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    tab <- truenames(x)
    old.other <- other(x)


    ## handle loc argument
    if(!is.null(loc)){
        loc <- as.character(loc)
        temp <- !loc %in% x@loc.fac
        if(any(temp)) { # si mauvais loci
            warning(paste("the following specified loci do not exist:", loc[temp]))
        }
        j <- x$loc.fac %in% loc
    } # end loc argument

    prevcall <- match.call()
    tab <- tab[i, j, ...,drop=FALSE]

    if(drop){
        allNb <- apply(tab, 2, sum, na.rm=TRUE) # allele absolute frequencies
        toKeep <- (allNb > 1e-10)
        tab <- tab[,toKeep, drop=FALSE]
    }

    res <- genpop(tab,prevcall=prevcall)

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

        res@other <- lapply(x@other, f1) # treat all elements

    } else {
        other(res) <- old.other
    } # end treatOther


    return(res)
})




##########################
# Method show for genind
##########################
setMethod ("show", "genind", function(object){
  x <- object
  cat("\n")
  cat("   #####################\n")
  cat("   ### Genind object ### \n")
  cat("   #####################")
  cat("\n- genotypes of individuals - \n")
  cat("\nS4 class: ", as.character(class(x)))

  cat("\n@call: ")
  print(x@call)

  p <- ncol(x@tab)
  len <- 7

  cat("\n@tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of genotypes\n" )

  cat("\n@ind.names: vector of ", length(x@ind.names), "individual names")
  cat("\n@loc.names: vector of ", length(x@loc.names), "locus names")

  if(!is.null(x@loc.nall)){
      cat("\n@loc.nall: number of alleles per locus")
  } else {
      cat("\n@loc.nall: NULL")
  }

  if(!is.null(x@loc.fac)){
      cat("\n@loc.fac: locus factor for the ", ncol(x@tab), "columns of @tab")
  } else {
      cat("\n@loc.fac: NULL")
  }

  if(!is.null(x@all.names)){
      cat("\n@all.names: list of ", length(x@all.names), "components yielding allele names for each locus")
  } else {
      cat("\n@all.names: NULL")
  }

  cat("\n@ploidy: ",x@ploidy)
  cat("\n@type: ",x@type)

  cat("\n\nOptionnal contents: ")
  cat("\n@pop: ", ifelse(is.null(x@pop), "- empty -", "factor giving the population of each individual"))
  cat("\n@pop.names: ", ifelse(is.null(x@pop.names), "- empty -", "factor giving the population of each individual"))

  cat("\n\n@other: ")
  if(!is.null(x@other)){
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  } else {
    cat("- empty -\n")
  }

  cat("\n")
}
) # end show method for genind




##########################
# Method show for genpop
##########################
setMethod ("show", "genpop", function(object){
  x <- object
  cat("\n")
  cat("       #####################\n")
  cat("       ### Genpop object ### \n")
  cat("       #####################")
  cat("\n- Alleles counts for populations - \n")
  cat("\nS4 class: ", as.character(class(x)))

  cat("\n@call: ")
  print(x@call)

  p <- ncol(x@tab)

  cat("\n@tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of alleles counts\n" )

  cat("\n@pop.names: vector of ", length(x@pop.names), "population names")
  cat("\n@loc.names: vector of ", length(x@loc.names), "locus names")

  if(!is.null(x@loc.nall)){
      cat("\n@loc.nall: number of alleles per locus")
  } else {
      cat("\n@loc.nall: NULL")
  }

  if(!is.null(x@loc.fac)){
      cat("\n@loc.fac: locus factor for the ", ncol(x@tab), "columns of @tab")
  } else {
      cat("\n@loc.fac: NULL")
  }

  if(!is.null(x@all.names)){
      cat("\n@all.names: list of ", length(x@all.names), "components yielding allele names for each locus")
  } else {
      cat("\n@all.names: NULL")
  }

  cat("\n@ploidy: ",x@ploidy)
  cat("\n@type: ",x@type)

  cat("\n\n@other: ")
  if(!is.null(x@other)){
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  } else {
    cat("- empty -\n")
  }

  cat("\n")

}
) # end show method for genpop





############################
# Method summary for genind
############################
if(!isGeneric("summary")){
    setGeneric("summary", function(object, ...) standardGeneric("summary"))
}
setMethod ("summary", signature(object="genind"), function(object, ...){
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
      cat("\n",listlab[1],res[[1]],"\n")
      for(i in 2:3){
          cat("\n",listlab[i],"\n")
          print(res[[i]])
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
  if(x@ploidy > 1){
      ## auxiliary function to compute observed heterozygosity
      temp <- seploc(x,truenames=FALSE,res.type="matrix")
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

      temp <- genind2genpop(x,pop=rep(1,nrow(x@tab)),quiet=TRUE)
      temp <- makefreq(temp,quiet=TRUE)$tab
      temp.names <- colnames(temp)
      temp <- as.vector(temp)
      names(temp) <- temp.names
      temp <- split(temp,x@loc.fac)
      ## temp is a list of alleles frequencies (one element per locus)

      res$Hexp <- unlist(lapply(temp,f2))
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
  cat("\n",listlab[1],res[[1]],"\n")
  for(i in 2:7){
    cat("\n",listlab[i],"\n")
    print(res[[i]])
  }

  return(invisible(res))
}) # end summary.genind





############################
# Method summary for genpop
############################
setMethod ("summary", signature(object="genpop"), function(object, ...){
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
      cat("\n",listlab[1],res[[1]],"\n")
      for(i in 2){
          cat("\n",listlab[i],"\n")
          print(res[[i]])
      }

      return(invisible(res))
  }


  ## codom case ##
  res$loc.nall <- x@loc.nall

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
  cat("\n",listlab[1],res[[1]],"\n")
  for(i in 2:4){
    cat("\n",listlab[i],"\n")
    print(res[[i]])
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

