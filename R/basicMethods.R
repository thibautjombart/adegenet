

#' Accessors for adegenet objects
#'
#' An accessor is a function that allows to interact with slots of an object in
#' a convenient way. Several accessors are available for \linkS4class{genind}
#' or \linkS4class{genpop} objects. The operator "\$" and "\$<-" are used to
#' access the slots, being equivalent to "@@" and "@@<-".\cr
#'
#' The operator "[" can be used to access components of the matrix slot "@@tab",
#' returning a \linkS4class{genind} or \linkS4class{genpop} object. This syntax
#' is the same as for a matrix; for instance:\cr - "obj[,]" returns "obj" \cr -
#' "obj[1:10,]" returns an object with only the first 10 genotypes (if "obj" is
#' a \linkS4class{genind}) or the first 10 populations (if "obj" is a
#' \linkS4class{genpop}) of "obj" \cr - "obj[1:10, 5:10]" returns an object
#' keeping the first 10 entities and the alleles 5 to 10.\cr -
#' "obj[loc=c("L1","L3")]" returns an object keeping only the loci specified in
#' the \code{loc} argument (using generic names, not true names; in this
#' example, only the first and the third locus would be retained)\cr -
#' "obj[1:3, drop=TRUE]" returns the first 3 genotypes/populations of "obj",
#' but retaining only alleles that are present in this subset (as opposed to
#' keeping all alleles of "obj", which is the default behavior).\cr
#'
#' The argument \code{treatOther} handles the treatment of objects in the
#' \code{@@other} slot (see details). The argument \code{drop} can be set to
#' TRUE to drop alleles that are no longer represented in the subset.
#'
#' The "[" operator can treat elements in the \code{@@other} slot as well. For
#' instance, if \code{obj@@other$xy} contains spatial coordinates, the
#' \code{obj[1:3,]@@other$xy} will contain the spatial coordinates of the
#' genotypes (or population) 1,2 and 3. This is handled through the argument
#' \code{treatOther}, a logical defaulting to TRUE. If set to FALSE, the
#' \code{@@other} returned unmodified.\cr
#'
#' Note that only matrix-like, vector-like and lists can be proceeded in
#' \code{@@other}. Other kind of objects will issue a warning an be returned as
#' they are, unless the argument \code{quiet} is left to TRUE, its default
#' value.\cr
#'
#' The \code{drop} argument can be set to TRUE to retain only alleles that are
#' present in the subset. To achieve better control of polymorphism of the
#' data, see \code{\link{isPoly}}.
#'
#' @name Accessors
#' @aliases $,genind-method $,genpop-method $<-,genind-method $<-,genpop-method
#' [,genind-method [,genpop-method nLoc nLoc,genind-method nLoc,genpop-method
#' nInd nInd,genind-method pop pop<- pop,genind-method pop<-,genind-method
#' locNames locNames,genind-method locNames,genpop-method locNames<-
#' locNames<-,genind-method locNames<-,genpop-method indNames
#' indNames,genind-method indNames<- indNames<-,genind-method ploidy
#' ploidy,genind-method ploidy,genpop-method ploidy<- ploidy<-,genind-method
#' ploidy<-,genpop-method alleles alleles,genind-method alleles,genpop-method
#' alleles<- alleles<-,genind-method alleles<-,genpop-method other
#' other,genind-method other,genpop-method other<- other<-,genind-method
#' other<-,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} or a \linkS4class{genpop} object.
#' @param withAlleles a logical indicating whether the result should be of the
#' form [locus name].[allele name], instead of [locus name].
#' @param \dots further arguments to be passed to other methods (currently not
#' used).
#' @return A \linkS4class{genind} or \linkS4class{genpop} object.
#' @section Methods: \describe{ \item{nInd}{returns the number of individuals
#' in the \code{genind} object} \item{nLoc}{returns the number of loci of the
#' object} \item{pop}{returns the population factor of the object, using true
#' (as opposed to generic) levels.} \item{pop<-}{replacement method for the
#' \code{@@pop} slot of an object. The content of \code{@@pop} and
#' \code{@@pop.names} is updated automatically.} \item{indNames}{returns the
#' true names of individuals.} \item{indNames<-}{sets the true names of
#' individuals using a vector of length \code{nInd(x)}.} \item{locNames}{returns the
#' true names of markers and/or alleles.} \item{locNames<-}{sets the true names
#' of markers using a vector of length \code{nLoc(x)}.} \item{ploidy}{returns the
#' ploidy of the data.} \item{ploidy<-}{sets the ploidy of the data using an
#' integer.} \item{alleles}{returns the alleles of each locus.}
#' \item{alleles<-}{sets the alleles of each locus using a list with one
#' character vector for each locus.} \item{other}{returns the content of the
#' \code{@@other} slot (misc. information); returns \code{NULL} if the slot is
#' empty or of length zero.} \item{other<-}{sets the content of the
#' \code{@@other} slot (misc. information); the provided value needs to be a
#' list; it not, provided value will be stored within a list.} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords manip
#' @examples
#'
#' data(nancycats)
#' nancycats
#' pop(nancycats) # get the populations
#' indNames(nancycats) # get the labels of individuals
#' locNames(nancycats) # get the labels of the loci
#' alleles(nancycats) # get the alleles
#'
#' # let's isolate populations 4 and 8
#' temp <- nancycats@@pop=="P04" | nancycats@@pop=="P08"
#' obj <- nancycats[temp,]
#' obj
#'
#' pop(obj)
#'
#' # let's isolate two markers, fca23 and fca90
#' locNames(nancycats)
#'
#' # they correspond to L2 and L7
#' nancycats$loc.fac
#' temp <- nancycats$loc.fac=="L2" | nancycats$loc.fac=="L7"
#' obj <- nancycats[,temp]
#' obj
#'
#' obj$loc.fac
#' locNames(obj)
#'
#' # or more simply
#' nancycats[loc=c("L2","L7")]
#' obj$loc.fac
#' locNames(obj)
#'
#' # using 'drop':
#' truenames(nancycats[1:2])$tab
#' truenames(nancycats[1:2, drop=TRUE])$tab
#'
#' # illustrate how 'other' slot is handled
#' colonies <- genind2genpop(nancycats)
#' colonies@@other$aChar <- "This will not be proceeded"
#' colonies123 <- colonies[1:3]
#' colonies
#' colonies@@other$xy
#'
#' # illustrate pop
#' obj <- nancycats[sample(1:100,10)]
#' obj$pop
#' obj$pop.names
#' pop(obj)
#' pop(obj) <- rep(c('b','a'), each=5)
#' obj$pop
#' obj$pop.names
#' pop(obj)
#'
#' # illustrate locNames
#' locNames(obj)
#' locNames(obj, withAlleles=TRUE)
#'
#'

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

  cat("\n\nOptional contents: ")
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

