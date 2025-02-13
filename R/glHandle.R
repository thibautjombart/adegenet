
# Function to subset raw vectors
.subsetbin <- function(x, i){
    # Take a raw vector, subset the bits and then convert to integers.
    xint   <- as.integer(rawToBits(x)[i])
    # Figure out how many zeroes are needed to pad the end.
    extra_bits <- length(xint) %% 8
    # https://github.com/thibautjombart/adegenet/issues/363
    # Fix to a bug caught by @maxcoulter. In the case 
    zeroes <- if (extra_bits > 0) 8 - extra_bits else 0
    # Convert the integer vector with zeroes on the end back into a raw vector.
    return(packBits(c(xint, rep(0L, zeroes))))
}

# old method for [] for SNPbin
.oldSNPbinset <- function(x, i){
    if (missing(i)) i <- TRUE
    temp <- .SNPbin2int(x) # data as integers with NAs
    x <- new("SNPbin", snp=temp[i], label=x@label, ploidy=x@ploidy)
    return(x)
}

# Zhian N. Kamvar
# Mon Aug 17 09:39:12 2015 ------------------------------
# 
# This function takes two steps:
#   1. Subset the missing positions
#   2. Subset the vectors of raw SNPs
# 
# Both steps are not exactly straighforward. Because the missing vector only
# represents the positions of missing data, it must be subset by value as
# opposed to position.
.SNPbinset <- function(x, i){
  if (missing(i)) i <- TRUE
  
  # Create a logical value indicating whether or not subsetting is necessary.
  we_take_all <- length(i) == 1 && is.logical(i) && i
  n.loc <- x@n.loc
  if (length(x@NA.posi) > 0){
    if (is.logical(i)){
      if (we_take_all){
        # Keep all of the data
        return(x)
      } else {
        # If the positons are logical, perhaps the best way to address this is
        # to match the TRUE positions to the NA.posi vector. Adding nomatch = 0 
        # avoids introducing NAs.
        namatches <- match(which(i), x@NA.posi, nomatch = 0)
        nas.kept  <- x@NA.posi[namatches]
      }
    } else if (is.character(i)){
      stop("Cannot subset a SNPbin object with a character vector", call. = FALSE)
    } else if (all(i < 0)){
      # For negative subscripts, find which ones they match and then
      # negate those. Luckily -0 is allowed.
      namatches <- match(abs(i), x@NA.posi, nomatch = 0)
      # Unfortunately, if nothing matches, then the default are zeroes. When you
      # subset a vector in R with only zero, you will get an empty vector. This
      # conditional makes sure that NA positions are retained.
      if (all(namatches == 0)){
        nas.kept  <- x@NA.posi
      } else {
        nas.kept  <- x@NA.posi[-namatches]
      }
      
    } else if (all(i > 0)){
      # Positive subscripts are much easier. First you find where the subscripts
      # match and then your subset with those positions. 
      namatches <- match(i, x@NA.posi, nomatch = 0)
      nas.kept  <- x@NA.posi[namatches]
    } else {
      stop("Cannot subset a SNPbin with mixed subscripts.", call. = FALSE)
    }
    # After we find out which missing positions we need to keep, we reset the 
    # missing positions to the subset data.
    if (length(nas.kept) > 0){
      old.posi  <- 1:n.loc
      x@NA.posi <- match(nas.kept, old.posi[i])
    } else {
      x@NA.posi <- nas.kept
    }
  }
  # Here we calculate the number of loci we will have left in the data.
  if (we_take_all){
    return(x)
  } else if (all(is.logical(i))){
    n.loc <- sum(i)
  } else if (any(i < 0)){
    n.loc <- n.loc - length(i)
  } else {
    n.loc <- length(i)
  }
  # Now we loop over all chromosomes and subset.
  x@snp   <- lapply(x@snp, .subsetbin, i)
  # Set the new value of the number of loci and return.
  x@n.loc <- n.loc
  return(x)
}

###############
## '[' operators
###############
## SNPbin

setMethod("[", signature(x="SNPbin", i="ANY"), function(x, i) {
    .SNPbinset(x, i)
}) # end [] for SNPbin




## genlight
setMethod("[", signature(x = "genlight", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., pop=NULL, treatOther=TRUE, quiet=TRUE, drop=FALSE) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    ori.n <- nInd(x)
    ori.p <- nLoc(x)

    ## recycle logicals if needed
    if(!is.null(i) && is.logical(i)) i <- rep(i, length=ori.n)
    if(!is.null(j) && is.logical(j)) j <- rep(j, length=ori.p)

    if (!is.null(pop) && !is.null(pop(x))){
      i <- .get_pop_inds(x, pop)
    }


    ## SUBSET INDIVIDUALS ##
    ## genotypes
    x@gen <- x@gen[i]

    ## ind names
    x@ind.names <- x@ind.names[i]

    ## ploidy
    if(!is.null(x@ploidy)) {
        ori.ploidy <- ploidy(x) <- ploidy(x)[i]
    } else {
        ori.ploidy <- NULL
    }

    ## pop
    if(!is.null(pop(x))) {
        ori.pop <- pop(x) <- factor(pop(x)[i])
    } else {
        ori.pop <- NULL
    }
    ## strata
    if(!is.null(x@strata)) {
        ori.strata <- x@strata <- x@strata[i, , drop = FALSE]
    } else {
        ori.strata <- NULL
    }


    ## HANDLE 'OTHER' SLOT ##
    nOther <- length(other(x))
    namesOther <- names(other(x))
    counter <- 0
    if(treatOther & !(is.logical(i) && all(i))){
        f1 <- function(obj,n=ori.n){
            counter <<- counter+1
            if(!is.null(dim(obj)) && nrow(obj)==ori.n) { # if the element is a matrix-like obj
                obj <- obj[i,,drop=FALSE]
            } else if(length(obj) == ori.n) { # if the element is not a matrix but has a length == n
                obj <- obj[i]
                if(is.factor(obj)) {obj <- factor(obj)}
            } else {if(!quiet) warning(paste("cannot treat the object",namesOther[counter]))}

            return(obj)
        } # end f1

        other(x) <- lapply(x@other, f1) # treat all elements

    } # end treatOther


    ## SUBSET LOCI ##

    ## handle ind.names, loc.names, chromosome, position, and alleles
    if (is.character(j)){
      j <- match(j, x@loc.names, nomatch = 0)
    }
    x@loc.names   <- x@loc.names[j]
    x@chromosome  <- chr(x)[j]
    x@position    <- position(x)[j]
    x@loc.all     <- alleles(x)[j]
    x@gen         <- lapply(x@gen, function(e) e[j])
    x@n.loc       <- x@gen[[1]]@n.loc

    return(x)
}) # end [] for genlight







######################
##
## c, cbind, rbind...
##
######################

################
## cbind SNPbin
################
##setMethod("cbind", signature("SNPbin"), function(..., deparse.level = 1) {
#' @export
cbind.SNPbin <- function(..., checkPloidy=TRUE){
    myList <- list(...)
    if(!all(sapply(myList, class)=="SNPbin")) stop("some objects are not SNPbin objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }


    if(checkPloidy && length(unique(sapply(myList, ploidy))) !=1 ) stop("objects have different ploidy levels")
    if(checkPloidy) {
        ori.ploidy <- ploidy(myList[[1]])
    } else {
        ori.ploidy <- NULL
    }
    x <- new("SNPbin", unlist(lapply(myList, as.integer)), ploidy=ori.ploidy)
    return(x)
} # end cbind.SNPbin
##})



#' @export
c.SNPbin <- function(...){
    return(cbind(...))
}




##################
## cbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
#' @export
cbind.genlight <- function(...){
      ## store arguments
    dots <- list(...)

    ## extract arguments which are genlight objects
    myList <- dots[sapply(dots, inherits, "genlight")]

    ## keep the rest in 'dots'
    dots <- dots[!sapply(dots, inherits, "genlight")]

    if(length(myList)==1 && is.list(myList[[1]])) myList <- myList[[1]]
    if(!all(sapply(myList, class)=="genlight")) stop("some objects are not genlight objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }

    ## different checks
    if(length(unique(sapply(myList, nInd))) > 1 ) stop("objects have different numbers of individuals")
    n.obj <- length(myList)
    n.ind <- nInd(myList[[1]])
    if(n.ind==0){
        warning("All objects are empty")
        return(NULL)
    }
    temp <- as.matrix(as.data.frame(lapply(myList, ploidy)))
    if(any(apply(temp,1,function(r) length(unique(r)))>1)) stop("non-consistent ploidy across datasets")
    ori.ploidy <- ploidy(myList[[1]])


    ## merge one individual at a time ##
    res <- list()
    for(i in 1:n.ind){
        res[[i]] <- Reduce(function(a,b) {cbind(a,b,checkPloidy=FALSE)}, lapply(myList, function(e) e@gen[[i]]) )
    }

    dots$gen <- res
    dots$Class <- "genlight"
    res <- do.call(new, dots)

    ## handle loc.names, alleles, etc. ##
    indNames(res) <- indNames(myList[[1]])
    locNames(res) <- unlist(lapply(myList, locNames))
    alleles(res) <- unlist(lapply(myList, alleles))
    pop(res) <- pop(myList[[1]])
    res@strata <- myList[[1]]@strata
    ploidy(res) <- ori.ploidy

    ## return object ##
    return(res)
} # end cbind.genlight
##})






##################
## rbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
#' @importFrom dplyr bind_rows
#' @export
rbind.genlight <- function(...){
    ## store arguments
    dots <- list(...)

    ## extract arguments which are genlight objects
    myList <- dots[sapply(dots, inherits, "genlight")]

    ## keep the rest in 'dots'
    dots <- dots[!sapply(dots, inherits, "genlight")]

    if(!all(sapply(myList, class)=="genlight")) stop("some objects are not genlight objects")

    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }

    if(length(unique(sapply(myList, nLoc))) !=1 ) stop("objects have different numbers of SNPs")

    ## build output
    dots$Class <- "genlight"
    dots$gen <- Reduce(c, lapply(myList, function(e) e@gen))
    res <- do.call(new, dots)
    locNames(res) <- locNames(myList[[1]])
    alleles(res)  <- alleles(myList[[1]])
    indNames(res) <- unlist(lapply(myList, indNames))
    pop(res)      <- factor(unlist(lapply(myList, pop)))

    # Hierarchies are tricky. Using dplyr's bind_rows.
    res <- .rbind_strata(myList, res)

    ## return object ##
    return(res)

} # end rbind.genlight






##########
## seppop
##########
setMethod("seppop", signature(x="genlight"), function(x, pop=NULL, treatOther=TRUE, keepNA = FALSE, quiet=TRUE, ...){
  .seppop_internal(
    x = x,
    pop = pop,
    treatOther = treatOther,
    keepNA = keepNA,
    quiet = quiet,
    ...
  )
}) # end seppop






##########
## seploc
##########
setMethod("seploc", signature(x="genlight"), function(x, n.block=NULL, block.size=NULL, random=FALSE,
                               parallel=FALSE, n.cores=NULL){
    ## CHECKS ##
    if(is.null(n.block) & is.null(block.size)) stop("n.block and block.size are both missing.")
    if(!is.null(n.block) & !is.null(block.size)) stop("n.block and block.size are both provided.")
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }


    ## GET BLOCK SIZE VECTOR ##
    P <- nLoc(x)

    ## n.block is given
    if(!is.null(n.block)){
        vec.blocksize <- rep(P %/% n.block, n.block)
        if(P %% n.block >0){
            vec.blocksize[1:(P %% n.block)] <- vec.blocksize[1:(P %% n.block)] + 1
        }

    }

    ## block.size is given
    if(!is.null(block.size)){
        vec.blocksize <- rep(block.size, P %/% block.size)
        if(P %% block.size >0){
             vec.blocksize <- c( vec.blocksize, P %% block.size)
        }
    }


    ## split data by blocks ##
    fac.block <- factor(rep(1:length(vec.blocksize), vec.blocksize))
    if(random){
        fac.block <- sample(fac.block)
    }

    if(parallel){
        if(random){
            res <- parallel::mclapply(levels(fac.block), function(lev) x[, sample(which(fac.block==lev))],
                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
        } else {
            res <- parallel::mclapply(levels(fac.block), function(lev) x[, which(fac.block==lev)],
                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
        }
    } else {
         if(random){
             res <- lapply(levels(fac.block), function(lev) x[, sample(which(fac.block==lev))])
         } else {
             res <- lapply(levels(fac.block), function(lev) x[, which(fac.block==lev)])
         }
    }

    ## return result ##
    names(res) <- paste("block", 1:length(res),sep=".")

    return(res)
}) # end seploc






###################
### TESTING
###################


## c, cbind, rbind ##
## a <- new("genlight", list(c(1,0,1), c(0,0,1,0)) )
## b <- new("genlight", list(c(1,0,1,1,1,1), c(1,0)) )
## locNames(a) <- letters[1:4]
## locNames(b) <- 1:6
## c <- cbind(a,b)
## identical(as.matrix(c),cbind(as.matrix(a), as.matrix(b))) # MUST BE TRUE
## identical(as.matrix(rbind(a,a)),rbind(as.matrix(a),as.matrix(a)))




## test subsetting with/without @other ##
## x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters, data.frame(2:4)))
## pop(x) <- c("pop1","pop1", "pop2")
