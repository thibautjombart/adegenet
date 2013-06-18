#################
# fstat function
#################
#
# Wrapper for fst estimator from hierfstat package
#
fstat <- function(x, pop=NULL, fstonly=FALSE){
    ## cat("\nSorry, hierfstat package has been disabled - this function will be restored in a future release.\n")
    ## return(invisible())
    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")
    if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
    checkType(x)

    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

    ## computations
    dat <- genind2hierfstat(x)[,-1]
    res <- varcomp.glob(levels=data.frame(pop), loci=dat)$F

    if(fstonly) {res <- res[1,1]}
    return(res)
}



###############
## pairwise.fst
###############
##
## pairwise fst sensu Nei (Ht - mean(Hs))/Ht
##
pairwise.fst <- function(x, pop=NULL, res.type=c("dist","matrix"), truenames=TRUE){
    ## MISC CHECKS ##
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!is.null(pop)){
        pop(x) <- pop
    }
    temp <- pop(x)
    if(is.null(temp)) stop("no grouping factor (pop) provided")
    if(length(levels(temp)) < 2){
        warning("There is only one pop - returning NULL")
        return(NULL)
    }

    res.type <- match.arg(res.type)


    ## COMPUTATIONS ##

    ## function to compute one Fst ##
    f1 <- function(pop1, pop2){ # pop1 and pop2 are genind obj. with a single pop each
        n1 <- nrow(pop1@tab)
        n2 <- nrow(pop2@tab)
        temp <- repool(pop1,pop2)
        b <- weighted.mean(Hs(temp), c(n1,n2)) # mean Hs is weighted for pop sizes
        pop(temp) <- NULL
        a <- Hs(temp)
        return((a-b)/a)
    }


    ## compute pairwise Fst for all pairs
    lx <- seppop(x,treatOther=FALSE)
    temp <- pop(x)
    levPop <- levels(temp)
    allPairs <- combn(1:length(levPop), 2)
    if(!is.matrix(allPairs)){
        allPairs <- matrix(allPairs,nrow=2)
    }

    vecRes <- numeric()
    for(i in 1:ncol(allPairs)){
        vecRes[i] <- f1(lx[[allPairs[1,i]]], lx[[allPairs[2,i]]])
    }


    squelres <- dist(1:length(levPop))
    res <- vecRes
    attributes(res) <- attributes(squelres)

    if(res.type=="matrix"){
        res <- as.matrix(res)
        if(truenames){
            lab <- x@pop.names
        } else {
            lab <- names(x@pop.names)
        }

        colnames(res) <- rownames(res) <- lab
    }

    return(res)
} # end of pairwise.fst
