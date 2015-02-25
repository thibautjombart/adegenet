##########################
## Function gstat.randtest
##########################
gstat.randtest <- function(x,pop=NULL, method=c("global","within","between"),
                           sup.pop=NULL, sub.pop=NULL, nsim=499){
    ##   cat("\nSorry, hierfstat package has been disabled - this function will be restored in a future release.\n")
    ## return(invisible())

    if(!is.genind(x)) stop("x is not a valid genind object")
    if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
    checkType(x)
    if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")

    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

    met <- tolower(method[1])
    if(met=="within" && is.null(sup.pop)) stop("Method 'within' chosen but 'sup.pop' is not provided.")
    if(met=="between" && is.null(sub.pop)) stop("Method 'between' chosen but 'sub.pop' is not provided.")

    ## make data for hierfstat
    X <- genind2hierfstat(x=x,pop=pop)

    ## compute obs gstat
    obs <- g.stats.glob(X)$g.stats

    pop <- X[,1]
    X <- X[,-1]

    ## simulations according one of the 3 different schemes
    ## note: for, lapply and sapply are all equivalent
    ## recursive functions would require options("expression") to be modified...
    sim <- vector(mode="numeric",length=nsim)

    if(met=="global"){

        sim <- sapply(1:nsim, function(i) g.stats.glob(cbind(sample(pop),X))$g.stats)

    } else if(met=="within"){

        if(length(sup.pop) != length(pop)) stop("pop and sup.pop do not have the same length.")
        sim <- sapply(1:nsim, function(i) g.stats.glob(cbind(pop,X[samp.within(sup.pop),]))$g.stats)

    } else if(met=="between"){

        if(length(sub.pop) != length(pop)) stop("pop and sub.pop do not have the same length.")
        sim <- sapply(1:nsim, function(i) g.stats.glob(cbind(pop,X[samp.between(sub.pop),]))$g.stats)

    } else {
        stop("Unknown method requested.")
    }

    prevcall <- match.call()

    res <- as.randtest(sim=sim, obs=obs, call=prevcall)

    return(res)

}
