####################
# scaleGen methods
####################
setGeneric("scaleGen", function(x,...){standardGeneric("scaleGen")})

setMethod("scaleGen", "genind", function(x, center=TRUE, scale=TRUE,
                                         missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    missing <- match.arg(missing)
    ## checkType(x)

    ## handle "missing" arg
    if(missing %in% c("0","mean")){
        x <- na.replace(x, method=missing, quiet=TRUE)
    }

    ## handle specific cases
    if(scale[1]){
        ## get allele freq
        temp <- apply(x$tab,2,mean,na.rm=TRUE)
        if(x@type=="codom"){
            ## coerce sum of alleles freq to one (in case of missing data)
            temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
            pbar <- unlist(temp)
        } else {
            pbar <- temp
        }

        scale <- sqrt(pbar*(1-pbar))
    }

    X <- x$tab
    ## handle truenames
    if(truenames){
        X <- truenames(x)
        if(is.list(X)) { X <- X$tab }
    }

    ## return result
    res <- scale(X, center=center, scale=scale)

    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})





setMethod("scaleGen", "genpop", function(x, center=TRUE, scale=TRUE,
                                         missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    missing <- match.arg(missing)
    ## checkType(x)

    ## make allele frequencies here
    if(x@type=="codom"){
        X <- makefreq(x,quiet=TRUE,missing=missing,truenames=truenames)$tab
    } else{
        X <- truenames(x) # keep binary data if type is PA
    }

    ## handle specific cases
    if(scale[1]){
        ## get allele freq
        temp <- apply(X,2,mean,na.rm=TRUE)
        if(x@type=="codom"){
            ## coerce sum of alleles freq to one (in case of missing data)
            temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
            pbar <- unlist(temp)
        } else{
            pbar <- temp
        }

        scale <- sqrt(pbar*(1-pbar))
    }

    ## return result

    res <- scale(X, center=center, scale=scale)

    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})
