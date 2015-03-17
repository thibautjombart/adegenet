

#' Compute the proportion of typed elements
#'
#' The generic function \code{propTyped} is devoted to investigating the
#' structure of missing data in adegenet objects.\cr
#'
#' Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
#' objects. They can return the proportion of available (i.e. non-missing) data
#' per individual/population, locus, or the combination of both in with case
#' the matrix indicates which entity (individual or population) was typed on
#' which locus.
#'
#' When \code{by} is set to "both", the result is a matrix of binary data with
#' entities in rows (individuals or populations) and markers in columns. The
#' values of the matrix are 1 for typed data, and 0 for NA.
#'
#' @name propTyped-methods
#' @aliases propTyped propTyped-methods propTyped,genind-method
#' propTyped,genpop-method
#' @docType methods
#' @param x a \linkS4class{genind} and \linkS4class{genpop} object
#' @param by a character being "ind","loc", or "both" for \linkS4class{genind}
#' object and "pop","loc", or "both" for \linkS4class{genpop} object. It
#' specifies whether proportion of typed data are provided by entity
#' ("ind"/"pop"), by locus ("loc") or both ("both"). See details.
#' @return A vector of proportion (when \code{by} equals "ind", "pop", or
#' "loc"), or a matrix of binary data (when \code{by} equals "both")
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods manip
#' @examples
#'
#' \dontrun{
#' data(nancycats)
#' propTyped(nancycats,by="loc")
#' propTyped(genind2genpop(nancycats),by="both")
#' }
#'

setGeneric("propTyped", function(x,...){
    standardGeneric("propTyped")
})



setMethod("propTyped","genind", function(x, by=c("ind","loc","both")){

    ## checkType(x)
    by <- match.arg(by)

    ## PA case ##
    if(x@type=="PA"){
        temp <- as.matrix(x)

        if(by=="ind"){
            res <- apply(temp,1,function(r) mean(is.na(r)))
            return(1-res)
        }

        if(by=="loc"){
            res <- apply(temp,2,function(c) mean(is.na(c)))
            return(1-res)
        }

        if(by=="both"){
            res <- !is.na(temp)
            return(res*1)
        }
    } # end PA case


    ## codom case ##

    ## auxil function f1
    f1 <- function(vec){
        if(any(is.na(vec))) return(0)
        else return(1)
    }

    ## temp is a list (one component / marker)
    ## with n values (0: not typed, 1: typed)
    kX <- seploc(x,res.type="matrix")
    temp <- lapply(kX, function(X) apply(X, 1, f1))

    ## by individual
    if(by=="ind"){
        temp <- as.data.frame(temp)
        res <- apply(temp,1,mean)
    }

    ## by locus
    if(by=="loc"){
        res <- unlist(lapply(temp,mean))
    }

    ## by individual and locus
    if(by=="both"){
        res <- as.matrix(as.data.frame(temp))
    }

    return(res)
})





setMethod("propTyped","genpop", function(x, by=c("pop","loc","both")){

    ## checkType(x)
    by <- match.arg(by)


    ## PA case ##
    if(x@type=="PA"){
        temp <- as.matrix(x)

        if(by=="pop"){
            res <- apply(temp,1,function(r) mean(is.na(r)))
            return(1-res)
        }

        if(by=="loc"){
            res <- apply(temp,2,function(c) mean(is.na(c)))
            return(1-res)
        }

        if(by=="both"){
            res <- !is.na(temp)
            return(res*1)
        }
    } # end PA case


    ## codom case ##


    ## auxil function f1
    f1 <- function(vec){
        if(any(is.na(vec))) return(0)
        else return(1)
    }

    ## weighted mean
    mean.w <- function(x,w=rep(1/length(x),length(x))){
        x <- x[!is.na(x)]
        w <- w[!is.na(x)]
        w <- w/sum(w)
        return(sum(x*w))
    }

    ## temp is a list (one component / marker)
    ## with n values (0: not typed, 1: typed)
    kX <- seploc(x,res.type="matrix")
    temp <- lapply(kX, function(X) apply(X, 1, f1))

    ## by individual
    if(by=="pop"){
        temp <- as.data.frame(temp)
        w <- unlist(lapply(kX, sum,na.rm=TRUE))
        res <- apply(temp,1,mean.w,w=w)
    }

    ## by locus
    if(by=="loc"){
        w <- apply(x@tab,1,sum,na.rm=TRUE)
        res <- unlist(lapply(temp,mean.w,w=w))
    }

    ## by individual and locus
    if(by=="both"){
        res <- as.matrix(as.data.frame(temp))
    }

    return(res)
})
