############
# propTyped
############
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
