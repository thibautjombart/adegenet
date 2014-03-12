
###############
## '[' operators
###############
## SNPbin
setMethod("[", signature(x="SNPbin", i="ANY"), function(x, i) {
    if (missing(i)) i <- TRUE
    temp <- .SNPbin2int(x) # data as integers with NAs
    x <- new("SNPbin", snp=temp[i], label=x@label, ploidy=x@ploidy)
    return(x)
}) # end [] for SNPbin




## genlight
setMethod("[", signature(x="genlight", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ..., treatOther=TRUE, quiet=TRUE, drop=FALSE) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    ori.n <- nInd(x)


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
    if(length(j)==1 && is.logical(j) && j){ # no need to subset SNPs
        return(x)
    } else { # need to subset SNPs
        old.other <- other(x)
        old.ind.names <- indNames(x)

        ## handle ind.names, loc.names, chromosome, position, and alleles
        new.loc.names <- locNames(x)[j]
        new.chr <- chr(x)[j]
        new.position <- position(x)[j]
        new.alleles <- alleles(x)[j]
        new.gen <- lapply(x@gen, function(e) e[j])
        ##x <- as.matrix(x)[, j, drop=FALSE] # maybe need to process one row at a time
        x <- new("genlight", gen=new.gen, pop=ori.pop, ploidy=ori.ploidy,
                 ind.names=old.ind.names, loc.names=new.loc.names,
                 chromosome=new.chr, position=new.position, alleles=new.alleles, other=old.other, parallel=FALSE,...)
    }

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



c.SNPbin <- function(...){
    return(cbind(...))
}




##################
## cbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
cbind.genlight <- function(...){
    myList <- list(...)
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

    res <- new("genlight",res,...)

    ## handle loc.names, alleles, etc. ##
    indNames(res) <- indNames(myList[[1]])
    locNames(res) <- unlist(lapply(myList, locNames))
    alleles(res) <- unlist(lapply(myList, alleles))
    pop(res) <- pop(myList[[1]])
    ploidy(res) <- ori.ploidy

    ## return object ##
    return(res)
} # end cbind.genlight
##})






##################
## rbind genlight
##################
##setMethod("cbind", signature(x="genlight"), function(..., deparse.level = 1) {
rbind.genlight <- function(...){
    myList <- list(...)
    if(!all(sapply(myList, class)=="genlight")) stop("some objects are not genlight objects")
    ## remove empty objects
    myList <- myList[sapply(myList,nLoc)>0 & sapply(myList,nInd)>0]
    if(length(myList)==0) {
        warning("All objects are empty")
        return(NULL)
    }

    if(length(unique(sapply(myList, nLoc))) !=1 ) stop("objects have different numbers of SNPs")


    ## build output
    res <- new("genlight", Reduce(c, lapply(myList, function(e) e@gen)), ...)
    locNames(res) <- locNames(myList[[1]])
    alleles(res) <- alleles(myList[[1]])
    indNames(res) <- unlist(lapply(myList, indNames))
    pop(res) <- factor(unlist(lapply(myList, pop)))

    ## return object ##
    return(res)

} # end rbind.genlight






##########
## seppop
##########
setMethod("seppop", signature(x="genlight"), function(x, pop=NULL, treatOther=TRUE, quiet=TRUE, ...){
    ## HANDLE POP ARGUMENT ##
    if(!is.null(pop)) {
        pop(x) <- pop
    }

    if(is.null(pop(x))) stop("pop not provided and pop(x) is NULL")

    ## PERFORM SUBSETTING ##
    kObj <- lapply(levels(pop(x)), function(lev) x[pop(x)==lev, , treatOther=treatOther, quiet=quiet, ...])
    names(kObj) <- levels(pop(x))

    return(kObj)
}) # end seppop






##########
## seploc
##########
setMethod("seploc", signature(x="genlight"), function(x, n.block=NULL, block.size=NULL, random=FALSE,
                               parallel=require(parallel), n.cores=NULL){
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
            res <- mclapply(levels(fac.block), function(lev) x[,sample(which(fac.block==lev))],
                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
        } else {
            res <- mclapply(levels(fac.block), function(lev) x[,which(fac.block==lev)],
                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
        }
    } else {
         if(random){
             res <- lapply(levels(fac.block), function(lev) x[,sample(which(fac.block==lev))])
         } else {
             res <- lapply(levels(fac.block), function(lev) x[,which(fac.block==lev)])
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
