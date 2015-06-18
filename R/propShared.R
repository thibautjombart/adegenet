## propShared computes the proportion of shared alleles
## in a genind object


######################
# Function propShared
######################
propShared <- function(obj){
    ## CHECK THAT THIS IS A VALID GENIND ##
    if(!inherits(obj,"genind")) stop("obj must be a genind object.")
    invisible(validObject(obj))


    ## GET MATRIX OF NB OF ALLELES ##
    x <- tab(obj)
    x[is.na(x)] <- 0L

    ## COMPUTE NB OF SHARED ALLELES ##
    n <- nInd(obj)
    resVec <- integer(n*(n-1)/2)

    res <- .C("nb_shared_all", as.integer(x), as.integer(resVec), as.integer(n), as.integer(ncol(obj$tab)), PACKAGE="adegenet")[[2]]
    attr(res,"Size") <- n
    attr(res,"Diag") <- FALSE
    attr(res,"Upper") <- FALSE
    class(res) <- "dist"
    res <- as.matrix(res)


    ## COMPUTE NB OF ALLELES TYPED IN COMMON ##
    tabNA <- propTyped(obj, by="both")
    tabTypCom <- tabNA %*% t(tabNA) * ploidy(obj)


    ## GET PROPORTIONS OF SHARED ALLELES ##
    res <- res/tabTypCom
    diag(res) <- 1L
    colnames(res) <-rownames(res) <- indNames(obj)
    return(res)
}

