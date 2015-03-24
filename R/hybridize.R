##
## Function hybridize takes two genind in inputs
## and generates hybrids individuals having one parent
## in both objects.
##

hybridize <- function(x1, x2, n, pop=NULL,
                      res.type=c("genind","df","STRUCTURE"),
                      file=NULL, quiet=FALSE, sep="/", hyb.label="h"){
    ## checks
    if(!is.genind(x1)) stop("x1 is not a valid genind object")
    if(!is.genind(x2)) stop("x2 is not a valid genind object")
    if(!all(ploidy(x1)==ploidy(x1)[1])) stop("varying ploidy (in x1) is not supported for this function")
    if(!all(ploidy(x2)==ploidy(x2)[1])) stop("varying ploidy (in x2) is not supported for this function")
    if(ploidy(x1)[1] %% 2 != 0) stop("not implemented for odd levels of ploidy")
    if(ploidy(x1)[1] != ploidy(x2)[1]) stop("x1 and x2 have different ploidy")
    checkType(x1)
    checkType(x2)

    n <- as.integer(n)
    ploidy <- ploidy(x1)[1]
    res.type <- match.arg(res.type)
    if(!all(locNames(x1)==locNames(x2))) stop("names of markers in x1 and x2 do not correspond")

    ## used variables
    n1 <- nInd(x1)
    n2 <- nInd(x2)
    k <- nLoc(x1)

    #### get frequencies for each locus
    y1 <- genind2genpop(x1,pop=factor(rep(1,n1)),quiet=TRUE)
    freq1 <- tab(y1, freq=TRUE) # get frequencies
    freq1 <- split(freq1, y1@loc.fac) # split by locus
    freq1 <- freq1[locNames(x1)] # ensure right order

    y2 <- genind2genpop(x2,pop=factor(rep(1,n2)),quiet=TRUE)
    freq2 <- tab(y2, freq=TRUE) # get frequencies
    freq2 <- split(freq2, y2@loc.fac) # split by locus
    freq2 <- freq1[locNames(x2)] # ensure right order

    #### sampling of gametes
    ## kX1 / kX2 are lists of tables of sampled gametes
    kX1 <- lapply(freq1, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX1) <- locNames(x1)
    for(i in 1:k) { colnames(kX1[[i]]) <- alleles(x1)[[i]]}
    kX2 <- lapply(freq2, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX2) <- locNames(x2)
    for(i in 1:k) { colnames(kX2[[i]]) <- alleles(x2)[[i]]}

    ## construction of zygotes ##
    ## individual gamete tables
    tab1 <- as.matrix(cbind.data.frame(kX1))
    tab2 <- as.matrix(cbind.data.frame(kX2))

    ## make empty matrix with all alleles in tab1 and tab2
    zyg.rownames <- .genlab(hyb.label,n)
    zyg.colnames <- sort(unique(c(colnames(tab1),colnames(tab2))))
    zyg <- matrix(0, nrow=n, ncol=length(zyg.colnames),
                  dimnames=list(zyg.rownames, zyg.colnames))

    ## add in the alleles
    zyg[, colnames(tab1)] <- zyg[, colnames(tab1)] + tab1
    zyg[, colnames(tab2)] <- zyg[, colnames(tab2)] + tab2
    zyg <- zyg
    zyg <- genind(zyg, type="codom", ploidy=ploidy)

    ## res.type=="STRUCTURE"
    if(res.type=="STRUCTURE"){
        ## res <- paste(gam1,gam2,sep=" ") # make df for the hybrids
        ## res <- as.data.frame(matrix(res,ncol=k))
        temp <- genind2df(repool(x1,x2,zyg), usepop=FALSE, sep=" ")
        res <- unlist(apply(temp,1,strsplit," "))
        res <- as.data.frame(matrix(res, nrow=nrow(temp), byrow=TRUE))
        colnames(res) <- rep(colnames(temp),each=ploidy)
        ## df1 <- genind2df(x1,sep=" ",usepop=FALSE) # make df with parents and hybrids
        ## df2 <- genind2df(x2,sep=" ",usepop=FALSE)
        ## res <- rbind.data.frame(df1,df2,res) # rbind the three df
        res[is.na(res)] <- "-9" # this is two missing alleles for STRUCTURE
        pop <- rep(1:3,c(nrow(x1@tab), nrow(x2@tab), n)) # make a pop identifier
        res <- cbind.data.frame(pop,res, stringsAsFactors = FALSE)
        names(res)[1] <- ""

        if(is.null(file)) {
            file <- gsub("[[:space:]]|:","-",date())
            file <- paste("hybrid",file,sep="_")
            file <- paste(file,"str",sep=".")
        }
        write.table(res, file=file,row.names = TRUE, col.names = TRUE, quote=FALSE)
        if(!quiet) cat("\nWrote results to file", file, "\n")

        return(invisible())
    }


    ## res.type=="df"
    if(res.type=="df"){
        ## res <- paste(gam1,gam2,sep=sep)
        ## res <- as.data.frame(matrix(res,ncol=k), stringsAsFactors=FALSE)
        ## names(res) <- x1@loc.names
        ## row.names(res) <- .genlab(hyb.label,n)

        res <- genind2df(zyg, sep=sep)
        return(res)
    }


    ## res.type=="genind"
    if(res.type=="genind"){
        ## res <- paste(gam1,gam2,sep="")
        ##         res <- as.data.frame(matrix(res,ncol=k), stringsAsFactors=FALSE)
        ##         names(res) <- x1@loc.names
        ##         row.names(res) <- .genlab(hyb.label,n)
        if(is.null(pop)){ # if pop is not provided, merge the two parent populations
            pop <- paste(deparse(substitute(x1)) , deparse(substitute(x2)), sep="-")
        }
        pop <- factor(rep(pop,n))

        res <- zyg
        pop(res) <- pop
        res@call <- match.call()

        return(res)
    }

} # end hybridize
