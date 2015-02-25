##
## Function hybridize takes two genind in inputs
## and generates hybrids individuals having one parent
## in both objects.
##



#' Simulated hybridization between two samples of populations
#' 
#' The function \code{hybridize} performs hybridization between two set of
#' genotypes stored in \linkS4class{genind} objects (referred as the "2
#' populations"). Allelic frequencies are derived for each population, and then
#' gametes are sampled following a multinomial distribution. \cr
#' 
#' The result consists in a set of 'n' genotypes, with different possible
#' outputs (see 'res.type' argument).
#' 
#' If the output is a STRUCTURE file, this file will have the following
#' caracteristics:\cr - file contains the genotypes of the parents, and then
#' the genotypes of hybrids\cr - the first column identifies genotypes\cr - the
#' second column identifies the population (1 and 2 for parents x1 and x2; 3
#' for hybrids)\cr - the first line contains the names of the markers\cr - one
#' row = one genotype (onerowperind will be true)\cr - missing values coded by
#' "-9" (the software's default)\cr
#' 
#' @param x1 a \linkS4class{genind} object
#' @param x2 a \linkS4class{genind} object
#' @param n an integer giving the number of hybrids requested
#' @param pop a character string giving naming the population of the created
#' hybrids. If NULL, will have the form "x1-x2"
#' @param res.type a character giving the type of output requested. Must be
#' "genind" (default), "df" (i.e. data.frame like in \code{\link{genind2df}}),
#' or "STRUCTURE" to generate a .str file readable by STRUCTURE (in which case
#' the 'file' must be supplied). See 'details' for STRUCTURE output.
#' @param file a character giving the name of the file to be written when
#' 'res.type' is "STRUCTURE"; if NULL, a the created file is of the form
#' "hybrids\_[the current date].str".
#' @param quiet a logical specifying whether the writing to a file (when
#' 'res.type' is "STRUCTURE") should be announced (FALSE, default) or not
#' (TRUE).
#' @param sep a character used to separate two alleles
#' @param hyb.label a character string used to construct the hybrids labels; by
#' default, "h", which gives labels: "h01", "h02", "h03",...
#' @return A \linkS4class{genind} object (by default), or a data.frame of
#' alleles (res.type="df"). No R output if res.type="STRUCTURE" (results
#' written to the specified file).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' ## Let's make some cattle hybrids
#' ##
#' data(microbov)
#' 
#' ## first, isolate each breed
#' temp <- seppop(microbov)
#' names(temp)
#' 
#' salers <- temp$Salers
#' zebu <- temp$Zebu
#' borgou <- temp$Borgou
#' somba <- temp$Somba
#' 
#' ## let's make some... Zeblers
#' zebler <- hybridize(salers, zebu, n=40)
#' 
#' ## and some Somgou
#' somgou <- hybridize(somba, borgou, n=40)
#' 
#' ## now let's merge all data into a single genind
#' newDat <- repool(microbov, zebler, somgou)
#' 
#' ## make a correspondance analysis
#' ## and see where hybrids are placed
#' X <- genind2genpop(newDat,missing="chi2",quiet=TRUE)
#' coa1 <- dudi.coa(as.data.frame(X$tab),scannf=FALSE,nf=3)
#' s.label(coa1$li,label=X$pop.names)
#' add.scatter.eig(coa1$eig,2,1,2)
#' 
#' }
#' 
#' @export hybridize
hybridize <- function(x1, x2, n, pop=NULL,
                      res.type=c("genind","df","STRUCTURE"), file=NULL, quiet=FALSE, sep="/", hyb.label="h"){
    ## checks
    if(!is.genind(x1)) stop("x1 is not a valid genind object")
    if(!is.genind(x2)) stop("x2 is not a valid genind object")
    if(x1@ploidy %% 2 != 0) stop("not implemented for odd levels of ploidy")
    if(x2@ploidy != x1@ploidy) stop("not implemented for genotypes with different ploidy levels")
    checkType(x1)
    checkType(x2)

    n <- as.integer(n)
    ploidy <- ploidy(x1)
    res.type <- match.arg(res.type)
    if(!all(locNames(x1)==locNames(x2))) stop("names of markers in x1 and x2 do not correspond")

    ## used variables
    n1 <- nInd(x1)
    n2 <- nInd(x2)
    k <- nLoc(x1)

    #### get frequencies for each locus
    y1 <- genind2genpop(x1,pop=factor(rep(1,n1)),missing="0",quiet=TRUE)
    freq1 <- makefreq(y1,quiet=TRUE)$tab
    freq1 <- split(freq1, y1@loc.fac)

    y2 <- genind2genpop(x2,pop=factor(rep(1,n2)),missing="0",quiet=TRUE)
    freq2 <- makefreq(y2,quiet=TRUE)$tab
    freq2 <- split(freq2, y2@loc.fac)

    #### sampling of gametes
    ## kX1 / kX2 are lists of tables of sampled gametes
    kX1 <- lapply(freq1, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX1) <- x1$loc.names
    for(i in 1:k) { colnames(kX1[[i]]) <- alleles(x1)[[i]]}
    kX2 <- lapply(freq2, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX2) <- x2$loc.names
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
    zyg <- zyg/ploidy
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
