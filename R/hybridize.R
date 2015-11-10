#'
#' Function hybridize takes two genind in inputs
#' and generates hybrids individuals having one parent
#' in both objects.
#'
#'The function \code{hybridize} performs hybridization between two set
#' of genotypes stored in \linkS4class{genind} objects (referred as the "2
#' populations"). Allelic frequencies are derived for each population,
#' and then gametes are sampled following a multinomial distribution. \cr
#'
#' The result consists in a set of 'n' genotypes, with different possible
#' outputs (see 'res.type' argument).
#'
#'
#' @export
#'
#' @param x1 a \linkS4class{genind} object
#' @param x2 a \linkS4class{genind} object
#' @param n an integer giving the number of hybrids requested
#' @param pop a character string giving naming the population of the
#'  created hybrids.
#' @param res.type a character giving the type of output requested. Must
#'  be "genind" (default), "df" (i.e. data.frame like in
#'  \code{\link{genind2df}}), or "STRUCTURE" to generate a .str file
#'  readable by STRUCTURE (in which case the 'file' must be supplied). See
#'  'details' for STRUCTURE output.
#' @param file a character giving the name of the file to be written
#'  when 'res.type' is "STRUCTURE"; if NULL, a the created file is of the
#'  form "hybrids\_[the current date].str".
#' @param quiet a logical specifying whether the writing to a file (when
#'    'res.type' is "STRUCTURE") should be announced (FALSE, default) or
#'    not (TRUE).
#' @param sep a character used to separate two alleles
#' @param hyb.label a character string used to construct the hybrids
#'  labels; by default, "h", which gives labels: "h01", "h02", "h03",...
#'
#' @return
#' A \linkS4class{genind} object (by default), or a data.frame of alleles
#' (res.type="df"). No R output if res.type="STRUCTURE" (results written
#' to the specified file).
#'
#' @details
#'   If the output is a STRUCTURE file, this file will have the following
#'  caracteristics:\cr
#'  - file contains the genotypes of the parents, and then the genotypes
#'  of hybrids\cr
#'  - the first column identifies genotypes\cr
#'  - the second column identifies the population (1 and 2 for parents x1 and x2;
#'  3 for hybrids)\cr
#'  - the first line contains the names of the markers\cr
#'  - one row = one genotype (onerowperind will be true)\cr
#'  - missing values coded by "-9" (the software's default)\cr
#'
#' @examples
#' \dontrun{
#' ## Let's make some cattle hybrids
#' data(microbov)
#'
#' ## first, isolate each breed
#' temp <- seppop(microbov)
#' names(temp)
#' salers <- temp$Salers
#' zebu <- temp$Zebu
#'
#' ## let's make some... Zeblers
#' zebler <- hybridize(salers, zebu, n=40,
#'                     pop="Zebler")
#'
#'
#' ## now let's merge all data into a single genind
#' newDat <- repool(microbov, zebler)
#'
#' ## make a correspondance analysis
#' ## and see where hybrids are placed
#' X <- genind2genpop(newDat, quiet=TRUE)
#' coa1 <- dudi.coa(tab(X),scannf=FALSE,nf=3)
#' s.label(coa1$li)
#' add.scatter.eig(coa1$eig,2,1,2)
#'
#' }
#'
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @seealso \code{\link{seploc}}, \code{\link{seppop}}, \code{\link{repool}}
#'
#'
#'
hybridize <- function(x1, x2, n, pop="hybrid",
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

    ## store a few variables
    n <- as.integer(n)
    ploidy <- ploidy(x1)[1]
    res.type <- match.arg(res.type)

    ## ensure different names for pop
    popNames(x1) <- "pop1"
    popNames(x2) <- "pop2"

    ## repool data
    x1x2 <- repool(x1, x2)
    x1 <- x1x2[pop=1]
    x2 <- x1x2[pop=2]

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
    freq2 <- freq2[locNames(x2)] # ensure right order
    
    #### sampling of gametes
    ## kX1 / kX2 are lists of tables of sampled gametes
    kX1 <- lapply(freq1, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX1) <- locNames(x1)
    vec.paste1<-NULL
    Vec.all1<-NULL
    for(i in 1:k) { 
      colnames(kX1[[i]]) <- alleles(x1)[[i]]
      ## Paste the alleles locus after locus
      vec.paste1<-c(vec.paste1, alleles(x1)[[i]])
      ## Paste the number of alleles, locus after locus
      Vec.all1<-c(Vec.all1, length(alleles(x1)[[i]]))
    }
    kX2 <- lapply(freq2, function(v) t(rmultinom(n,ploidy/2,v)))
    names(kX2) <- locNames(x2)
    vec.paste2<-NULL
    Vec.all2<-NULL
    for(i in 1:k) { 
      colnames(kX2[[i]]) <- alleles(x2)[[i]]
      vec.paste2<-c(vec.paste2, alleles(x2)[[i]])
      Vec.all2<-c(Vec.all2, length(alleles(x2)[[i]]))
    }
    
    
    ## construction of zygotes ##
    ## individual gamete tables
    tab1 <- as.matrix(cbind.data.frame(kX1))
    ## Force the names of the columns for tab1 and tab2 with the pattern "locNames.allele"
    colnames(tab1)<-paste(rep(locNames(x1), Vec.all1), ".",vec.paste1, sep = "")
    tab2 <- as.matrix(cbind.data.frame(kX2))
    colnames(tab2)<-paste(rep(locNames(x2), Vec.all2), ".",vec.paste2, sep = "")
    

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
        ## if(is.null(pop)){ # if pop is not provided, merge the two parent populations
        ##     pop <- paste(deparse(substitute(x1)) , deparse(substitute(x2)), sep="-")
        ## }
        pop <- factor(rep(pop,n))

        res <- zyg
        pop(res) <- pop
        res@call <- match.call()

        return(res)
    }

} # end hybridize
