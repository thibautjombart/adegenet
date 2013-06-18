######################################
##
## The code below implements import
## from alignement data.
##
######################################



################
# DNAbin2genind
################
DNAbin2genind <- function(x, pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100){

    ## MISC CHECKS ##
    if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
    if(is.list(x)) {
        x <- as.matrix(x)
    }

    if(is.null(colnames(x))) {
        colnames(x) <- 1:ncol(x)
    }


    ## FUNCTION TO PROCESS ONE LOCUS ##
    ## INPUTS:
    ## locus is a column of a DNAbin matrix
    ## posi is the index of this column
    ## OUTPUTS:
    ## returns NULL if no polymorphism
    ## returns a disjonctive table with named columns otherwise
    ## column names are given as [position.allele]
    processLocus <- function(locus, posi){
        vec <- as.character(locus)
        vec[!vec %in% exp.char] <- NA
        N <- sum(!is.na(vec)) # N: number of sequences
        if(N==0 || sum(table(vec)/N >= polyThres )<2) return(NULL) # escape if untyped locus or no SNPs
        alleles <- unique(na.omit(vec))
        ## vec[!vec %in% alleles] <- NA # this would replace rare alleles by NAs
        out <- sapply(alleles, function(e) 1*(vec==e))
        colnames(out) <- paste(posi, alleles, sep=".")
        return(out)
    }


    ## PROCESS ALL LOCI ##
    ## get disjonctive matrix ##
    ## system.time(res@tab <- Reduce(cbind, lapply(1:ncol(x), function(i) processLocus(x[,i], i)))) # works, but Reduce is real slow
    temp <- lapply(1:ncol(x), function(i) processLocus(x[,i], i)) # process all loci, return a list
    col.names <- unlist(sapply(temp, colnames))
    temp <-  as.matrix(data.frame(temp[!sapply(temp, is.null)])) # remove NULL slots, list -> matrix
    if(is.null(temp) || ncol(temp)==0){
        cat("\nNo polymorphism detected - returning NULL.\n")
        return(NULL)
    }

    ## sort out col/row names ##
    colnames(temp) <- col.names # restore correct names
    rownames(temp) <- rownames(x)

    ## create genind output ##
    res <- genind(temp, ploidy=1, type="codom")
    rm(temp) # remove temp
    gc() # collect garbage
    res$call <- match.call()

    return(res)
} # end DNAbin2genind







####################
## alignment2genind
####################
alignment2genind <- function(x, pop=NULL, exp.char=c("a","t","g","c"), na.char="-", polyThres=1/100){

    ## misc checks
    if(!require(seqinr)) stop("The package seqinr is required.")
    if(!inherits(x,"alignment")) stop("x is not a alignment object")
    N <- length(x$seq)
    if(!is.null(x$nam) && length(x$nam)!=N) stop("Inconsistent names in x (length of x$nam and x$seq do not match). ")


    ## check that na.char does not overide specified exp.char
    if(!is.null(na.char) && na.char %in% exp.char){
        na.char <- na.char[!na.char %in% exp.char]
        if(length(na.char)==0) na.char <- NULL
    }


    ## convert alignment to matrix of characters
    mat <- sapply(x$seq, s2c, USE.NAMES=FALSE)
    if(nrow(mat)!=x$nb){
        mat <- t(mat)
    }

    rownames(mat) <- x$nam

    if(is.null(colnames(x))) {
        colnames(mat) <- 1:ncol(mat)
    }

    ## replace NAs
    if(is.null(na.char)){
        if(is.null(exp.char)) stop("both exp.char and na.char are NULL")
        temp <- paste(exp.char, collapse="", sep="")
        if(any(exp.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must begin the regexp
        }
        temp <- paste("[^", temp, "]", sep="") # anything but the expected is NA
        mat <- gsub(temp,NA,mat)
    } else {
        temp <- paste(na.char, collapse="", sep="")
        if(any(na.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must start the regexp
        }
        temp <- paste("[", temp, "]", sep="")
        mat <- gsub(temp,NA,mat)
    }

    ## keep only columns with polymorphism (i.e., SNPs)
    isPoly <- function(vec){
        N <- sum(!is.na(vec)) # N: number of sequences
        temp <- table(vec)/N
        if(sum(temp > polyThres) >= 2) return(TRUE)
        return(FALSE)
    }

    toKeep <- apply(mat, 2, isPoly)
    if(sum(toKeep)==0) stop("No polymorphic site detected")

    mat <- mat[,toKeep, drop=FALSE]

    ## build output
    res <- df2genind(mat, pop=pop, ploidy=1, ncode=1, type="codom")
    res$call <- match.call()

    if(!is.null(x$com)){
        res@other$com <- x$com
    }

    return(res)
} # end alignment2genind







## #################
## ## findMutations
## #################

## ## GENERIC
## findMutations <- function(...){
##     UseMethod("findMutations")
## }

## ## METHOD FOR DNABIN
## findMutations.DNAbin <- function(x, from=NULL, to=NULL, ...){
##     ## CHECKS ##
##     if(!require(ape)) stop("the ape package is needed")
##     if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
##     x <- as.matrix(x)

##     ## function to pull out mutations from sequence a to b ##
##     NUCL <- c('a','t','g','c')
##     f1 <- function(a,b){
##         seqa <- as.character(x[a,])
##         seqb <- as.character(x[b,])
##         temp <- which(seqa != seqb)
##         ori <- seqa[temp]
##         mut <- seqb[temp]
##         names(ori) <- names(mut) <- temp
##         toRemove <- !ori %in% NUCL | !mut %in% NUCL
##         ori <- ori[!toRemove]
##         mut <- mut[!toRemove]
##         if(all(toRemove)) return(NULL)
##         res <- data.frame(ori,mut)
##         names(res) <- rownames(x)[c(a,b)]
##         res$short <- paste(row.names(res),":",res[,1],"->",res[,2],sep="")
##         return(res)
##     }

##     ## GET LIST OF PAIRS TO COMPARE ##
##     ## handle NULL
##     if(is.null(from)) from <- 1:nrow(x)
##     if(is.null(to)) to <- 1:nrow(x)

##     ## get pairs
##     pairs <- expand.grid(from, to)

##     ## remove unwanted comparisons
##     pairs <- pairs[pairs[,1]!=pairs[,2],,drop=FALSE]

##     ## GET NUMBER OF MUTATIONS ##
##     out <- lapply(1:nrow(pairs), function(i) f1(pairs[i,1], pairs[i,2]))
##     names(out) <- paste(rownames(x)[pairs[,1]], rownames(x)[pairs[,2]],sep="->")

##     return(out)

## } # end findMutations







## ##################
## ## graphMutations
## ##################

## ## GENERIC
## graphMutations <- function(...){
##     UseMethod("graphMutations")
## }

## ## METHOD FOR DNABIN
## graphMutations.DNAbin <- function(x, from=NULL, to=NULL, plot=TRUE, edge.curved=TRUE, ...){
##     if(!require(igraph)) stop("igraph is required")

##     ## GET MUTATIONS ##
##     x <- findMutations(x, from=from, to=to)

##     ## GET GRAPH ##
##     from <- gsub("->.*","",names(x))
##     to <- gsub(".*->","",names(x))
##     vnames <- sort(unique(c(from,to)))
##     dat <- data.frame(from,to,stringsAsFactors=FALSE)
##     out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(vnames, label=vnames))

##     ## SET ANNOTATIONS FOR THE BRANCHES ##
##     annot <- unlist(lapply(x, function(e) paste(e$short, collapse="\n")))
##     E(out)$label <- annot
##     E(out)$curved <- edge.curved

##     ## PLOT / RETURN ##
##     if(plot) plot(out, ...)

##     return(out)
## } # end graphMutations














## ###############
## ## transiProb
## ###############
## ##
## ## proba/distance based on transition prob from one sequence to another
## ## time is taken into account
## ## output: matrix with term proba(rowIdx to colIdx)
## ##
## transiProb <- function(x, mu, dates, result=c("prob","dist")){
##     ## MISC CHECKS ##
##     if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
##     if(!require(ape)) stop("The package ape is required.")
##     result <- match.arg(result)

##     ## COMPUTATIONS ##

##     ## get numbers of differing nucleotides between sequences
##     seq.length <- ncol(as.matrix(x))
##     D <- as.matrix(dist.dna(x, model="raw")) * seq.length
##     ## if(sum(D-round(D)) > 1e-10){ # make sure we've got integers there
##     ##         warning("Number of nucleotides are not all integers")
##     ##     }
##     D <- round(D)

##     ## compute matrix T (time between sequences)
##     if(inherits(dates,"POSIXct")){ # dates in POSIXct format
##         temp <- outer(dates, dates, difftime, unit="days")
##         T <- -matrix(as.numeric(temp),ncol=length(dates))
##     } else { # dates are numeric
##         T <- -outer(dates, dates, "-")
##     }

##     ## spot negative times
##     toSetToNull <- T < 1e-15

##     ## compute proba(no change @ a site) term
##     mu <- mu/365 # express mu per day
##     p1 <- exp(-T*mu) + (1-exp(-T*mu))/4
##     p1[toSetToNull] <- 0
##     res <- dbinom(D, size=seq.length, prob=(1-p1))

##     ## PROCESS/RETURN RESULT
##     if(result=="prob"){ # return probabilities
##         res[toSetToNull] <- 0
##         diag(res) <- 1
##     } else { # return d = -log(proba)
##         res <- -log(res)
##         res[toSetToNull] <- 1e15
##         diag(res) <- 0
##     }

##     return(res)
## } # end transiProb
