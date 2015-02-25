###################################################################
## Fonctions designed to import files from other softwares
## into genind objects
##
## currently supported formats are :
## .gtx (GENETIX)
## .dat (Fstat)
## .gen (Genepop)
## .stru (STRUCTURE)
##
## Thibaut Jombart, avril 2006
## t.jombart@imperial.ac.uk
##
##################################################################



#####################
# Function df2genind
#####################
df2genind <- function(X, sep=NULL, ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL, missing=NA, ploidy=2, type=c("codom","PA")){

    if(is.data.frame(X)) X <- as.matrix(X)
    if (!inherits(X, "matrix")) stop ("X is not a matrix")

    res <- list()
    type <- match.arg(type)
    ## checkType(type)


    ## type-independent stuff ##
    n <- nrow(X)
    nloc <- ncol(X)
    ploidy <- as.integer(ploidy)
    if(ploidy < 1L) stop("ploidy cannot be less than 1")

    if(is.null(ind.names)) {ind.names <- rownames(X)}
    if(is.null(loc.names)) {loc.names <- colnames(X)}

    ## pop optionnelle
    if(!is.null(pop)){
        if(length(pop)!= n) stop("length of factor pop differs from nrow(X)")
        pop <- as.factor(pop)
    }


    ## PA case ##
    if(toupper(type)=="PA"){
        ## preliminary stuff
        mode(X) <- "numeric"
        rownames(X) <- ind.names
        colnames(X) <- loc.names

        ## Erase entierely non-typed loci
        temp <- apply(X,2,function(c) all(is.na(c)))
        if(any(temp)){
            X <- X[,!temp]
            warning("entirely non-type marker(s) deleted")
        }

        ## Erase entierely non-type individuals
        temp <- apply(X,1,function(r) all(is.na(r)))
        if(any(temp)){
            X <- X[!temp,,drop=FALSE]
            pop <- pop[!temp]
            warning("entirely non-type individual(s) deleted")
        }

        ## erase non-polymorphic loci
        temp <- apply(X, 2, function(loc) length(unique(loc[!is.na(loc)]))==1)
        if(any(temp)){
            X <- X[,!temp,drop=FALSE]
            warning("non-polymorphic marker(s) deleted")
        }

        prevcall <- match.call()

        res <- genind( tab=X, pop=pop, prevcall=prevcall, ploidy=ploidy, type="PA")

        return(res)
    } # end type PA


    ## codom case ##

    ## make sure X is in character mode
    mode(X) <- "character"


    ## find or check the number of coding characters, 'ncode'
    if(is.null(sep)){
        if(!is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            if(ncode <  max(temp) ) stop("some character strings exceed the provided ncode.")
        }
        if(is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            ncode <- max(temp)
        }
        if((ncode %% ploidy)>0) stop(paste(ploidy,"alleles cannot be coded by a total of",
                                           ncode,"characters", sep=" "))
    }

    ## ERASE ENTIRELY NON-TYPE LOCI AND INDIVIDUALS
    tempX <- X
    if(!is.null(sep)) tempX <- gsub(sep,"",X)
    ## turn NANANA, 00000, ... into NA
    tempX <- gsub("^0*$",NA,tempX)
    tempX <- gsub("(NA)+",NA,tempX)

    ## Erase entierely non-typed loci
    temp <- apply(tempX,2,function(c) all(is.na(c)))
    if(any(temp)){
        X <- X[,!temp,drop=FALSE]
        tempX <- tempX[,!temp,drop=FALSE]
        loc.names <- loc.names[!temp]
        nloc <- ncol(X)
        warning("entirely non-type marker(s) deleted")
    }

    ## Erase entierely non-type individuals
    temp <- apply(tempX,1,function(r) all(is.na(r)))
    if(any(temp)){
        X <- X[!temp,,drop=FALSE]
        tempX <- tempX[!temp,,drop=FALSE]
        pop <- pop[!temp]
        ind.names <- ind.names[!temp]
        n <- nrow(X)
        warning("entirely non-type individual(s) deleted")
    }

    n <- nrow(X)
    ## SET NAs IN X
    X[is.na(tempX)] <- NA

    # ind.names <- rownames(X) this erases the real labels
    # note: if X is kept as a matrix, duplicate row names are no problem


    ## function to fill a matrix of char 'M' with the required
    ## number of zero, targetN being the total number of char required
    fillWithZero <- function(M, targetN){
        naIdx <- is.na(M)
        keepCheck <- any(nchar(M) < targetN)
        while(keepCheck){
            mat0 <- matrix("", ncol=ncol(M), nrow=nrow(M))
            mat0[nchar(M) < targetN] <- "0"
            M <-  matrix(paste(mat0, M, sep=""), nrow=nrow(mat0))
            keepCheck <- any(nchar(M) < targetN)
        }

        ## restore NA (otherwise we're left with "NA")
        M[naIdx] <- NA
        return(M)
    }

    ## CHECK STRING LENGTH IF NO SEPARATOR PROVIDED
    if(is.null(sep) | ploidy==as.integer(1)){
        ##     ## now check all strings and make sure they all have 'ncode' characters
        ##         ## NA are temporarily coded as "00", "000" or "000000" to fit the check
        ##         keepCheck <- any(nchar(X) < ncode)
        ##         missAll <- paste(rep("0",ncode/ploidy),collapse="")
        ##         missTyp <- paste(rep("0",ncode),collapse="")
        ##         X[is.na(X)] <- missTyp

        ##         while(keepCheck){
        ##             mat0 <- matrix("", ncol=ncol(X), nrow=nrow(X))
        ##             mat0[nchar(X) < ncode] <- "0"
        ##             X <-  matrix(paste(mat0, X, sep=""), nrow=nrow(mat0))
        ##             keepCheck <- any(nchar(X) < ncode)
        ##         }

        X <- fillWithZero(X,targetN=ncode)

        ## now split X by allele
        splitX <- list()
        for(i in 1:ploidy){
            splitX[[i]] <- substr(X,1,ncode/ploidy)
            X <- sub(paste("^.{",ncode/ploidy,"}",sep=""),"",X)
        }

    } # END CHECK STRING LENGTH WITHOUT SEP


    ## CHECK STRING LENGTH WITH SEPARATOR PROVIDED
    if(!is.null(sep)){
        if(ploidy > 1){
            temp <- t(as.matrix(as.data.frame(strsplit(X,sep))))
            splitX <- list()
            for(i in 1:ncol(temp)){
                splitX[[i]] <- matrix(temp[,i], nrow=n)
            } # each matrix of splitX contains typing for 1 allele
        } else {
            splitX <- list()
            splitX[[1]] <- X
        }

        ## get the right ncode
        temp <- unlist(splitX)
        temp <- temp[!is.na(temp)]
        ncode <- max(nchar(temp))*ploidy
        splitX <- lapply(splitX, function(Y) fillWithZero(Y,targetN=ncode/ploidy))
    } # END CHECK STRING LENGTH WITH SEP


    ## AT THIS STAGE, splitX IS A LIST OF MATRICES,
    ## EACH GIVING TYPING FOR AN ALLELE

    ## fetch all possible alleles per locus
    loc.all <- list()
    for(i in 1:nloc){
        temp <- unlist(lapply(splitX,function(e) e[,i]))
        loc.all[[i]] <- sort(unique(temp[!is.na(temp)]))
    }

    names(loc.all) <- loc.names
    ## loc.all is a list whose element are vectors of sorted possible alleles at a locus
    temp <- lapply(1:nloc, function(i) matrix(0,nrow=n,ncol=length(loc.all[[i]]),
       dimnames=list(NULL,loc.all[[i]])) )

    names(temp) <- loc.names
    # note: keep rownames as NULL in case of duplicates
    ## temp is a list whose elements are one matrix (indiv x alleles) for each marker

    ## now tables in 'temp' are filled up
    findall <- function(cha,loc.all){
        if(is.na(cha)) return(NULL)
        return(which(cha==loc.all))
    }

    for(k in 1:ploidy){
        for(i in 1:n){
            for(j in 1:nloc){
                allIdx <- findall(splitX[[k]][i,j],loc.all[[j]])
                temp[[j]][i,allIdx] <- temp[[j]][i,allIdx] + 1
                if(is.null(allIdx)) {temp[[j]][i,] <- NA}
            }
        }
    }

    ## beware: colnames are wrong when there is only one allele in a locus
    ## right colnames are first generated
    nall <- unlist(lapply(temp,ncol))
    loc.rep <- rep(names(nall),nall)
    col.lab <- paste(loc.rep,unlist(loc.all,use.names=FALSE),sep=".")

    ## mat <- as.matrix(cbind.data.frame(temp)) # ! does not work for huge numbers of alleles
    mat <- matrix(unlist(temp), nrow=nrow(temp[[1]]))
    mat <- mat/ploidy
    colnames(mat) <- col.lab
    rownames(mat) <- ind.names

    if(!is.na(missing)){
      if(missing==0) {mat[is.na(mat)] <- 0}
      if(toupper(missing)=="MEAN") {
        moy <- apply(mat,2,function(c) mean(c,na.rm=TRUE))
        for(j in 1:ncol(mat)) {mat[,j][is.na(mat[,j])] <- moy[j]}
      }
    }

    prevcall <- match.call()

    res <- genind( tab=mat, pop=pop, prevcall=prevcall, ploidy=ploidy, type=type)

    return(res)
} # end df2genind





########################################
# Function read.genetix
# code based on previous ade4 functions
########################################
read.genetix <- function(file=NULL,missing=NA,quiet=FALSE) {
    if(!quiet) cat("\n Converting data from GENETIX to a genind object... \n")


    ## read from file
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed

    if(toupper(.readExt(file)) != "GTX") stop("File extension .gtx expected")
    ## retrieve first infos
    nloc <- as.integer(scan(file,nlines=1,what="character",quiet=TRUE)[1])
    npop <- as.integer(scan(file,nlines=1,skip=1,what="character",quiet=TRUE)[1])
    txt <- scan(file,skip=2,what="character",sep="\n",quiet=TRUE)
    txt <- gsub("\t"," ",txt)
    ## check that nloc is consistent with actual nloc (bug-report 1.2-2.02)
    temp <- temp <- .rmspaces(txt[length(txt)])
    nlocbis <- length(unlist(strsplit(temp, "[[:space:]]+")))-1
    if(nloc != nlocbis) {
        warning(paste("\n== Genetix file error == \n",
                      "Indicated number of locus (", nloc, ")\n",
                      "does not match actual number (", nlocbis, ").\n",
                      "Using ", nlocbis, " as number of locus.\n",
                      "Please check your file.", sep=""))
        nloc <- nlocbis
    }
    loc.names <- txt[seq(1,by=2,length=nloc)]
    txt <- txt[-(1:(nloc*2))]

    ## retrieve populations infos
    pop.names <- vector(mode="character",length=npop)
    pop.nind <- vector(mode="integer",length=npop)
    index <- 1
    temp <- vector(mode="integer",length=npop)
    for(i in 1:npop){
        pop.names[i] <- txt[index]
        pop.nind[i] <- as.numeric(txt[index+1])
        temp[i] <- index
        index <- index + pop.nind[i] + 2
    }
    pop.names <- .rmspaces(pop.names)

    ## retrieve genotypes infos
    txt <- txt[-c(temp,temp+1)]
    txt <- .rmspaces(txt)
    txt <- sapply(1:length(txt),function(i) unlist(strsplit(txt[i],"([[:space:]]+)|([[:blank:]]+)")) )
    X <- t(txt)
    if(ncol(X) == (nloc+1)){
        rownames(X) <- X[,1]
        X <- X[,-1]
    } else{
        rownames(X) <- 1:nrow(X)
    }

    colnames(X) <- loc.names

    ## make a factor "pop" if there is more than one population
    pop <- factor(rep(pop.names,pop.nind))

    ## pass X to df2genind
    res <- df2genind(X=X, ncode=6, pop=pop, missing=missing, ploidy=2)
    res@call <- match.call()

    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end read.genetix





######################
# Function read.fstat
######################
read.fstat <- function(file,missing=NA,quiet=FALSE){
    ##if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "DAT") stop("File extension .dat expected")

    if(!quiet) cat("\n Converting data from a FSTAT .dat file to a genind object... \n\n")

    call <- match.call()
    txt <- scan(file,what="character",sep="\n",quiet=TRUE)
    txt <- gsub("\t"," ",txt)

                                        # read first infos
    info <- unlist(strsplit(txt[1],"([[:space:]]+)"))
                                        # npop <- as.numeric(info[1]) ## no longer used
    nloc <- as.numeric(info[2])

    loc.names <- txt[2:(nloc+1)]

                                        # build genotype matrix
    txt <- txt[-(1:(nloc+1))]
    txt <- .rmspaces(txt)
    txt <- sapply(1:length(txt),function(i) unlist(strsplit(txt[i],"([[:space:]]+)|([[:blank:]]+)")) )
    X <- t(txt)
    pop <- factor(X[,1])
    if(length(levels(pop)) == 1 ) pop <- NULL
    X <- X[,-1]

    colnames(X) <- loc.names
    rownames(X) <- 1:nrow(X)

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)
                                        # beware : fstat files do not yield ind names
    res@ind.names <- rep("",length(res@ind.names))
    names(res@ind.names) <- rownames(res@tab)
    res@call <- call

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.fstat





##########################
# Function read.genepop
##########################
read.genepop <- function(file,missing=NA,quiet=FALSE){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "GEN") stop("File extension .gen expected")

    if(!quiet) cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")

    prevcall <- match.call()

    txt <- scan(file,sep="\n",what="character",quiet=TRUE)
    if(!quiet) cat("\nFile description: ",txt[1], "\n")
    txt <- txt[-1]
    txt <- gsub("\t", " ", txt)

  # two cases for locus names:
  # 1) all on the same row, separated by ","
  # 2) one per row
  # ! spaces and tab allowed
  # a bug was reported by S. Devillard, occuring
  # when the two cases occur together,
  # that is:
  # loc1,
  # loc2,
  # ...

  ### former version
  #1
  #if(length(grep(",",txt[1])) > 0){
  #  loc.names <- unlist(strsplit(txt[1],","))
  #  loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
  #  loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)
  #  nloc <- length(loc.names)

  #  txt <- txt[-1]
  #} else { #2
  #  nloc <- min(grep("POP",toupper(txt)))-1
  #  loc.names <- txt[1:nloc]
  #  loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
  #  loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)

  #  txt <- txt[-(1:nloc)]
  #}

    ## new strategy (shorter): isolate the 'locus names' part and then parse it.
    locinfo.idx <- 1:(min(grep("POP",toupper(txt)))-1)
    locinfo <- txt[locinfo.idx]
    locinfo <- paste(locinfo,collapse=",")
    loc.names <- unlist(strsplit(locinfo,"([,]|[\n])+"))
    loc.names <- .rmspaces(loc.names)
    nloc <- length(loc.names)
    txt <- txt[-locinfo.idx]

    ## locus names have been retreived

    ## build the pop factor
    ## and correct the genotypes splited on more than 1 line
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
    npop <- length(pop.idx)
    ## correction for splited genotype
    ## isolated by the absence of comma on a line not containing "pop"
    nocomma <- which(! (1:length(txt)) %in% grep(",",txt))
    splited <- nocomma[which(! nocomma %in% pop.idx)]
    if(length(splited)>0){
        for(i in sort(splited,decreasing=TRUE)){
            txt[i-1] <- paste(txt[i-1],txt[i],sep=" ")
        }
        txt <- txt[-splited]
    }
    ## end correction

    ## reevaluate pop index
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))

    txt[length(txt)+1] <- "POP"
    nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt)))-1
    pop <- factor(rep(1:npop,nind.bypop))

    txt <- txt[-c(pop.idx,length(txt))]

    temp <- sapply(1:length(txt),function(i) strsplit(txt[i],","))
    ## temp is a list with nind elements, first being ind. name and 2nd, genotype

    ind.names <- sapply(temp,function(e) e[1])
    ind.names <- .rmspaces(ind.names)
    ## individuals' name are now clean

    vec.genot <- sapply(temp,function(e) e[2])
    vec.genot <- .rmspaces(vec.genot)

    ## X is a individual x locus genotypes matrix
    X <- matrix(unlist(strsplit(vec.genot,"[[:space:]]+")),ncol=nloc,byrow=TRUE)

    rownames(X) <- ind.names
    colnames(X) <- loc.names

 ##  # correct X to fulfill the genetix format
##   f1 <- function(char){
##     paste("00", substr(char,1,1), "00", substr(char,2,2), sep="")
##   }

##   f2 <- function(char){
##     paste("0", substr(char,1,2), "0", substr(char,3,4), sep="")
##   }

##   if(all(nchar(X)==2)) {X <- apply(X,c(1,2),f1)}
##   if(all(nchar(X)==4)) {X <- apply(X,c(1,2),f2)}

    ## give right pop names
    ## beware: genepop takes the name of the last individual of a sample as this sample's name
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)
    res@call <- prevcall

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.genepop





############################
# Function read.structure
############################
read.structure <- function(file, n.ind=NULL, n.loc=NULL,  onerowperind=NULL, col.lab=NULL, col.pop=NULL, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, missing=NA, ask=TRUE, quiet=FALSE){

    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(!toupper(.readExt(file)) %in% c("STR","STRU")) stop("File extension .stru expected")

    ## set defaults for non optional arguments without default values
    if(!ask){
        if(is.null(col.lab)) col.lab <- as.integer(0)
        if(is.null(col.pop)) col.pop <- as.integer(0)
        if(is.null(row.marknames)) row.marknames <- as.integer(0)
    }

    ## required questions
    if(is.null(n.ind)){
        cat("\n How many genotypes are there? ")
        n.ind <- as.integer(readLines(n = 1))
    }

    if(is.null(n.loc)){
        cat("\n How many markers are there? ")
        n.loc <- as.integer(readLines(n = 1))
    }

    if(is.null(col.lab)){
        cat("\n Which column contains labels for genotypes ('0' if absent)? ")
        col.lab <- as.integer(readLines(n = 1))
    }

    if(is.null(col.pop)){
        cat("\n Which column contains the population factor ('0' if absent)? ")
        col.pop <- as.integer(readLines(n = 1))
    }

    if(is.null(col.others) & ask){
        cat("\n Which other optional columns should be read (press 'return' when done)? ")
        col.others <- scan(quiet=TRUE)
        if(length(col.others) == 0)  col.others <- NULL
    }

    if(is.null(row.marknames)){
        cat("\n Which row contains the marker names ('0' if absent)? ")
        row.marknames <- as.integer(readLines(n = 1))
    }

    if(is.null(onerowperind)){
        cat("\n Are genotypes coded by a single row (y/n)? ")
        onerowperind <- toupper(readLines(n = 1))
        if(onerowperind == "Y") {
            onerowperind <- TRUE
        } else {
            onerowperind <- FALSE
        }
    }

    if(is.null(NA.char)){
        cat("\n What is the code for missing data (default is '-9')? ")
        NA.char <- as.character(readLines(n = 1))
    }

    ## message to console
    if(!quiet) cat("\n Converting data from a STRUCTURE .stru file to a genind object... \n\n")

    ## read the file
    txt <- scan(file,sep="\n",what="character",quiet=TRUE)

    ## remove empty lines and spaces/tabs at the end of a line
    temp <- grep("^[[:space:]]*$",txt)
    if(length(temp) > 0) {
        txt <- txt[-temp]
    }

    txt <- gsub("([[:blank:]]+)$","",txt)

    ## isolate each useful component of the file
    ## matrix of data
    if(onerowperind) {
        n <- n.ind
        p <- 2*n.loc
    } else{
        n <- 2*n.ind
        p <- n.loc
    }

    lastline <- length(txt)
    mat <- txt[(lastline-n+1):lastline]
    mat <- t(as.data.frame(strsplit(mat,"[[:blank:]]+")))
    rownames(mat) <- 1:n
    gen <- mat[, (ncol(mat)-p+1):ncol(mat)]


    ## markers names
    if(row.marknames != 0) {
        loc.names <- .rmspaces(txt[row.marknames])
        loc.names <- unlist(strsplit(loc.names,"[[:blank:]]+"))
    } else {
        loc.names <- .genlab("L",n.loc)
    }

    ## genotypes labels
    if(col.lab !=0) {
        ind.names <- mat[, col.lab]
    } else {
        ind.names <- .genlab("",n.ind)
    }

    ## population factor
    if(col.pop !=0) {
        pop <- factor(mat[, col.pop])
    } else {
        pop <- NULL
    }

    ## other variables
    if(!is.null(col.others)){
        X.other <- mat[,col.others]
    }

    ## transformations if onerowperind is FALSE
    if(!onerowperind) {
        temp <- seq(1,n,by=2)
        ind.names <- ind.names[temp]
        if(length(ind.names) < n.ind) warning("Duplicated identifier for genotypes")
        pop <- pop[temp]
        if(exists("X.other")) X.other <- X.other[temp]

        ## make sur that all strings in gen have the same number of characters
        ncode <- max(nchar(gen))
        keepCheck <- any(nchar(gen) < ncode)

        while(keepCheck){
            mat0 <- matrix("", ncol=ncol(gen), nrow=nrow(gen))
            mat0[nchar(gen) < ncode] <- "0"
            gen <-  matrix(paste(mat0, gen, sep=""), nrow=nrow(mat0))
            keepCheck <- any(nchar(gen) < ncode)
        }

        ## reorder matrix of genotypes
        X <- t(sapply(temp, function(i) paste(gen[i,],gen[i+1,],sep="") ))

    } else { # else of "if(!onerowperind)"
        temp <- seq(1,p-1,by=2)
        X <- paste(gen[,temp] , gen[,temp+1], sep="")
        X <- matrix(X, nrow=n.ind)
    }

    ## replace missing values by NAs
    X <- gsub(NA.char,NA,X)
    rownames(X) <- ind.names
    colnames(X) <- loc.names

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)

    res@call <- match.call()

    if(exists("X.other")) {res@other <- list(X=X.other)}

    return(res)

}




#########################
# Function import2genind
#########################
import2genind <- function(file,missing=NA,quiet=FALSE, ...){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    ext <- .readExt(file)
    ext <- toupper(ext)

    if(ext == "GTX")
        return(read.genetix(file,missing=missing,quiet=quiet))

    if(ext == "DAT")
        return(read.fstat(file,missing=missing,quiet=quiet))

    if(ext == "GEN")
        return(read.genepop(file,missing=missing,quiet=quiet))

    if(ext %in% c("STR","STRU"))
        return(read.structure(file,missing=missing,quiet=quiet, ...))

    ## evaluated only if extension is not supported
    cat("\n File format (",ext,") not supported.\n")
    cat("\nSupported formats are:\nGENETIX (.gtx) \nFSTAT (.dat) \nGenepop (.gen)\n \nSTRUCTURE (.str)\n")

    return(invisible())
}







#######################
# Function read.snp
#######################
read.snp <- function(file, quiet=FALSE, chunkSize=1000,
                     parallel=require("parallel"), n.cores=NULL, ...){
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(ext != "SNP") warning("wrong file extension - '.snp' expected")
    if(!quiet) cat("\n Reading biallelic SNP data file into a genlight object... \n\n")
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }

    call <- match.call()


    ## HANDLE THE COMMENTS ##
    if(!quiet) cat("\n Reading comments... \n")

    count <- 0L
    i <- 0
    while(count < 2L){
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=i, nmax=1, blank.lines.skip=FALSE)
        if(length(grep(">>>>", txt)>0)){
            count <- count + 1L
        }
        i <- i+1
        if(count==0L && i>10){
            warning("No comment section at the beginning of the file. Format may be wrong.")
            i <- 0
            break
        }
    }

    lines.to.skip <- i


    ## READ GENERAL DATA (>>) ##
    if(!quiet) cat("\n Reading general information... \n")

    misc.info <- list()

    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=1)

    while(length(grep(">>", txt))>0){
        itemName <- gsub(">>","", txt)
        itemName <- gsub("(^[[:space:]]+)|([[:space:]]+$)", "", itemName)
        misc.info[itemName] <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip + 1, nmax=1)
        lines.to.skip <-lines.to.skip + 2
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=1)
    }

    ## transform each character string into a vector
    misc.info <- sapply(misc.info, function(e) unlist(strsplit(e,"[[:space:]]+")))

    ## READ GENOTYPE DATA ##
    ## one genotype is read/converted at a time to spare RAM
    if(!quiet) cat("\n Reading",ifelse(is.null(misc.info$population),"",length(misc.info$population)), "genotypes... \n")


    res <- list() # this will be a list of SNPbin objects

    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunkSize*2)

    ID.INDIV <- grep(">", txt)

    COUNT <- 0 # used to count the nb reads

    while(length(ID.INDIV)>0){
        COUNT <- COUNT + 1
        if(!quiet) {
            if(COUNT %% 5 == 0){
                cat(length(res)+length(ID.INDIV))
            } else {
                cat(".")
            }
        }

        ind.lab <- gsub(">","", txt[ID.INDIV])
        ind.lab <- gsub("(^[[:space:]]+)|([[:space:]]+$)", "", ind.lab)
        temp <- strsplit(txt[ID.INDIV+1], "")
        temp <- lapply(temp, function(e) suppressWarnings(as.integer(e)))
        if(parallel){
            res <- c(res, mclapply(temp, function(e) new("SNPbin", e),
                                   mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) )
        } else {
            res <- c(res, lapply(temp, function(e) new("SNPbin", e)) )
        }
        names(res)[(length(res)-length(ID.INDIV)+1):length(res)] <- ind.lab
        lines.to.skip <-lines.to.skip + length(txt)
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunkSize*2)
        ID.INDIV <- grep(">", txt)
    }


    ## CHECK CONSISTENCY ##
    if(!quiet) cat("\n Checking consistency... \n")

    n.loc <- unique(sapply(res, nLoc))
    n.ind <- length(res)
    other <- list()

    if(length(n.loc)>1) {
        print(n.loc)
        warning("!!! Differing numbers of loci detected between individuals !!!")
    }

    if(!is.null(misc.info$position) && length(misc.info$position)!=n.loc) {
        other$position <- misc.info$position
        misc.info$position <- NULL
        warning("vector of positions of the SNPs does not match the number of SNPs - storing this information in @other")
    }
    if(!is.null(misc.info$allele) && length(misc.info$allele)!=n.loc) {
        other$allele <- misc.info$allele
        misc.info$allele <- NULL
        warning("vector of alleles of the SNPs does not match the number of SNPs - storing this information in @other")
    }
    if(!is.null(misc.info$chromosome) && length(misc.info$chromosome)!=n.loc) {
        other$chromosome <- misc.info$chromosome
        misc.info$chromosome <- NULL
        warning("vector of chromosomes of the SNPs does not match the number of SNPs - storing this information in @other")
    }
    if(!is.null(misc.info$population) && length(misc.info$population)!=n.ind) {
        other$population <- misc.info$population
        misc.info$population <- NULL
        warning("vector of population of the individuals does not match the number of individuals - storing this information in @other")
    }
    if(!is.null(misc.info$ploidy) && length(misc.info$ploidy)>1 && length(misc.info$ploidy)!=n.ind) {
        other$ploidy <- misc.info$ploidy
        misc.info$ploidy <- NULL
        warning("vector of ploidy of the individuals has more than one value but does not match the number of individuals - storing this information in @other")
    }


    ## BUILD OUTPUT ##
    if(!quiet) cat("\n Building final object... \n")

    ind.names <- names(res)
    if(!is.null(misc.info$chromosome)){
        other <- list(chromosome = misc.info$chromosome)
    }

    res <- new("genlight", gen=res, ind.names=ind.names, position=misc.info$position, loc.all=misc.info$allele, ploidy=misc.info$ploidy, pop=misc.info$population, other=other, parallel=parallel)

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.snp







####################
## extract.PLINKmap
####################
extract.PLINKmap <- function(file, x=NULL){
    ## CHECK EXTENSION ##
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(ext != "MAP") warning("wrong map.file extension - '.map' expected")


    ## READ FILE ##
    ## find nb of columns
    txt <- scan(file,what="character",sep="\n",quiet=TRUE,  nlines=1)
    nb.col <- length( unlist(strsplit(txt,"[[:blank:]]+")))

    ## read file
    txt <- scan(file,what="character",sep="\t",quiet=TRUE)
    txt <- matrix(txt, ncol=4, byrow=TRUE)


    ## EXTRACT INFO AND RETURN OBJECT ##
    ## return a genlight
    if(!is.null(x)){
        ## match data
        ord <- match(locNames(x), txt[,2]) # check that it is the 2nd column
        if(!inherits(x, "genlight")) stop("x is not a genlight object")
        other(x)$chromosome <- factor(txt[ord,1])
        other(x)$position <- as.integer(txt[ord,4])

        return(x)
    }

    ## return a list
    res <- list(chromosome=factor(txt[ord,1]), position=as.integer(txt[ord,4]))

    return(res)
} # end extract.PLINKmap







########################
## Function read.PLINK
########################
read.PLINK <- function(file, map.file=NULL, quiet=FALSE, chunkSize=1000,
                       parallel=require("parallel"), n.cores=NULL, ...){
    ## HANDLE ARGUMENTS ##
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(ext != "RAW") warning("wrong file extension - '.raw' expected")
    if(!quiet) cat("\n Reading PLINK raw format into a genlight object... \n\n")
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }


    ## READ NAMES OF LOCI ##
    if(!quiet) cat("\n Reading loci information... \n")

    loc.names <- scan(file,what="character",sep=" ",quiet=TRUE,  nlines=1, blank.lines.skip=FALSE)
    n.loc <- length(loc.names) - 6
    misc.info <- lapply(1:6,function(i) NULL)
    names(misc.info) <- loc.names[1:6]
    loc.names <- loc.names[7:length(loc.names)]
    loc.names <- gsub("_[1-9]$","",loc.names)

    ## READ GENOTYPES ##
    if(!quiet) cat("\n Reading and converting genotypes... \n")

    res <- list() # this will be a list of SNPbin objects

    ## initialize reading
    lines.to.skip <- 1
    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunkSize)
    txt <- lapply(txt, function(e) unlist(strsplit(e,"[[:blank:]]+") ))

    COUNT <- 0 # used to count the nb reads

    while(length(txt)>0){
        COUNT <- COUNT + 1
        if(!quiet) {
            if(COUNT %% 5 == 0){
                cat(length(res)+length(txt))
            } else {
                cat(".")
            }
        }


        ## handle misc info
        temp <- lapply(txt, function(e) e[1:6])
        for(i in 1:6){
            misc.info[[i]] <- c(misc.info[[i]], unlist(lapply(temp, function(e) e[[i]])) )
        }


        ## build SNPbin objects
        txt <- lapply(txt, function(e) suppressWarnings(as.integer(e[-(1:6)])))

        if(parallel){
            res <- c(res, mclapply(txt, function(e) new("SNPbin", snp=e, ploidy=2),
                                   mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) )
        } else {
            res <- c(res, lapply(txt, function(e) new("SNPbin", snp=e, ploidy=2)) )
        }

        lines.to.skip <-lines.to.skip + length(txt)

        ## read lines
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=chunkSize)
        txt <- lapply(txt, function(e) unlist(strsplit(e,"[[:blank:]]+") ))
    }


    ## MAKE A FEW CHECKS ##
    if(!all(sapply(res, nLoc)==n.loc)) stop(paste("some individuals do not have",n.loc,"SNPs."))


    ## BUILD FINAL OBJECT ##
    if(!quiet) cat("\n Building final object... \n")

    res <- new("genlight",res, ploidy=2, parallel=parallel)
    indNames(res) <- misc.info$IID
    pop(res) <- misc.info$FID
    locNames(res) <- loc.names
    misc.info <- misc.info[c("SEX", "PHENOTYPE", "PAT","MAT")]
    names(misc.info) <- tolower(names(misc.info))
    misc.info$sex[misc.info$sex==1] <- "m"
    misc.info$sex[misc.info$sex==2] <- "f"
    misc.info$sex <- factor(misc.info$sex)
    misc.info$phenotype[misc.info$phenotype==1] <- "control"
    misc.info$phenotype[misc.info$phenotype==2] <- "case"
    misc.info$phenotype <- factor(misc.info$phenotype)
    other(res) <- misc.info


    ## HANDLE MAP FILE INFO ##
    if(!is.null(map.file)){
        res <- extract.PLINKmap(map.file, res)
    }


    ## RETURN OUTPUT ##
    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end read.PLINK







###########################
## Function fasta2genlight
###########################
fasta2genlight <- function(file, quiet=FALSE, chunkSize=1000, saveNbAlleles=FALSE,
                       parallel=require("parallel"), n.cores=NULL, ...){
    ## HANDLE ARGUMENTS ##
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(!ext %in% c("FASTA", "FA", "FAS")) warning("wrong file extension - '.fasta', '.fa' or '.fas' expected")
    if(!quiet) cat("\n Converting FASTA alignment into a genlight object... \n\n")
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }


    ## PRIOR CHECKS ##
    ## find nb of lines per genome
    lines.to.skip <- 0
    txt <- scan(file,what="character",sep="\n",quiet=TRUE, nmax=1)

    while(length(grep("^>.+", txt))<2){
        lines.to.skip <- lines.to.skip + 1
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, nmax=lines.to.skip)
    }

    LINES.PER.IND <- lines.to.skip-1


    ## find length of a genome
    NLOC <- sum(nchar(txt[2:LINES.PER.IND]))


    ## SCAN ALL POSITIONS AND IDENTIFY SNPs ##
    if(!quiet) cat("\n Looking for polymorphic positions... \n")

    ## read all genomes by chunks
    ## initialize
    lines.to.skip <- 0
    IND.LAB <- NULL
    POOL <- as.list(rep("-", NLOC))
    COUNT <- 0 # used to count the nb reads

    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=LINES.PER.IND*chunkSize)

    ## read and process chunks
    while(length(txt)>0){
        COUNT <- COUNT + 1
        if(!quiet) {
            for(i in 1:(COUNT*chunkSize)) cat(".")
        }


        nb.ind <- length(grep("^>", txt))
        IND.LAB <- c(IND.LAB, sub(">","",txt[grep("^>", txt)])) # find individuals' labels
        txt <- split(txt, rep(1:nb.ind, each=LINES.PER.IND)) # split per individuals
        if(parallel){
            txt <- mclapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split=""),
                            mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) # each genome -> one vector
        } else {
            txt <- lapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split="")) # each genome -> one vector
        }

        ## POOL contains all alleles of each position
        temp <- as.list(apply(matrix(unlist(txt), byrow=TRUE, nrow=length(txt)),2,unique)) # alleles current genomes
        POOL <- mapply(function(x,y) unique(c(x,y)), POOL, temp, SIMPLIFY=FALSE) # update global pool

        lines.to.skip <- lines.to.skip + nb.ind*LINES.PER.IND
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=LINES.PER.IND*chunkSize)
    }


    ## analyse pool of alleles
    letterOK <- c("a","g","c","t","A","G","C","T")
    POOL <- lapply(POOL, function(e) e[e %in% letterOK]) # keep only proper letters
    ## POOL <- lapply(POOL, setdiff, "-")
    nb.alleles <- sapply(POOL, length)
    snp.posi <- nb.alleles==2
    if(all(!snp.posi)){
        warning("No polymorphism in the alignment - returning empty object")
        return(new("genlight"))
    }
    sec.all <- unlist(lapply(POOL[snp.posi], function(e) e[2]))



    ## RE-READ DATA, CONVERT SNPs TO GENLIGHT ##
    if(!quiet) cat("\n Extracting SNPs from the alignment... \n")

    ## initialize
    lines.to.skip <- 0
    COUNT <- 0 # used to count the nb reads
    res <- list()

    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=LINES.PER.IND*chunkSize)

    ## returns a vector of nb of second alleles or NAs
    f1 <- function(vec){
        out <- as.integer(vec==sec.all)
        out[!vec %in% letterOK] <- NA
        return(out)
    }

    ## read and process chunks
    while(length(txt)>0){
        COUNT <- COUNT + 1
        if(!quiet) {
            for(i in 1:(COUNT*chunkSize)) cat(".")
        }


        ## read SNPs
        nb.ind <- length(grep("^>", txt))
        txt <- split(txt, rep(1:nb.ind, each=LINES.PER.IND)) # split per individuals
        if(parallel){
            txt <- mclapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split="")[[1]][snp.posi],
                                        mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) # each genome -> one SNP vector
        } else {
            txt <- lapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split="")[[1]][snp.posi]) # each genome -> one SNP vector
        }

        ## convert to genlight
        ##res <- c(res, lapply(txt, function(e) new("SNPbin", as.integer(e==sec.all))))
        res <- c(res, lapply(txt, function(e) new("SNPbin", f1(e))))

        lines.to.skip <- lines.to.skip + nb.ind*LINES.PER.IND
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, nmax=LINES.PER.IND*chunkSize)
    }



    ## BUILD FINAL OBJECT ##
    if(!quiet) cat("\n Building final object... \n")

    res <- new("genlight",res, ploidy=1, parallel=parallel)
    indNames(res) <- IND.LAB
    alleles(res) <- sapply(POOL[snp.posi], paste, collapse="/")
    position(res) <- which(snp.posi)
    if(saveNbAlleles) other(res) <- list(nb.all.per.loc=nb.alleles)


    ## RETURN OUTPUT ##
    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end fasta2genlight








###########################
## Function fasta2DNAbin
###########################
fasta2DNAbin <- function(file, quiet=FALSE, chunkSize=10, snpOnly=FALSE){

    ## HANDLE ARGUMENTS ##
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(!ext %in% c("FASTA", "FA", "FAS")) warning("wrong file extension - '.fasta', '.fa' or '.fas' expected")
    if(!quiet) cat("\n Converting FASTA alignment into a DNAbin object... \n\n")


    ## PRIOR CHECKS ##
    ## find nb of lines per genome ##

    ## find length of a single line of sequence
    if(!quiet) cat("\n Finding the size of a single genome... \n\n")
    lines.to.skip <- 0
    txt <- scan(file,what="character",sep="\n",quiet=TRUE, nmax=1)

    while(length(grep("^>.+", txt))==1){
        lines.to.skip <- lines.to.skip + 1
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, n=1)
    }

    char.per.line <- nchar(txt)
    nbOfLinesToRead <- max(round(1e6/char.per.line),1)#

    ## read file to find genome size, 1e6 characters at a time
    nblocks <- 1
    txt <- scan(file,what="character",sep="\n",quiet=TRUE, n=nbOfLinesToRead*nblocks)

    while(length(grep("^>.+", txt))<2){
        nblocks <- nblocks+1
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, n=nbOfLinesToRead*nblocks)
    }

    ## this is the nb of lines for one genome
    ## including the first line of annotation
    LINES.PER.IND <- diff(grep("^>.+", txt))[1]

    ## this is the length of a genome single
    GENOMESIZE <- sum(nchar(txt[2:LINES.PER.IND]))
    if(!quiet) cat("\n genome size is:", format(GENOMESIZE, big.mark=","), "nucleotides \n")
    if(!quiet) cat("\n(",format(LINES.PER.IND, big.mark=","), " lines per genome )\n")


    ## START READING / CONVERTING GENOMES ##
    if(!quiet) cat("\n Importing sequences... \n")

    ## read all genomes by chunks
    ## initialize
    lines.to.skip <- 0
    IND.LAB <- NULL
    COUNT <- 0 # used to count the nb reads

    txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, n=LINES.PER.IND*chunkSize)
    res <- raw()

    ## read and process chunks
    while(length(txt)>0){
        ## clean memory ##
        invisible(gc())

        ## progression... ##
        COUNT <- COUNT + 1
        if(!quiet) {
            for(i in 1:(COUNT*chunkSize)) cat(".")
        }

        ## process txt ##
        nb.ind <- length(grep("^>", txt))
        IND.LAB <- c(IND.LAB, sub(">","",txt[grep("^>", txt)])) # find individuals' labels
        txt <- split(txt, rep(1:nb.ind, each=LINES.PER.IND)) # split per individuals
        txt <- lapply(txt, function(e) unlist(strsplit(tolower(paste(e[-1], collapse="")), split=""))) # each genome -> one vector

        ## convert character vectors to DNAbin output
        res <- c(res, unlist(lapply(txt, as.DNAbin)))

        ## ## POOL contains all alleles of each position
        ## temp <- as.list(apply(matrix(unlist(txt), byrow=TRUE, nrow=length(txt)),2,unique)) # alleles current genomes
        ## POOL <- mapply(function(x,y) unique(c(x,y)), POOL, temp, SIMPLIFY=FALSE) # update global pool

        ## scan file further ##
        lines.to.skip <- lines.to.skip + nb.ind*LINES.PER.IND
        txt <- scan(file,what="character",sep="\n",quiet=TRUE, skip=lines.to.skip, n=LINES.PER.IND*chunkSize)
    }

    ## GET FINAL OBJECT ##
    if(!quiet) cat("\n Forming final object... \n")

    ## form matrix ##
    res <- matrix(res, nrow=length(IND.LAB), byrow=TRUE)
    class(res) <- "DNAbin"
    rownames(res) <- IND.LAB

    ## extract snps if needed ##
    if(snpOnly){
        if(!quiet) cat("\n Extracting SNPs... \n")
        snp.posi <- seg.sites(res)
        if(length(snp.posi)==0) warning("no polymorphic site in the sequences")
        res <- res[,seg.sites(res),drop=FALSE]
        colnames(res) <- snp.posi
    }

    ## RETURN OUTPUT ##
    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end fasta2DNAbin
