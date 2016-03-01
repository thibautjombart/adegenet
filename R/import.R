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
## Revised March 2015
## t.jombart@imperial.ac.uk
##
##################################################################


######################
## Function df2genind
######################
#' Convert a data.frame of allele data to a genind object.
#' 
#' The function \code{df2genind} converts a data.frame (or a matrix) into a 
#' \linkS4class{genind} object. The data.frame must meet the following 
#' requirements:
#' \itemize{
#' \item genotypes are in row (one row per genotype)
#' \item markers/loci are in columns
#' \item each element is a string of characters coding alleles, ideally
#' separated by a character string (argument \code{sep}); if no separator is
#' used, the number of characters coding alleles must be indicated (argument
#' \code{ncode}).}
#' 
#' See \code{\link{genind2df}} to convert \linkS4class{genind} objects back to
#' such a data.frame.
#' 
#' === Details for the \code{sep} argument ===\cr this character is directly 
#' used in reguar expressions like \code{gsub}, and thus require some characters
#' to be preceeded by double backslashes. For instance, "/" works but "|" must
#' be coded as "\\|".
#' 
#' @aliases df2genind
#' @param X a matrix or a data.frame containing allelle data only (see 
#'   decription)
#' @param sep a character string separating alleles. See details.
#' @param ncode an optional integer giving the number of characters used for 
#'   coding one genotype at one locus. If not provided, this is determined from 
#'   data.
#' @param ind.names optinal, a vector giving the individuals names; if NULL,
#'   taken from rownames of X. If factor or numeric, vector is converted to
#'   character.
#' @param loc.names an optional character vector giving the markers names; if 
#'   NULL, taken from colnames of X.
#' @param pop an optional factor giving the population of each individual.
#' @param NA.char a character string corresponding to missing allele (to be
#'   treated as NA)
#' @param ploidy an integer indicating the degree of ploidy of the genotypes.
#' @param type a character string indicating the type of marker: 'codom' stands 
#'   for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for 
#'   'presence/absence' markers (e.g. AFLP, RAPD).
#' @param strata an optional data frame that defines population stratifications 
#'   for your samples. This is especially useful if you have a hierarchical or 
#'   factorial sampling design.
#' @param hierarchy a hierarchical formula that explicitely defines hierarchical
#'   levels in your strata. see \code{\link{hierarchy}} for details.
#'   
#' @return an object of the class \linkS4class{genind} for \code{df2genind}; a 
#'   matrix of biallelic genotypes for \code{genind2df}
#'   
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Zhian N. Kamvar 
#'   \email{kamvarz@@science.oregonstate.edu}
#'   
#' @seealso \code{\link{genind2df}}, \code{\link{import2genind}}, 
#'   \code{\link{read.genetix}}, \code{\link{read.fstat}}, 
#'   \code{\link{read.structure}}
#'   
#' @keywords manip
#' @examples
#' 
#' ## simple example
#' df <- data.frame(locusA=c("11","11","12","32"),
#' locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
#' row.names(df) <- .genlab("genotype",4)
#' df
#' 
#' obj <- df2genind(df, ploidy=2, ncode=1)
#' obj
#' tab(obj)
#' 
#' 
#' ## converting a genind as data.frame
#' genind2df(obj)
#' genind2df(obj, sep="/")
#' 
#' @export
#' 
df2genind <- function(X, sep=NULL, ncode=NULL, ind.names=NULL, loc.names=NULL,
                      pop=NULL, NA.char="", ploidy=2, type=c("codom","PA"),
                      strata = NULL, hierarchy = NULL){

    ## CHECKS ##
    if(is.data.frame(X)) X <- as.matrix(X)
    if (!inherits(X, "matrix")) stop ("X is not a matrix")
    res <- list()
    type <- match.arg(type)
    if (is.null(sep) && is.null(ncode) && any(ploidy > 1)){
        stop("Not enough information to convert data: please indicate the separator (sep=...) or the number of characters coding an allele (ncode=...)")
    }
    if(length(NA.char)>1) {
        warning("NA.char has several values; only the first one will be considered")
        NA.char <- NA.char[1]
    }

    # If by any chance provided ind.names are of class int/factor, they are (silently?) converted to characters.
    if (!is.character(ind.names) & !is.null(ind.names)) {
      ind.names <- as.character(ind.names)
    }


    ## TYPE-INDEPENDENT STUFF ##
    ## misc variables
    n <- nrow(X)
    nloc <- ncol(X)
    if (length(ploidy) < n){
      if (length(ploidy) == 1){
        ploidy <- rep(as.integer(ploidy), length=n)
      } else {
        undefined <- length(ploidy)/n
        msg <- paste0("\nPloidy is undefined for ",
                      undefined*100,
                      "% of data.\n",
                      "This must be a single integer indicating the ploidy of",
                      "the entire data set or vector of integers the same",
                      "length as the number of samples.")
        stop(msg)
      }

    }
    if(any(ploidy < 1L)) stop("ploidy cannot be less than 1")

    ## check individual labels
    if(is.null(ind.names)) ind.names <- rownames(X)
    if(is.null(ind.names)) ind.names <- .genlab("",n)
    if(any(duplicated(ind.names))){
        warning("duplicate labels detected for some individuals; using generic labels")
        ind.names <- .genlab("",n)
    }
    rownames(X) <- ind.names

    ## check locus labels
    if(is.null(loc.names)) loc.names <- colnames(X)
    if(is.null(loc.names)) loc.names <- .genlab("loc",nloc)
    if(any(duplicated(loc.names))){
        warning("duplicate labels detected for some loci; using generic labels")
        loc.names <- .genlab("loc",nloc)
    }
    if(length(grep("[.]", loc.names))>0L){
        warning("character '.' detected in names of loci; replacing with '_'")
        gsub("[.]","_", loc.names)
    }
    colnames(X) <- loc.names


    ## pop argument
    if(!is.null(pop)){
        if(length(pop)!= n) stop("length of factor pop differs from nrow(X)")
        pop <- as.factor(pop)
    }
    
    ## check alleles for periods
    if (length(grep("[.]", X)) > 0L){
      if (is.null(sep) || sep != "_"){
        warning("character '.' detected in names of loci; replacing with '_'")
        replacement <- "_"
      } else {
        warning("character '.' detected in names of loci; replacing with 'p'")
        replacement <- "p"
      }
      X <- apply(X, 2, function(i) gsub("[.]", replacement, i))
    }

    ## PRESENCE/ABSENCE MARKERS ##
    if(toupper(type)=="PA"){
        ## preliminary stuff
        rownames(X) <- ind.names
        colnames(X) <- loc.names

        ## Erase entierely non-typed loci
        temp <- colSums(is.na(X))==nrow(X)
        if(any(temp)){
            X <- X[,!temp]
            warning("entirely non-type marker(s) deleted")
        }

        ## Erase entierely non-type individuals
        temp <- rowSums(is.na(X))==ncol(X)
        if(any(temp)){
            X <- X[!temp,,drop=FALSE]
            if(!is.null(pop)) pop <- pop[!temp]
            ploidy <- ploidy[!temp]
            ind.names <- ind.names[!temp]
            warning("entirely non-type individual(s) deleted")
        }

        ## erase non-polymorphic loci
        temp <- apply(X, 2, function(loc) length(unique(loc[!is.na(loc)]))==1)
        if(any(temp)){
            X <- X[,!temp,drop=FALSE]
            loc.names <- loc.names[!temp]
            nloc <- ncol(X)
            warning("non-polymorphic marker(s) deleted")
        }

        prevcall <- match.call()

        res <- genind(tab=X, pop=pop, prevcall=prevcall, ploidy=ploidy,
                      type = "PA", strata = strata, hierarchy = hierarchy)

        return(res)
    } # end type PA


    ## CODOMINANT MARKERS ##
    ## make sure X is in character mode
    mode(X) <- "character"


    ## HANDLE MISSING SEPARATORS
    if(is.null(sep) && any(ploidy>1)){
        ## check that ncode is provided
        if(is.null(ncode)) stop("please indicate either the separator (sep) or the number of characters coding an allele (ncode).")

        ## add "/" as separator
        X <- gsub(paste("([[:alnum:]]{",ncode,"})",sep=""), "\\1/", X)
        X <- sub("/$","",X)
        sep <- "/"
    }

    ## HANDLE NAs
    ## find all strings which are in fact NAs
    NA.list <- unlist(lapply(unique(ploidy), function(nrep) paste(rep(NA.char, nrep), collapse=sep)))
    NA.list <- unique(c(NA.list, NA.char))

    ## replace NAs
    X[X %in% NA.list] <- NA

    ## erase entirely non-type loci
    toRemove <- which(colSums(is.na(X))==nrow(X))
    if(length(toRemove) > 0){
        X <- X[,-toRemove, drop = FALSE]
        loc.names <- loc.names[-toRemove]
        warning("entirely non-type marker(s) deleted")
    }


    ## erase entierely non-type individuals
    toRemove <- which(rowSums(is.na(X))==ncol(X))
    if(length(toRemove) > 0){
        X <- X[-toRemove, , drop = FALSE]
        ind.names <- rownames(X)
        ploidy <- ploidy[-toRemove]
        if(!is.null(pop)) pop <- pop[-toRemove]
        warning("entirely non-type individual(s) deleted")
    }


    ## TRANSLATE DATA INTO ALLELE COUNTS ##
    ## get dimensions of X
    nloc <- ncol(X)
    nind <- nrow(X)

    ## unfold data for each cell of the table
    if (any(ploidy > 1)){
        allele.data <- strsplit(X, sep)
        n.items <- sapply(allele.data, length)
        locus.data <- rep(rep(loc.names, each=nind), n.items)
        ind.data <- rep(rep(ind.names,ncol(X)), n.items)
        allele.data <- unlist(allele.data)
    } else {
        n.items     <- rep(1, length(X))
        locus.data  <- rep(rep(loc.names, each=nind), n.items)
        ind.data    <- rep(rep(ind.names, ncol(X)), n.items)
        allele.data <- unlist(X)
    }


    ## identify NAs
    NA.posi <- which(is.na(allele.data))
    NA.ind <- ind.data[NA.posi]
    NA.locus <- locus.data[NA.posi]

    ## remove NAs
    if(length(NA.posi)>0){
        allele.data <- allele.data[-NA.posi]
        locus.data <- locus.data[-NA.posi]
        ind.data <- ind.data[-NA.posi]
    }

    ## get matrix of allele counts
    allele.data <- paste(locus.data, allele.data, sep=".")
    allele.data <- factor(allele.data, levels=unique(allele.data))
    out         <- table(ind.data, allele.data)
    out         <- out[ind.names, , drop = FALSE] # table sorts alphabetically. This resets.

    ## force type 'matrix'
    class(out) <- NULL
    dimnames(out) <- list(rownames(out), colnames(out))

    ## restore NAs
    ## 
    ## Thanks to Klaus Schliep for the proposed speedup:
    ## 
    # if (length(NA.posi) > 0) {
    #     out.colnames <- colnames(out)
    #     NA.row <- match(NA.ind, rownames(out))
    #     loc <- paste0(NA.locus, "\\.")
    #     uloc <- unique(loc)
    #     loc.list <- lapply(uloc, grep, out.colnames)
    #     NA.col <- match(loc, uloc)
    #     out[cbind(rep(NA.row, unlist(lapply(loc.list, length))[NA.col]), unlist(loc.list[NA.col]))] <- NA
    #  }  
    ## This one is modified from above to make everything more explicit. 
    if (length(NA.posi) > 0) {
      out.colnames <- colnames(out)
      NA.row <- match(NA.ind, rownames(out))
      loc <- paste0(NA.locus, "\\.")
      uloc <- unique(loc)
      loc.list <- lapply(uloc, grep, out.colnames)
      NA.col <- match(loc, uloc)
      
      # Coordinates for missing rows
      missing.ind <- vapply(loc.list, length, integer(1))[NA.col]
      missing.ind <- rep(NA.row, missing.ind)
      # Coordinates for missing columns
      missing.loc <- unlist(loc.list[NA.col], use.names = FALSE)
      
      missing_coordinates <- matrix(0L, nrow = length(missing.ind), ncol = 2L)
      missing_coordinates[, 1] <- missing.ind
      missing_coordinates[, 2] <- missing.loc
      
      out[missing_coordinates] <- NA
    }


    ## call upon genind constructor
    prevcall <- match.call()
    out <- genind(tab=out, pop=pop, prevcall=prevcall, ploidy=ploidy, type=type,
                  strata = strata, hierarchy = hierarchy)

    return(out)
} # end df2genind







########################################
## Function read.genetix
## code based on previous ade4 functions
########################################
#'
#' Reading data from GENETIX
#'
#' The function \code{read.genetix} reads GENETIX data files (.gtx) and convert
#' them into a \linkS4class{genind} object.
#'
#' Note: \code{read.genetix} is meant for DIPLOID DATA ONLY. Haploid data with
#' the GENETIX format can be read into R using \code{read.table} or
#' \code{read.csv} after removing headers and 'POP' lines, and then converted
#' using \code{\link{df2genind}}.
#'
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}
#' @references Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F.
#' (1996-2004) GENETIX 4.05, logiciel sous Windows TM pour la genetique des
#' populations. Laboratoire Genome, Populations, Interactions, CNRS UMR 5000,
#' Universite de Montpellier II, Montpellier (France). \cr
#' @keywords manip
#' @examples
#'
#' obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))
#' obj
#'
#' @export read.genetix
read.genetix <- function(file=NULL,quiet=FALSE) {
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
    res <- df2genind(X=X, ncode=3, pop=pop, ploidy=2, NA.char="000")
    res@call <- match.call()

    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end read.genetix





######################
## Function read.fstat
######################
#' Reading data from Fstat
#'
#' The function \code{read.fstat} reads Fstat data files (.dat) and convert
#' them into a \linkS4class{genind} object.
#'
#' Note: \code{read.fstat} is meant for DIPLOID DATA ONLY. Haploid data with
#' the Hierfstat format can be read into R using \code{read.table} or
#' \code{read.csv} after removing headers and 'POP' lines, and then converted
#' using \code{\link{df2genind}}.
#'
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.genetix}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}
#' @references Fstat (version 2.9.3). Software by Jerome Goudet.
#' http://www2.unil.ch/popgen/softwares/fstat.htm\cr
#' @keywords manip
#' @examples
#'
#' obj <- read.fstat(system.file("files/nancycats.dat",package="adegenet"))
#' obj
#'
#' @export read.fstat
read.fstat <- function(file, quiet=FALSE){
    ##if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "DAT") stop("File extension .dat expected")

    if(!quiet) cat("\n Converting data from a FSTAT .dat file to a genind object... \n\n")

    call <- match.call()
    txt <- scan(file,what="character",sep="\n",quiet=TRUE)
    txt <- gsub("\t"," ",txt)

    ## read length of allele
    ncode <- as.integer(unlist(strsplit(txt[1], " "))[4])
    NA.char <- paste(rep("0",ncode),collapse="")

    ## read first infos
    info <- unlist(strsplit(txt[1],"([[:space:]]+)"))
    ## npop <- as.numeric(info[1]) ## no longer used
    nloc <- as.numeric(info[2])

    loc.names <- txt[2:(nloc+1)]

    ## build genotype matrix
    txt <- txt[-(1:(nloc+1))]
    txt <- .rmspaces(txt)
    txt <- sapply(1:length(txt),function(i) unlist(strsplit(txt[i],"([[:space:]]+)|([[:blank:]]+)")) )
    X <- t(txt)
    pop <- factor(X[,1])
    if(length(levels(pop)) == 1 ) pop <- NULL
    X <- X[,-1]

    colnames(X) <- loc.names
    rownames(X) <- 1:nrow(X)

    ## replace all possible missing data coding by NA.char
    allNAs <- sapply(1:8, function(i) paste(rep("0",i),collapse=""))
    X[X %in% allNAs] <- NA.char

    ## call df2genind
    res <- df2genind(X=X,pop=pop, ploidy=2, ncode=ncode, NA.char=NA.char)
    res@call <- call

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.fstat





##########################
## Function read.genepop
##########################
#' Reading data from Genepop
#'
#' The function \code{read.genepop} reads Genepop data files (.gen) and convert
#' them into a \linkS4class{genind} object.
#'
#' Note: \code{read.genepop} is meant for DIPLOID DATA ONLY. Haploid data with
#' the Genepop format can be read into R using \code{read.table} or
#' \code{read.csv} after removing headers and 'POP' lines, and then converted
#' using \code{\link{df2genind}}.
#'
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param ncode an integer indicating the number of characters used to code an allele.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genetix}}
#' @references Raymond M. & Rousset F, (1995). GENEPOP (version 1.2):
#' population genetics software for exact tests and ecumenicism. \emph{J.
#' Heredity}, \bold{86}:248-249 \cr
#' @keywords manip
#' @examples
#'
#' obj <- read.genepop(system.file("files/nancycats.gen",package="adegenet"))
#' obj
#'
#' @export read.genepop
read.genepop <- function(file, ncode=2L, quiet=FALSE){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "GEN") stop("File extension .gen expected")

    if(!quiet) cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")

    prevcall <- match.call()

    txt <- scan(file,sep="\n",what="character",quiet=TRUE)
    if(!quiet) cat("\nFile description: ",txt[1], "\n")
    txt <- txt[-1]
    txt <- gsub("\t", " ", txt)
    NA.char <- paste(rep("0",ncode), collapse="")

    ## two cases for locus names:
    ## 1) all on the same row, separated by ","
    ## 2) one per row
    ## ! spaces and tab allowed
    ## a bug was reported by S. Devillard, occuring
    ## when the two cases occur together,
    ## that is:
    ## loc1,
    ## loc2,
    ## ...


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

    # If there are any duplicate names, make them unique and issue a warning. Else
    # use existing individual names.
    if (any(duplicated(ind.names))) {
      rownames(X) <- .genlab("", nrow(X))
      warning("Duplicate individual names detected. Coercing them to be unique.")
    } else {
      rownames(X) <- ind.names
    }
    colnames(X) <- loc.names

    ## give right pop names
    ## beware: genepop takes the name of the last individual of a sample as this sample's name
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names

    ## check that data are consistent with NCODE and ploidy=2
    if(!all(unique(nchar(X))==(ncode*2))) stop(paste("some alleles are not encoded with", ncode,
                                                     "characters\nCheck 'ncode' argument"))

    res <- df2genind(X=X,pop=pop, ploidy=2, ncode=ncode, NA.char=NA.char)
    res@call <- prevcall

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.genepop





############################
## Function read.structure
############################
#' Reading data from STRUCTURE
#'
#' The function \code{read.structure} reads STRUCTURE data files (.str ou
#' .stru) and convert them into a \linkS4class{genind} object. By default, this
#' function is interactive and asks a few questions about data content. This
#' can be disabled (for optional questions) by turning the 'ask' argument to
#' FALSE. However, one has to know the number of genotypes, of markers and if
#' genotypes are coded on a single or on two rows before importing data.
#'
#' Note: \code{read.structure} is meant for DIPLOID DATA ONLY. Haploid data
#' with the STRUCTURE format can easily be read into R using \code{read.table}
#' or \code{read.csv} and then converted using \code{\link{df2genind}}.
#'
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param n.ind an integer giving the number of genotypes (or 'individuals') in
#' the dataset
#' @param n.loc an integer giving the number of markers in the dataset
#' @param onerowperind a STRUCTURE coding option: are genotypes coded on a
#' single row (TRUE), or on two rows (FALSE, default)
#' @param col.lab an integer giving the index of the column containing labels
#' of genotypes. '0' if absent.
#' @param col.pop an integer giving the index of the column containing
#' population to which genotypes belong. '0' if absent.
#' @param col.others an vector of integers giving the indexes of the columns
#' containing other informations to be read. Will be available in @@other of the
#' created object.
#' @param row.marknames an integer giving the index of the row containing the
#' names of the markers. '0' if absent.
#' @param NA.char the character string coding missing data. "-9" by default.
#' Note that in any case, series of zero (like "000") are interpreted as NA
#' too.
#' @param pop an optional factor giving the population of each individual.
#' @param sep a character string used as separator between alleles.
#' @param ask a logical specifying if the function should ask for optional
#' informations about the dataset (TRUE, default), or try to be as quiet as
#' possible (FALSE).
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.fstat}}, \code{\link{read.genetix}},
#' \code{\link{read.genepop}}
#' @references Pritchard, J.; Stephens, M. & Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. \emph{Genetics},
#' \bold{155}: 945-959
#' @keywords manip
#' @examples
#'
#' obj <- read.structure(system.file("files/nancycats.str",package="adegenet"),
#'   onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE)
#'
#' obj
#'
#' @export read.structure
read.structure <- function(file, n.ind=NULL, n.loc=NULL,  onerowperind=NULL,
                           col.lab=NULL, col.pop=NULL, col.others=NULL,
                           row.marknames=NULL, NA.char="-9", pop=NULL,
                           sep=NULL, ask=TRUE, quiet=FALSE){

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

    } else { # if onerowperind
        temp <- seq(1,p-1,by=2)
        X <- paste(gen[,temp] , gen[,temp+1], sep="/")
        X <- matrix(X, nrow=n.ind)
        sep <- "/"
    }

    ## replace missing values by NAs
    X <- gsub(NA.char,NA,X)
    rownames(X) <- ind.names
    colnames(X) <- loc.names

    res <- df2genind(X=X,pop=pop, ploidy=2,sep=sep,ncode=ncode)

    res@call <- match.call()

    if(exists("X.other")) {res@other <- list(X=X.other)}

    return(res)

}




#########################
## Function import2genind
#########################
#'
#' Importing data from several softwares to a genind object
#'
#' Their are several ways to import genotype data to a \linkS4class{genind}
#' object: i) from a data.frame with a given format (see
#' \code{\link{df2genind}}), ii) from a file with a recognized extension, or
#' iii) from an alignement of sequences (see \code{\link{DNAbin2genind}}).\cr
#'
#' The function \code{import2genind} detects the extension of the file given in
#' argument and seeks for an appropriate import function to create a
#' \code{genind} object.\cr Current recognized formats are :\cr - GENETIX files
#' (.gtx) \cr - Genepop files (.gen) \cr - Fstat files (.dat) \cr - STRUCTURE
#' files (.str or .stru) \cr
#'
#' Beware: same data in different formats are not expected to produce exactly
#' the same \code{genind} objects.\cr For instance, conversions made by GENETIX
#' to Fstat may change the the sorting of the genotypes; GENETIX stores
#' individual names whereas Fstat does not; Genepop chooses a sample's name
#' from the name of its last genotype; etc.
#'
#' @aliases import2genind
#' @param file a character string giving the path to the file to convert, with
#' the appropriate extension.
#' @param quiet logical stating whether a conversion message must be printed
#' (TRUE,default) or not (FALSE).
#' @param \dots other arguments passed to the appropriate 'read' function
#' (currently passed to \code{read.structure})
#' @return an object of the class \code{genind}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}
#' @references Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F.
#' (1996-2004) GENETIX 4.05, logiciel sous Windows TM pour la genetique des
#' populations. Laboratoire Genome, Populations, Interactions, CNRS UMR 5000,
#' Universite de Montpellier II, Montpellier (France). \cr
#'
#' Pritchard, J.; Stephens, M. & Donnelly, P. (2000) Inference of population
#' structure using multilocus genotype data. \emph{Genetics}, \bold{155}:
#' 945-959
#'
#' Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population genetics
#' software for exact tests and ecumenicism. \emph{J. Heredity},
#' \bold{86}:248-249 \cr
#'
#' Fstat (version 2.9.3). Software by Jerome Goudet.
#' http://www2.unil.ch/popgen/softwares/fstat.htm\cr
#'
#' Excoffier L. & Heckel G.(2006) Computer programs for population genetics
#' data analysis: a survival guide \emph{Nature}, \bold{7}: 745-758
#' @keywords manip
#' @examples
#'
#' import2genind(system.file("files/nancycats.gtx",
#' package="adegenet"))
#'
#' import2genind(system.file("files/nancycats.dat",
#' package="adegenet"))
#'
#' import2genind(system.file("files/nancycats.gen",
#' package="adegenet"))
#'
#' import2genind(system.file("files/nancycats.str",
#' package="adegenet"), onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE)
#'
import2genind <- function(file, quiet=FALSE, ...){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    ext <- .readExt(file)
    ext <- toupper(ext)

    if(ext == "GTX")
        return(read.genetix(file,quiet=quiet))

    if(ext == "DAT")
        return(read.fstat(file, quiet=quiet))

    if(ext == "GEN")
        return(read.genepop(file, quiet=quiet, ...))

    if(ext %in% c("STR","STRU"))
        return(read.structure(file, quiet=quiet, ...))

    ## evaluated only if extension is not supported
    cat("\n File format (",ext,") not supported.\n")
    cat("\nSupported formats are:\nGENETIX (.gtx) \nFSTAT (.dat) \nGenepop (.gen)\n \nSTRUCTURE (.str)\n")

    return(invisible())
}







#######################
## Function read.snp
#######################
#' Reading Single Nucleotide Polymorphism data
#'
#' The function \code{read.snp} reads a SNP data file with extension '.snp' and
#' converts it into a \linkS4class{genlight} object. This format is devoted to
#' handle biallelic SNP only, but can accommodate massive datasets such as
#' complete genomes with considerably less memory than other formats.
#'
#' The function reads data by chunks of a few genomes (minimum 1, no maximum)
#' at a time, which allows one to read massive datasets with negligible RAM
#' requirements (albeit at a cost of computational time). The argument
#' \code{chunkSize} indicates the number of genomes read at a time. Increasing
#' this value decreases the computational time required to read data in, while
#' increasing memory requirements.
#'
#' A description of the .snp format is provided in an example file distributed
#' with adegenet (see example below).
#'
#' === The .snp format ===
#'
#' Details of the .snp format can be found in the example file distributed with
#' adegenet (see below), or on the adegenet website (type \code{adegenetWeb()}
#' in R).
#'
#' @param file a character string giving the path to the file to convert, with
#' the extension ".snp".
#' @param quiet logical stating whether a conversion messages should be printed
#' (TRUE,default) or not (FALSE).
#' @param chunkSize an integer indicating the number of genomes to be read at a
#' time; larger values require more RAM but decrease the time needed to read
#' the data.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed (see details).
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @param \dots other arguments to be passed to other functions - currently not
#' used.
#' @return an object of the class \code{"\linkS4class{genlight}"}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{?genlight} for a description of the class
#' \code{"\linkS4class{genlight}"}.
#'
#' - \code{\link{read.PLINK}}: read SNPs in PLINK's '.raw' format.
#'
#' - \code{\link{fasta2genlight}}: extract SNPs from alignments with fasta
#' format.
#'
#' - \code{\link{df2genind}}: convert any multiallelic markers into adegenet
#' \code{"\linkS4class{genlight}"}.
#'
#' - \code{\link{import2genind}}: read multiallelic markers from various
#' software into adegenet.\cr
#' @keywords manip
#' @examples
#'
#' \dontrun{
#' ## show the example file ##
#' ## this is the path to the file:
#' system.file("files/exampleSnpDat.snp",package="adegenet")
#'
#' ## show its content:
#' file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))
#'
#'
#' ## read the file
#' obj <-
#' read.snp(system.file("files/exampleSnpDat.snp",package="adegenet"), chunk=2)
#' obj
#' as.matrix(obj)
#' ploidy(obj)
#' alleles(obj)
#' locNames(obj)
#' }
#'
#' @export read.snp
#'
read.snp <- function(file, quiet=FALSE, chunkSize=1000,
                     parallel=FALSE, n.cores=NULL, ...){
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(ext != "SNP") warning("wrong file extension - '.snp' expected")
    if(!quiet) cat("\n Reading biallelic SNP data file into a genlight object... \n\n")
    ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
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
            res <- c(res, parallel::mclapply(temp, function(e) new("SNPbin", e),
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
#' @export
#' @rdname read.PLINK
#' @aliases extract.PLINKmap
#'
#' @param x an optional object of the class \code{"\linkS4class{genlight}"}, in which
#' the information read is stored; if provided, information is matched against
#' the names of the loci in \code{x}, as returned by \code{locNames(x)}; if not
#' provided, a list of two components is returned, containing chromosome and
#' position information.
#'
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
#' Reading PLINK Single Nucleotide Polymorphism data
#'
#' The function \code{read.PLINK} reads a data file exported by the PLINK
#' software with extension '.raw' and converts it into a \code{"\linkS4class{genlight}"}
#' object. Optionally, information about SNPs can be read from a ".map" file,
#' either by specifying the argument \code{map.file} in \code{read.PLINK}, or
#' using \code{extract.PLINKmap} to add information to an existing
#' \code{"\linkS4class{genlight}"} object.
#'
#' The function reads data by chunks of several genomes (minimum 1, no maximum)
#' at a time, which allows one to read massive datasets with negligible RAM
#' requirements (albeit at a cost of computational time). The argument
#' \code{chunkSize} indicates the number of genomes read at a time. Increasing
#' this value decreases the computational time required to read data in, while
#' increasing memory requirements.
#'
#' See details for the documentation about how to export data using PLINK to
#' the '.raw' format.
#'
#' === Exporting data from PLINK ===
#'
#' Data need to be exported from PLINK using the option "--recodeA" (and NOT
#' "--recodeAD"). The PLINK command should therefore look like: \code{plink
#' --file data --recodeA}. For more information on this topic, please look at
#' this webpage: \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}
#'
#' @aliases read.PLINK read.plink
#' @param file for \code{read.PLINK} a character string giving the path to the
#' file to convert, with the extension ".raw"; for \code{extract.PLINKmap}, a
#' character string giving the path to a file with extension ".map".
#' @param map.file an optional character string indicating the path to a ".map"
#' file, which contains information about the SNPs (chromosome, position). If
#' provided, this information is processed by \code{extract.PLINKmap} and
#' stored in the \code{@@other} slot.
#' @param quiet logical stating whether a conversion messages should be printed
#' (TRUE,default) or not (FALSE).
#' @param chunkSize an integer indicating the number of genomes to be read at a
#' time; larger values require more RAM but decrease the time needed to read
#' the data.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE, default), or not (FALSE);
#' requires the package \code{parallel} to be installed (see details).
#' @param n.cores if \code{parallel} is TRUE, the number of cores to be used in
#' the computations; if NULL, then the maximum number of cores available on the
#' computer is used.
#' @param \dots other arguments to be passed to other functions - currently not
#' used.
#'
#' @return - read.PLINK: an object of the class \code{"\linkS4class{genlight}"}
#'
#' - extract.PLINKmap: if a \code{"\linkS4class{genlight}"} is provided as argument
#' \code{x}, this object incorporating the new information about SNPs in the
#' \code{@@other} slot (with new components 'chromosome' and 'position');
#' otherwise, a list with two components containing chromosome and position
#' information.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso - \code{?genlight} for a description of the class
#' \code{"\linkS4class{genlight}"}.
#'
#' - \code{\link{read.snp}}: read SNPs in adegenet's '.snp' format.
#'
#' - \code{\link{fasta2genlight}}: extract SNPs from alignments with fasta
#' format.
#'
#' - other import function in adegenet: \code{\link{import2genind}},
#' \code{\link{df2genind}}, \code{\link{read.genetix}}
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genepop}}.
#'
#' - another function \code{read.plink} is available in the package
#' \code{snpMatrix}.
#' @keywords manip
#' @export
#' @rdname read.PLINK
read.PLINK <- function(file, map.file=NULL, quiet=FALSE, chunkSize=1000,
                       parallel=require("parallel"), n.cores=NULL, ...){
    ## HANDLE ARGUMENTS ##
    ext <- .readExt(file)
    ext <- toupper(ext)
    if(ext != "RAW") warning("wrong file extension - '.raw' expected")
    if(!quiet) cat("\n Reading PLINK raw format into a genlight object... \n\n")
    ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
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
            res <- c(res, parallel::mclapply(txt, function(e) new("SNPbin", snp=e, ploidy=2L),
                                   mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE) )
        } else {
            res <- c(res, lapply(txt, function(e) new("SNPbin", snp=e, ploidy=2L)) )
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

    res <- new("genlight",res, ploidy=2L, parallel=parallel)
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
    ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
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
            txt <- parallel::mclapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split=""),
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
            txt <- parallel::mclapply(txt, function(e) strsplit(paste(e[-1], collapse=""), split="")[[1]][snp.posi],
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
