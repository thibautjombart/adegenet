

###############
##
##   CLASSES
##
###############

###############
## SNPbin class
###############
setClass("SNPbin", representation(snp = "list",
                                  n.loc = "integer",
                                  NA.posi = "integer",
                                  label = "charOrNULL",
                                  ploidy = "integer"),
         prototype(snp = list(), n.loc = 0L, label = NULL, ploidy = 1L))




###############
## genlight class
###############
setClass("genlight", representation(gen = "list",
                                    n.loc = "integer",
                                    ind.names = "charOrNULL",
                                    loc.names = "charOrNULL",
                                    loc.all = "charOrNULL",
                                    chromosome = "factorOrNULL",
                                    position = "intOrNULL",
                                    ploidy = "intOrNULL",
                                    pop = "factorOrNULL",
                                    other = "list"),
         prototype(gen = list(), n.loc = 0L, ind.names = NULL, loc.names = NULL, loc.all = NULL,
                   chromosome = NULL, position = NULL, ploidy=NULL, pop=NULL, other=list()))








#####################
##
##   CONSTRUCTORS
##
#####################

####################
## SNPbin constructor
####################
setMethod("initialize", "SNPbin", function(.Object, ...) {
    x <- .Object
    input <- list(...)
    if(length(input)==1) names(input) <- "snp"
    if(length(input)>1 && ! "snp" %in% names(input)) names(input)[1] <- "snp"

    ## handle snp data ##
    if(!is.null(input$snp) && length(input$snp)>0 && !all(is.na(input$snp))){
        ## a vector of raw is provided
        if(is.raw(input$snp)){
            x@snp <-list(input$snp)
        }

        ## a list of raw vectors is provided
        if(is.list(input$snp)){
            if(all(sapply(input$snp, class)=="raw")){
                x@snp <- input$snp
            }
        }

        ## a numeric/integer vector is provided
        ## conversion from a vector of 0/1 (integers)
        if(is.numeric(input$snp) | is.integer(input$snp)){
            input$snp <- as.integer(input$snp)
            ## determine ploidy
            if(is.null(input$ploidy)){
                input$ploidy <- max(input$snp, na.rm=TRUE)
                if(input$ploidy==0) input$ploidy <- 1
            }
            input$ploidy <- as.integer(input$ploidy)
            if(input$ploidy<1) stop("Ploidy is less than 1")

            ## check values in the vector
            if(any(input$snp<0 | input$snp>input$ploidy, na.rm=TRUE)){
                stop("Values of SNPs < 0 or > ploidy")
            }


            ## handle ploidy (may have to split info into binary vectors)
            x@snp <- list()
            i <- max(input$snp, na.rm=TRUE) # determine max nb of alleles
            if(i > 1){ # haplotype can be 0/1/2/...
                j <- 0 # index for the length of the list @snp
                while(i > 0){
                    j <- j+1 # update length of the result
                    temp <- as.integer(input$snp==i)
                    x@snp[[j]] <- .bin2raw(temp)$snp # make a vector of 1
                    input$snp <- input$snp-temp # deflate data (subtract the recoded alleles)
                    i <- max(input$snp, na.rm=TRUE) # update the max nb of alleles
                }
            } else { # haplotype is only 0/1/NA
                x@snp[[1]] <- .bin2raw(input$snp)$snp
            }
            x@n.loc <- length(input$snp)
            x@NA.posi <- which(is.na(input$snp))
            x@ploidy <- input$ploidy
            return(x)
        }
    }

    ## handle full-NA data
    if(!is.null(input$snp) && all(is.na(input$snp))){
        x@snp <- list()
        x@n.loc <- length(input$snp)
        x@snp[[1]] <- .bin2raw(rep(0L, length(input$snp)))$snp
        x@NA.posi <- 1:length(input$snp)
        if(!is.null(input$ploidy)){
            x@ploidy <- input$ploidy
        } else {
            x@ploidy <- as.integer(NA)
        }
        return(x)
    }


    ## handle n.loc ##
    if(!is.null(input$n.loc)){
        x@n.loc <- as.integer(input$n.loc)
    } else {
        if(!is.null(input$snp)){
            warning("number of SNPs (n.loc) not provided to the genlight constructor - using the maximum number given data coding.")
            x@n.loc <- as.integer(length(x@snp)*8)
        } else {
            x@n.loc <- 0L
        }
    }


    ## handle NA.posi ##
    if(!is.null(input$NA.posi)){
        x@NA.posi <- as.integer(input$NA.posi)
    }


    ## handle ploidy ##
    if(!is.null(input$ploidy)){
        x@ploidy <- as.integer(input$ploidy)
    }


    return(x)
}) # end SNPbin constructor







########################
## genlight constructor
########################
setMethod("initialize", "genlight", function(.Object, ..., parallel=require("parallel"), n.cores=NULL) {
    if(parallel && !require(parallel)) stop("parallel package requested but not installed")
    if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
    }

    x <- .Object
    input <- list(...)
    if(length(input)==1 && is.null(names(input))) names(input) <- "gen"
    if(length(input)>1 && ! "gen" %in% names(input)) names(input)[1] <- "gen"


    ## HANDLE INPUT$GEN ##
    if(!is.null(input$gen)){
        ## input$gen is a list of SNPbin ##
        if(is.list(input$gen) && all(sapply(input$gen, class)=="SNPbin")){
            ## check nb of loci in each SNPbin
            if(length(unique(sapply(input$gen, nLoc)))>1) {
                warning("SNPbin objects have different numbers of loci")
                input$gen <- lapply(input$gen, as.integer)
            } else { # all seems fine
                x@gen <- input$gen
                if(is.null(input$ind.names)){
                    input$ind.names <- names(input$gen)
                }
            }
        }


        ## input$gen is a matrix or a data.frame
        if((is.matrix(input$gen) & !inherits(input$gen,"snp.matrix")) | is.data.frame(input$gen)){
            if(is.null(input$ind.names)){
                input$ind.names <- rownames(input$gen)
            }
            if(is.null(input$loc.names)){
                input$loc.names <- colnames(input$gen)
                if(is.data.frame(input$gen)){ # do not use names if these are the default names of a data.frame
                    if(identical(colnames(input$gen), paste("V", 1:ncol(input$gen), sep=""))){
                        input$loc.names <- NULL
                    }
                }
            }
            ##input$gen <- lapply(1:nrow(input$gen), function(i) as.integer(input$gen[i,]))
            if(parallel){
                x@gen <- mclapply(1:nrow(input$gen), function(i) new("SNPbin", as.integer(input$gen[i,])),
                                  mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
            } else {
                x@gen <- lapply(1:nrow(input$gen), function(i) new("SNPbin", as.integer(input$gen[i,])) )
            }
        }


        ## input$gen is a list of integers/numeric ##
        if(is.list(input$gen) && !is.data.frame(input$gen) && all(sapply(input$gen, class) %in% c("integer","numeric"))){
            ## check length consistency
            lengthvec <- sapply(input$gen, length)

            ## complete with NA is necessary
            if(length(unique(lengthvec))>1) {
                warning("Genotypes have variable length; completing shorter ones with NAs.")
                for(i in 1:length(input$gen)){
                    input$gen[[i]] <- c(input$gen[[i]], rep(NA, max(lengthvec)-length(input$gen[[i]])))
                }
            }

            ## name individuals if needed
            if(is.null(input$ind.names)){
                input$ind.names <- names(input$gen)
            }

            ## create SNPbin list
            if(parallel){
                x@gen <- mclapply(input$gen, function(e) new("SNPbin",e), mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
            } else {
                x@gen <- lapply(input$gen, function(e) new("SNPbin",e))
            }
        }


        ## input$gen is a snp.matrix object ##
        if(inherits(input$gen,"snp.matrix")){
            if(!require(snpMatrix)){
                cat("\nThe package snp.matrix is needed for this conversion.")
                cat("\nTo install it, type:")
                cat("\n  source(\"http://bioconductor.org/biocLite.R\")")
                cat("\n  biocLite(\"snpMatrix\")\n")
                x@gen <- NULL
            } else {

                ## function to convert one indiv
                f1 <- function(x){
                    res <- as.integer(x)
                    res[res==0] <- NA
                    res <- res-1
                    return(new("SNPbin", as.integer(res), ploidy=2))
                }

                ## create SNPbin list
                if(parallel){
                    x@gen <- mclapply(1:nrow(input$gen), function(i) f1(input$gen[i,]), mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE)
                } else {
                    x@gen <- lapply(1:nrow(input$gen), function(i) f1(input$gen[i,]))
                }

                ## handle names
                if(is.null(input$ind.names)) {input$ind.names <- rownames(input$gen)}
                if(is.null(input$loc.names)) {input$loc.names <- colnames(input$gen)}
            }
        }
    }


    if(length(x@gen) > 0) { # if non-emtpy object
        ## HANDLE INPUT$IND.NAMES ##
        if(!is.null(input$ind.names)){
            input$ind.names <- as.character(input$ind.names)

            ## check length consistency
            if(length(input$ind.names) != nInd(x)){
                warning("Inconsistent length for ind.names - storing this argument in @other.")
                if(is.null(input$other)) {
                    input$other <- list(ind.names.wrong.length=input$ind.names)
                } else {
                    input$other$ind.names.wrong.length <- input$ind.names
                }

            } else {
                ## assign value to the output object
                x@ind.names <- input$ind.names
                ## ## name list and each SNPbin ## THIS DUPLICATES THE INFORMATION
                ## names(x@gen) <- input$ind.names
                ## for(i in 1:length(x@gen)){
                ##     x@gen[[i]]@label <- input$ind.names[i]
                ## }

            }


        }


        ## HANDLE INPUT$N.LOC ##
        if(!is.null(input$n.loc)){ # n.loc is provided
            input$n.loc <- as.integer(input$n.loc)

            ## check length consistency
            if(input$n.loc != nLoc(x@gen[[1]])) {
                warning("Inconsistent number of loci (n.loc) - ignoring this argument.")
            } else {
                x@n.loc <- input$n.loc
            }
        } else { # n.loc is not provided
            x@n.loc <- nLoc(x@gen[[1]])
        }


        ## HANDLE INPUT$PLOIDY ##
        ## note: if not provided, @ploidy is NULL (saves some space)
        if(!is.null(input$ploidy)){ # ploidy is provided
            input$ploidy <- as.integer(input$ploidy)
            input$ploidy <- rep(input$ploidy, length=length(x@gen))
            x@ploidy <- input$ploidy
        }


        ## HANDLE INPUT$LOC.NAMES ##
        if(!is.null(input$loc.names) && length(input$loc.names)>0){ # ploidy is provided
            input$loc.names <- as.character(input$loc.names)

            ## check length consistency
            if(length(input$loc.names) != x@n.loc){ # if problem, store in @other
                warning("Inconsistent length for loc.names - storing this argument in @other.")
                if(is.null(input$other)) {
                    input$other <- list(loc.names.wrong.length=input$loc.names)
                } else {
                    input$other$loc.names.wrong.length <- input$loc.names
                }
            } else {
                x@loc.names <- input$loc.names
            }
        }


        ## HANDLE INPUT$LOC.ALL ##
        if(!is.null(input$loc.all) && length(input$loc.all)>0){ # ploidy is provided
            input$loc.all <- as.character(input$loc.all)

            ## check length consistency
            if(length(input$loc.all) != x@n.loc){
                warning("Inconsistent length for loc.all - storing this argument in @other.")
                if(is.null(input$other)) {
                    input$other <- list(loc.all.wrong.length=input$loc.all)
                } else {
                    input$other$loc.all.wrong.length <- input$loc.all
                }
            } else {
                ## check string consistency (format is e.g. "a/t")
                temp <- grep("^[[:alpha:]]{1}/[[:alpha:]]{1}$", input$loc.all)
                if(any(! 1:nLoc(x@gen[[1]]) %in% temp)){ # if problem, store in @other
                    ## input$loc.all <- gsub("[[:space:]]","", input$loc.all)
                    warning("Miss-formed strings in loc.all (must be e.g. 'c/g') - storing this argument in @other.")
                    if(is.null(input$other)) {
                        input$other <- list(loc.all.misformed=input$loc.all)
                    } else {
                        input$other$loc.all.misformed <- input$loc.all
                    }
                } else {
                    x@loc.all <- input$loc.all
                }
            }
        }


        ## HANDLE CHROMOSOME ##
        if(!is.null(input$chromosome)){
            if(length(input$chromosome) != x@n.loc) { # if wrong length, store in @other
                warning("chromosome argument has inconsistent length - storing this argument in @other")
                if(is.null(input$other)) {
                    input$other <- list(chromosome.wrong.length=input$chromosome)
                } else {
                    input$other$chromosome.wrong.length <- input$chromosome
                }
            } else {
                x@chromosome <- factor(input$chromosome)
            }
        }


        ## HANDLE POSITION ##
        if(!is.null(input$position)){
            if(length(input$position) != x@n.loc) { # if wrong length, store in @other
                warning("position argument has inconsistent length - storing this argument in @other")
                if(is.null(input$other)) {
                    input$other <- list(position.wrong.length=input$position)
                } else {
                    input$other$position.wrong.length <- input$position
                }
            } else {
                x@position <- as.integer(input$position)
            }
        }


        ## HANDLE INPUT$POP ##
        if(!is.null(input$pop)){
            ## check length consistency
            if(length(input$pop) != nInd(x)){
                warning("Inconsistent length for pop - ignoring this argument.")
                if(is.null(input$other)) {
                    input$other <- list(pop.wrong.length=input$pop)
                } else {
                    input$other$pop.wrong.length <- input$pop
                }
            } else {
                x@pop <- factor(input$pop)
            }
        }


    } # end if non-empty @gen


    ## HANDLE INPUT$OTHER ##
    if(!is.null(input$other)){
        x@other <- input$other
    }


    ## RETURN OBJECT ##
    names(x@gen) <- NULL # do not store ind.names twice
    return(x)
}) # end genlight constructor










################################
##
##   METHODS AND ACCESSORS
##
################################

###############
## show SNPbin
###############
setMethod ("show", "SNPbin", function(object){
    cat(" === S4 class SNPbin ===")
    if(!is.null(object@label)) {
        cat("\n", object@label)
    }
    cat("\n", nLoc(object), "SNPs coded as bits")
    cat("\n Ploidy:", object@ploidy)
    temp <- round(length(object@NA.posi)/nLoc(object) *100,2)
    cat("\n ", length(object@NA.posi), " (", temp," %) missing data\n", sep="")
}) # end show method





###############
## show genlight
###############
setMethod ("show", "genlight", function(object){
    cat(" === S4 class genlight ===")
    cat("\n", nInd(object), "genotypes, ", nLoc(object),  "binary SNPs")
    temp <- unique(ploidy(object))
    if(!is.null(temp)){
        if(length(temp)==1){
            cat("\n Ploidy:", temp)
        } else {
            temp <- summary(ploidy(object))
            cat("\n Ploidy statistics (min/median/max):", temp[1], "/", temp[3], "/", temp[6])
        }
    }
    temp <- sapply(object@gen, function(e) length(e@NA.posi))
    if(length(temp>1)){
        cat("\n ", sum(temp), " (", round(sum(temp)/(nInd(object)*nLoc(object)),2)," %) missing data", sep="")
    }

    if(!is.null(pop(object))){
        cat("\n @pop: individual membership for", length(levels(pop(object))), "populations")
    }

    if(!is.null(chr(object))){
        cat("\n @chromosome: chromosome of the SNPs")
    }

    if(!is.null(position(object))){
        cat("\n @position: position of the SNPs")
    }

    if(!is.null(alleles(object))){
        cat("\n @alleles: alleles of the SNPs")
    }

    if(!is.null(object@loc.names)){
        cat("\n @loc.names: labels of the SNPs")
    }

    if(!is.null(other(object))){
        cat("\n @other: ")
        cat("a list containing: ")
        cat(ifelse(is.null(names(other(object))), paste(length(other(object)),"unnamed elements"),
                   paste(names(other(object)), collapse= "  ")), "\n")
    }

    cat("\n")
}) # end show method






############
## accessors
############

## nLoc
setMethod("nLoc","SNPbin", function(x,...){
    return(x@n.loc)
})

setMethod("nLoc","genlight", function(x,...){
    return(x@n.loc)
})


## nInd
setMethod("nInd","genlight", function(x,...){
    return(length(x@gen))
})


## $
setMethod("$","SNPbin",function(x,name) {
    return(slot(x,name))
})

setMethod("$","genlight",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","SNPbin",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})

setMethod("$<-","genlight",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})


## names
setMethod("names", signature(x = "SNPbin"), function(x){
    return(slotNames(x))
})

setMethod("names", signature(x = "genlight"), function(x){
    return(slotNames(x))
})


## ploidy
setMethod("ploidy","SNPbin", function(x,...){
    return(x@ploidy)
})

setMethod("ploidy","genlight", function(x,...){
    if(nInd(x)>0){
        if(!is.null(x@ploidy)){
            res <- x@ploidy
        } else {
            res <- sapply(x@gen, function(e) e@ploidy)
        }
        names(res) <- x@ind.names
        return(res)
    } else {
        return(NULL)
    }
})


setReplaceMethod("ploidy","SNPbin",function(x,value) {
    if(is.null(value)){
        slot(x, "ploidy", check=TRUE) <- value
        return(x)
    }
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value)>1) warning("Several ploidy numbers provided; using only the first integer")
    slot(x,"ploidy",check=TRUE) <- value[1]
    return(x)
})

setReplaceMethod("ploidy","genlight",function(x,value) {
    if(is.null(value)){
        slot(x, "ploidy", check=TRUE) <- value
        return(x)
    }
    value <- as.integer(value)
    if(any(value)<1) stop("Negative or null values provided")
    if(any(is.na(value))) stop("NA values provided")
    if(length(value) == 1) value <- rep(value, length=nInd(x))
    if(length(value) != nInd(x)) stop("Length of the provided vector does not match nInd(x)")
    slot(x,"ploidy",check=TRUE) <- value
    return(x)
})


## locNames
setMethod("locNames","genlight", function(x,...){
    ## if loc names provided, return them
    if(!is.null(x@loc.names)) return(x@loc.names)

    ## otherwise, look for position / alleles
    if(!is.null(res <- position(x))){
        if(!is.null(alleles(x))){
            res <- paste(res, alleles(x), sep=".")
        }
        return(res)
    }
})


setReplaceMethod("locNames","genlight",function(x,value) {
    if(is.null(value)){
        slot(x, "loc.names", check=TRUE) <- value
        return(x)
    }
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    slot(x,"loc.names",check=TRUE) <- as.character(value)
    return(x)
})


## indNames
setMethod("indNames","genlight", function(x,...){
    if(length(x@ind.names)==0) return(NULL)
    return(x@ind.names)
})


setReplaceMethod("indNames","genlight",function(x,value) {
    if(is.null(value)){
        slot(x, "ind.names", check=TRUE) <- value
        return(x)
    }
    value <- as.character(value)
    if(length(value) != nInd(x)) stop("Vector length does no match number of individuals")
    slot(x,"ind.names",check=TRUE) <- value
    return(x)
})


## alleles
setMethod("alleles","genlight", function(x,...){
    if(length(x@loc.all)==0) return(NULL)
    return(x@loc.all)
})

setReplaceMethod("alleles","genlight", function(x, value){
    if(is.null(value)){
        slot(x, "loc.all", check=TRUE) <- value
        return(x)
    }
    value <- as.character(value)
    if(length(value)!=nLoc(x)) stop("replacement vector must be of length nLoc(x)")
    temp <- grep("^[[:alpha:]]{1}/[[:alpha:]]{1}$", value)
    if(any(! 1:nLoc(x) %in% temp)) stop("Miss-formed strings in replacement (must be e.g. 'c/g')")
    slot(x, "loc.all", check=TRUE) <- value
    return(x)
})


## chromosome
setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
setGeneric("chromosome<-", function(x, value) standardGeneric("chromosome<-"))
setGeneric("chr", function(x,...) standardGeneric("chr"))
setGeneric("chr<-", function(x, value) standardGeneric("chr<-"))

setMethod("chromosome","genlight", function(x,...){
    if(length(x@chromosome)==0) return(NULL)
    return(x@chromosome)
})

setMethod("chr","genlight", function(x,...){
    return(chromosome(x))
})


setReplaceMethod("chromosome","genlight",function(x,value) {
    if(is.null(value)){
        slot(x, "chromosome", check=TRUE) <- value
        return(x)
    }
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    slot(x,"chromosome",check=TRUE) <- factor(value)
    return(x)
})

setReplaceMethod("chr","genlight",function(x,value) {
    chromosome(x) <- value
    return(x)
})



## position
setGeneric("position", function(x, ...) standardGeneric("position"))
setGeneric("position<-", function(x, value) standardGeneric("position<-"))

setMethod("position","genlight", function(x,...){
    if(length(x@position)==0) return(NULL)
    return(x@position)
})


setReplaceMethod("position","genlight",function(x,value) {
    if(is.null(value)){
        slot(x, "position", check=TRUE) <- value
        return(x)
    }
    if(length(value) != nLoc(x)) stop("Vector length does no match number of loci")
    slot(x,"position",check=TRUE) <- as.integer(value)
    return(x)
})



## NA.posi
setGeneric("NA.posi", function(x, ...) standardGeneric("NA.posi"))

setMethod("NA.posi","SNPbin", function(x,...){
    return(x@NA.posi)
})

setMethod("NA.posi","genlight", function(x,...){
    res <- lapply(x@gen, function(e) e@NA.posi)
    names(res) <- indNames(x)
    return(res)
})


## pop
setMethod("pop","genlight", function(x){
    return(x@pop)
})


setReplaceMethod("pop","genlight",function(x,value) {
    if(is.null(value) | length(value)==0){
        slot(x, "pop", check=TRUE) <- NULL
        return(x)
    }
    if(length(value) != nInd(x)) stop("Vector length does no match number of individuals")
    slot(x,"pop", check=TRUE) <- factor(value)
    return(x)
})



## other
setMethod("other","genlight", function(x,...){
    if(length(x@other)==0) return(NULL)
    return(x@other)
})


setReplaceMethod("other","genlight",function(x,value) {
    if( !is.null(value) && (!is.list(value) | is.data.frame(value)) ) {
        value <- list(value)
    }
    slot(x,"other",check=TRUE) <- value
    return(x)
})








###################
##
##   CONVERSIONS
##
###################

############
## .bin2raw
###########
## each byte takes a value on [0,255]

## function to code multiple SNPs on a byte
## 8 combinations of SNPs can be coded onto a single byte (0->255)
.bin2raw <- function(vecSnp){
    ## handle missing data
    NAposi <- which(is.na(vecSnp))
    if(length(NAposi)>0){
        vecSnp[is.na(vecSnp)] <- 0L
    }


    nbBytes <- length(vecSnp) %/% 8
    if(length(vecSnp) %% 8 > 0) {nbBytes <- nbBytes +1}
    ori.length <- length(vecSnp)
    new.length <- 8*nbBytes
    vecSnp <- c(vecSnp, rep(0, new.length-ori.length)) # fill the end with 0 of necessary


    ## map info to bytes (0:255)
    vecSnp <- as.integer(vecSnp)
    ##vecRaw <- integer(nbBytes) # no longer needed - sending raw type directly
    vecRaw <- raw(nbBytes)

    vecRaw <- .C("binIntToBytes", vecSnp, length(vecSnp), vecRaw, nbBytes, PACKAGE="adegenet")[[3]]
    ## vecraw <- sapply(seq(1, by=8, length=nbBytes), function(i) which(apply(SNPCOMB,1, function(e) all(temp[i:(i+7)]==e))) ) # old R version

    ## return result
    res <- list(snp=vecRaw, n.loc=as.integer(ori.length), NA.posi=as.integer(NAposi))
    return(res)
} # end .bin2raw






###########
## .raw2bin
###########
## convert vector of raw to 0/1 integers
.raw2bin <- function(x){
    if(!is.raw(x)) stop("x is not of class raw")
    ## SNPCOMB <- as.matrix(expand.grid(rep(list(c(0,1)), 8)))
    ## colnames(SNPCOMB) <- NULL
    ## res <- unlist(lapply(as.integer(x), function(i) SNPCOMB[i+1,]))
    res <- .C("bytesToBinInt", x, length(x), integer(length(x)*8), PACKAGE="adegenet")[[3]]
    return(res)
} # end .raw2bin





#############
## .SNPbin2int
#############
## convert SNPbin to integers (0/1/2...)
.SNPbin2int <- function(x){
    ##res <- lapply(x@snp, .raw2bin)
    resSize <- length(x@snp[[1]])*8
    res <- .C("bytesToInt", unlist(x@snp), length(x@snp[[1]]), length(x@snp), integer(resSize), as.integer(resSize), PACKAGE="adegenet")[[4]][1:nLoc(x)]
    ##res <- lapply(res, function(e) e[1:x@n.loc])
    ##res <- as.integer(Reduce("+", res))
    if(length(x@NA.posi)>0){
        res[x@NA.posi] <- NA
    }
    return(res)
} # end .SNPbin2int






#############
## as methods
#############
## KLUDGE - needed for as.matrix.genlight to be dispatched correctly (R-2.12.1)
setGeneric("as.matrix")


## SNPbin/genlight -> other
setAs("SNPbin", "integer", def=function(from){
    res <- .SNPbin2int(from)
    return(res)
})


as.integer.SNPbin <- function(x, ...){
    return(as(x, "integer"))
}


setAs("genlight", "matrix", def=function(from){
    res <- unlist(lapply(from@gen, as.integer))
    res <- matrix(res, ncol=nLoc(from), nrow=nInd(from), byrow=TRUE)
    colnames(res) <- locNames(from)
    ## if(!is.null(alleles(from))){
    ##     colnames(res) <- paste(locNames(from), alleles(from), sep=ifelse(is.null(locNames(from)), "", "."))
    ## }
    rownames(res) <- indNames(from)
    return(res)
})


as.matrix.genlight <- function(x, ...){
    return(as(x, "matrix"))
}


setAs("genlight", "data.frame", def=function(from){
    return(as.data.frame(as.matrix(from)))
})


as.data.frame.genlight <- function(x, ...){
    return(as(x, "data.frame"))
}


setAs("genlight", "list", def=function(from){
    res <- lapply(from@gen, as.integer)
    names(res) <- indNames(from)
    return(res)
})


as.list.genlight <- function(x, ...){
    return(as(x, "list"))
}




## other -> SNPbin/genlight
setGeneric("as.SNPbin", function(x, ...) standardGeneric("as.SNPbin"))
setGeneric("as.genlight", function(x, ...) standardGeneric("as.genlight"))

setAs("integer", "SNPbin", def=function(from){
    res <- new("SNPbin", from)
    return(res)
})

setAs("numeric", "SNPbin", def=function(from){
    res <- new("SNPbin", from)
    return(res)
})


setMethod("as.SNPbin", "integer", function(x, ...) as(x, "SNPbin"))
setMethod("as.SNPbin", "numeric", function(x, ...) as(x, "SNPbin"))


setAs("matrix", "genlight", def=function(from){
    return(new("genlight", from))
})


setAs("data.frame", "genlight", def=function(from){
    return(new("genlight", from))
})


setAs("list", "genlight", def=function(from){
    return(new("genlight", from))
})


## setAs("snp.matrix", "genlight", def=function(from){
##     return(new("genlight", from))
## })




setMethod("as.genlight", "matrix", function(x, ...) as(x, "genlight"))
setMethod("as.genlight", "data.frame", function(x, ...) as(x, "genlight"))
setMethod("as.genlight", "list", function(x, ...) as(x, "genlight"))
## setMethod("as.genlight", "snp.matrix", function(x, ...) as(x, "genlight"))






################################
## testing SNPbin
##
##
## library(adegenet)
## dat <- c(1,0,0,1,0,NA,1,0,0,0,0,1)
## x <- new("SNPbin",dat)
## as.integer(x)


## HAPLOID DATA - NO NA
## dat <- sample(c(0L,1L), 1e6, replace=TRUE)
## x <- new("SNPbin", dat)
## identical(as(x, "integer"),dat) # SHOULD NORMALLY BE TRUE
## all(as(x, "integer") == dat, na.rm=TRUE) # MUST BE TRUE
## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION


## HAPLOID DATA - WITH NAs
## dat <- sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE)
## x <- new("SNPbin", dat)
## identical(as(x, "integer"),dat) # SHOULD NORMALLY BE TRUE
## all(as(x, "integer") == dat, na.rm=TRUE) # MUST BE TRUE
## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION



 ## DIPLOID DATA
## dat <- sample(c(0:2,NA), 1e6, prob=c(.4, .4, .195 ,.005), replace=TRUE)
## x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION



 ## POLYLOID DATA
## dat <- sample(c(0:5,NA), 1e6, prob=c(rep(.995/6,6), 0.005), replace=TRUE)
##  x <- new("SNPbin", dat)

## identical(as(x, "integer"),dat) # MUST BE TRUE

## object.size(dat)/object.size(x) # EFFICIENCY OF CONVERSION




################################
## testing genlight
##
##


## SIMPLE TESTS
## dat <- c(1,0,0,1,0,0,1,1,1,0,1)
## x <- new("SNPbin",dat)$snp[[1]]
## as.integer(x)==dat


## library(adegenet)
## dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
## x <- new("genlight", dat)
## x
## as.list(x)
## as.matrix(x)

## identical(x, new("genlight", as.list(x))) # round trip - list - MUST BE TRUE
## identical(x, new("genlight", as.matrix(x))) # round trip - matrix - MUST BE TRUE
## identical(x, new("genlight", as.data.frame(x))) # round trip - data.frame - MUST BE TRUE

## ## test subsetting
## identical(as.list(x[c(1,3)]), as.list(x)[c(1,3)]) # MUST BE TRUE
## identical(x, x[]) # MUST BE TRUE
## all.equal(t(as.matrix(as.data.frame(dat)))[,1:3], as.matrix(x[,1:3])) # MUST BE TRUE




## ## ## BIG SCALE TEST - HAPLOID DATA WITH NA
## library(adegenet)
## dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE))
## names(dat) <- paste("indiv", 1:length(dat))
## print(object.size(dat), unit="aut")

## system.time(x <- new("genlight", dat)) # conversion + time taken
## print(object.size(x), unit="au")
## object.size(dat)/object.size(x) # conversion efficiency


## ## time taken by subsetting (quite long, +- 35sec)
## system.time(y <- x[1:10, 1:5e5])



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
