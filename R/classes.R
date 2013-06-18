########################################################################
# adegenet classes definitions. All classes are S4.
#
# Thibaut Jombart, November 2007
# t.jombart@imperial.ac.uk
########################################################################

###############################
# Two classes of R object are
# defined :
# gen - common part to genind and genpop
# genind - genotypes of individuals
# genpop - allelic frequencies of populations
###############################


###############################################################
###############################################################
# AUXILIARY FUNCTIONS
###############################################################
###############################################################




###############################################################
###############################################################
# CLASSES DEFINITION
###############################################################
###############################################################

##.initAdegenetClasses <- function(){


####################
# Unions of classes
####################
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))
setClassUnion("charOrNULL", c("character","NULL"))
setClassUnion("callOrNULL", c("call","NULL"))
setClassUnion("intOrNum", c("integer","numeric","NULL"))
setClassUnion("intOrNULL", c("integer","NULL"))


####################
# virtual class gen
####################
.gen.valid <- function(object){
  # this function tests only the consistency
  # of the length of each component
  p <- ncol(object@tab)
  k <- length(unique(object@loc.names))


  if(!is.null(object@loc.fac)){
      if(length(object@loc.fac) != p) {
          cat("\ninvalid length for loc.fac\n")
          return(FALSE)
      }

      if(length(levels(object@loc.fac)) != k) {
          cat("\ninvalid number of levels in loc.fac\n")
          return(FALSE)
      }
  }

  if(!is.null(object@loc.nall)){
      if(length(object@loc.nall) != k) {
          cat("\ninvalid length in loc.nall\n")
          return(FALSE)
      }
  }

  temp <- table(object@loc.names[object@loc.names!=""])
  if(any(temp>1)) {
      warning("\nduplicate names in loc.names:\n")
      print(temp[temp>1])
  }

  if(!is.null(object@all.names)){
      if(length(unlist(object@all.names)) != p) {
          cat("\ninvalid length in all.names\n")
          return(FALSE)
      }
  }

  return(TRUE)

}# end .gen.valid


setClass("gen", representation(tab = "matrix",
                               loc.names = "character",
                               loc.fac = "factorOrNULL",
                               loc.nall = "intOrNum",
                               all.names = "listOrNULL",
                               call = "callOrNULL",
                               "VIRTUAL"),
         prototype(tab=matrix(ncol=0,nrow=0), loc.nall=integer(0), call=NULL))

setValidity("gen", .gen.valid)





########################
# virtual class indInfo
########################
setClass("indInfo", representation(ind.names = "character",
                                   pop = "factorOrNULL",
                                   pop.names = "charOrNULL",
                                   ploidy = "integer",
                                   type = "character",
                                   other = "listOrNULL", "VIRTUAL"),
         prototype(pop=NULL, pop.names = NULL, type = "codom", ploidy = as.integer(2), other = NULL))





###############
# Class genind
###############
.genind.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)

    if(length(object@ind.names) != nrow(object@tab)) {
        cat("\ninvalid length in ind.names\n")
        return(FALSE)
    }

    temp <- table(object@ind.names[object@ind.names!=""])
    if(any(temp>1)) {
        warning("\nduplicate names in ind.names:\n")
        print(temp[temp>1])
    }

    if(!is.null(object@pop)){ # check pop

        if(length(object@pop) != nrow(object@tab)) {
            cat("\npop is given but has invalid length\n")
            return(FALSE)
        }

        if(is.null(object@pop.names)) {
            cat("\npop is provided without pop.names")
        }


        if(length(object@pop.names) != length(levels(object@pop))) {
            cat("\npop.names has invalid length\n")
            return(FALSE)
        }

        temp <- table(object@pop.names[object@pop.names!=""])
        if(any(temp>1)) {
            warning("\nduplicate names in pop.names:\n")
            print(temp[temp>1])
        }

    } # end check pop

    ## check ploidy
    if(object@ploidy < as.integer(1)){
        cat("\nploidy inferior to 1\n")
        return(FALSE)
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        cat("\nunknown type of marker\n")
        return(FALSE)
    }


    return(TRUE)
} #end .genind.valid

setClass("genind", contains=c("gen", "indInfo"))
setValidity("genind", .genind.valid)



########################
# virtual class popInfo
########################
setClass("popInfo", representation(pop.names = "character", ploidy = "integer",
                                   type = "character", other = "listOrNULL", "VIRTUAL"),
         prototype(type = "codom", ploidy = as.integer(2), other = NULL))



###############
# Class genpop
###############
.genpop.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)
    if(length(object@pop.names) != nrow(object@tab)) {
        cat("\ninvalid length in pop.names\n")
        return(FALSE)
    }

    temp <- table(object@pop.names[object@pop.names!=""])
    if(any(temp>1)) {
        warning("\nduplicate names in pop.names:\n")
        print(temp[temp>1])
    }

     ## check ploidy
    if(object@ploidy < as.integer(1)){
        cat("\nploidy inferior to 1\n")
        return(FALSE)
    }

    ## check type of marker
    if(!object@type %in% c("codom","PA") ){
        cat("\nunknown type of marker\n")
        return(FALSE)
    }

    return(TRUE)
} #end .genpop.valid

setClass("genpop", contains=c("gen", "popInfo"))
setValidity("genpop", .genpop.valid)







###############################################################
###############################################################
# MAIN CLASS METHODS
###############################################################
###############################################################



#################
# Function names
#################
setMethod("names", signature(x = "genind"), function(x){
    return(slotNames(x))
})

setMethod("names", signature(x = "genpop"), function(x){
    return(slotNames(x))
})





##################
# Function genind
##################
## constructor of a genind object
genind <- function(tab,pop=NULL,prevcall=NULL,ploidy=2,type=c("codom","PA")){
    ## handle arguments
    X <- as.matrix(tab)
    if(is.null(colnames(X))) stop("tab columns have no name.")
    if(is.null(rownames(X))) {rownames(X) <- 1:nrow(X)}

    type <- match.arg(type)
    ploidy <- as.integer(ploidy)
    nind <- nrow(X)


    ## HANDLE LABELS ##

    ## loc names is not type-dependent
    temp <- colnames(X)
    ## temp <- gsub("[.].*$","",temp)
    temp <- gsub("[.][^.]*$", "", temp)
    temp <- .rmspaces(temp)
    loc.names <- unique(temp)
    nloc <- length(loc.names)
    loc.codes <- .genlab("L",nloc)
    names(loc.names) <- loc.codes

    ## ind names is not type-dependent either
    ind.codes <- .genlab("", nind)
    ind.names <- .rmspaces(rownames(X))
    names(ind.names) <- ind.codes
    rownames(X) <- ind.codes


    if(type=="codom"){
        ## loc.nall
        loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]
        loc.nall <- as.integer(loc.nall)
        names(loc.nall) <- loc.codes

        ## loc.fac
        loc.fac <- rep(loc.codes,loc.nall)

        ## alleles name
        temp <- colnames(X)
        temp <- gsub("^.*[.]","",temp)
        temp <- .rmspaces(temp)
        all.names <- split(temp,loc.fac)
        all.codes <- lapply(all.names,function(e) .genlab("",length(e)))
        for(i in 1:length(all.names)){
            names(all.names[[i]]) <- all.codes[[i]]
        }

        colnames(X) <- paste(loc.fac,unlist(all.codes),sep=".")
        loc.fac <- as.factor(loc.fac)
    } else { # end if type=="codom" <=> if type=="PA"
        colnames(X) <- loc.codes
        loc.fac <- NULL
        all.names <- NULL
        loc.nall <- NULL
    }

    ## Ideally I should use an 'initialize' method here
    res <- new("genind")
    res@tab <- X
    res@ind.names <- ind.names
    res@loc.names <- loc.names
    res@loc.nall <- loc.nall
    res@loc.fac <- loc.fac
    res@all.names <- all.names

    ## populations name (optional)
    ## beware, keep levels of pop sorted in
    ## there order of appearance
    if(!is.null(pop)) {
        # convert pop to a factor if it is not
        if(!is.factor(pop)) {pop <- factor(pop)}
        pop.lab <- .genlab("P",length(levels(pop)) )
        # put pop levels in appearance order
        pop <- as.character(pop)
        pop <- factor(pop, levels=unique(pop))
        temp <- pop
        # now levels are correctly ordered
        levels(pop) <- pop.lab
        res@pop <- pop
        pop.names <- as.character(levels(temp))
        names(pop.names) <- as.character(levels(res@pop))
        res@pop.names <- pop.names
    }

    ## ploidy
    plo <- as.integer(ploidy)
    if(plo < as.integer(1)) stop("ploidy inferior to 1")
    res@ploidy <- plo

    ## type of marker
    res@type <- as.character(type)

    if(is.null(prevcall)) {prevcall <- match.call()}
    res@call <- prevcall

    return(res)

} # end genind

######################
# alias for as.genind
######################
as.genind <- genind



##################
# Function genpop
##################
genpop <- function(tab,prevcall=NULL,ploidy=as.integer(2),type=c("codom","PA")){

    ## handle args
    X <- as.matrix(tab)
    if(is.null(colnames(X))) stop("tab columns have no name.")
    if(is.null(rownames(X))) {rownames(X) <- 1:nrow(X)}

    type <- match.arg(type)
    ploidy <- as.integer(ploidy)
    npop <- nrow(X)


    ## HANDLE LABELS ##

    ## loc names is not type-dependent
    temp <- colnames(X)
    ## temp <- gsub("[.].*$","",temp)
    temp <- gsub("[.][^.]*$", "", temp)
    temp <- .rmspaces(temp)
    loc.names <- unique(temp)
    nloc <- length(loc.names)
    loc.codes <- .genlab("L",nloc)
    names(loc.names) <- loc.codes

    ## pop names is not type-dependent either
    pop.codes <- .genlab("", npop)
    pop.names <- .rmspaces(rownames(X))
    names(pop.names) <- pop.codes
    rownames(X) <- pop.codes

    ## type-dependent stuff
    if(type=="codom"){
        ## loc.nall
        loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]
        loc.nall <- as.integer(loc.nall)
        names(loc.nall) <- loc.codes

        ## loc.fac
        loc.fac <- rep(loc.codes,loc.nall)

        ## alleles name
        temp <- colnames(X)
        temp <- gsub("^.*[.]","",temp)
        temp <- .rmspaces(temp)
        all.names <- split(temp,loc.fac)
        all.codes <- lapply(all.names,function(e) .genlab("",length(e)))
        for(i in 1:length(all.names)){
            names(all.names[[i]]) <- all.codes[[i]]
        }

        rownames(X) <- pop.codes
        colnames(X) <- paste(loc.fac,unlist(all.codes),sep=".")
        loc.fac <- as.factor(loc.fac)
    } else { # end if type=="codom" <=> if type=="PA"
        colnames(X) <- loc.codes
        loc.fac <- NULL
        all.names <- NULL
        loc.nall <- NULL
    }

    res <- new("genpop")

    res@tab <- X
    res@pop.names <- pop.names
    res@loc.names <- loc.names
    res@loc.nall <- loc.nall
    res@loc.fac <- loc.fac
    res@all.names <- all.names
    res@ploidy <- ploidy
    res@type <- as.character(type)

    if(is.null(prevcall)) {prevcall <- match.call()}
    res@call <- prevcall

    return(res)

} # end genpop



######################
# alias for as.genpop
######################
as.genpop <- genpop

