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
# CLASSES DEFINITION
###############################################################
###############################################################

#' Virtual classes for adegenet
#'
#' These virtual classes are only for internal use in adegenet
#'
#'
#' @name virtualClasses
#' @aliases indInfo-class popInfo-class gen-class callOrNULL-class
#' charOrNULL-class factorOrNULL-class intOrNum-class listOrNULL-class
#' intOrNULL-class
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords classes
NULL

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




#' adegenet formal class (S4) for individual genotypes
#'
#' The S4 class \code{genind} is used to store individual genotypes.\cr It
#' contains several components described in the 'slots' section).\cr The
#' \code{summary} of a \code{genind} object invisibly returns a list of
#' component.  The function \code{.valid.genind} is for internal use.  The
#' function \code{genind} creates a genind object from a valid table of alleles
#' corresponding to the \code{@@@tab} slot.  Note that as in other S4 classes,
#' slots are accessed using @@ instead of \$.
#'
#'
#' @aliases genind-class print,genind-method show,genind-method
#' names,genind-method summary,genind-method .valid.genind
#' @section Slots: \describe{ \item{list("tab")}{matrix of genotypes (in rows)
#' for all alleles (in columns). The table differs depending on the
#' \code{@@type} slot:\cr - 'codom': values are frequencies ; '0' if the
#' genotype does not have the corresponding allele, '1' for an homozygote and
#' 0.5 for an heterozygte.\cr - 'PA': values are presence/absence of
#' alleles.\cr In all cases, rows and columns are given generic
#' names.}\item{:}{matrix of genotypes (in rows) for all alleles (in columns).
#' The table differs depending on the \code{@@type} slot:\cr - 'codom': values
#' are frequencies ; '0' if the genotype does not have the corresponding
#' allele, '1' for an homozygote and 0.5 for an heterozygte.\cr - 'PA': values
#' are presence/absence of alleles.\cr In all cases, rows and columns are given
#' generic names.} \item{list("loc.names")}{character vector containing the
#' real names of the loci}\item{:}{character vector containing the real names
#' of the loci} \item{list("loc.fac")}{locus factor for the columns of
#' \code{tab}}\item{:}{locus factor for the columns of \code{tab}}
#' \item{list("loc.nall")}{integer vector giving the number of alleles per
#' locus}\item{:}{integer vector giving the number of alleles per locus}
#' \item{list("all.names")}{list having one component per locus, each
#' containing a character vector of alleles names}\item{:}{list having one
#' component per locus, each containing a character vector of alleles names}
#' \item{list("call")}{the matched call}\item{:}{the matched call}
#' \item{list("ind.names")}{character vector containing the real names of the
#' individuals. Note that as Fstat does not store these names, objects
#' converted from .dat files will contain empty
#' \code{ind.names}.}\item{:}{character vector containing the real names of the
#' individuals. Note that as Fstat does not store these names, objects
#' converted from .dat files will contain empty \code{ind.names}.}
#' \item{list("ploidy")}{ an integer indicating the degree of ploidy of the
#' genotypes. Beware: 2 is not an integer, but as.integer(2) is.}\item{:}{ an
#' integer indicating the degree of ploidy of the genotypes. Beware: 2 is not
#' an integer, but as.integer(2) is.} \item{list("type")}{ a character string
#' indicating the type of marker: 'codom' stands for 'codominant' (e.g.
#' microstallites, allozymes); 'PA' stands for 'presence/absence' (e.g.
#' AFLP).}\item{:}{ a character string indicating the type of marker: 'codom'
#' stands for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).} \item{list("pop")}{(optional) factor giving
#' the population of each individual}\item{:}{(optional) factor giving the
#' population of each individual} \item{list("pop.names")}{(optional) vector
#' giving the real names of the populations}\item{:}{(optional) vector giving
#' the real names of the populations} \item{list("other")}{(optional) a list
#' containing other information}\item{:}{(optional) a list containing other
#' information} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{as.genind}}, \code{\link{is.genind}},
#' \code{\link{genind2genpop}}, \code{\link{genpop}},
#' \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.genepop}}, \code{\link{read.fstat}},
#' \code{\link{na.replace}}\cr
#'
#' Related classes:\cr - \linkS4class{genpop} for storing data per
#' populations\cr
#'
#' - \linkS4class{genlight} for an efficient storage of binary SNPs
#' genotypes\cr
#' @keywords classes manip multivariate
#' @examples
#'
#' showClass("genind")
#'
#' obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"),missing="mean")
#' obj
#' validObject(obj)
#' summary(obj)
#'
#' \dontrun{
#' # test inter-colonies structuration
#' if(require(hierfstat)){
#' gtest <- gstat.randtest(obj,nsim=99)
#' gtest
#' plot(gtest)
#' }
#'
#' # perform a between-class PCA
#' pca1 <- dudi.pca(obj@@tab,scannf=FALSE,scale=FALSE)
#' pcabet1 <- between(pca1,obj@@pop,scannf=FALSE)
#' pcabet1
#'
#' s.class(pcabet1$ls,obj@@pop,sub="Inter-class PCA",possub="topleft",csub=2)
#' add.scatter.eig(pcabet1$eig,2,xax=1,yax=2)
#'
#' }
#'
NULL
setClass("genind", contains=c("gen", "indInfo"))
setValidity("genind", .genind.valid)






#' adegenet formal class (S4) for allele counts in populations
#'
#' An object of class \code{genpop} contain alleles counts for several loci.\cr
#' It contains several components (see 'slots' section).\cr Such object is
#' obtained using \code{genind2genpop} which converts individuals genotypes of
#' known population into a \code{genpop} object.  Note that the function
#' \code{summary} of a \code{genpop} object returns a list of components.  Note
#' that as in other S4 classes, slots are accessed using @@ instead of \$.
#'
#'
#' @aliases genpop-class dist,genpop,ANY,ANY,ANY,missing-method
#' names,genpop-method show,genpop-method summary,genpop-method
#' @section Slots: \describe{ \item{list("tab")}{matrix of alleles counts for
#' each combinaison of population -in rows- and alleles -in columns-. Rows and
#' columns are given generic names.}\item{:}{matrix of alleles counts for each
#' combinaison of population -in rows- and alleles -in columns-. Rows and
#' columns are given generic names.} \item{list("loc.names")}{character vector
#' containing the real names of the loci}\item{:}{character vector containing
#' the real names of the loci} \item{list("loc.fac")}{locus factor for the
#' columns of \code{tab}}\item{:}{locus factor for the columns of \code{tab}}
#' \item{list("loc.nall")}{integer vector giving the number of alleles per
#' locus}\item{:}{integer vector giving the number of alleles per locus}
#' \item{list("all.names")}{list having one component per locus, each
#' containing a character vector of alleles names}\item{:}{list having one
#' component per locus, each containing a character vector of alleles names}
#' \item{list("call")}{the matched call}\item{:}{the matched call}
#' \item{list("pop.names")}{character vector containing the real names of the
#' populations}\item{:}{character vector containing the real names of the
#' populations} \item{list("ploidy")}{ an integer indicating the degree of
#' ploidy of the genotypes. Beware: 2 is not an integer, but as.integer(2)
#' is.}\item{:}{ an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.} \item{list("type")}{ a
#' character string indicating the type of marker: 'codom' stands for
#' 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).}\item{:}{ a character string indicating the
#' type of marker: 'codom' stands for 'codominant' (e.g. microstallites,
#' allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
#' \item{list("other")}{(optional) a list containing other
#' information}\item{:}{(optional) a list containing other information} }
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{as.genpop}},
#' \code{\link{is.genpop}},\code{\link{makefreq}}, \code{\link{genind}},
#' \code{\link{import2genind}}, \code{\link{read.genetix}},
#' \code{\link{read.genepop}}, \code{\link{read.fstat}},
#' \code{\link{na.replace}}
#' @keywords classes manip multivariate
#' @examples
#'
#' obj1 <- import2genind(system.file("files/nancycats.gen",
#' package="adegenet"))
#' obj1
#'
#'
#' obj2 <- genind2genpop(obj1)
#' obj2
#'
#' \dontrun{
#' data(microsatt)
#' # use as.genpop to convert convenient count tab to genpop
#' obj3 <- as.genpop(microsatt$tab)
#' obj3
#'
#' all(obj3@@tab==microsatt$tab)
#' all(obj3@@pop.names==rownames(microsatt$tab))
#' # it worked
#'
#' # perform a correspondance analysis
#' obj4 <- genind2genpop(obj1,missing="chi2")
#' ca1 <- dudi.coa(as.data.frame(obj4@@tab),scannf=FALSE)
#' s.label(ca1$li,sub="Correspondance Analysis",csub=2)
#' add.scatter.eig(ca1$eig,2,xax=1,yax=2,posi="top")
#'
#' }
#'

setClass("popInfo", representation(pop.names = "character", ploidy = "integer",
                                   type = "character", other = "listOrNULL", "VIRTUAL"),
         prototype(type = "codom", ploidy = as.integer(2), other = NULL))


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






#' genind constructor
#'
#' Constructor for \linkS4class{genind} objects.\cr The function \code{genind}
#' creates a \linkS4class{genind} object from a matrix of allelic frequency
#' where genotypes are in rows and alleles in columns. This table must have
#' correct names for rows and columns.\cr
#'
#' The function \code{as.genind} is an alias for \code{genind} function.\cr
#'
#' \code{is.genind} tests if an object is a valid genind object.\cr
#'
#' Note: to get the manpage about \linkS4class{genind}, please type 'class ?
#' genind'.
#'
#'
#' @aliases genind-methods genind as.genind is.genind
#' @param tab A table corresponding to the @@tab slot of a genind object, with
#' individuals in rows and alleles in columns.  Its content depends on
#' \code{type} (type of marker).\cr - 'codom': table contains allele
#' frequencies (numeric values summing to 1).\cr - 'PA': table contains binary
#' values, which indicate presence(1)/absence(0) of alleles.\cr
#' @param pop a factor giving the population of each genotype in 'x'
#' @param prevcall call of an object
#' @param ploidy an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.
#' @param type a character string indicating the type of marker: 'codom' stands
#' for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP).
#' @param x an object
#' @return For \code{genind} and \code{as.genind}, a genind object. For
#' \code{is.genind}, a logical.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\linkS4class{genind} class}, and \code{\link{import2genind}}
#' for importing from various types of file.\cr
#'
#' Related classes:\cr - \linkS4class{genpop} for storing data per
#' populations\cr
#'
#' - \linkS4class{genlight} for an efficient storage of binary SNPs
#' genotypes\cr
#' @keywords manip
#' @examples
#'
#' data(nancycats)
#' nancycats@@loc.names
#'
#' # isolate one marker, fca23
#' obj <- seploc(nancycats)$"fca23"
#' obj
#'
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




#' Converting genind/genpop objects to other classes
#'
#' These S3 and S4 methods are used to coerce \linkS4class{genind} and
#' \linkS4class{genpop} objects to matrix-like objects. In most cases, this is
#' equivalent to calling the \code{@@tab} slot. An exception to this is the
#' convertion to \code{\link[ade4]{ktab}} objects used in the ade4 package as
#' inputs for K-tables methods (e.g. Multiple Coinertia Analysis).\cr
#'
#'
#' @name as methods in adegenet
#' @aliases as-method as,genind,data.frame-method as,genpop,data.frame-method
#' as,genind,matrix-method as,genpop,matrix-method as,genind,genpop-method
#' ktab-class as,genind,ktab-method as,genpop,ktab-method
#' coerce,genind,data.frame-method coerce,genpop,data.frame-method
#' coerce,genind,matrix-method coerce,genpop,matrix-method
#' coerce,genind,genpop-method coerce,genind,ktab-method
#' coerce,genpop,ktab-method as.data.frame.genind as.data.frame.genpop
#' as.matrix.genind as.matrix.genpop as.genpop.genind as.ktab.genind
#' as.ktab.genpop
#' @docType methods
#' @section Usage: \code{as(object, Class)}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' data(microbov)
#' x <- na.replace(microbov,method="0")
#' as(x[1:3],"data.frame")
#'
#' ## dudi functions attempt to convert their first argument
#' ## to a data.frame; so they can be used on genind/genpop objects.
#' ## perform a PCA
#' pca1 <- dudi.pca(x, scale=FALSE, scannf=FALSE)
#' pca1
#'
#' x <- genind2genpop(microbov,miss="chi2")
#' x <- as(x,"ktab")
#' class(x)
#' ## perform a STATIS analysis
#' statis1 <- statis(x, scannf=FALSE)
#' statis1
#' plot(statis1)
#'
#' }
#'
as.genind <- genind




#' genpop constructor
#'
#' Constructor for \linkS4class{genpop} objects.\cr The function \code{genpop}
#' creates a \linkS4class{genpop} object from a matrix of alleles counts where
#' genotypes are in rows and alleles in columns. This table must have correct
#' names for rows and columns.\cr
#'
#' The function \code{as.genpop} is an alias for \code{genpop} function.\cr
#'
#' \code{is.genpop} tests if an object is a valid genpop object.\cr
#'
#' Note: to get the manpage about \linkS4class{genpop}, please type 'class ?
#' genpop'.
#'
#'
#' @aliases genpop-methods genpop as.genpop is.genpop
#' @param tab a pop x alleles matrix which terms are numbers of alleles, i.e.
#' like in a genpop object
#' @param prevcall call of an object
#' @param ploidy an integer indicating the degree of ploidy of the genotypes.
#' Beware: 2 is not an integer, but as.integer(2) is.
#' @param type a character string indicating the type of marker: 'codom' stands
#' for 'codominant' (e.g. microstallites, allozymes); 'PA' stands for
#' 'presence/absence' (e.g. AFLP, RAPD).
#' @param x an object
#' @return For \code{genpop} and \code{as.genpop}, a genpop object. For
#' \code{is.genpop}, a logical.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\linkS4class{genpop} class}, and \code{\link{genind2genpop}}
#' for conversion from a genind to a genpop object.
#' @keywords manip
#' @examples
#'
#' data(nancycats)
#' obj <- genind2genpop(nancycats)
#'
#' # isolate one locus, fca77
#' obj <- seploc(obj)$"fca77"
#' obj
#'
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

