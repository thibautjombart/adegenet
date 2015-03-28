
#' genind constructor
#'
#' The function \code{new} has a method for building \linkS4class{genind} objects.
#' See the class description of \linkS4class{genind} for more information on this data structure.
#' The functions \code{genind} and \code{as.genind} are aliases for \code{new("genind", ...)}.
#'
#' Most users do not need using the constructor, but merely to convert raw allele data using \code{\link{df2genind}} and related functions.
#'
#' @export
#' @docType methods
#'   
#' @aliases initialize,genind-methods
#' @aliases genind
#' @aliases as.genind
#'   
#' @rdname new.genind
#'   
#' @param .Object prototyped object (generated automatically when calling 'new')
#' @param tab A matrix of integers corresponding to the @@tab slot of a genind
#'   object, with individuals in rows and alleles in columns, and containing
#'   either allele counts (if type="codom") or allele presence/absence (if
#'   type="PA")
#' @param pop an optional factor with one value per row in \code{tab} indicating
#'   the population of each individual
#' @param prevcall an optional call to be stored in the object
#' @param ploidy an integer vector indicating the ploidy of the individual; each
#'   individual can have a different value; if only one value is provided, it is
#'   recycled to generate a vector of the right length.
#' @param type a character string indicating the type of marker: codominant
#'   ("codom") or presence/absence ("PA")
#' @param strata a data frame containing population hierarchies or
#'   stratifications in columns. This must be the same length as the number of
#'   individuals in the data set.
#' @param hierarchy a hierarchical formula defining the columns of the strata
#'   slot that are hierarchical. Defaults to NULL.
#' @param ... further arguments passed to other methods (currently not used)
#'
#' @return a \linkS4class{genind} object
#'
#' @seealso the description of the \linkS4class{genind} class; \code{\link{df2genind}}
#'
setMethod("initialize", "genind", function(.Object, tab, pop=NULL, prevcall=NULL, ploidy=2L, type=c("codom","PA"), strata = NULL, hierarchy = NULL, ...){
   ## HANDLE ARGUMENTS ##
    out <- .Object
    if(is.null(colnames(tab))) stop("tab columns have no name.")
    if(is.null(rownames(tab))) {rownames(tab) <- 1:nrow(tab)}

    ## force matrix & integer
    old.rownames <- rownames(tab)
    old.colnames <- colnames(tab)
    old.dim <- dim(tab)
    if(typeof(tab)!="integer"){
        tab <- as.integer(tab)
        dim(tab) <- old.dim
        rownames(tab) <- old.rownames
        colnames(tab) <- old.colnames
    }
    type <- match.arg(type)
    nind <- nrow(tab)
    ploidy <- as.integer(ploidy)
    ploidy <- rep(ploidy, length=nind)

    ## HANDLE LABELS ##
    ## loc names is not type-dependent
    temp <- gsub("[.][^.]*$", "", old.colnames)
    temp <- .rmspaces(temp)
    loc.names <- unique(temp)
    nloc <- length(loc.names)

    ## ind names is not type-dependent either
    ## only use generic label if no name or duplicates
    if(is.null(rownames(tab))) {
        rownames(tab) <- .genlab("", nind)
    }
    ind.names <- rownames(tab)
    if(length(unique(ind.names))!=length(ind.names)) {
        warning("duplicate labels detected for some individuals; using generic labels")
        rownames(tab) <- ind.names <- .genlab("", nind)
    }

    if (!is.null(strata)){
      # Make sure that the hierarchies are factors.
      strata <- data.frame(lapply(strata, function(f) factor(f, unique(f))))
      rownames(strata) <- rownames(tab)  
    } 
    
    if (!is.null(strata) && !is.null(hierarchy)){
      if (is.language(hierarchy)){
        the_names <- all.vars(hierarchy)
        if (all(the_names %in% names(strata))){
          ## TODO: CHECK HIERARCHY HERE
        } else {
          warning("hierarchy names do not match names of strata. Setting slot to NULL")
          hierarchy <- NULL
        }
      } else {
        warning("hierarchy must be a formula. Setting slot to NULL.")
        hierarchy <- NULL
      }
    }


    if(type=="codom"){
        ## loc.nall
        loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]
        loc.nall <- as.integer(loc.nall)
        names(loc.nall) <- loc.names

        ## loc.fac
        loc.fac <- factor(rep(loc.names,loc.nall),levels=loc.names)

        ## alleles name
        temp <- colnames(tab)
        temp <- gsub("^.*[.]","",temp)
        temp <- .rmspaces(temp)
        all.names <- split(temp,loc.fac)
        all.names <- all.names[loc.names]

    } else { # end if type=="codom" <=> if type=="PA"
        loc.fac   <- NULL
        all.names <- NULL
        loc.nall  <- NULL
    }

    ## Ideally I should use an 'initialize' method here
    out@tab       <- tab
    out@ind.names <- ind.names
    out@loc.names <- loc.names
    out@loc.nall  <- loc.nall
    out@loc.fac   <- loc.fac
    out@all.names <- all.names
    out@strata    <- strata
    out@hierarchy <- hierarchy
    ## populations name (optional)
    ## beware, keep levels of pop sorted in
    ## there order of appearance
    if(!is.null(pop)) {
        # convert pop to a factor if it is not
        if(!is.factor(pop)) {pop <- factor(pop)}
        out@pop <- pop
        out@pop.names <- levels(pop)
    }

    ## ploidy
    plo <- as.integer(ploidy)
    if(any(plo < 1L)) stop("ploidy inferior to 1")
    out@ploidy <- plo

    ## type of marker
    out@type <- as.character(type)

    if(is.null(prevcall)) {prevcall <- match.call()}
    out@call <- prevcall

    return(out)
})


#' @export
#' @rdname new.genind
genind <- function(...){
    out <- new("genind", ...)
    return(out)
} # end genind

#' @export
#' @rdname new.genind
as.genind <- function(...){
    out <- new("genind", ...)
    return(out)
} # end genind






#' genpop constructor
#'
#' The function \code{new} has a method for building \linkS4class{genpop} objects.
#' See the class description of \linkS4class{genpop} for more information on this data structure.
#' The functions \code{genpop} and \code{as.genpop} are aliases for \code{new("genpop", ...)}.
#'
#' Most users do not need using the constructor, but merely to convert raw allele data using \code{\link{genind2genpop}}.
#'
#' @export
#' @docType methods
#'
#' @aliases initialize,genpop-methods
#' @aliases genpop
#' @aliases as.genpop
#'
#' @rdname new.genpop
#'
#' @param .Object prototyped object (generated automatically when calling 'new')
#' @param tab A matrix of integers corresponding to the @@tab slot of a genpop object, with individuals in rows and alleles in columns, and containing either allele counts
#' @param prevcall an optional call to be stored in the object
#' @param ploidy an integer vector indicating the ploidy of the individual; each individual can have a different value; if only one value is provided, it is recycled to generate a vector of the right length.
#' @param type a character string indicating the type of marker: codominant ("codom") or presence/absence ("PA")
#' @param ... further arguments passed to other methods (currently not used)
#'
#' @return a \linkS4class{genpop} object
#'
#' @seealso the description of the \linkS4class{genpop} class; \code{\link{df2genind}} and related functions for reading raw allele data
#'
##################
# Function genpop
##################
setMethod("initialize", "genpop", function(.Object, tab, prevcall=NULL, ploidy=2L, type=c("codom","PA"), ...){
    ## HANDLE ARGS ##
    out <- .Object
    if(is.null(colnames(tab))) stop("tab columns have no name.")
    if(is.null(rownames(tab))) {rownames(tab) <- 1:nrow(tab)}

    ## force matrix & integer
    old.rownames <- rownames(tab)
    old.colnames <- colnames(tab)
    old.dim <- dim(tab)
    if(typeof(tab)!="integer"){
        tab <- as.integer(tab)
        dim(tab) <- old.dim
        rownames(tab) <- old.rownames
        colnames(tab) <- old.colnames
    }
    type <- match.arg(type)
    ploidy <- as.integer(ploidy)
    npop <- nrow(tab)


    ## HANDLE LABELS ##

    ## loc names is not type-dependent
    temp <- gsub("[.][^.]*$", "", old.colnames)
    temp <- .rmspaces(temp)
    loc.names <- unique(temp)
    nloc <- length(loc.names)

    ## pop names is not type-dependent either
    ## only use generic label if no name or duplicates
    if(is.null(rownames(tab))) {
        rownames(tab) <- .genlab("", npop)
    }
    pop.names <- rownames(tab)
    if(length(unique(pop.names))!=length(pop.names)) {
        warning("duplicate labels detected for some populations; using generic labels")
        rownames(tab) <- pop.names <- .genlab("", npop)
    }

    if(type=="codom"){
        ## loc.nall
        loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]
        loc.nall <- as.integer(loc.nall)
        names(loc.nall) <- loc.names

        ## loc.fac
        loc.fac <- rep(loc.names,loc.nall)

        ## alleles name
        temp <- colnames(tab)
        temp <- gsub("^.*[.]","",temp)
        temp <- .rmspaces(temp)
        all.names <- split(temp,loc.fac)
        all.names <- all.names[loc.names]
        loc.fac <- as.factor(loc.fac)

    } else { # end if type=="codom" <=> if type=="PA"
        loc.fac <- NULL
        all.names <- NULL
        loc.nall <- NULL
    }

    ## build final output
    out@tab <- tab
    out@pop.names <- pop.names
    out@loc.names <- loc.names
    out@loc.nall <- loc.nall
    out@loc.fac <- loc.fac
    out@all.names <- all.names
    out@ploidy <- ploidy
    out@type <- as.character(type)

    if(is.null(prevcall)) {prevcall <- match.call()}
    out@call <- prevcall

    return(out)
})



#' @export
#' @rdname new.genpop
genpop <- function(...){
    out <- new("genpop", ...)
    return(out)
} # end genpop

#' @export
#' @rdname new.genpop
as.genpop <- function(...){
    out <- new("genpop", ...)
    return(out)
} # end genpop
