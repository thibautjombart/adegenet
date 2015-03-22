########################################################################
# hierarchy methods definitions. 
#
# Zhian Kamvar, March 2015
# kamvarz@science.oregonstate.edu
########################################################################
#==============================================================================#
#==============================================================================#
# Workhorse functions for hierarchy methods. These avoids needless copying and
# pasting of code to dispatch the methods to two unrelated classes 
# (genind and genlight).
#==============================================================================#
#==============================================================================#
.getHier <- function(x, formula = NULL, combine = TRUE){
  if (is.null(x@hierarchy)) return(NULL)
  if (is.null(formula)) return(x@hierarchy)
  vars <- all.vars(formula)
  if (any(!vars %in% names(x@hierarchy))){
    stop(.hier_incompatible_warning(vars, x@hierarchy))
  }
  if (combine){
    hier <- .make_hierarchy(formula, x@hierarchy)
  } else {
    hier <- x@hierarchy[all.vars(formula)]
  }
  invisible(return(hier))
}

.setHier <- function(x, value){
  if (!inherits(value, "data.frame")){
    stop(paste(substitute(value), "is not a data frame"))
  }
  if (nrow(value) != nInd(x)){
    stop("Number of rows in data frame not equal to number of individuals in object.")
  }
  value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
  x@hierarchy <- value
  return(x)
}

.nameHier <- function(x, value){
  if (is.null(x@hierarchy)){
    warning("Cannot name an empty hierarchy")
    return(x)
  }
  if (is.language(value)){
    value <- all.vars(value)
  }
  if (!is.vector(value) | length(value) != length(x@hierarchy)){
    stop(paste("Hierarchy, needs a vector argument of length", length(x@hierarchy)))
  }
  names(x@hierarchy) <- value
  return(x)
}

.splitHier <- function(x, value, sep = "_"){
  if (is.null(x@hierarchy)){
    warning("Cannot split an empty hierarchy")
    return(x)
  }
  if (is.language(value)){
    # valterms <- attr(terms(value), "term.labels")
    # valterms <- valterms[length(valterms)]
    # valterms <- gsub(":", sep, valterms)
    value    <- all.vars(value)
  } else {
    stop("value must be a formula.")
  }
  if (length(value) < 1){
    stop("value must have more than one hierarchical level.")
  }
  hierarchy  <- x@hierarchy
  if (length(hierarchy) > 1){
    warning("Hierarchy must be length 1. Taking the first column.")
    hierarchy <- hierarchy[1]
  }
  seps     <- gregexpr(sep, hierarchy[[1]])
  sepmatch <- vapply(seps, function(val) all(as.integer(val) > 0), logical(1))
  seps     <- vapply(seps, length, numeric(1))
  all_seps_match <- all(sepmatch)
  given_seps     <- length(value) - 1 
  if (!all_seps_match | all(seps != given_seps)){
    seps <- ifelse(all_seps_match, seps[1], 0) + 1
    msg1 <- paste("\n  Data has", seps, ifelse(seps == 1, "level", "levels"),
                  "of hierarchy with the separator", sep, ".")
    msg2 <- paste("Here is the fist column of the data:", hierarchy[1, ])
    stop(paste(msg1, "\n ", msg2))
  }
  x@hierarchy <- colsplit(as.character(hierarchy[[1]]), pattern = sep, value)
  x@hierarchy <- data.frame(lapply(x@hierarchy, function(f) factor(f, levels = unique(f))))
  # names(hierarchy) <- value
  # x@hierarchy      <- hierarchy
  return(x) 
}

.addHier <- function(x, value, name = "NEW"){
  if (is.null(x@hierarchy)){
    hierarchy <- data.frame(vector(mode = "character", length = nInd(x)))
    wasNULL <- TRUE
  } else {
    hierarchy  <- x@hierarchy  
    wasNULL <- FALSE    
  }

  if ((is.vector(value) | is.factor(value)) & length(value) == nrow(hierarchy)){
    value <- factor(value, levels = unique(value))
    NEW <- data.frame(value)
    names(NEW) <- name
    hierarchy <- cbind(hierarchy, NEW)
  } else if (is.data.frame(value) && nrow(value) == nrow(hierarchy)){
    value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
    hierarchy <- cbind(hierarchy, value)
  } else {
    stop("value must be a vector or data frame.")
  }
  if (wasNULL){
    hierarchy <- hierarchy[-1, , drop = FALSE]
  }
  x@hierarchy <- hierarchy
  return(x) 
}

.setPop <- function(x, formula = NULL){
  if (is.null(x@hierarchy)){
    warning("Cannot set the population from an empty hierarchy")
    return(x)
  }
  if (is.null(formula) | !is.language(formula)){
    stop(paste(substitute(formula), "must be a valid formula object."))
  }
  vars <- all.vars(formula)
  if (!all(vars %in% names(x@hierarchy))){
    stop(.hier_incompatible_warning(vars, x@hierarchy))
  }
  pop(x) <- .make_hierarchy(formula, x@hierarchy)[[length(vars)]]
  return(x)
}

#==============================================================================#
#==============================================================================#
# Formal methods definitions and documentation.
#==============================================================================#
#==============================================================================#

#==============================================================================#
#' Access and manipulate the population hierarchy for genind or genlight objects.
#' 
#' The following methods allow the user to quickly change the hierarchy or
#' population of a genind or genlight object. 
#' 
#' @export 
#' @rdname hierarchy-methods
#' @aliases gethierarchy,genind-method gethierarchy,genlight-method
#' @param x a genind or genlight object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param combine if \code{TRUE}, the levels will be combined according to the
#' formula argument. If it is \code{FALSE}, the levels will not be combined.
#' @docType methods
#==============================================================================#
gethierarchy <- function(x, formula = NULL, combine = TRUE){
  standardGeneric("gethierarchy")
} 

#' @export
setGeneric("gethierarchy")

setMethod(
  f = "gethierarchy",
  signature(x = "genind"),
  definition = function(x, formula = NULL, combine = TRUE){
    .getHier(x, formula = formula, combine = combine)
  })

#==============================================================================#
#' @export
#' @rdname hierarchy-methods
#' @aliases sethierarchy<-,genind-method sethierarchy<-,genlight-method
#' @param value a data frame OR vector OR formula (see details).
#' @docType methods
#'   
#' @details \subsection{Function Specifics}{ \itemize{ \item
#' \strong{gethierarchy()} - This will retrieve the data from the
#' \emph{hierarchy} slot in the \linkS4class{genind} object. You have the
#' option to choose specific heirarchical levels using a formula (see below) and
#' you can choose to combine the hierarchical levels (default) \item
#' \strong{sethierarchy()} - Set or reset the hierarchical levels in your
#' \linkS4class{genind} object. \item \strong{namehierarchy()} - Rename the
#' hierarchical levels. \item \strong{splithierarchy()} - It is often
#' difficult to import files with several levels of hierarchy as most data
#' formats do not allow unlimited population levels. This is circumvented by
#' collapsing all hierarchical levels into a single population factor with a
#' common separator for each observation. This function will then split those
#' hierarchies for you, but it works best on a hierarchy that only has a single
#' column in it. See the rootrot example below. \item \strong{addhierarchy()} -
#' Add levels to your population hierarchy. If you have extra hierarchical
#' levels you want to add to your population hierarchy, you can use this method
#' to do so. You can input a data frame or a vector, but if you put in a vector,
#' you have the option to name it. }}
#' 
#' \subsection{Argument Specifics}{
#' 
#' These functions allow the user to seamlessly assign the hierarchical levels
#' of their \code{\linkS4class{genind}} object. Note that there are two ways
#' of performing all methods (except for \code{gethierarchy()}). They
#' essentially do the same thing except that the assignment method (the one with
#' the "\code{<-}") will modify the object in place whereas the non-assignment 
#' method will not modify the original object. Due to convention, everything 
#' right of the assignment is termed \code{value}. To avoid confusion, here is a
#' guide to the inputs: \itemize{ \item \strong{sethierarchy()} This will be a 
#' \code{\link{data.frame}} that defines the hierarchy for each individual in 
#' the rows. \item \strong{namehierarchy()} This will be either a 
#' \code{\link{vector}} or a \code{\link{formula}} that will define the names. 
#' \item \strong{splithierarchy()} This will be a \code{\link{formula}} argument
#' with the same number of levels as the hierarchy you wish to split. \item 
#' \strong{addhierarchy()} This will be a \code{\link{vector}} or 
#' \code{\link{data.frame}} with the same length as the number of individuals in
#' your data. }}
#' 
#' \subsection{Details on Formulas}{
#' 
#' The preferred use of these functions is with a \code{\link{formula}} object. 
#' Specifically, a hierarchical formula argument is used to assign the levels of
#' the hierarchy. An example of a hierarchical formula would be:\cr 
#' \code{~Country/City/Neighborhood}\cr or \cr \code{~Country + Country:City + 
#' Country:City:Neighborhood}\cr of course, the first method is slightly easier 
#' to read. It is important to use hiearchical formulas when specifying 
#' hierarchies as other types of formulas (eg. 
#' \code{~Country*City*Neighborhood}) might give spurious results.}
#' 
#' @seealso \code{\link{setpop}} \code{\link{genind}}
#'   \code{\link{as.genind}}
#'   
#' @author Zhian N. Kamvar
#' @examples
#' # let's look at the microbov data set:
#' data(microbov)
#' microbov
#' 
#' # We see that we have three vectors of different names in the 'other' slot. 
#' ?microbov
#' # These are Country, Breed, and Species
#' names(other(microbov))
#' 
#' # Let's set the hierarchy
#' sethierarchy(microbov) <- data.frame(other(microbov))
#' microbov
#' 
#' # And change the names so we know what they are
#' namehierarchy(microbov) <- ~Country/Breed/Species
#' 
#' # let's see what the hierarchy looks like by Species and Breed:
#' head(gethierarchy(microbov, ~Breed/Species))
#' 
#==============================================================================#
sethierarchy <- function(x, value){
  standardGeneric("sethierarchy")
} 

#' @export
setGeneric("sethierarchy")

setMethod(
  f = "sethierarchy",
  signature(x = "genind"),
  definition = function(x, value){
    .setHier(x, value)
  })

setMethod(
  f = "sethierarchy",
  signature(x = "genlight"),
  definition = function(x, value){
    .setHier(x, value)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases sethierarchy,genind-method sethierarchy,genlight-method
#' @docType methods
#==============================================================================#
"sethierarchy<-" <- function(x, value){
  standardGeneric("sethierarchy<-")
}  

#' @export
setGeneric("sethierarchy<-")

setMethod(
  f = "sethierarchy<-",
  signature(x = "genind"),
  definition = function(x, value){
    return(sethierarchy(x, value))
  })

setMethod(
  f = "sethierarchy<-",
  signature(x = "genlight"),
  definition = function(x, value){
    return(sethierarchy(x, value))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy,genind-method namehierarchy,genlight-method
#' @docType methods
#==============================================================================#
namehierarchy <- function(x, value){
  standardGeneric("namehierarchy")
}  

#' @export
setGeneric("namehierarchy")

setMethod(
  f = "namehierarchy",
  signature(x = "genind"),
  definition = function(x, value){
    .nameHier(x, value)
  })

setMethod(
  f = "namehierarchy",
  signature(x = "genlight"),
  definition = function(x, value){
    .nameHier(x, value)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy<-,genind-method namehierarchy<-,genlight-method
#' @docType methods
#==============================================================================#
"namehierarchy<-" <- function(x, value){
  standardGeneric("namehierarchy<-")
}  

#' @export
setGeneric("namehierarchy<-")

setMethod(
  f = "namehierarchy<-",
  signature(x = "genind"),
  definition = function(x, value){
    return(namehierarchy(x, value))
  })

setMethod(
  f = "namehierarchy<-",
  signature(x = "genlight"),
  definition = function(x, value){
    return(namehierarchy(x, value))
  })
#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy,genind-method splithierarchy,genlight-method
#' @docType methods
#' @param sep a \code{character} indicating the character used to separate
#' hierarchical levels. This defaults to "_".
#' @importFrom reshape2 colsplit
#==============================================================================#
splithierarchy <- function(x, value, sep = "_"){
  standardGeneric("splithierarchy")
}  

#' @export
setGeneric("splithierarchy")

setMethod(
  f = "splithierarchy",
  signature(x = "genind"),
  definition = function(x, value, sep = "_"){
    .splitHier(x, value, sep = sep) 
  })

setMethod(
  f = "splithierarchy",
  signature(x = "genlight"),
  definition = function(x, value, sep = "_"){
    .splitHier(x, value, sep = sep) 
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy<-,genind-method splithierarchy<-,genlight-method
#' @docType methods
#==============================================================================#
"splithierarchy<-" <- function(x, sep = "_", value){
  standardGeneric("splithierarchy<-")
}  

#' @export
setGeneric("splithierarchy<-")

setMethod(
  f = "splithierarchy<-",
  signature(x = "genind"),
  definition = function(x, sep = "_", value){
    return(splithierarchy(x, value, sep))
  })

setMethod(
  f = "splithierarchy<-",
  signature(x = "genlight"),
  definition = function(x, sep = "_", value){
    return(splithierarchy(x, value, sep))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy,genind-method addhierarchy,genlight-method
#' @param name an optional name argument for use with addhierarchy if supplying
#'   a vector. Defaults to "NEW".
#' @docType methods
#==============================================================================#
addhierarchy <- function(x, value, name = "NEW"){
  standardGeneric("addhierarchy")
}  

#' @export
setGeneric("addhierarchy")

setMethod(
  f = "addhierarchy",
  signature(x = "genind"),
  definition = function(x, value, name = "NEW"){
    .addHier(x, value, name = name)
  })

setMethod(
  f = "addhierarchy",
  signature(x = "genlight"),
  definition = function(x, value, name = "NEW"){
    .addHier(x, value, name = name)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy<-,genind-method addhierarchy<-,genlight-method
#' @docType methods
#==============================================================================#
"addhierarchy<-" <- function(x, name = "NEW", value){
  standardGeneric("addhierarchy<-")
}  

#' @export
setGeneric("addhierarchy<-")

setMethod(
  f = "addhierarchy<-",
  signature(x = "genind"),
  definition = function(x, name = "NEW", value){
    return(addhierarchy(x, value, name))
  })

setMethod(
  f = "addhierarchy<-",
  signature(x = "genlight"),
  definition = function(x, name = "NEW", value){
    return(addhierarchy(x, value, name))
  })

#==============================================================================#
#' Manipulate the population factor of genind objects.
#' 
#' The following methods allow the user to quickly change the population of a 
#' genind object. 
#' 
#' @export 
#' @rdname population-methods
#' @param x a genind or genlight object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param value same as formula
#' @aliases setpop,genind-method setpop,genlight-method
#' @docType methods 
#' @author Zhian N. Kamvar
#' @examples
#' 
#' data(microbov)
#' 
#' sethierarchy(microbov) <- data.frame(other(microbov))
#' 
#' # Currently set on just 
#' head(pop(microbov)) 
#' 
#' # setting the hierarchy to both Pop and Subpop
#' setpop(microbov) <- ~coun/breed 
#' head(pop(microbov))
#' 
#' \dontrun{
#' 
#' # Can be used to create objects as well.
#' microbov.old <- setpop(microbov, ~spe) 
#' head(pop(microbov.old))
#' }
#==============================================================================#
setpop <- function(x, formula = NULL) standardGeneric("setpop")

#' @export
setGeneric("setpop")

setMethod(
  f = "setpop",
  signature(x = "genind"),
  definition = function(x, formula = NULL){
    .setPop(x, formula = formula)
  })

setMethod(
  f = "setpop",
  signature(x = "genlight"),
  definition = function(x, formula = NULL){
    .setPop(x, formula = formula)
  })
#==============================================================================#
#' @export
#' @rdname population-methods
#' @aliases setpop<-,genind-method setpop<-,genlight-method
#' @docType methods
#==============================================================================#
"setpop<-" <- function(x, value) standardGeneric("setpop<-")

#' @export
setGeneric("setpop<-")

setMethod(
  f = "setpop<-",
  signature(x = "genind"),
  definition = function(x, value){
    return(setpop(x, value))
  })

setMethod(
  f = "setpop<-",
  signature(x = "genlight"),
  definition = function(x, value){
    return(setpop(x, value))
  })

#==============================================================================#
#==============================================================================#
# internal functions utilized
#==============================================================================#
#==============================================================================#

#==============================================================================#
# A function for creating a population hierarchy using a formula and data frame
# 
# hier = a nested formula such as ~ A/B/C where C is nested within B, which is
# nested within A.
#
# df = a data frame containing columns corresponding to the variables in hier.
#
# example:
# df <- data.frame(list(a = letters, b = LETTERS, c = 1:26))
# newdf <- .make_hierarchy(~ a/b/c, df)
# df[names(newdf)] <- newdf # Add new columns.
#
# Public functions utilizing this function:
#
# # setpop, gethierarchy
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
.make_hierarchy <- function(hier, df, expand_label = FALSE){
  newlevs <- attr(terms(hier), "term.labels")
  levs <- all.vars(hier)
  if (length(levs) > 1){
    newlevs <- gsub(":", "_", newlevs)
  }
  if (!all(levs %in% names(df))){
    stop(.hier_incompatible_warning(levs, df))
  }
  newdf <- df[levs[1]]
  if (!expand_label){
    newlevs <- levs
  }
  lapply(1:length(levs), function(x) newdf[[newlevs[x]]] <<- as.factor(.pop_combiner(df, levs[1:x])))
  return(newdf)
}

#==============================================================================#
# This will be used to join heirarchical population vectors for the purposes of
# maintaining hierarchy. 
# Public functions utilizing this function:
# # 
#
# Internal functions utilizing this function:
# # .make_hierarchy
#==============================================================================#

.pop_combiner <- function(df, hier=c(1), sep="_"){
  if(!is.list(df)){
    warning("df must be a data frame or a list")
    return(df)
  }
  else{
    if(length(hier)==1){
      return(df[[hier]])
    }
    else{
      comb <- vector(length=length(df[[hier[1]]]))
      comb <- df[[hier[1]]]
      lapply(hier[-1], function(x) comb <<- paste(comb, df[[x]], sep=sep))
      return(comb)
    }
  }
}

#==============================================================================#
# A function that will quit the function if a level in the hierarchy is not
# present in the given data frame.
#
# Public functions utilizing this function:
# # setpop gethierarchy poppr.amova
#
# Internal functions utilizing this function:
# # .make_hierarchy make_ade_df
#==============================================================================#
.hier_incompatible_warning <- function(levs, df){
  msg <- paste("One or more levels in the given hierarchy is not present", 
               "in the data frame.",
               "\nHierarchy:\t", paste(levs, collapse = ", "), "\nData:\t\t", 
               paste(names(df), collapse = ", "))
  return(msg)
}