########################################################################
# strata methods definitions. 
#
# Zhian Kamvar, March 2015
# kamvarz@science.oregonstate.edu
########################################################################
#==============================================================================#
#==============================================================================#
# Workhorse functions for strata methods. These avoids needless copying and
# pasting of code to dispatch the methods to two unrelated classes 
# (genind and genlight).
#==============================================================================#
#==============================================================================#
.getStrata <- function(x, formula = NULL, combine = TRUE, call = match.call()){
  if (is.null(x@strata)) return(NULL)
  if (is.null(formula)) return(x@strata)
  vars <- all.vars(formula)
  if (any(!vars %in% names(x@strata))){
    stop(.strata_incompatible_warning(vars, x@strata), call. = FALSE)
  }
  if (combine){
    strata <- .make_strata(formula, x@strata)
  } else {
    strata <- x@strata[all.vars(formula)]
  }
  invisible(return(strata))
}

.setStrata <- function(x, value, call = match.call()){
  if (is.null(value)){
    x@strata <- value
    return(x)
  }
  if (!inherits(value, "data.frame")){
    callval <- as.character(call["value"])
    stop(paste(callval, "is not a data frame"), call. = FALSE)
  }
  if (nrow(value) != nInd(x)){
    stop("Number of rows in data frame not equal to number of individuals in object.", call. = FALSE)
  }
  value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
  x@strata <- value
  return(x)
}

.nameStrata <- function(x, value, call = match.call()){
  if (missing(value)){
    return(names(x@strata))
  }
  if (is.null(x@strata)){
    warning("Cannot name an empty strata", call. = FALSE)
    return(x)
  }
  if (is.language(value)){
    value <- all.vars(value)
  }
  if (!is.vector(value) | length(value) != length(x@strata)){
    stop(paste("nameStrata needs a vector argument of length", length(x@strata)), call. = FALSE)
  }
  names(x@strata) <- value
  return(x)
}

.splitStrata <- function(x, value, sep = "_", call = match.call()){
  if (is.null(x@strata)){
    warning("Cannot split empty strata", call. = FALSE)
    return(x)
  }
  if (is.language(value)){
    # valterms <- attr(terms(value), "term.labels")
    # valterms <- valterms[length(valterms)]
    # valterms <- gsub(":", sep, valterms)
    value    <- all.vars(value)
  } else {
    stop("Value must be a formula.", call. = FALSE)
  }
  if (length(value) < 1){
    stop("Value must have more than one level.", call. = FALSE)
  }
  strata  <- x@strata
  if (length(strata) > 1){
    warning("Strata must be length 1. Taking the first column.", call. = FALSE)
    strata <- strata[1]
  }
  seps     <- gregexpr(sep, strata[[1]])
  sepmatch <- vapply(seps, function(val) all(as.integer(val) > 0), logical(1))
  seps     <- vapply(seps, length, numeric(1))
  all_seps_match <- all(sepmatch)
  given_seps     <- length(value) - 1 
  if (!all_seps_match | all(seps != given_seps)){
    seps <- ifelse(all_seps_match, seps[1], 0) + 1
    msg1 <- paste("\n  Data has", seps, ifelse(seps == 1, "level", "levels"),
                  "of strata with the separator", sep, ".")
    msg2 <- paste("Here is the fist column of the data:", strata[1, ])
    stop(paste(msg1, "\n ", msg2), call. = FALSE)
  }
  x@strata <- colsplit(as.character(strata[[1]]), pattern = sep, value)
  x@strata <- data.frame(lapply(x@strata, function(f) factor(f, levels = unique(f))))
  # names(strata) <- value
  # x@strata      <- strata
  return(x) 
}

.addStrata <- function(x, value, name = "NEW", call = match.call()){
  if (is.null(x@strata)){
    strata <- data.frame(vector(mode = "character", length = nInd(x)))
    wasNULL <- TRUE
  } else {
    strata  <- x@strata  
    wasNULL <- FALSE    
  }

  if ((is.vector(value) | is.factor(value)) & length(value) == nrow(strata)){
    value <- factor(value, levels = unique(value))
    NEW <- data.frame(value)
    names(NEW) <- name
    strata <- cbind(strata, NEW)
  } else if (is.data.frame(value) && nrow(value) == nrow(strata)){
    value <- data.frame(lapply(value, function(f) factor(f, unique(f))))
    strata <- cbind(strata, value)
  } else if (is.data.frame(value)){
    callval <- as.character(call["value"])
    msg     <- .not_enough_rows_warning(callval, nrow(value), nrow(strata))
    stop(msg, call. = FALSE)
  } else {
    stop("value must be a vector or data frame.", call. = FALSE)
  }
  if (wasNULL){
    strata <- strata[-1, , drop = FALSE]
  }
  x@strata <- strata
  return(x) 
}

.setPop <- function(x, formula = NULL, call = match.call()){
  if (is.null(x@strata)){
    warning("Cannot set the population from an empty strata", call. = FALSE)
    return(x)
  }
  if (is.null(formula) | !is.language(formula)){
    callform <- as.character(call["formula"])
    stop(paste(callform, "must be a valid formula object."), call. = FALSE)
  }
  vars <- all.vars(formula)
  if (!all(vars %in% names(x@strata))){
    stop(.strata_incompatible_warning(vars, x@strata), call. = FALSE)
  }
  pop(x) <- .make_strata(formula, x@strata)[[length(vars)]]
  return(x)
}

#==============================================================================#
#==============================================================================#
# Formal methods definitions and documentation.
#==============================================================================#
#==============================================================================#

#==============================================================================#
#' Access and manipulate the population strata for genind or genlight objects.
#' 
#' The following methods allow the user to quickly change the strata of a genind
#' or genlight object.
#' 
#' @export 
#' @rdname strata-methods
#' @aliases strata,genind-method strata,genlight-method
#' @param x a genind or genlight object
#' @param formula a nested formula indicating the order of the population
#' strata.
#' @param combine if \code{TRUE} (default), the levels will be combined according to the
#' formula argument. If it is \code{FALSE}, the levels will not be combined.
#' @param value a data frame OR vector OR formula (see details).
#' @docType methods
#'   
#' @details \subsection{Function Specifics}{ \itemize{ \item \strong{strata()} -
#'   Use this function to view or define population stratification of a 
#'   \linkS4class{genind} or \linkS4class{genlight} object. \item 
#'   \strong{nameStrata()} - View or rename the different levels of strata. 
#'   \item \strong{splitStrata()} - Split strata that are combined with a common
#'   separator. This function should only be used once during a workflow. 
#'   \itemize{ \item \emph{Rationale:} It is often difficult to import files 
#'   with several levels of strata as most data formats do not allow unlimited 
#'   population levels. This is circumvented by collapsing all population strata
#'   into a single population factor with a common separator for each 
#'   observation. } \item \strong{addStrata()} - Add levels to your population 
#'   strata. This is ideal for adding groups defined by 
#'   \code{\link{find.clusters}}. You can input a data frame or a vector, but if
#'   you put in a vector, you have the option to name it. }}
#'   
#'   \subsection{Argument Specifics}{
#'   
#'   These functions allow the user to seamlessly carry all possible population 
#'   stratification with their \linkS4class{genind} or \linkS4class{genlight} 
#'   object. Note that there are two ways of performing all methods: \itemize{ 
#'   \item modifying: \code{strata(myData) <- myStrata} \item preserving: 
#'   \code{myNewData <- strata(myData, value = myStrata)} } They essentially do 
#'   the same thing except that the modifying assignment method (the one with 
#'   the "\code{<-}") will modify the object in place whereas the non-assignment
#'   method will preserve the original object (unless you overwrite it). Due to 
#'   convention, everything right of the assignment is termed \code{value}. To 
#'   avoid confusion, here is a guide to the argument \strong{\code{value}} for 
#'   each function: \itemize{ \item \strong{strata()} \code{value = }a 
#'   \code{\link{data.frame}} that defines the strata for each individual in the
#'   rows. \item \strong{nameStrata()} \code{value = }a \code{\link{vector}} or 
#'   a \code{\link{formula}} that will define the names. \item 
#'   \strong{splitStrata()} \code{value = }a \code{\link{formula}} argument with
#'   the same number of levels as the strata you wish to split. \item 
#'   \strong{addStrata()} \code{value = }a \code{\link{vector}} or 
#'   \code{\link{data.frame}} with the same length as the number of individuals 
#'   in your data. }}
#'   
#'   \subsection{Details on Formulas}{
#'   
#'   The preferred use of these functions is with a \code{\link{formula}} 
#'   object. Specifically, a hierarchical formula argument is used to assign the
#'   levels of the strata. An example of a hierarchical formula would
#'   be:\tabular{r}{ \code{~Country/City/Neighborhood}} This convention was
#'   chosen as it becomes easier to type and makes intuitive sense when defining
#'   a hierarchy. Note: it is important to use hiearchical
#'   formulas when specifying hierarchies as other types of formulas (eg. 
#'   \code{~Country*City*Neighborhood}) will give incorrect results.}
#'   
#' @seealso \code{\link{setPop}} \code{\link{genind}}
#'   \code{\link{as.genind}}
#'   
#' @author Zhian N. Kamvar
#' @examples
#' # let's look at the microbov data set:
#' data(microbov)
#' microbov
#' 
#' # We see that we have three vectors of different names in the 'other' slot. 
#' # ?microbov
#' # These are Country, Breed, and Species
#' names(other(microbov))
#' 
#' # Let's set the strata
#' strata(microbov) <- data.frame(other(microbov))
#' microbov
#' 
#' # And change the names so we know what they are
#' nameStrata(microbov) <- ~Country/Breed/Species
#' 
#' \dontrun{
#' # let's see what the strata looks like by Species and Breed:
#' head(strata(microbov, ~Breed/Species))
#' 
#' # If we didn't want the last column combined with the first, we can set
#' # combine = FALSE
#' head(strata(microbov, ~Breed/Species, combine = FALSE))
#' 
#' #### USING splitStrata ####
#' 
#' # For the sake of example, we'll imagine that we have imported our data set
#' # with all of the stratifications combined. 
#' setPop(microbov) <- ~Country/Breed/Species
#' strata(microbov) <- NULL
#' 
#' # This is what our data would look like after import.
#' microbov
#' 
#' # To set our strata here, we need to use the functions strata and splitStrata
#' strata(microbov) <- data.frame(x = pop(microbov))
#' microbov # shows us that we have "one" level of stratification
#' head(strata(microbov)) # all strata are separated by "_"
#' 
#' splitStrata(microbov) <- ~Country/Breed/Species
#' microbov # Now we have all of our strata named and split
#' head(strata(microbov)) # all strata are appropriately named and split.
#' }
#==============================================================================#
strata <- function(x, formula = NULL, combine = TRUE, value){
  standardGeneric("strata")
} 

#' @export
setGeneric("strata")

setMethod(
  f = "strata",
  signature(x = "genind"),
  definition = function(x, formula = NULL, combine = TRUE, value){
    theCall <- match.call()
    if (missing(value)){
      .getStrata(x, formula = formula, combine = combine, theCall)  
    } else {
      .setStrata(x, value, theCall)
    }
  })

setMethod(
  f = "strata",
  signature(x = "genlight"),
  definition = function(x, formula = NULL, combine = TRUE, value){
    theCall <- match.call()
    if (missing(value)){
      .getStrata(x, formula = formula, combine = combine, theCall)  
    } else {
      .setStrata(x, value, theCall)
    }
    
  })


#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases strata<-,genind-method strata<-,genlight-method
#' @docType methods
#==============================================================================#
"strata<-" <- function(x, value){
  standardGeneric("strata<-")
}  

#' @export
setGeneric("strata<-")

setMethod(
  f = "strata<-",
  signature(x = "genind"),
  definition = function(x, value){
    theCall <- match.call()
    return(.setStrata(x, value, theCall))
  })

setMethod(
  f = "strata<-",
  signature(x = "genlight"),
  definition = function(x, value){
    theCall <- match.call()
    return(.setStrata(x, value, theCall))
  })

#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases nameStrata,genind-method nameStrata,genlight-method
#' @docType methods
#==============================================================================#
nameStrata <- function(x, value){
  standardGeneric("nameStrata")
}  

#' @export
setGeneric("nameStrata")

setMethod(
  f = "nameStrata",
  signature(x = "genind"),
  definition = function(x, value){
    .nameStrata(x, value)
  })

setMethod(
  f = "nameStrata",
  signature(x = "genlight"),
  definition = function(x, value){
    .nameStrata(x, value)
  })

#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases nameStrata<-,genind-method nameStrata<-,genlight-method
#' @docType methods
#==============================================================================#
"nameStrata<-" <- function(x, value){
  standardGeneric("nameStrata<-")
}  

#' @export
setGeneric("nameStrata<-")

setMethod(
  f = "nameStrata<-",
  signature(x = "genind"),
  definition = function(x, value){
    return(nameStrata(x, value))
  })

setMethod(
  f = "nameStrata<-",
  signature(x = "genlight"),
  definition = function(x, value){
    return(nameStrata(x, value))
  })
#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases splitStrata,genind-method splitStrata,genlight-method
#' @docType methods
#' @param sep a \code{character} indicating the character used to separate
#' hierarchical levels. This defaults to "_".
#' @importFrom reshape2 colsplit
#==============================================================================#
splitStrata <- function(x, value, sep = "_"){
  standardGeneric("splitStrata")
}  

#' @export
setGeneric("splitStrata")

setMethod(
  f = "splitStrata",
  signature(x = "genind"),
  definition = function(x, value, sep = "_"){
    .splitStrata(x, value, sep = sep) 
  })

setMethod(
  f = "splitStrata",
  signature(x = "genlight"),
  definition = function(x, value, sep = "_"){
    .splitStrata(x, value, sep = sep) 
  })

#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases splitStrata<-,genind-method splitStrata<-,genlight-method
#' @docType methods
#==============================================================================#
"splitStrata<-" <- function(x, sep = "_", value){
  standardGeneric("splitStrata<-")
}  

#' @export
setGeneric("splitStrata<-")

setMethod(
  f = "splitStrata<-",
  signature(x = "genind"),
  definition = function(x, sep = "_", value){
    return(splitStrata(x, value, sep))
  })

setMethod(
  f = "splitStrata<-",
  signature(x = "genlight"),
  definition = function(x, sep = "_", value){
    return(splitStrata(x, value, sep))
  })

#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases addStrata,genind-method addStrata,genlight-method
#' @param name an optional name argument for use with addStrata if supplying
#'   a vector. Defaults to "NEW".
#' @docType methods
#==============================================================================#
addStrata <- function(x, value, name = "NEW"){
  standardGeneric("addStrata")
}  

#' @export
setGeneric("addStrata")

setMethod(
  f = "addStrata",
  signature(x = "genind"),
  definition = function(x, value, name = "NEW"){
    theCall <- match.call()
    .addStrata(x, value, name = name, theCall)
  })

setMethod(
  f = "addStrata",
  signature(x = "genlight"),
  definition = function(x, value, name = "NEW"){
    theCall <- match.call()
    .addStrata(x, value, name = name, theCall)
  })

#==============================================================================#
#' @export 
#' @rdname strata-methods
#' @aliases addStrata<-,genind-method addStrata<-,genlight-method
#' @docType methods
#==============================================================================#
"addStrata<-" <- function(x, name = "NEW", value){
  standardGeneric("addStrata<-")
}  

#' @export
setGeneric("addStrata<-")

setMethod(
  f = "addStrata<-",
  signature(x = "genind"),
  definition = function(x, name = "NEW", value){
    theCall <- match.call()
    return(.addStrata(x, value, name = name, theCall))
  })

setMethod(
  f = "addStrata<-",
  signature(x = "genlight"),
  definition = function(x, name = "NEW", value){
    theCall <- match.call()
    return(.addStrata(x, value, name = name, theCall))
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
#' strata.
#' @param value same as formula
#' @aliases setPop,genind-method setPop,genlight-method
#' @docType methods 
#' @author Zhian N. Kamvar
#' @examples
#' 
#' data(microbov)
#' 
#' strata(microbov) <- data.frame(other(microbov))
#' 
#' # Currently set on just 
#' head(pop(microbov)) 
#' 
#' # setting the strata to both Pop and Subpop
#' setPop(microbov) <- ~coun/breed 
#' head(pop(microbov))
#' 
#' \dontrun{
#' 
#' # Can be used to create objects as well.
#' microbov.old <- setPop(microbov, ~spe) 
#' head(pop(microbov.old))
#' }
#==============================================================================#
setPop <- function(x, formula = NULL) standardGeneric("setPop")

#' @export
setGeneric("setPop")

setMethod(
  f = "setPop",
  signature(x = "genind"),
  definition = function(x, formula = NULL){
    .setPop(x, formula = formula)
  })

setMethod(
  f = "setPop",
  signature(x = "genlight"),
  definition = function(x, formula = NULL){
    .setPop(x, formula = formula)
  })
#==============================================================================#
#' @export
#' @rdname population-methods
#' @aliases setPop<-,genind-method setPop<-,genlight-method
#' @docType methods
#==============================================================================#
"setPop<-" <- function(x, value) standardGeneric("setPop<-")

#' @export
setGeneric("setPop<-")

setMethod(
  f = "setPop<-",
  signature(x = "genind"),
  definition = function(x, value){
    return(setPop(x, value))
  })

setMethod(
  f = "setPop<-",
  signature(x = "genlight"),
  definition = function(x, value){
    return(setPop(x, value))
  })

#==============================================================================#
#==============================================================================#
# internal functions utilized
#==============================================================================#
#==============================================================================#

#==============================================================================#
# A function for creating a population strata using a formula and data frame
# 
# hier = a nested formula such as ~ A/B/C where C is nested within B, which is
# nested within A.
#
# df = a data frame containing columns corresponding to the variables in hier.
#
# example:
# df <- data.frame(list(a = letters, b = LETTERS, c = 1:26))
# newdf <- .make_strata(~ a/b/c, df)
# df[names(newdf)] <- newdf # Add new columns.
#
# Public functions utilizing this function:
#
# # setPop, strata
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
.make_strata <- function(strata, df, expand_label = FALSE){
  newlevs <- attr(terms(strata), "term.labels")
  levs <- all.vars(strata)
  if (length(levs) > 1){
    newlevs <- gsub(":", "_", newlevs)
  }
  if (!all(levs %in% names(df))){
    stop(.strata_incompatible_warning(levs, df), call. = FALSE)
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
# maintaining strata. 
# Public functions utilizing this function:
# # 
#
# Internal functions utilizing this function:
# # .make_strata
#==============================================================================#

.pop_combiner <- function(df, strata=c(1), sep="_"){
  if(!is.list(df)){
    warning("df must be a data frame or a list", call. = FALSE)
    return(df)
  }
  else{
    if(length(strata)==1){
      return(df[[strata]])
    }
    else{
      comb <- vector(length=length(df[[strata[1]]]))
      comb <- df[[strata[1]]]
      lapply(strata[-1], function(x) comb <<- paste(comb, df[[x]], sep=sep))
      return(comb)
    }
  }
}

#==============================================================================#
# A function that will quit the function if a level in the strata is not
# present in the given data frame.
#
# Public functions utilizing this function:
# # setPop strata poppr.amova
#
# Internal functions utilizing this function:
# # .make_strata make_ade_df
#==============================================================================#
.strata_incompatible_warning <- function(levs, df){
  msg <- paste("One or more levels in the given strata is not present", 
               "in the data frame.\n", paste(rep("-", 78), collapse = ""),
               "\nstrata:\t", paste(levs, collapse = ", "), "\nData:\t", 
               paste(names(df), collapse = ", "))
  return(msg)
}

#==============================================================================#
# A function that will send a message if the number of rows in a data frame does
# not match what is expected.
#
# Public functions utilizing this function:
# # addStrata
#
# Internal functions utilizing this function:
# # .addStrata
#==============================================================================#
.not_enough_rows_warning <- function(dfname, dfrow, nind){
  msg <- paste("The data frame or vector does not have enough rows\n", 
               paste(rep("-", 78), collapse = ""),
               "\n", dfname, ":\t", dfrow, "\nData :\t", nind)
  return(msg)
}
