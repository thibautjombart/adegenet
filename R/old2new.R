##################
# Methods old2new
##################
#' Convert objects with obsolete classes into new objects 
#'  
#' The genind and genlight objects have changed in Adegenet version 2.0. They 
#' have each gained strata and hierarchy slots. What's more is that the genind 
#' objects have been optimized for storage and now store the tab slot as
#' integers instead of numerics. This function will convert old genind or
#' genlight objects to new ones seamlessly.
#'  
#' @name old2new 
#' @rdname old2new
#' @aliases old2new old2new,ANY-method old2new-methods old2new,genind-method 
#' old2new,genpop-method old2new,genlight-method
#' @docType methods 
#' @param object a genind or genlight object from version 1.4 or earlier.
#' 
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}\cr
#' Zhian N. Kamvar \email{kamvarz@@science.oregonstate.edu}
#' @keywords methods classes manip 
#' @export
old2new <- function(object) standardGeneric("old2new")

#' @export
setGeneric("old2new")

setMethod("old2new", "genind", function(object){
  object@hierarchy <- NULL
  object@strata    <- NULL
  if (object@type == "codom"){
    
    if (!is.integer(object@tab) && length(object@ploidy == 1)){
      xtab   <- object@tab
      newtab <- as.integer(object@tab * object@ploidy)
      object@tab <- matrix(newtab, nrow = nrow(xtab), ncol = ncol(xtab),
                           dimnames = dimnames(xtab))
      object@ploidy <- rep(object@ploidy, nrow(object@tab))  
    }
    
    names(object@loc.names) <- NULL
    rownames(object@tab)    <- object@ind.names
    names(object@all.names) <- object@loc.names
    object@all.names     <- lapply(object@all.names, setNames, NULL)
    colnames(object@tab) <- unlist(lapply(object@loc.names, 
                                          function(i) paste(i, object@all.names[[i]], 
                                                            sep = ".")), 
                                   use.names = FALSE)
    levels(object@loc.fac) <- object@loc.names
    names(object@loc.n.all) <- object@loc.names
  }
  
  if (!is.null(object@pop)){
    names(object@pop.names) <- NULL
    levels(object@pop) <- object@pop.names
  }
  object@call <- match.call()
  return(object)
  
#   res@tab <- as.matrix(x$tab)
#   res@ind.names <- as.character(x$ind.names)
#   res@loc.names <- as.character(x$loc.names)
#   res@loc.n.all <- as.integer(x$loc.n.all)
#   res@loc.fac <- as.factor(x$loc.fac)
#   res@all.names <- as.list(x$all.names)
#   if(!is.null(x$pop)) {
#       res@pop <- as.factor(x$pop)
#       theoLength <- theoLength + 1
#   }
#   if(!is.null(x$pop.names)) {
#       res@pop.names <- as.character(x$pop.names)
#       theoLength <- theoLength + 1
#   }
#   res@call <- match.call()
#   res@ploidy <- as.integer(2)
#   res@type <- "codom"
# 
#   if(length(object) > theoLength) warning("optional content else than pop and pop.names was not converted")

  # return(res)
})

setMethod("old2new", "genlight", function(object){
  object@strata    <- NULL
  object@hierarchy <- NULL
  return(object)
})

setMethod("old2new", "genpop", function(object){
  
  ## names(object@pop.names) <- NULL
  ## names(object@loc.names) <- NULL
  
  rownames(object@tab) <- object@pop.names
  colnames(object@tab) <- unlist(lapply(object@loc.names, 
                                        function(i) paste(i, object@all.names[[i]], 
                                                          sep = ".")), 
                                 use.names = FALSE)
  names(object@all.names) <- object@loc.names
  object@all.names     <- lapply(object@all.names, setNames, NULL)
  object@call <- match.call()
  return(object)
  #   x <- object
#   res <- new("genpop")
# 
#   res@tab <- as.matrix(x$tab)
#   res@pop.names <- as.character(x$pop.names)
#   res@loc.names <- as.character(x$loc.names)
#   res@loc.n.all <- as.integer(x$loc.n.all)
#   res@loc.fac <- as.factor(x$loc.fac)
#   res@all.names <- as.list(x$all.names)
#   res@ploidy <- as.integer(2)
#   res@type <- "codom"
# 
#   res@call <- match.call()
# 
#   if(length(object)>7) warning("optional content was not converted")
# 
#   return(res)
})
