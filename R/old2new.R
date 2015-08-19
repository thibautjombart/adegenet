##################
# Methods old2new
##################
.get_slot <- function(object, theslot = "strata"){
  if (.hasSlot(object, theslot)) slot(object, theslot)
  else NULL
}
#' Convert objects with obsolete classes into new objects
#'
#' The genind and genlight objects have changed in Adegenet version 2.0. They
#' have each gained strata and hierarchy slots. What's more is that the genind
#' objects have been optimized for storage and now store the tab slot as
#' integers instead of numerics. This function will convert old genind or
#' genlight objects to new ones seamlessly.
#'
#' @rdname old2new
#' @param object a genind or genlight object from version 1.4 or earlier.
#' @param donor a new object to place all the data into.
#' @aliases old2new
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}\cr
#' Zhian N. Kamvar \email{kamvarz@@science.oregonstate.edu}
#' @keywords manip
#' @export
old2new_genind <- function(object, donor = new("genind")){
  if (missing(donor) || !is(donor, "genind")){
    donor <- new("genind")
  } 
  donor_tab       <- slot(object, "tab")
  donor_loc.fac   <- slot(object, "loc.fac")
  donor_all.names <- slot(object, "all.names")
  donor_ind.names <- slot(object, "ind.names")
  donor_loc.names <- slot(object, "loc.names")
  donor_loc.n.all <- slot(object, "loc.nall")
  donor_all.names <- slot(object, "all.names")
  donor_ploidy    <- slot(object, "ploidy")
  donor_type      <- slot(object, "type")
  donor_other     <- slot(object, "other")
  donor_call      <- slot(object, "call")
  donor_pop       <- slot(object, "pop")
  donor_pop.names <- slot(object, "pop.names")
  
  if (donor_type == "codom"){
    
    if (!is.integer(donor_tab) && length(donor_ploidy == 1)){
      xtab   <- donor_tab
      newtab <- as.integer(donor_tab * donor_ploidy)
      donor_tab <- matrix(newtab, nrow = nrow(xtab), ncol = ncol(xtab),
                          dimnames = dimnames(xtab))
      donor_ploidy <- rep(donor_ploidy, nrow(donor_tab))
    }
    
    names(donor_loc.names) <- NULL
    rownames(donor_tab)    <- donor_ind.names
    names(donor_all.names) <- donor_loc.names
    donor_all.names     <- lapply(donor_all.names, setNames, NULL)
    colnames(donor_tab) <- unlist(lapply(donor_loc.names,
                                         function(i) paste(i, donor_all.names[[i]],
                                                           sep = ".")),
                                  use.names = FALSE)
    levels(donor_loc.fac) <- donor_loc.names
    names(donor_loc.n.all) <- donor_loc.names
  }
  
  if (!is.null(donor_pop)){
    levels(donor_pop) <- donor_pop.names
  }
  donor_call <- match.call()
  donor@tab       <- donor_tab       
  donor@loc.fac   <- donor_loc.fac   
  donor@all.names <- donor_all.names 
  donor@loc.n.all <- donor_loc.n.all 
  donor@all.names <- donor_all.names 
  donor@ploidy    <- donor_ploidy    
  donor@type      <- donor_type      
  donor@other     <- donor_other     
  donor@call      <- donor_call      
  donor@pop       <- donor_pop       
  donor@strata    <- .get_slot(object, "strata")
  donor@hierarchy <- .get_slot(object, "hierarchy")
  return(donor)
}



#' @rdname old2new
#' @export
old2new_genlight <- function(object, donor = new("genlight")){
  if (missing(donor) || !is(donor, "genlight")){
    donor <- new("genlight")
  } 
  slot(donor, "gen") <- slot(object, "gen") 
  slot(donor, "n.loc") <- slot(object, "n.loc") 
  slot(donor, "ind.names") <- slot(object, "ind.names") 
  slot(donor, "loc.names") <- slot(object, "loc.names") 
  slot(donor, "loc.all") <- slot(object, "loc.all") 
  slot(donor, "chromosome") <- slot(object, "chromosome") 
  slot(donor, "position") <- slot(object, "position") 
  slot(donor, "ploidy") <- slot(object, "ploidy") 
  slot(donor, "pop") <- slot(object, "pop")
  slot(donor, "other") <- slot(object, "other") 
  donor@strata    <- .get_slot(object, "strata")
  donor@hierarchy <- .get_slot(object, "hierarchy")
  return(donor)
}

#' @rdname old2new
#' @export
old2new_genpop <- function(object, donor = new("genpop")){

  ## names(donor_pop.names) <- NULL
  ## names(donor_loc.names) <- NULL
  donor_tab       <- object@tab
  donor_loc.fac   <- object@loc.fac
  donor_all.names <- object@all.names
  donor_loc.n.all <- object@loc.nall
  donor_all.names <- object@all.names
  donor_ploidy    <- object@ploidy
  donor_type      <- object@type
  donor_other     <- object@other
  donor_call      <- object@call
  
  donor_pop.names <- object@pop.names
  donor_loc.names <- object@loc.names
  
  rownames(donor_tab) <- donor_pop.names
  colnames(donor_tab) <- unlist(lapply(donor_loc.names,
                                        function(i) paste(i, donor_all.names[[i]],
                                                          sep = ".")),
                                 use.names = FALSE)
  names(donor_all.names) <- donor_loc.names
  donor_all.names     <- lapply(donor_all.names, setNames, NULL)
  donor_call <- match.call()
  
  donor@tab       <- donor_tab
  donor@loc.fac   <- donor_loc.fac
  donor@all.names <- donor_all.names
  donor@loc.n.all <- donor_loc.n.all
  donor@all.names <- donor_all.names
  donor@ploidy    <- donor_ploidy
  donor@type      <- donor_type
  donor@other     <- donor_other
  donor@call      <- donor_call
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
}
