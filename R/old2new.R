
#' Convert objects with obsolete classe into new objects
#'
#' Adegenet classes changed from S3 to S4 types starting from version 1.1-0.
#' \code{old2new} has two methods for genind and genpop objects, so that old
#' adegenet objects can be retrieved and used in recent versions.
#'
#' Optional content but \code{$pop} and \code{$pop.names} will not be
#' converted. These are to be coerced into a list and set in the \code{@@other}
#' slot of the new object.
#'
#' @name old2new
#' @aliases old2new old2new,ANY-method old2new-methods old2new,genind-method
#' old2new,genpop-method
#' @docType methods
#' @param object a genind or genpop object in S3 version, i.e. prior
#' adegenet\_1.1-0
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @keywords methods classes manip

setGeneric("old2new",  function(object) standardGeneric("old2new"))

setMethod("old2new", "genind", function(object){
  x <- object
  res <- new("genind")
  theoLength <- 7

  res@tab <- as.matrix(x$tab)
  res@ind.names <- as.character(x$ind.names)
  res@loc.names <- as.character(x$loc.names)
  res@loc.nall <- as.integer(x$loc.nall)
  res@loc.fac <- as.factor(x$loc.fac)
  res@all.names <- as.list(x$all.names)
  if(!is.null(x$pop)) {
      res@pop <- as.factor(x$pop)
      theoLength <- theoLength + 1
  }
  if(!is.null(x$pop.names)) {
      res@pop.names <- as.character(x$pop.names)
      theoLength <- theoLength + 1
  }
  res@call <- match.call()
  res@ploidy <- as.integer(2)
  res@type <- "codom"

  if(length(object) > theoLength) warning("optional content else than pop and pop.names was not converted")

  return(res)
})


setMethod("old2new", "genpop", function(object){
  x <- object
  res <- new("genpop")

  res@tab <- as.matrix(x$tab)
  res@pop.names <- as.character(x$pop.names)
  res@loc.names <- as.character(x$loc.names)
  res@loc.nall <- as.integer(x$loc.nall)
  res@loc.fac <- as.factor(x$loc.fac)
  res@all.names <- as.list(x$all.names)
  res@ploidy <- as.integer(2)
  res@type <- "codom"

  res@call <- match.call()

  if(length(object)>7) warning("optional content was not converted")

  return(res)
})
