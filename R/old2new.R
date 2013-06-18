##################
# Methods old2new
##################
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
