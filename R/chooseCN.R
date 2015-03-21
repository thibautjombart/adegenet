#####################
# Function chooseCN
#####################


#' Function to choose a connection network
#'
#' The function \code{chooseCN} is a simple interface to build a connection
#' network (CN) from xy coordinates. The user chooses from 6 types of graph and
#' one additional weighting scheme.  \code{chooseCN} calls functions from
#' appropriate packages, handles non-unique coordinates and returns a
#' connection network either with classe \code{nb} or \code{listw}. For graph
#' types 1-4, duplicated locations are not accepted and will issue an error.
#'
#' There are 7 kinds of graphs proposed: \cr Delaunay triangulation (type 1)\cr
#' Gabriel graph (type 2)\cr Relative neighbours (type 3)\cr Minimum spanning
#' tree (type 4)\cr Neighbourhood by distance (type 5)\cr K nearests neighbours
#' (type 6)\cr Inverse distances (type 7)\cr
#'
#' The last option (type=7) is not a true neighbouring graph: all sites are
#' neighbours, but the spatial weights are directly proportional to the
#' inversed spatial distances.\cr Also not that in this case, the output of the
#' function is always a \code{listw} object, even if \code{nb} was
#' requested.\cr
#'
#' The choice of the connection network has been discuted on the adegenet
#' forum. Please search the archives from adegenet website (section 'contact')
#' using 'graph' as keyword.
#'
#' @param xy an matrix or data.frame with two columns for x and y coordinates.
#' @param ask a logical stating whether graph should be chosen interactively
#' (TRUE,default) or not (FALSE). Set to FALSE if \code{type} is provided.
#' @param type an integer giving the type of graph (see details).
#' @param result.type a character giving the class of the returned object.
#' Either "nb" (default) or "listw", both from \code{spdep} package. See
#' details.
#' @param d1 the minimum distance between any two neighbours. Used if
#' \code{type=5.}
#' @param d2 the maximum distance between any two neighbours. Used if
#' \code{type=5}. Can also be a character: "dmin" for the minimum distance so
#' that each site has at least one connection, or "dmax" to have all sites
#' connected (despite the later has no sense).
#' @param k the number of neighbours per point. Used if \code{type=6}.
#' @param a the exponent of the inverse distance matrix. Used if \code{type=7}.
#' @param dmin the minimum distance between any two distinct points. Used to
#' avoid infinite spatial proximities (defined as the inversed spatial
#' distances). Used if \code{type=7}.
#' @param plot.nb a logical stating whether the resulting graph should be
#' plotted (TRUE, default) or not (FALSE).
#' @param edit.nb a logical stating whether the resulting graph should be
#' edited manually for corrections (TRUE) or not (FALSE, default).
#' @return Returns a connection network having the class \code{nb} or
#' \code{listw}. The xy coordinates are passed as attribute to the created
#' object.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{spca}}
#' @keywords spatial utilities
#' @examples
#'
#' \dontrun{
#' data(nancycats)
#'
#' par(mfrow=c(2,2))
#' cn1 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=1)
#' cn2 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=2)
#' cn3 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=3)
#' cn4 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=4)
#' par(mfrow=c(1,1))
#' }
#'
#' @export chooseCN
#' @importFrom spdep "tri2nb" "gabrielneigh" "graph2nb" "relativeneigh" "dnearneigh" "knearneigh" "knn2nb" "nb2listw" "mat2listw" "listw2mat" "lag.listw" "card"
#' @import ade4
#'
chooseCN <- function(xy,ask=TRUE, type=NULL, result.type="nb", d1=NULL, d2=NULL, k=NULL,
                     a=NULL, dmin=NULL, plot.nb=TRUE, edit.nb=FALSE){

  if(is.data.frame(xy)) xy <- as.matrix(xy)
  if(ncol(xy) != 2) stop("xy does not have two columns.")
  if(any(is.na(xy))) stop("NA entries in xy.")
  result.type <- tolower(result.type)
   if(is.null(type) & !ask) stop("Non-interactive mode but no graph chosen; please provide a value for 'type' argument.")

  ## if(!require(spdep, quietly=TRUE)) stop("spdep library is required.")

  res <- list()

  if(!is.null(d2)){
      if(d2=="dmin"){
          tempmat <- as.matrix(dist(xy))
          d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
          d2min <- d2min * 1.0001 # to avoid exact number problem
          d2 <- d2min
      } else if(d2=="dmax"){
          d2max <- max(dist(xy))
          d2max <- d2max * 1.0001 # to avoid exact number problem
          d2 <- d2max
      }
  } # end handle d2

  d1.first <- d1
  d2.first <- d2
  k.first <- k

  ## handle type argument
  if(!is.null(type)){
      type <- as.integer(type)
      if(type < 1 |type > 7) stop("type must be between 1 and 7")
      ask <- FALSE
  }

  ## check for uniqueness of coordinates
  if(any(xyTable(xy)$number>1)){ # if duplicate coords
      DUPLICATE.XY <- TRUE
  } else {
      DUPLICATE.XY <- FALSE
  }


  ## if(is.null(type) & !ask) { type <- 1 }

  ### begin large while ###
  chooseAgain <- TRUE
  while(chooseAgain){
    # re-initialisation of some variables
    d1 <- d1.first
    d2 <- d2.first
    k <- k.first

  ## read type from console
    if(ask){
      temp <- TRUE
      while(temp){
        cat("\nChoose a connection network:\n")
        cat("\t Delaunay triangulation (type 1)\n")
        cat("\t Gabriel graph (type 2)\n")
        cat("\t Relative neighbours (type 3)\n")
        cat("\t Minimum spanning tree (type 4)\n")
        cat("\t Neighbourhood by distance (type 5)\n")
        cat("\t K nearest neighbours (type 6)\n")
        cat("\t Inverse distances (type 7)\n")
        cat("Answer: ")

        type <- as.integer(readLines(n = 1))
        temp <- type < 1 |type > 7
        if(temp) cat("\nWrong answer\n")

        if(type %in% 1:4 & DUPLICATE.XY){
            cat("\n\n== PROBLEM DETECTED ==")
            cat("\nDuplicate locations detected\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).\n")
            temp <- TRUE
        }

      } # end while
    }
    ##

    ## warning about duplicate xy coords
    if(type %in% 1:4 & DUPLICATE.XY){
        stop("Duplicate locations detected and incompatible with graph type 1-4.\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).")
    }

    ## graph types
    ## type 1: Delaunay
    if(type==1){
      ## if(!require(tripack, quietly=TRUE)) stop("tripack library is required.")
      cn <- tri2nb(xy)
    }

    # type 2: Gabriel
    if(type==2){
      cn <- gabrielneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }

    ## type 3: Relative neighbours
    if(type==3){
      cn <- relativeneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }

    ## type 4: Minimum spanning tree
    if(type==4){
      cn <- ade4::mstree(dist(xy)) # there is also a spdep::mstree
      cn <- neig2nb(cn)
    }

    ## type 5: Neighbourhood by distance
    if(type==5){
      if(is.null(d1) |is.null(d2)){
        tempmat <- as.matrix(dist(xy))
        d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
        d2min <- d2min * 1.0001 # to avoid exact number problem
        d2max <- max(dist(xy))
        d2max <- d2max * 1.0001 # to avoid exact number problem
        dig <- options("digits")
        options("digits=5")
        cat("\n Enter minimum distance: ")
        d1 <- as.numeric(readLines(n = 1))
        cat("\n Enter maximum distance \n(dmin=", d2min, ", dmax=", d2max, "): ")
        d2 <- readLines(n = 1)
        ## handle character
        if(d2=="dmin") {
            d2 <- d2min
        } else if(d2=="dmax") {
            d2 <- d2max
        } else {
            d2 <- as.numeric(d2)
        }
        ## restore initial digit option
        options(dig)
      }
    # avoid that a point is its neighbour
      dmin <- mean(dist(xy))/100000
      if(d1<dmin) d1 <- dmin
      if(d2<d1) stop("d2 < d1")
      cn <- dnearneigh(x=xy, d1=d1, d2=d2)
    }

    ## type 6: K nearests
    if(type==6){
      if(is.null(k)) {
        cat("\n Enter the number of neighbours: ")
        k <- as.numeric(readLines(n = 1))
      }
      cn <- knearneigh(x=xy, k=k)
      cn <- knn2nb(cn, sym=TRUE)
    }

    ## type 7: inverse distances
    if(type==7){
        if(is.null(a)) {
            cat("\n Enter the exponent: ")
            a <- as.numeric(readLines(n = 1))
        }
        cn <- as.matrix(dist(xy))
        if(is.null(dmin)) {
            cat("\n Enter the minimum distance \n(range = 0 -", max(cn),"): ")
            dmin <- as.numeric(readLines(n = 1))
        }
        if(a<1) { a <- 1 }
        thres <- mean(cn)/1e8
        if(dmin > thres) dmin <- thres
        cn[cn < dmin] <- dmin
        cn <- 1/(cn^a)
        diag(cn) <- 0
        cn <- prop.table(cn,1)
        plot.nb <- FALSE
        edit.nb <- FALSE
        result.type <- "listw"
    } # end type 7

    ## end graph types

    if(ask & plot.nb) {
      plot(cn,xy)
      cat("\nKeep this graph (y/n)? ")
    ans <- tolower(readLines(n=1))
      if(ans=="n") {chooseAgain <- TRUE} else {chooseAgain <- FALSE}
    }
    else if(plot.nb){
      plot(cn,xy)
      chooseAgain <- FALSE
    }
  else {chooseAgain <- FALSE}

  }
### end large while

  if(edit.nb) {cn <- edit(cn,xy)}

  if(result.type == "listw") {
      if(type!=7) {
          cn <- nb2listw(cn, style="W", zero.policy=TRUE)
      } else {
          cn <- mat2listw(cn)
          cn$style <- "W"
      }
  }

  res <- cn

  attr(res,"xy") <- xy

  return(res)

} # end chooseCN

