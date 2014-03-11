\encoding{UTF-8}
\name{chooseCN}
\alias{chooseCN}
\title{Function to choose a connection network}
\description{The function \code{chooseCN} is a simple interface to
  build a connection network (CN) from xy coordinates. The user chooses
  from 6 types of graph and one additional weighting scheme.
  \code{chooseCN} calls functions from appropriate packages, handles
  non-unique coordinates and returns a connection network either with
  classe \code{nb} or \code{listw}. For graph types 1-4, duplicated
  locations are not accepted and will issue an error.
}
\usage{
chooseCN(xy, ask = TRUE, type = NULL, result.type = "nb", d1 = NULL, 
    d2 = NULL, k = NULL,  a=NULL, dmin=NULL, plot.nb = TRUE, edit.nb = FALSE) 
}
\arguments{
  \item{xy}{an matrix or data.frame with two columns for x and y coordinates.}
  \item{ask}{a logical stating whether graph should be chosen
    interactively (TRUE,default) or not (FALSE). Set to FALSE if \code{type} is
  provided.}
  \item{type}{an integer giving the type of graph (see details).}
  \item{result.type}{a character giving the class of the returned
    object. Either "nb" (default) or "listw", both from \code{spdep}
    package. See details.}
  \item{d1}{the minimum distance between any two neighbours. Used if
    \code{type=5.}}
  \item{d2}{the maximum distance between any two neighbours. Used if
    \code{type=5}. Can also be a character: "dmin" for the minimum
    distance so that each site has at least one connection, or "dmax" to
  have all sites connected (despite the later has no sense).}
  \item{k}{the number of neighbours per point. Used if
    \code{type=6}.}
  \item{a}{the exponent of the inverse distance matrix. Used if
    \code{type=7}.}
  \item{dmin}{the minimum distance between any two distinct points. Used
    to avoid infinite spatial proximities (defined as the inversed
    spatial distances). Used if \code{type=7}.}
  \item{plot.nb}{a logical stating whether the resulting graph should be
    plotted (TRUE, default) or not  (FALSE).}
  \item{edit.nb}{a logical stating whether the resulting graph should be
    edited manually for corrections (TRUE) or not  (FALSE, default).}  
}
\details{There are 7 kinds of graphs proposed: \cr
  Delaunay triangulation (type 1)\cr
  Gabriel graph (type 2)\cr
  Relative neighbours (type 3)\cr
  Minimum spanning tree (type 4)\cr
  Neighbourhood by distance (type 5)\cr
  K nearests neighbours (type 6)\cr
  Inverse distances (type 7)\cr

  The last option (type=7) is not a true neighbouring graph: all sites are
  neighbours, but the spatial weights are directly proportional to the
  inversed spatial distances.\cr
  Also not that in this case, the output of the function is always a
  \code{listw} object, even if \code{nb} was requested.\cr

  The choice of the connection network has been discuted on the adegenet
  forum. Please search the archives from adegenet website (section
  'contact') using 'graph' as keyword. 
}
\value{Returns a connection network having the class \code{nb} or
  \code{listw}. The xy coordinates are passed as attribute to the
  created object.
}
\seealso{\code{\link{spca}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(nancycats)
if(require(spdep)){

par(mfrow=c(2,2))
cn1 <- chooseCN(nancycats@other$xy,ask=FALSE,type=1)
cn2 <- chooseCN(nancycats@other$xy,ask=FALSE,type=2)
cn3 <- chooseCN(nancycats@other$xy,ask=FALSE,type=3)
cn4 <- chooseCN(nancycats@other$xy,ask=FALSE,type=4)
par(mfrow=c(1,1))
}
}
}
\keyword{spatial}
\keyword{utilities}