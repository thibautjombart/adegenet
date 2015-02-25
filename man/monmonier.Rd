\encoding{UTF-8}
\name{monmonier}
\alias{monmonier}
\alias{optimize.monmonier}
\alias{plot.monmonier}
\alias{print.monmonier}

\title{Boundary detection using Monmonier algorithm}
\description{
  The Monmonier's algorithm detects boundaries among vertices of a
  valuated graph. This is achieved by finding the path
  exhibiting the largest distances between connected vertices.\cr

  The highest distance between two connected vertices (i.e. neighbours) is
  found, giving the starting point of the path. Then, the algorithm
  seeks the highest distance between immediate neighbours, and so on
  until a threshold value is attained.
  This threshold can be chosen from the plot of sorted distances between
  connected vertices: a boundary will likely result in an abrupt decrease
  of these values.\cr
 
  When several paths are looked for, the previous paths are taken into
  account, and cannot be either crossed or redrawn. Monmonier's
  algorithm can be used to assess the boundaries between patches of
  homogeneous observations.\cr
    
  Although Monmonier algorithm was initially designed for Voronoi
  tesselation, this implementation generalizes this algorithm to different
  connection networks. The \code{optimize.monmonier} function produces a
  \code{monmonier} object by trying several starting points, and
  returning the best boundary (i.e. largest sum of local distances). This
  is designed to avoid the algorithm to be trapped by a single strong
  local difference inside an homogeneous patch.
}

\usage{
monmonier(xy, dist, cn, threshold=NULL, bd.length=NULL, nrun=1,
skip.local.diff=rep(0,nrun),scanthres=is.null(threshold), allowLoop=TRUE)

optimize.monmonier(xy, dist, cn, ntry=10, bd.length=NULL, return.best=TRUE,
display.graph=TRUE, threshold=NULL, scanthres=is.null(threshold), allowLoop=TRUE)

\method{plot}{monmonier}(x, variable=NULL,
displayed.runs=1:x$nrun, add.arrows=TRUE,
col='blue', lty=1, bwd=4, clegend=1, csize=0.7,
method=c('squaresize','greylevel'), sub='', csub=1, possub='topleft',
cneig=1, pixmap=NULL, contour=NULL, area=NULL, add.plot=FALSE, \dots)

\method{print}{monmonier}(x, \dots)
}
\arguments{
  \item{xy}{a matrix yielding the spatial coordinates of the objects, with two columns respectively giving X and Y}
  \item{dist}{an object of class \code{dist}, giving the distances between the objects}
  \item{cn}{a connection network of class \code{nb} (package \code{spdep})}
  \item{threshold}{a number giving the minimal distance between two
    neighbours crossed by the path; by default, this is the third quartile of all the distances between neighbours}
  \item{bd.length}{an optional integer giving the requested length of
    the boundaries (in number of local differences)}
  \item{nrun}{is a integer giving the number of runs of the algorithm, that is, the number of paths to search, being one by default}
  \item{skip.local.diff}{is a vector of integers, whose length is the number of paths (\code{nrun}); each integer gives the number of starting point to skip, to avoid being stuck in a local difference between two neighbours into an homogeneous patch; none are skipped by default}
  \item{scanthres}{a logical stating whether the threshold sould be
    chosen from the barplot of sorted distances between neighbours}
  \item{allowLoop}{a logical specifying whether the boundary can loop
    (TRUE, default) or not (FALSE)}
  \item{ntry}{an integer giving the number of different starting points
    tried.}
  \item{return.best}{a logical stating whether the best monmonier object
    should be returned (TRUE, default) or not (FALSE)}
   \item{display.graph}{a logical whether the scores of each try should
    be plotted (TRUE, default) or not}
  \item{x}{a monmonier object}
  \item{variable}{a variable to be plotted using \code{s.value} (package
  \code{ade4})}
  \item{displayed.runs}{an integer vector giving the rank of the paths to represent}
  \item{add.arrows}{a logical, stating whether arrows should indicate the direction of the path (TRUE) or not (FALSE, used by default)}
  \item{col}{a characters vector giving the colors to be used for each boundary; recycled is needed; 'blue' is used by default}
  \item{lty}{a characters vector giving the type of line to be used for each boundary; 1 is used by default}
  \item{bwd}{a number giving the boundary width factor, applying to every segments of the paths; 4 is used by default}
  \item{clegend}{like in \code{s.value}, the size factor of the legend if a variable is represented}
  \item{csize}{like in \code{s.value}, the size factor of the squares used to represent a variable}
  \item{method}{like in \code{s.value}, a character giving the method to be used to represent the variable, either 'squaresize' (by default) or 'greylevel'}
  \item{sub}{a string of characters giving the subtitle of the plot}
  \item{csub}{the size factor of the subtitle}
  \item{possub}{the position of the subtitle; available choices are 'topleft' (by default), 'topright', 'bottomleft', and 'bottomright'}
  \item{cneig}{the size factor of the connection network}
  \item{pixmap}{an object of the class \code{pixmap} displayed in the map background}
  \item{contour}{ a data frame with 4 columns to plot the contour of the map: each row gives a segment (x1,y1,x2,y2)}
  \item{area}{a data frame of class 'area' to plot a set of surface units in contour}
  \item{add.plot}{a logical stating whether the plot should be added to the current one (TRUE), or displayed in a new window (FALSE, by default)}
  \item{\dots}{further arguments passed to other methods}
}
\details{
The function \code{monmonier} returns a list of the class \code{monmonier}, which contains the general informations about the algorithm, and about each run.
When displayed, the width of the boundaries reflects their 'strength'.
Let a segment MN be part of the path, M being the middle of AB, N of CD.
Then the boundary width for MN is proportionnal to (d(AB)+d(CD))/2.\cr

As there is no perfect method to display graphically a quantitative
variable (see for instance the differences between the two methods of
\code{s.value}), the boundaries provided by this algorithm seem
sometimes more reliable than the boundaries our eyes perceive (or miss).\cr

}

\value{
Returns an object of class \code{monmonier}, which contains the following elements :
	\item{run1 (run2, \dots)}{for each run, a list containing a dataframe giving the path coordinates, and a vector of the distances between neighbours of the path}
	\item{nrun}{the number of runs performed, i.e. the number of boundaries in the monmonier object}
	\item{threshold}{the threshold value, minimal distance between
	  neighbours accounted for by the algorithm}
	\item{xy}{the matrix of spatial coordinates}
	\item{cn}{the connection network of class \code{nb}}
	\item{call}{the call of the function}
}
\references{
Monmonier, M. (1973) Maximum-difference barriers: an alternative numerical regionalization method. \emph{Geographic Analysis}, \bold{3}, 245--261.

Manni, F., Guerard, E. and Heyer, E. (2004) Geographic patterns of (genetic, morphologic, linguistic) variation: how barriers can be detected by "Monmonier's algorithm". \emph{Human Biology}, \bold{76}, 173--190
}

\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}

\seealso{\code{\link{spca}},\code{\link[spdep]{edit.nb}}}

\examples{
if(require(spdep)){

### non-interactive example

# est-west separation
load(system.file("files/mondata1.rda",package="adegenet"))
cn1 <- chooseCN(mondata1$xy,type=2,ask=FALSE)
mon1 <- monmonier(mondata1$xy,dist(mondata1$x1),cn1,threshold=2)
plot(mon1,mondata1$x1)
plot(mon1,mondata1$x1,met="greylevel",add.arr=FALSE,col="red",bwd=6,lty=2)

# square in the middle
load(system.file("files/mondata2.rda",package="adegenet"))
cn2 <- chooseCN(mondata2$xy,type=1,ask=FALSE)
mon2 <- monmonier(mondata2$xy,dist(mondata2$x2),cn2,threshold=2)
plot(mon2,mondata2$x2,method="greylevel",add.arr=FALSE,bwd=6,col="red",csize=.5)

### genetic data example
\dontrun{
data(sim2pop)

if(require(hierfstat)){
## try and find the Fst
fstat(sim2pop,fst=TRUE)
# Fst = 0.038
}

## run monmonier algorithm

# build connection network
gab <- chooseCN(sim2pop@other$xy,ask=FALSE,type=2)

# filter random noise 
pca1 <- dudi.pca(sim2pop@tab,scale=FALSE, scannf=FALSE, nf=1)

# run the algorithm
mon1 <- monmonier(sim2pop@other$xy,dist(pca1$l1[,1]),gab,scanthres=FALSE)

# graphical display 
plot(mon1,var=pca1$l1[,1])
temp <- sim2pop@pop
levels(temp) <- c(17,19)
temp <- as.numeric(as.character(temp))
plot(mon1)
points(sim2pop@other$xy,pch=temp,cex=2)
legend("topright",leg=c("Pop A", "Pop B"),pch=c(17,19))


### interactive example

# north-south separation
xy <- matrix(runif(120,0,10), ncol=2)
x1 <- rnorm(60)
x1[xy[,2] > 5] <- x1[xy[,2] > 5]+3
cn1 <- chooseCN(xy,type=1,ask=FALSE)
mon1 <- optimize.monmonier(xy,dist(x1)^2,cn1,ntry=10)

# graphics
plot(mon1,x1,met="greylevel",csize=.6)

# island in the middle
x2 <- rnorm(60)
sel <- (xy[,1]>3.5 & xy[,2]>3.5 & xy[,1]<6.5 & xy[,2]<6.5)
x2[sel] <- x2[sel]+4
cn2 <- chooseCN(xy,type=1,ask=FALSE)
mon2 <- optimize.monmonier(xy,dist(x2)^2,cn2,ntry=10)

# graphics
plot(mon2,x2,method="greylevel",add.arr=FALSE,bwd=6,col="red",csize=.5)
}
}
}
\keyword{multivariate}
\keyword{spatial}
