\encoding{UTF-8}
\name{spcaIllus}
\alias{spcaIllus}
\docType{data}
\title{Simulated data illustrating the sPCA}
\description{
  Datasets illustrating the spatial Principal Component Analysis
  (Jombart et al. 2009).
  These data were simulated using various models using Easypop (2.0.1).
  Spatial coordinates were defined so that different spatial patterns
  existed in the data. The \code{spca-illus} is a list containing the
  following \linkS4class{genind} or \linkS4class{genpop} objects:\cr
  - dat2A: 2 patches \cr
  - dat2B: cline between two pop \cr
  - dat2C: repulsion among individuals from the same gene pool \cr
  - dat3: cline and repulsion \cr
  - dat4: patches and local alternance \cr

  See "source" for a reference providing simulation details.
}
\usage{data(spcaIllus)}
\format{
  \code{spcaIllus} is list of 5 components being either genind or genpop objects.
}
\source{
   Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
Revealing cryptic spatial patterns in genetic variability by a new
multivariate method. \emph{Heredity}, \bold{101}, 92--103.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\seealso{ \code{\link{spca}} }
\references{
 Jombart, T., Devillard, S., Dufour, A.-B. and Pontier, D.
Revealing cryptic spatial patterns in genetic variability by a new
multivariate method. \emph{Heredity}, \bold{101}, 92--103.

 
  Balloux F (2001) Easypop (version 1.7): a computer program for oppulation
  genetics simulations \emph{Journal of Heredity}, \bold{92}: 301-302
}
\examples{
if(require(spdep)){

data(spcaIllus)
attach(spcaIllus)
opar <- par(no.readonly=TRUE)
## comparison PCA vs sPCA

# PCA
pca2A <- dudi.pca(dat2A$tab,center=TRUE,scale=FALSE,scannf=FALSE)
pca2B <- dudi.pca(dat2B$tab,center=TRUE,scale=FALSE,scannf=FALSE)
pca2C <- dudi.pca(dat2C$tab,center=TRUE,scale=FALSE,scannf=FALSE)
pca3 <- dudi.pca(dat3$tab,center=TRUE,scale=FALSE,scannf=FALSE,nf=2)
pca4 <- dudi.pca(dat4$tab,center=TRUE,scale=FALSE,scannf=FALSE,nf=2)

# sPCA
spca2A <-spca(dat2A,xy=dat2A$other$xy,ask=FALSE,type=1,
plot=FALSE,scannf=FALSE,nfposi=1,nfnega=0)

spca2B <- spca(dat2B,xy=dat2B$other$xy,ask=FALSE,type=1,
plot=FALSE,scannf=FALSE,nfposi=1,nfnega=0)

spca2C <- spca(dat2C,xy=dat2C$other$xy,ask=FALSE,
type=1,plot=FALSE,scannf=FALSE,nfposi=0,nfnega=1)

spca3 <- spca(dat3,xy=dat3$other$xy,ask=FALSE,
type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)

spca4 <- spca(dat4,xy=dat4$other$xy,ask=FALSE,
type=1,plot=FALSE,scannf=FALSE,nfposi=1,nfnega=1)

# an auxiliary function for graphics
plotaux <- function(x,analysis,axis=1,lab=NULL,...){
neig <- NULL
if(inherits(analysis,"spca")) neig <- nb2neig(analysis$lw$neighbours)
xrange <- range(x$other$xy[,1])
xlim <- xrange + c(-diff(xrange)*.1 , diff(xrange)*.45)
yrange <- range(x$other$xy[,2])
ylim <- yrange + c(-diff(yrange)*.45 , diff(yrange)*.1)

s.value(x$other$xy,analysis$li[,axis],include.ori=FALSE,addaxes=FALSE,
cgrid=0,grid=FALSE,neig=neig,cleg=0,xlim=xlim,ylim=ylim,...)

par(mar=rep(.1,4))
if(is.null(lab)) lab = gsub("[P]","",x$pop)
text(x$other$xy, lab=lab, col="blue", cex=1.2, font=2)
add.scatter({barplot(analysis$eig,col="grey");box();
title("Eigenvalues",line=-1)},posi="bottomright",ratio=.3)
}

# plots
plotaux(dat2A,pca2A,sub="dat2A - PCA",pos="bottomleft",csub=2)
plotaux(dat2A,spca2A,sub="dat2A - sPCA glob1",pos="bottomleft",csub=2)

plotaux(dat2B,pca2B,sub="dat2B - PCA",pos="bottomleft",csub=2)
plotaux(dat2B,spca2B,sub="dat2B - sPCA glob1",pos="bottomleft",csub=2)

plotaux(dat2C,pca2C,sub="dat2C - PCA",pos="bottomleft",csub=2)
plotaux(dat2C,spca2C,sub="dat2C - sPCA loc1",pos="bottomleft",csub=2,axis=2)

par(mfrow=c(2,2))
plotaux(dat3,pca3,sub="dat3 - PCA axis1",pos="bottomleft",csub=2)
plotaux(dat3,spca3,sub="dat3 - sPCA glob1",pos="bottomleft",csub=2)
plotaux(dat3,pca3,sub="dat3 - PCA axis2",pos="bottomleft",csub=2,axis=2)
plotaux(dat3,spca3,sub="dat3 - sPCA loc1",pos="bottomleft",csub=2,axis=2)

plotaux(dat4,pca4,lab=dat4$other$sup.pop,sub="dat4 - PCA axis1",
pos="bottomleft",csub=2)
plotaux(dat4,spca4,lab=dat4$other$sup.pop,sub="dat4 - sPCA glob1",
pos="bottomleft",csub=2)
plotaux(dat4,pca4,lab=dat4$other$sup.pop,sub="dat4 - PCA axis2",
pos="bottomleft",csub=2,axis=2)
plotaux(dat4,spca4,lab=dat4$other$sup.pop,sub="dat4 - sPCA loc1",
pos="bottomleft",csub=2,axis=2)

# color plot
par(opar)
colorplot(spca3, cex=4, main="colorplot sPCA dat3")
text(spca3$xy[,1], spca3$xy[,2], dat3$pop)

colorplot(spca4, cex=4, main="colorplot sPCA dat4")
text(spca4$xy[,1], spca4$xy[,2], dat4$other$sup.pop)

# detach data
detach(spcaIllus)
}
}
\keyword{datasets}
\keyword{spatial}

