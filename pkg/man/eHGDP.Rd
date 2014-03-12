\encoding{UTF-8}
\name{eHGDP}
\alias{eHGDP}
\docType{data}
\title{Extended HGDP-CEPH dataset}
\description{
This dataset consists of 1350 individuals from native Human populations
distributed worldwide typed at 678 microsatellite loci. The original
HGDP-CEPH panel [1-3] has been extended by several native American
populations [4]. This dataset was used to illustrate the Discriminant
Analysis of Principal Components (DAPC, [5]).
}
\usage{data(eHGDP)}
\format{
    \code{eHGDP} is a genind object with a data frame named \code{popInfo} as supplementary
    component (\code{eHGDP@other$popInfo)}, which contains the following
    variables:
    \describe{
        \item{Population: }{a character vector indicating populations.}
        \item{Region: }{a character vector indicating the geographic region of each population.}
        \item{Label: }{a character vector indicating the correspondance
	  with population labels used in the genind object (i.e., as
	  output by \code{pop(eHGDP)}).}
	\item{Latitude,Longitude: }{geographic coordinates of the
	  populations, indicated as north and east degrees.}
    }
}
\source{
  Original panel by Human Genome Diversity Project (HGDP) and Centre
  d'Etude du Polymorphisme Humain (CEPH). See reference [4] for Native
  American populations.

  This copy of the dataset was prepared by Francois Balloux.
}
\references{
[1] Rosenberg NA, Pritchard JK, Weber JL, Cann HM, Kidd KK, et
al. (2002) Genetic structure of human populations. \emph{Science} 298:
2381-2385.

[2] Ramachandran S, Deshpande O, Roseman CC, Rosenberg NA, Feldman MW,
et al. (2005) Support from the relationship of genetic and geographic
distance in human populations for a serial founder effect originating in
Africa. \emph{Proc Natl Acad Sci U S A} 102: 15942-15947.

[3] Cann HM, de Toma C, Cazes L, Legrand MF, Morel V, et al. (2002) A
human genome diversity cell line panel.  \emph{Science} 296: 261-262.

[4] Wang S, Lewis CM, Jakobsson M, Ramachandran S, Ray N, et al. (2007)
Genetic Variation and Population Structure in Native
Americans. \emph{PLoS Genetics} 3: e185.

[5] Jombart, T., Devillard, S. and Balloux, F.  Discriminant analysis of
principal components: a new method for the analysis of genetically
structured populations. Submitted to \emph{BMC genetics}.
}
\examples{
\dontrun{
## LOAD DATA
data(eHGDP)
eHGDP


## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(eHGDP, all.contrib=TRUE, scale=FALSE,
n.pca=200, n.da=80) # takes 2 minutes
dapc1

## (see ?dapc for details about the output)



## SCREEPLOT OF EIGENVALUES
barplot(dapc1$eig, main="eHGDP - DAPC eigenvalues",
col=c("red","green","blue", rep("grey", 1000)))



## SCATTERPLOTS
## (!) Note: colors may be inverted with respect to [5]
## as signs of principal components are arbitrary
## and change from one computer to another
##
## axes 1-2
s.label(dapc1$grp.coord[,1:2], clab=0, sub="Axes 1-2")
par(xpd=T)
colorplot(dapc1$grp.coord[,1:2], dapc1$grp.coord, cex=3, add=TRUE)
add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)

## axes 2-3
s.label(dapc1$grp.coord[,2:3], clab=0, sub="Axes 2-3")
par(xpd=T)
colorplot(dapc1$grp.coord[,2:3], dapc1$grp.coord, cex=3, add=TRUE)
add.scatter.eig(dapc1$eig,10,1,2, posi="bottomright", ratio=.3, csub=1.25)



## MAP DAPC1 RESULTS
if(require(maps)){

xy <- cbind(eHGDP$other$popInfo$Longitude, eHGDP$other$popInfo$Latitude)

par(mar=rep(.1,4))
map(fill=TRUE, col="lightgrey")
colorplot(xy, -dapc1$grp.coord, cex=3, add=TRUE, trans=FALSE)
}



## LOOK FOR OTHER CLUSTERS
## to reproduce results of the reference paper, use :
## grp <- find.clusters(eHGDP, max.n=50, n.pca=200, scale=FALSE)
## and then
## plot(grp$Kstat, type="b", col="blue")

grp <- find.clusters(eHGDP, max.n=30, n.pca=200,
scale=FALSE, n.clust=4) # takes about 2 minutes
names(grp)

## (see ?find.clusters for details about the output)



## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc2 <- dapc(eHGDP, pop=grp$grp, all.contrib=TRUE,
scale=FALSE, n.pca=200, n.da=80) # takes around a 1 minute
dapc2


## PRODUCE SCATTERPLOT
scatter(dapc2) # axes 1-2
scatter(dapc2,2,3) # axes 2-3


## MAP DAPC2 RESULTS
if(require(maps)){
xy <- cbind(eHGDP$other$popInfo$Longitude,
eHGDP$other$popInfo$Latitude)

myCoords <- apply(dapc2$ind.coord, 2, tapply, pop(eHGDP), mean)

par(mar=rep(.1,4))
map(fill=TRUE, col="lightgrey")
colorplot(xy, myCoords, cex=3, add=TRUE, trans=FALSE)
}

}
}
\keyword{datasets}
