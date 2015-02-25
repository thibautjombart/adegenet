\encoding{UTF-8}
\name{sim2pop}
\alias{sim2pop}
\docType{data}
\title{Simulated genotypes of two georeferenced populations}
\description{
  This simple data set was obtained by sampling two populations evolving
  in a island model, simulated using Easypop (2.0.1). See \code{source}
  for simulation details. Sample sizes were respectively 100 and 30
  genotypes. The genotypes were given spatial coordinates so that both
  populations were spatially differentiated.
}
\usage{data(sim2pop)}
\format{
  \code{sim2pop} is a genind object with a matrix of xy coordinates as
  supplementary component.
}
\source{
  Easypop version 2.0.1 was run with the following parameters:\cr
  - two diploid populations, one sex, random mating\cr
  - 1000 individuals per population\cr
  - proportion of migration: 0.002\cr
  - 20 loci\cr
  - mutation rate: 0.0001 (KAM model)\cr
  - maximum of 50 allelic states\cr
  - 1000 generations (last one taken)\cr

}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\references{
  Balloux F (2001) Easypop (version 1.7): a computer program for oppulation
  genetics simulations \emph{Journal of Heredity}, \bold{92}: 301-302
}
\examples{
\dontrun{
data(sim2pop)

if(require(hierfstat)){
## try and find the Fst
temp <- genind2hierfstat(sim2pop)
varcomp.glob(temp[,1],temp[,-1])
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
temp <- sim2pop@pop
levels(temp) <- c(17,19)
temp <- as.numeric(as.character(temp))
plot(mon1)
points(sim2pop@other$xy,pch=temp,cex=2)
legend("topright",leg=c("Pop A", "Pop B"),pch=c(17,19))
}
}
\keyword{datasets}
