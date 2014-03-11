\name{scaleGen-methods}
\docType{methods}
\alias{scaleGen}
\alias{scaleGen-methods}
\alias{scaleGen,genind-method}
\alias{scaleGen,genpop-method}
\title{ Compute scaled allele frequencies }
\description{
  The generic function \code{scaleGen} is an analogue to the
  \code{scale} function, but is designed with further arguments giving
  scaling options.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
  objects.
  Both return data.frames of scaled allele frequencies.
}
\usage{
\S4method{scaleGen}{genind}(x, center=TRUE, scale=TRUE,
         missing=c("NA","0","mean"),truenames=TRUE)
\S4method{scaleGen}{genpop}(x, center=TRUE, scale=TRUE,
         missing=c("NA","0","mean"),truenames=TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{center}{a logical stating whether alleles frequencies should be
    centred to mean zero (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{scale}{a logical stating whether alleles frequencies should be
    scaled (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
 \item{missing}{a character giving the treatment for missing values. Can be "NA", "0" or "mean"}
 }
 \value{
   A matrix of scaled allele frequencies with genotypes
  (\linkS4class{genind}) or populations in (\linkS4class{genpop}) in
  rows and alleles in columns.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## load data
data(microbov)
obj <- genind2genpop(microbov)

## compare different scaling
X1 <- scaleGen(obj)
X2 <- scaleGen(obj,met="bin")

## compute PCAs
pcaObj <- dudi.pca(obj,scale=FALSE,scannf=FALSE) # pca with no scaling
pcaX1 <- dudi.pca(X1,scale=FALSE,scannf=FALSE,nf=100) # pca with usual scaling
pcaX2 <- dudi.pca(X2,scale=FALSE,scannf=FALSE,nf=100) # pca with scaling for binomial variance

## get the loadings of alleles for the two scalings
U1 <- pcaX1$c1
U2 <- pcaX2$c1


## find an optimal plane to compare loadings
## use a procustean rotation of loadings tables
pro1 <- procuste(U1,U2,nf=2)

## graphics
par(mfrow=c(2,2))
# eigenvalues
barplot(pcaObj$eig,main="Eigenvalues\n no scaling")
barplot(pcaX1$eig,main="Eigenvalues\n usual scaling")
barplot(pcaX2$eig,main="Eigenvalues\n 'binomial' scaling")
# differences between loadings of alleles
s.match(pro1$scor1,pro1$scor2,clab=0,sub="usual -> binom (procustean rotation)")

}
}
\keyword{methods}
\keyword{manip}
