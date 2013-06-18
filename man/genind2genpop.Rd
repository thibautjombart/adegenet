\encoding{UTF-8}
\name{genind2genpop}
\alias{genind2genpop}
\title{Conversion from a genind to a genpop object}
\description{The function \code{genind2genpop} converts genotypes data
  (genind) into alleles counts per population (genpop).}
\usage{
genind2genpop(x,pop=NULL,missing=c("NA","0","chi2"),quiet=FALSE,
              process.other=FALSE, other.action=mean)
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    note provided, seeked in x@pop, but if given, the argument
    prevails on x@pop.}
  \item{missing}{can be "NA", "0", or "chi2". See details for more information.}
  \item{quiet}{logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
  \item{process.other}{a logical indicating whether the \code{@other}
    slot should be processed (see details).}
  \item{other.action}{a function to be used when processing the
    \code{@other} slot. By default, 'mean' is used.}
}
\details{
  === 'missing' argument ===\cr
  The values of the 'missing' argument in \code{genind2genpop}
  have the following effects:\cr
  - "NA": if all genotypes of a population for a given allele are
  missing, count value will be NA\cr
  
  - "0": if all genotypes of a population for a given allele are
  missing, count value will be 0\cr
  
  - "chi2": if all genotypes of a population for a given allele are
  missing, count value will be that of a theoretical count in of a
  Chi-squared test. This is obtained by the product of the margins
  sums divided by the total number of alleles.\cr
  
  === processing the \code{@other} slot ===\cr
  Essentially, \code{genind2genpop} is about aggregating data per
  population. The function can do the same for all numeric items in the
  \code{@other} slot provided they have the same length (for vectors) or
  the same number of rows (matrix-like objects) as the number of
  genotypes. When the case is encountered and if \code{process.other} is
  TRUE, then these objects are processed using the function defined in
  \code{other.action} per population. For instance, spatial coordinates
  of genotypes would be averaged to obtain population coordinates.
}
\value{
  A genpop object. The component @other in 'x' is passed to the created
  genpop object.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{
  \linkS4class{genind}, \linkS4class{genpop}, \code{\link{na.replace}}
}
\examples{
## simple conversion
data(nancycats)
nancycats
catpop <- genind2genpop(nancycats)
catpop
summary(catpop)

## processing the @other slot
data(sim2pop)
sim2pop$other$foo <- letters
sim2pop
dim(sim2pop$other$xy) # matches the number of genotypes
sim2pop$other$foo # does not match the number of genotypes

obj <- genind2genpop(sim2pop, process.other=TRUE)
obj$other # the new xy is the populations' centre

pch <- as.numeric(pop(sim2pop))
col <- pop(sim2pop)
levels(col) <- c("blue","red")
col <- as.character(col)
plot(sim2pop$other$xy, pch=pch, col=col)
text(obj$other$xy, lab=row.names(obj$other$xy), col=c("blue","red"), cex=2, font=2)

}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}
