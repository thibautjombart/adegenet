\encoding{UTF-8}
\name{nancycats}
\alias{nancycats}
\docType{data}
\title{Microsatellites genotypes of 237 cats from 17 colonies of Nancy (France)}
\description{
This data set gives the genotypes of 237 cats (Felis catus L.) for 9
microsatellites markers. The individuals are divided into 17 colonies
whose spatial coordinates are also provided.
}
\usage{data(nancycats)}
\format{
    \code{nancycats} is a genind object with spatial coordinates of the
    colonies as a supplementary components (@xy). Beware: these
    coordinates are given for the true names (stored in @pop.names) and
    not for the generic names (used in @pop).
}
\source{
Dominique Pontier (UMR CNRS 5558, University Lyon1, France)
}
\references{
  Devillard, S.; Jombart, T. & Pontier, D. Disentangling spatial and
  genetic structure of stray cat (Felis catus L.) colonies in urban
  habitat using: not all colonies are equal. submitted to \emph{Molecular Ecology}
}
\examples{
\dontrun{
data(nancycats)
nancycats

# summary's results are stored in x
x <- summary(nancycats)

# some useful graphics
barplot(x$loc.nall,ylab="Alleles numbers",main="Alleles numbers
per locus")

plot(x$pop.eff,x$pop.nall,type="n",xlab="Sample size",ylab="Number of alleles")
text(x$pop.eff,y=x$pop.nall,lab=names(x$pop.nall))

par(las=3)
barplot(table(nancycats@pop),ylab="Number of genotypes",main="Number of genotypes per colony")

# are cats structured among colonies ?
if(require(hierfstat)){

gtest <- gstat.randtest(nancycats,nsim=99)
gtest
plot(gtest)


dat <- genind2hierfstat(nancycats)

Fstat <- varcomp.glob(dat$pop,dat[,-1])
Fstat
}
}
}
\keyword{datasets}
