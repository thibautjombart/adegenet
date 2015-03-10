\encoding{UTF-8}
\name{SequencesToGenind}
\alias{DNAbin2genind}
\alias{alignment2genind}
\title{ Importing data from an alignement of sequences to a genind object}
\description{
  These functions take an alignement of sequences and translate SNPs
  into a \linkS4class{genind} object. Note that only polymorphic loci
  are retained.\cr

  Currently, accepted sequence formats are:\cr
  - DNAbin (ape package): function DNAbin2genind\cr
  - alignment (seqinr package): function alignment2genind\cr
}
\usage{
DNAbin2genind(x, pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100)

alignment2genind(x, pop=NULL, exp.char=c("a","t","g","c"), na.char="-",
                 polyThres=1/100)

}
\arguments{
 \item{x}{an object containing aligned sequences.}
 \item{pop}{an optional factor giving the population to which each sequence belongs.}
 \item{exp.char}{a vector of single character providing expected values;
  all other characters will be turned to NA.}
 \item{na.char}{a vector of single characters providing values that should be
   considered as NA. If not NULL, this is used instead of \code{exp.char}.}
 \item{polyThres}{the minimum frequency of a minor allele for a locus to
  be considered as polymorphic (defaults to 0.01).}
}

\value{an object of the class \linkS4class{genind}}

\seealso{\code{\link{import2genind}}, \code{\link{read.genetix}},
  \code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genepop}}, \code{\link[ape]{DNAbin}},
  \code{\link[seqinr]{as.alignment}}.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(woodmouse)
x <- DNAbin2genind(woodmouse)
x
genind2df(x)
}

if(require(seqinr)){
mase.res   <- read.alignment(file=system.file("sequences/test.mase",package="seqinr"),
format = "mase")
mase.res
x <- alignment2genind(mase.res)
x
locNames(x) # list of polymorphic sites
genind2df(x)

## look at Euclidean distances
D <- dist(truenames(x))
D

## summarise with a PCoA
pco1 <- dudi.pco(D, scannf=FALSE,nf=2)
scatter(pco1, posi="bottomright")
title("Principal Coordinate Analysis\n-based on proteic distances-")

}
}
\keyword{manip}
