\encoding{UTF-8}
\name{export}
%% \alias{genind2genotype}
\alias{genind2hierfstat}
\title{ Conversion functions from adegenet to other R packages}
\description{ Thes functions provide conversion to other packages. They
  include:
  \itemize{
    \item \code{genind2hierfstat}: converts \linkS4class{genind} into
    \code{hierfstat}'s original format.
  }
}
\usage{
%% genind2genotype(x,pop=NULL,res.type=c("matrix","list"))
genind2hierfstat(x,pop=NULL)
}
\arguments{
  \item{x}{ a \code{genind} object.}
  \item{pop}{a factor giving the population of each individual. If NULL,
  it is seeked in x\$pop. If NULL again, all individuals are assumed
  from the same population.}
%% \item{res.type}{a character (if a vector, only the first element is
%%   retained), indicating the type of result returned.}
 }
 \value{
   \code{genind2hierfstat} returns a data frame where individuals are
   in rows. The first columns is a population factor (but stored as
   integer); each other column is a locus. Genotypes are coded as
   integers (e.g., 44 is an homozygote 4/4, 56 is an heterozygote
   5/6).\cr
}
\references{
  Jerome Goudet (2005). HIERFSTAT, a package for R to compute and test
  hierarchical F-statistics. \emph{Molecular Ecology}, \bold{5}:184-186 \cr
  
  Fstat (version 2.9.3). Software by Jerome
  Goudet. http://www2.unil.ch/popgen/softwares/fstat.htm\cr
}
\seealso{\code{\link{import2genind}}}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\keyword{manip}