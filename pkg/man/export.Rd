\encoding{UTF-8}
\name{export}
\alias{genind2genotype}
\alias{genind2hierfstat}
\title{ Conversion functions from adegenet to other R packages}
\description{ The function \code{genind2genotype} and
  \code{genind2hierfstat} convert a \code{genind} object into,
  respectively, a list of genotypes (class \code{genotypes}, package
  \code{genetics}), and a data.frame to be used by the functions of the
  package \code{hierfstat}.
}
\usage{
genind2genotype(x,pop=NULL,res.type=c("matrix","list"))
genind2hierfstat(x,pop=NULL)
}
\arguments{
  \item{x}{ a \code{genind} object.}
  \item{pop}{a factor giving the population of each individual. If NULL,
  it is seeked in x\$pop. If NULL again, all individuals are assumed
  from the same population.}
\item{res.type}{a character (if a vector, only the first element is
  retained), indicating the type of result returned.}
 }
 \value{The function \code{genind2genotype} converts a \code{genind}
   object into \code{genotypes} (package \code{genetics}).\cr
   If res.type is set to "matrix" (default), the returned value is a
   individuals x locus matrix whose columns have the class
   \code{genotype}. Such data can be used by \code{LDheatmap} package to
   compute linkage disequilibrium.\cr
   
   If res.type is set to "list", the returned value is a list of
   \code{genotypes} sorted first by locus and then by population.)\cr

   \code{genind2hierfstat} returns a data frame where individuals are
   in rows. The first columns is a population factor (but stored as
   integer); each other column is a locus. Genotypes are coded as
   integers (e.g., 44 is an homozygote 4/4, 56 is an heterozygote
   5/6).\cr
}
\references{
  Gregory Warnes and Friedrich Leisch (2007). genetics: Population
  Genetics. R package version 1.2.1.

  Jerome Goudet (2005). HIERFSTAT, a package for R to compute and test
  hierarchical F-statistics. \emph{Molecular Ecology}, \bold{5}:184-186 \cr
  
  Fstat (version 2.9.3). Software by Jerome
  Goudet. http://www2.unil.ch/popgen/softwares/fstat.htm\cr
}
\seealso{\code{\link{import2genind}}}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\keyword{manip}