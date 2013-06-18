\encoding{UTF-8}
\name{read.genepop}
\alias{read.genepop}
\title{ Reading data from Genepop }
\description{
  The function \code{read.genepop} reads Genepop data files (.gen) and
  convert them into a \linkS4class{genind} object.

  Note: \code{read.genepop} is meant for DIPLOID DATA ONLY. Haploid
  data with the Genepop format can be read into R using
  \code{read.table} or \code{read.csv} after removing headers and 'POP'
  lines, and then converted using \code{\link{df2genind}}.
}
\usage{
read.genepop(file,missing=NA,quiet=FALSE)
}
\arguments{
   \item{file}{ a character string giving the path to the file to
    convert, with the appropriate extension.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{quiet}{ logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
}
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA.\cr
  
  - 0: allelic frequencies are set to 0  on all alleles of the concerned
  locus. Recommended for a PCA on compositionnal data.\cr
  
  - "mean": missing values are replaced by the mean frequency of the
  corresponding allele, computed on the whole set of
  individuals. Recommended for a centred PCA.\cr
}
\value{an object of the class \code{genind}}
\references{
  Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population
  genetics software for exact tests and ecumenicism. \emph{J. Heredity},
  \bold{86}:248-249 \cr
}
\seealso{
  \code{\link{import2genind}}, \code{\link{df2genind}},
  \code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genetix}}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
obj <- read.genepop(system.file("files/nancycats.gen",package="adegenet"))
obj
}
\keyword{manip}
