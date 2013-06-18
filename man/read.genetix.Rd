\encoding{UTF-8}
\name{read.genetix}
\alias{read.genetix}
\title{ Reading data from GENETIX}
\description{
  The function \code{read.genetix} reads GENETIX data files (.gtx) and
  convert them into a \linkS4class{genind} object.

  Note: \code{read.genetix} is meant for DIPLOID DATA ONLY. Haploid
  data with the GENETIX format can be read into R using
  \code{read.table} or \code{read.csv} after removing headers and 'POP'
  lines, and then converted using \code{\link{df2genind}}.
}
\usage{
read.genetix(file=NULL,missing=NA,quiet=FALSE)
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
  Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F. (1996-2004)
  GENETIX 4.05, logiciel sous Windows TM pour la genetique des
  populations. Laboratoire Genome, Populations, Interactions, CNRS UMR
  5000, Université de Montpellier II, Montpellier (France). \cr
}
\seealso{
\code{\link{import2genind}}, \code{\link{df2genind}},
  \code{\link{read.fstat}}, \code{\link{read.structure}}, \code{\link{read.genepop}}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))
obj
}
\keyword{manip}
