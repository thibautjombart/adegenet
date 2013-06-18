\encoding{UTF-8}
\name{import}
\alias{import2genind}
\alias{import2genind}
\title{ Importing data from several softwares to a genind object}
\description{
  Their are several ways to import genotype data to a
  \linkS4class{genind} object: i) from a data.frame with a given format
  (see \code{\link{df2genind}}), ii) from a file with a recognized
  extension, or iii) from an alignement of sequences (see \code{\link{DNAbin2genind}}).\cr
  
  The function \code{import2genind} detects the extension of
  the file given in argument and seeks for an appropriate import
  function to create a \code{genind} object.\cr
  Current recognized formats are :\cr
  - GENETIX files (.gtx) \cr
  - Genepop files (.gen) \cr
  - Fstat files (.dat) \cr
  - STRUCTURE files (.str or .stru) \cr
}
\usage{
import2genind(file,missing=NA,quiet=FALSE, \dots)
}
\arguments{
   \item{file}{ a character string giving the path to the file to
    convert, with the appropriate extension.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{quiet}{ logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
  \item{\dots}{other arguments passed to the appropriate 'read' function
  (currently passed to \code{read.structure})}
}
   
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA.\cr
  
  - 0: allelic frequencies are set to 0  on all alleles of the concerned
  locus. Recommended for a PCA on compositionnal data.\cr
  
  - "mean": missing values are replaced by the mean frequency of the
  corresponding allele, computed on the whole set of
  individuals. Recommended for a centred PCA.\cr

  Beware: same data in different formats are not expected to produce 
  exactly the same \code{genind} objects.\cr
  For instance, conversions made by GENETIX to Fstat may change the
  the sorting of the genotypes; GENETIX stores individual names whereas
  Fstat does not; Genepop chooses a sample's name from the name of its last
  genotype; etc.
}
\value{an object of the class \code{genind}}
\references{
  Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F. (1996-2004)
  GENETIX 4.05, logiciel sous Windows TM pour la génétique des
  populations. Laboratoire Génome, Populations, Interactions, CNRS UMR
  5000, Université de Montpellier II, Montpellier (France). \cr

  Pritchard, J.; Stephens, M. & Donnelly, P. (2000) Inference of population structure using multilocus genotype data. \emph{Genetics}, \bold{155}: 945-959
 
  Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population
  genetics software for exact tests and ecumenicism. \emph{J. Heredity},
  \bold{86}:248-249 \cr
  
  Fstat (version 2.9.3). Software by Jerome
  Goudet. http://www2.unil.ch/popgen/softwares/fstat.htm\cr

  Excoffier L. & Heckel G.(2006) Computer programs for population genetics data analysis: a survival guide \emph{Nature}, \bold{7}: 745-758
}
\seealso{\code{\link{import2genind}}, \code{\link{read.genetix}},
  \code{\link{read.fstat}}, \code{\link{read.structure}}, \code{\link{read.genepop}}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
import2genind(system.file("files/nancycats.gtx",
package="adegenet"))

import2genind(system.file("files/nancycats.dat",
package="adegenet"))

import2genind(system.file("files/nancycats.gen",
package="adegenet"))

import2genind(system.file("files/nancycats.str",
package="adegenet"), onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE)
}
\keyword{manip}
