\encoding{UTF-8}
\name{read.structure}
\alias{read.structure}
\title{ Reading data from STRUCTURE}
\description{
   The function \code{read.structure} reads STRUCTURE data files (.str
   ou .stru) and convert them into a \linkS4class{genind} object. By
   default, this function is interactive and asks a few questions about
   data content. This can be disabled (for optional questions) by
   turning the 'ask' argument to FALSE. However, one has to know the
   number of genotypes, of markers and if genotypes are coded on a
   single or on two rows before importing data.

   Note: \code{read.structure} is meant for DIPLOID DATA ONLY. Haploid
   data with the STRUCTURE format can easily be read into R using
   \code{read.table} or \code{read.csv} and then converted using
   \code{\link{df2genind}}.
 }
\usage{
read.structure(file, n.ind=NULL, n.loc=NULL, onerowperind=NULL,
               col.lab=NULL, col.pop=NULL, col.others=NULL, row.marknames=NULL,
               NA.char="-9", pop=NULL, missing=NA, ask=TRUE, quiet=FALSE)
}
\arguments{
  \item{file}{ a character string giving the path to the file to
    convert, with the appropriate extension.}
  \item{n.ind}{an integer giving the number of genotypes (or
    'individuals') in the dataset}
  \item{n.loc}{an integer giving the number of markers in the dataset}
  \item{onerowperind}{a STRUCTURE coding option: are genotypes coded on
    a single row (TRUE), or on two rows (FALSE, default)}
  \item{col.lab}{an integer giving the index of the column containing
    labels of genotypes. '0' if absent.}
  \item{col.pop}{an integer giving the index of the column containing
    population to which genotypes belong. '0' if absent.}
  \item{col.others}{an vector of integers giving the indexes of the
    columns containing other informations to be read. Will be available
    in @other of the created object.}
 \item{row.marknames}{an integer giving the index of the row containing
    the names of the markers. '0' if absent.}
  \item{NA.char}{the character string coding missing data. "-9" by
    default. Note that in any case, series of zero (like "000") are
    interpreted as NA too.}
  \item{pop}{an optional factor giving the population of each individual.}
  \item{ask}{a logical specifying if the function should ask for
    optional informations about the dataset (TRUE, default), or try to
    be as quiet as possible (FALSE).}
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
  Pritchard, J.; Stephens, M. & Donnelly, P. (2000) Inference of
  population structure using multilocus genotype data. \emph{Genetics},
  \bold{155}: 945-959
}
\seealso{
  \code{\link{import2genind}}, \code{\link{df2genind}},
  \code{\link{read.fstat}}, \code{\link{read.genetix}},
  \code{\link{read.genepop}}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
obj <- read.structure(system.file("files/nancycats.str",package="adegenet"),
  onerowperind=FALSE, n.ind=237, n.loc=9, col.lab=1, col.pop=2, ask=FALSE)

obj
}
\keyword{manip}
