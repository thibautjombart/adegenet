\encoding{UTF-8}
\name{read.snp}
\alias{read.snp}
\title{ Reading Single Nucleotide Polymorphism data}
\description{
  The function \code{read.snp} reads a SNP data file with extension '.snp' and
  converts it into a \linkS4class{genlight} object. This format is
  devoted to handle biallelic SNP only, but can accomodate massive
  datasets such as complete genomes with considerably less memory than
  other formats.

  The function reads data by chunks of a few genomes (minimum 1, no
  maximum) at a time, which allows one to read massive datasets with
  negligible RAM requirements (albeit at a cost of computational
  time). The argument \code{chunkSize} indicates the number of genomes
  read at a time. Increasing this value decreases the computational time
  required to read data in, while increasing memory requirements.

  A description of the .snp format is provided in an example file
  distributed with adegenet (see example below).
}
\usage{
read.snp(file, quiet=FALSE, chunkSize = 1000, parallel = require("parallel"),
    n.cores = NULL, \dots)
}
\arguments{
   \item{file}{ a character string giving the path to the file to
    convert, with the extension ".snp".}
  \item{quiet}{ logical stating whether a conversion messages should be
    printed (TRUE,default) or not (FALSE).}
  \item{chunkSize}{an integer indicating the number of genomes to be
  read at a time; larger values require more RAM but decrease the time
  needed to read the data.}
  \item{parallel}{a logical indicating whether multiple cores -if
    available- should be used for the computations (TRUE, default), or
    not (FALSE); requires the package \code{parallel} to be installed
    (see details).}
  \item{n.cores}{if \code{parallel} is TRUE, the number of cores to be
    used in the computations; if NULL, then the maximum number of cores
    available on the computer is used.}
  \item{\dots}{other arguments to be passed to other functions -
  currently not used.}
}
\details{
  === The .snp format ===
  
  Details of the .snp format can be found in the example file
  distributed with adegenet (see below), or on the adegenet website
  (type \code{adegenetWeb()} in R).
}
\value{an object of the class \linkS4class{genlight}}
\seealso{
  - \code{?genlight} for a description of the class \linkS4class{genlight}.

  - \code{\link{read.PLINK}}: read SNPs in PLINK's '.raw' format.
  
  - \code{\link{fasta2genlight}}: extract SNPs from alignments with
  fasta format.
  
  - \code{\link{df2genind}}: convert any multiallelic markers into
  adegenet \linkS4class{genind}.
  
  - \code{\link{import2genind}}: read multiallelic markers from various
  software into adegenet.\cr
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## show the example file ##
## this is the path to the file:
system.file("files/exampleSnpDat.snp",package="adegenet")

## show its content:
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))


## read the file
obj <-
read.snp(system.file("files/exampleSnpDat.snp",package="adegenet"), chunk=2)
obj
as.matrix(obj)
ploidy(obj)
alleles(obj)
locNames(obj)
}
}
\keyword{manip}
