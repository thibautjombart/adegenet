\encoding{UTF-8}
\name{read.PLINK}
\alias{read.PLINK}
\alias{read.plink}
\alias{extract.PLINKmap}
\title{ Reading PLINK Single Nucleotide Polymorphism data}
\description{
  The function \code{read.PLINK} reads a data file exported by the PLINK
  software with extension '.raw' and converts it into a
  \linkS4class{genlight} object. Optionally, information about SNPs can
  be read from a ".map" file, either by specifying the argument
  \code{map.file} in \code{read.PLINK}, or using \code{extract.PLINKmap}
  to add information to an existing \linkS4class{genlight} object.

  The function reads data by chunks of several genomes (minimum 1, no
  maximum) at a time, which allows one to read massive datasets with
  negligible RAM requirements (albeit at a cost of computational
  time). The argument \code{chunkSize} indicates the number of genomes
  read at a time. Increasing this value decreases the computational time
  required to read data in, while increasing memory requirements.

  See details for the documentation about how to export data using PLINK
  to the '.raw' format.
}
\usage{
read.PLINK(file, map.file=NULL, quiet=FALSE, chunkSize=1000,
           parallel=require("parallel"), n.cores=NULL, \dots)

extract.PLINKmap(file, x=NULL)
}
\arguments{
  \item{file}{for \code{read.PLINK} a character string giving the path
    to the file to convert, with the extension ".raw"; for
    \code{extract.PLINKmap}, a character string giving the path to a
    file with extension ".map".}
  \item{map.file}{an optional character string indicating the path to a
    ".map" file, which contains information about the SNPs (chromosome,
    position). If provided, this information is processed by
    \code{extract.PLINKmap} and stored in the \code{@other} slot.}
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
  \item{x}{an optional object of the class \linkS4class{genlight}, in
    which the information read is stored; if provided, information is
    matched against the names of the loci in \code{x}, as returned by
    \code{locNames(x)}; if not provided, a list of two components is
    returned, containing chromosome and position information.}
}
\details{
  === Exporting data from PLINK ===
  
  Data need to be exported from PLINK using the option "--recodeA" (and
  NOT "--recodeAD"). The PLINK command should therefore look like:
  \code{plink --file data --recodeA}. For more information on this
  topic, please look at this webpage:
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}

}
\value{
  - read.PLINK: an object of the class \linkS4class{genlight}

  - extract.PLINKmap: if a \linkS4class{genlight} is provided as
  argument \code{x}, this object incorporating the new information about
  SNPs in the \code{@other} slot (with new components 'chromosome' and
  'position'); otherwise, a list with two components containing
  chromosome and position information.
}
\seealso{
  - \code{?genlight} for a description of the class \linkS4class{genlight}.

  - \code{\link{read.snp}}: read SNPs in adegenet's '.snp' format.

  - \code{\link{fasta2genlight}}: extract SNPs from alignments with
  fasta format.

  - other import function in adegenet: \code{\link{import2genind}},
  \code{\link{df2genind}}, \code{\link{read.genetix}}
  \code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genepop}}.

  - another function \code{read.plink} is available in the package
  \code{snpMatrix}.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\keyword{manip}
