\encoding{UTF-8}
\docType{methods}
\name{seploc}
\alias{seploc}
\alias{seploc-methods}
\alias{seploc,ANY-method}
\alias{seploc,genind-method}
\alias{seploc,genpop-method}
\alias{seploc,genlight-method}
\title{ Separate data per locus}
\description{
  The function \code{seploc} splits an object (\linkS4class{genind},
  \linkS4class{genpop} or \linkS4class{genlight}) by marker. For
  \linkS4class{genind} and \linkS4class{genpop} objects, the method
  returns a list of objects whose components each correspond to a
  marker. For \linkS4class{genlight} objects, the methods returns blocks
  of SNPs.
}
\usage{
\S4method{seploc}{genind}(x,truenames=TRUE,res.type=c("genind","matrix"))
\S4method{seploc}{genpop}(x,truenames=TRUE,res.type=c("genpop","matrix"))
\S4method{seploc}{genlight}(x, n.block=NULL, block.size=NULL, random=FALSE,
       parallel=require(parallel), n.cores=NULL)
}
\arguments{
    \item{x}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
    \item{truenames}{a logical indicating whether true names should be
      used (TRUE, default) instead of generic labels (FALSE).}
    \item{res.type}{a character indicating the type of returned results,
      a genind or genpop object (default) or a matrix of data
      corresponding to the 'tab' slot.}
    \item{n.block}{an integer indicating the number of blocks of SNPs to
      be returned.}
    \item{block.size}{an integer indicating the size (in number of SNPs)
      of the blocks to be returned.}
    \item{random}{should blocks be formed of contiguous SNPs, or should
      they be made or randomly chosen SNPs.}
    \item{parallel}{a logical indicating whether multiple cores -if
      available- should be used for the computations (TRUE, default), or
      not (FALSE); requires the package \code{parallel} to be installed.}
    \item{n.cores}{if \code{parallel} is TRUE, the number of cores to be
      used in the computations; if NULL, then the maximum number of cores
      available on the computer is used.}
}
\value{The function \code{seploc} returns an list of objects of the
  same class as the initial object, or a list of matrices similar to
  x\$tab.\cr
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{\code{\link{seppop}}, \code{\link{repool}}}
\examples{
\dontrun{
## example on genind objects
data(microbov)

# separate all markers
obj <- seploc(microbov)
names(obj)

obj$INRA5


## example on genlight objects
x <- glSim(100, 1000, 0, ploidy=2) # simulate data
x <- x[,order(glSum(x))] # reorder loci by frequency of 2nd allele
glPlot(x, main="All data") # plot data
foo <- seploc(x, n.block=3) # form 3 blocks
foo
glPlot(foo[[1]], main="1st block") # plot 1st block
glPlot(foo[[2]], main="2nd block") # plot 2nd block
glPlot(foo[[3]], main="3rd block") # plot 3rd block

foo <- seploc(x, block.size=600, random=TRUE) # split data, randomize loci
foo # note the different block sizes
glPlot(foo[[1]])
}
}
\keyword{manip}