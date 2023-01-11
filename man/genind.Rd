\encoding{UTF-8}
\name{genind class}
\alias{genind-class}
\alias{print,genind-method}
\alias{show,genind-method}
\alias{names,genind-method}
\alias{summary,genind-method}
\alias{print,genindSummary-method}
\alias{print.genindSummary}
\alias{is.genind}
\alias{.valid.genind}
\title{adegenet formal class (S4) for individual genotypes}
\description{The S4 class \code{genind} is used to store individual genotypes.\cr
  It contains several components described in the 'slots' section).\cr
  The \code{summary} of a \code{genind} object invisibly returns a list of component.
  The function \code{.valid.genind} is for internal use.
  The function \code{genind} creates a genind object from a valid table
  of alleles corresponding to the \code{@tab} slot.
  Note that as in other S4 classes, slots are accessed using @ instead
  of $.
}
\section{Slots}{
  \describe{
    \item{\code{tab}:}{ \strong{(accessor: \code{\link{tab}})}
      matrix integers containing genotypes data for
      individuals (in rows) for all alleles (in columns). The table
      differs depending on the \code{@type} slot:\cr
      - 'codom': values
      are numbers of alleles, summing up to the
      individuals' ploidies.\cr
      - 'PA': values are presence/absence of alleles.\cr
      In all cases, rows and columns are given generic names.}
    \item{\code{loc.fac}:}{\strong{(accessor: \code{\link{locFac}})} locus 
    factor for the columns of \code{tab}}
    \item{\code{loc.n.all}:}{\strong{(accessor: \code{\link{nAll}})} integer 
    vector giving the number of observed alleles per locus (see note)}
    \item{\code{all.names}:}{\strong{(accessor: \code{\link{alleles}})} list 
    having one component per locus, each containing a character vector of allele
    names}
    \item{\code{ploidy}:}{\strong{(accessor: \code{\link{ploidy}})} an integer 
    vector indicating the degree of ploidy of the genotypes. Beware: 2 is not 
    an integer, but 2L or as.integer(2) is.}
    \item{\code{type}:}{ a character string indicating the type of
      marker: 'codom' stands for 'codominant' (e.g. microstallites,
      allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
    \item{\code{call}:}{the matched call}
    \item{\code{strata}:}{\strong{(accessor: \code{\link{strata}})} (optional) 
    data frame giving levels of population stratification for each individual}
    \item{\code{hierarchy}:}{\strong{(accessor: \code{\link{hier}})} (optional, 
    currently unused) a hierarchical \code{\link{formula}} defining the 
    hierarchical levels in the \code{@@strata} slot.}
    \item{\code{pop}:}{\strong{(accessor: \code{\link{pop}})} (optional) factor 
    giving the population of each individual}
    \item{\code{other}:}{\strong{(accessor: \code{\link{other}})} (optional) a 
    list containing other information}
  }
}
\section{Note:}{
  The \code{loc.n.all} slot will reflect the number of columns per locus that 
  contain at least one observation. This means that the sum of the this vector 
  will not necessarily equal the number of columns in the data unless you use 
  \code{drop = TRUE} when subsetting. 
}
\section{Extends}{
Class \code{"\linkS4class{gen}"}, directly.
Class \code{"\linkS4class{indInfo}"}, directly.
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "genind")}: give the names of the
      components of a genind object}
    \item{print}{\code{signature(x = "genind")}: prints a genind object}
    \item{show}{\code{signature(object = "genind")}: shows a genind
      object (same as print)}
    \item{summary}{\code{signature(object = "genind")}: summarizes a
      genind object, invisibly returning its content or suppress printing of auxiliary information by specifying \code{verbose = FALSE}}
  }
}
\seealso{\code{\link{as.genind}}, \code{\link{genind2genpop}},
  \code{\link{genpop}}, \code{\link{import2genind}},
  \code{\link{read.genetix}}, \code{\link{read.genepop}},
  \code{\link{read.fstat}}\cr

  Related classes:\cr
  - \linkS4class{genpop} for storing data per populations\cr

  - \linkS4class{genlight} for an efficient storage of binary SNPs genotypes\cr
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
showClass("genind")

obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))
obj
validObject(obj)
summary(obj)

\dontrun{
# test inter-colonies structuration
if(require(hierfstat)){
gtest <- gstat.randtest(obj,nsim=99)
gtest
plot(gtest)
}

# perform a between-class PCA
pca1 <- dudi.pca(scaleGen(obj, NA.method="mean"),scannf=FALSE,scale=FALSE)
pcabet1 <- between(pca1,obj@pop,scannf=FALSE)
pcabet1

s.class(pcabet1$ls,obj@pop,sub="Inter-class PCA",possub="topleft",csub=2)
add.scatter.eig(pcabet1$eig,2,xax=1,yax=2)

}
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}

