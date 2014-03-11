\encoding{UTF-8}
\name{genind class}
\alias{genind-class}
\alias{print,genind-method}
\alias{show,genind-method}
\alias{names,genind-method}
\alias{summary,genind-method}
\alias{.valid.genind}
\title{adegenet formal class (S4) for individual genotypes}
\description{The S4 class \code{genind} is used to store individual genotypes.\cr
  It contains several components described in the 'slots' section).\cr
  The \code{summary} of a \code{genind} object invisibly returns a list of component.
  The function \code{.valid.genind} is for internal use.
  The function \code{genind} creates a genind object from a valid table
  of alleles corresponding to the \code{@tab} slot.
  Note that as in other S4 classes, slots are accessed using @ instead
  of \$.
}
\section{Slots}{
  \describe{
    \item{\code{tab}:}{matrix of genotypes (in rows) for all alleles (in
      columns). The table differs depending on the \code{@type} slot:\cr
      - 'codom': values are frequencies ; '0' if the genotype does not have
      the corresponding allele, '1' for an homozygote and 0.5 for an
      heterozygte.\cr
      - 'PA': values are presence/absence of alleles.\cr
      In all cases, rows and columns are given generic names.}
    \item{\code{loc.names}:}{character vector containing the real names of the loci}
    \item{\code{loc.fac}:}{locus factor for the columns of \code{tab}}
    \item{\code{loc.nall}:}{integer vector giving the number of alleles per locus}
    \item{\code{all.names}:}{list having one component per locus, each containing a character vector of alleles names}
    \item{\code{call}:}{the matched call}
    \item{\code{ind.names}:}{character vector containing the real names of the
      individuals. Note that as Fstat does not store these names, objects
      converted from .dat files will contain empty \code{ind.names}.}
    \item{\code{ploidy}:}{ an integer indicating the degree of ploidy of
      the genotypes. Beware: 2 is not an integer, but as.integer(2) is.}
    \item{\code{type}:}{ a character string indicating the type of
      marker: 'codom' stands for 'codominant' (e.g. microstallites,
      allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
    \item{\code{pop}:}{(optional) factor giving the population of each individual}
    \item{\code{pop.names}:}{(optional) vector giving the real names of the
      populations}
    \item{\code{other}:}{(optional) a list containing other information}
  }
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
      genind object, invisibly returning its content}
  }
}
\seealso{\code{\link{as.genind}}, \code{\link{is.genind}}, \code{\link{genind2genpop}},
  \code{\link{genpop}}, \code{\link{import2genind}},
  \code{\link{read.genetix}}, \code{\link{read.genepop}},
  \code{\link{read.fstat}}, \code{\link{na.replace}}\cr

  Related classes:\cr
  - \linkS4class{genpop} for storing data per populations\cr

  - \linkS4class{genlight} for an efficient storage of binary SNPs genotypes\cr
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
showClass("genind")

obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"),missing="mean")
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
pca1 <- dudi.pca(obj@tab,scannf=FALSE,scale=FALSE)
pcabet1 <- between(pca1,obj@pop,scannf=FALSE)
pcabet1

s.class(pcabet1$ls,obj@pop,sub="Inter-class PCA",possub="topleft",csub=2)
add.scatter.eig(pcabet1$eig,2,xax=1,yax=2)

}
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}

