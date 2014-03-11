\encoding{UTF-8}
\name{genpop class}
\alias{genpop-class}
\alias{dist,genpop,ANY,ANY,ANY,missing-method}
\alias{names,genpop-method}
\alias{show,genpop-method}
\alias{summary,genpop-method}
\title{adegenet formal class (S4) for allele counts in populations}
\description{An object of class \code{genpop} contain alleles counts
  for several loci.\cr
  It contains several components (see 'slots' section).\cr
  Such object is obtained using \code{genind2genpop} which converts
  individuals genotypes of known population into a \code{genpop} object.
  Note that the function \code{summary} of a \code{genpop} object
  returns a list of components.
  Note that as in other S4 classes, slots are accessed using @ instead
  of \$.
}
\section{Slots}{
         \describe{
    \item{\code{tab}:}{matrix of alleles counts for each combinaison of population
    -in rows- and alleles -in columns-. Rows and columns are given generic names.}
    \item{\code{loc.names}:}{character vector containing the real names of the loci}
    \item{\code{loc.fac}:}{locus factor for the columns of \code{tab}}
    \item{\code{loc.nall}:}{integer vector giving the number of alleles per locus}
    \item{\code{all.names}:}{list having one component per locus, each containing a character vector of alleles names}
    \item{\code{call}:}{the matched call}
    \item{\code{pop.names}:}{character vector containing the real names of the populations}
    \item{\code{ploidy}:}{ an integer indicating the degree of ploidy of
      the genotypes. Beware: 2 is not an integer, but as.integer(2) is.}
    \item{\code{type}:}{ a character string indicating the type of
      marker: 'codom' stands for 'codominant' (e.g. microstallites,
      allozymes); 'PA' stands for 'presence/absence' (e.g. AFLP).}
    \item{\code{other}:}{(optional) a list containing other information}
  }
}
\section{Extends}{
Class \code{"\linkS4class{gen}"}, directly.
Class \code{"\linkS4class{popInfo}"}, directly.
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "genpop")}: give the names of the
      components of a genpop object}
    \item{print}{\code{signature(x = "genpop")}: prints a genpop object}
    \item{show}{\code{signature(object = "genpop")}: shows a genpop
      object (same as print)}
    \item{summary}{\code{signature(object = "genpop")}: summarizes a
      genpop object, invisibly returning its content}
  }
}
\seealso{\code{\link{as.genpop}}, \code{\link{is.genpop}},\code{\link{makefreq}}, \code{\link{genind}}, \code{\link{import2genind}}, \code{\link{read.genetix}}, \code{\link{read.genepop}}, \code{\link{read.fstat}}, \code{\link{na.replace}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
obj1 <- import2genind(system.file("files/nancycats.gen",
package="adegenet"))
obj1


obj2 <- genind2genpop(obj1)
obj2

\dontrun{
data(microsatt)
# use as.genpop to convert convenient count tab to genpop
obj3 <- as.genpop(microsatt$tab)
obj3

all(obj3@tab==microsatt$tab)
all(obj3@pop.names==rownames(microsatt$tab))
# it worked

# perform a correspondance analysis
obj4 <- genind2genpop(obj1,missing="chi2")
ca1 <- dudi.coa(as.data.frame(obj4@tab),scannf=FALSE)
s.label(ca1$li,sub="Correspondance Analysis",csub=2)
add.scatter.eig(ca1$eig,2,xax=1,yax=2,posi="top")

}
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}

