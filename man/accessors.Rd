\encoding{UTF-8}
\docType{methods}
\name{Accessors}
\alias{$,genind-method}
\alias{$,genpop-method}
\alias{$<-,genind-method}
\alias{$<-,genpop-method}
\alias{[,genind-method}
\alias{[,genind,ANY,ANY,ANY-method}
\alias{[,genpop-method}
\alias{[,genpop,ANY,ANY,ANY-method}
\alias{nLoc}
\alias{nLoc,gen-method}
\alias{nLoc,genind-method}
\alias{nLoc,genpop-method}
\alias{locFac}
\alias{locFac,gen-method}
\alias{locFac,genind-method}
\alias{locFac,genpop-method}
\alias{nAll}
\alias{nAll,gen-method}
\alias{nAll,genind-method}
\alias{nAll,genpop-method}
\alias{nPop}
\alias{nPop,genind-method}
\alias{nPop,genpop-method}
\alias{nInd}
\alias{nInd,genind-method}
\alias{pop}
\alias{pop<-}
\alias{pop,genind-method}
\alias{pop<-,gen-method}
\alias{pop<-,genind-method}
\alias{locNames}
\alias{locNames,gen-method}
\alias{locNames,genind-method}
\alias{locNames,genpop-method}
\alias{locNames<-}
\alias{locNames<-,gen-method}
\alias{locNames<-,genind-method}
\alias{locNames<-,genpop-method}
\alias{popNames}
\alias{popNames,genind-method}
\alias{popNames,genpop-method}
\alias{popNames<-}
\alias{popNames<-,genind-method}
\alias{popNames<-,genpop-method}
\alias{indNames}
\alias{indNames,genind-method}
\alias{indNames<-}
\alias{indNames<-,genind-method}
\alias{ploidy}
\alias{ploidy,genind-method}
\alias{ploidy,genpop-method}
\alias{ploidy<-}
\alias{ploidy<-,genind-method}
\alias{ploidy<-,genpop-method}
\alias{alleles}
\alias{alleles,gen-method}
\alias{alleles,genind-method}
\alias{alleles,genpop-method}
\alias{alleles<-}
\alias{alleles<-,gen-method}
\alias{alleles<-,genind-method}
\alias{alleles<-,genpop-method}
\alias{other}
\alias{other,gen-method}
\alias{other,genind-method}
\alias{other,genpop-method}
\alias{other<-}
\alias{other<-,gen-method}
\alias{other<-,genind-method}
\alias{other<-,genpop-method}

\title{Accessors for adegenet objects}
\description{
  An accessor is a function that allows to interact with slots of an
  object in a convenient way. Several accessors are available for \linkS4class{genind} or
  \linkS4class{genpop} objects. The operator "\$" and "\$<-" are used to
  access the slots, being equivalent to "@" and "@<-".\cr

  The operator "[" is a flexible way to subset data by individuals,
  populations, alleles, and loci. When using a matrix-like syntax,
  subsetting will apply to the dimensios of the @tab slot. In addition,
  specific arguments \code{loc} and \code{pop} can be used to indicate
  subsets of loci and populations. The argument \code{drop} is a logical
  indicating if alleles becoming non-polymorphic in a new dataset should
  be removed (default: FALSE). Examples:\cr
  \itemize{
  \item "obj[i,j]" returns "obj" with a subset 'i' of individuals and 'j' of
  alleles.\cr
  
  \item "obj[1:10,]" returns an object with only the first 10 genotypes (if
  "obj" is a \linkS4class{genind}) or the first 10 populations (if
  "obj" is a \linkS4class{genpop}) \cr
  
  \item "obj[1:10, 5:10]" returns an object keeping the first 10 entities and
  the alleles 5 to 10.\cr
  
  \item "obj[loc=c(1,3)]" returns an object keeping only the 1st and 3rd
  loci, using \code{locNames(obj)} as reference; logicals, or named loci
  also work; this overrides other subsetting of alleles.\cr
  
 \item "obj[pop=2:4]" returns an object keeping only individuals from the
  populations 2, 3 and 4, using \code{popNames(obj)} as reference;
  logicals, or named populations also work; this overrides other
  subsetting of individuals.\cr

  \item "obj[i=1:2, drop=TRUE]" returns an object keeping only the first
  two individuals (or populations), dropping the alleles no longer
  present in the data.\cr
}
  The argument \code{treatOther} handles the treatment of objects in the
  \code{@other} slot (see details). The argument \code{drop} can be set
  to TRUE to drop alleles that are no longer represented in the subset.
}
\section{Methods}{
  \describe{
    \item{nInd}{returns the number of individuals in the \code{genind} object}
    \item{nLoc}{returns the number of loci}
    \item{nAll}{returns the number of observed alleles in each locus}
    \item{nPop}{returns the number of populations}
    \item{pop}{returns a factor assigning individuals to populations.}
    \item{pop<-}{replacement method for the \code{@pop} slot of an
      object.}
    \item{popNames}{returns the names of populations.}
    \item{popNames<-}{sets the names of populations using a vector of
      length \code{nPop(x)}.}
    \item{indNames}{returns the names of individuals.}
    \item{indNames<-}{sets the names of individuals using a vector of
      length \code{nInd(x)}.}
    \item{locNames}{returns the names of markers and/or alleles.}
    \item{locNames<-}{sets the names of markers using a vector of
      length \code{nLoc(x)}.}
    \item{locFac}{returns a factor that defines which locus each column of the \code{@tab} slot belongs to}
    \item{ploidy}{returns the ploidy of the data.}
    \item{ploidy<-}{sets the ploidy of the data using an integer.}
    \item{alleles}{returns the alleles of each locus.}
    \item{alleles<-}{sets the alleles of each locus using a list with
      one character vector for each locus.}
    \item{other}{returns the content of the \code{@other} slot
      (misc. information); returns \code{NULL} if the slot is onlyObserved or of
      length zero.}
    \item{other<-}{sets the content of the \code{@other} slot
      (misc. information); the provided value needs to be a list; it
      not, provided value will be stored within a list.}
  }
}
\usage{
nInd(x, \dots)
nLoc(x, \dots)
nAll(x, onlyObserved = FALSE, \dots)
nPop(x, \dots)
pop(x)
indNames(x, \dots)
\S4method{indNames}{genind}(x, \dots)
locNames(x, \dots)
\S4method{locNames}{genind}(x, withAlleles=FALSE, \dots)
\S4method{locNames}{genpop}(x, withAlleles=FALSE, \dots)
popNames(x, \dots)
\S4method{popNames}{genind}(x, \dots)
popNames(x, \dots)
\S4method{popNames}{genpop}(x, \dots)
ploidy(x, \dots)
\S4method{ploidy}{genind}(x, \dots)
\S4method{ploidy}{genpop}(x, \dots)
\S4method{other}{genind}(x, \dots)
\S4method{other}{genpop}(x, \dots)
}
\arguments{
  \item{x}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
  \item{onlyObserved}{a logical indicating whether the allele count should
    also include the alleles with onlyObserved columns in the matrix. Defaults
    to \code{FALSE}, which will report only the observed alleles in the
    given population. \code{onlyObserved = TRUE} will be the equivalent of
    \code{table(locFac(x))}, but faster.}
  \item{withAlleles}{a logical indicating whether the result should be
    of the form [locus name].[allele name], instead of [locus name].}
  \item{\dots}{further arguments to be passed to other methods
  (currently not used).}
}
\value{
  A \linkS4class{genind} or \linkS4class{genpop} object.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\details{
  The "[" operator can treat elements in the \code{@other} slot as
  well. For instance, if \code{obj@other$xy} contains spatial
  coordinates, the \code{obj[1:3, ]@other$xy} will contain the spatial
  coordinates of the genotypes (or population) 1,2 and 3. This is
  handled through the argument \code{treatOther}, a logical defaulting
  to TRUE. If set to FALSE, the \code{@other} returned unmodified.\cr
  
  Note that only matrix-like, vector-like and lists can be proceeded in
  \code{@other}. Other kind of objects will issue a warning an be
  returned as they are, unless the argument \code{quiet} is left to
  TRUE, its default value.\cr

  The \code{drop} argument can be set to TRUE to retain only alleles
  that are present in the subset. To achieve better control of
  polymorphism of the data, see \code{\link{isPoly}}.
  
  \code{nAll()} reflects the number of columns per locus present in the current
  gen object. If \code{onlyObserved = TRUE}, then the number of columns with at
  least one non-missing allele is shown. 
}
\examples{
data(nancycats)
nancycats
pop(nancycats) # get the populations
indNames(nancycats) # get the labels of individuals
locNames(nancycats) # get the labels of the loci
alleles(nancycats)  # get the alleles
nAll(nancycats)     # count the number of alleles

head(tab(nancycats)) # get allele counts

# get allele frequencies, replace NAs
head(tab(nancycats, freq = TRUE, NA.method = "mean")) 

# let's isolate populations 4 and 8
popNames(nancycats)
obj <- nancycats[pop=c(4, 8)]
obj
popNames(obj)
pop(obj)
nAll(obj, onlyObserved = TRUE) # count number of alleles among these two populations
nAll(obj) # count number of columns in the data
all(nAll(obj, onlyObserved = TRUE) == lengths(alleles(obj))) # will be FALSE since drop = FALSE
all(nAll(obj) == lengths(alleles(obj))) # will be FALSE since drop = FALSE

# let's isolate two markers, fca23 and fca90
locNames(nancycats)
obj <- nancycats[loc=c("fca23","fca90")]
obj
locNames(obj)

# illustrate pop
obj <- nancycats[sample(1:100, 10)]
pop(obj)
pop(obj) <- rep(c('b', 'a'), each = 5)
pop(obj)

# illustrate locNames
locNames(obj)
locNames(obj, withAlleles = TRUE)
locNames(obj)[1] <- "newLocus"
locNames(obj)
locNames(obj, withAlleles=TRUE)

# illustrate how 'other' slot is handled
data(sim2pop)
nInd(sim2pop)
other(sim2pop[1:6]) # xy is subsetted automatically
other(sim2pop[1:6, treatOther=FALSE]) # xy is left as is

}
\keyword{manip}
