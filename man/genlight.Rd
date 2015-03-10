\name{genlight-class}
\docType{class}
\alias{genlight}
\alias{genlight-class}
\alias{[,genlight-method}
\alias{[,genlight,ANY,ANY-method}
\alias{initialize,genlight-method}
\alias{show,genlight-method}
\alias{nLoc,genlight-method}
\alias{nInd,genlight-method}
\alias{$,genlight-method}
\alias{$<-,genlight-method}
\alias{names,genlight-method}
\alias{ploidy,genlight-method}
\alias{ploidy<-,genlight-method}
\alias{locNames,genlight-method}
\alias{locNames<-,genlight-method}
\alias{indNames,genlight-method}
\alias{indNames<-,genlight-method}
\alias{alleles,genlight-method}
\alias{alleles<-,genlight-method}
\alias{chromosome}
\alias{chromosome<-}
\alias{chromosome,genlight-method}
\alias{chromosome<-,genlight-method}
\alias{chr}
\alias{chr<-}
\alias{chr,genlight-method}
\alias{chr<-,genlight-method}
\alias{position}
\alias{position<-}
\alias{position,genlight-method}
\alias{position<-,genlight-method}
\alias{pop,genlight-method}
\alias{pop<-,genlight-method}
\alias{NA.posi}
\alias{NA.posi,genlight-method}
\alias{other,genlight-method}
\alias{other<-,genlight-method}
\alias{as.matrix.genlight}
\alias{as.data.frame.genlight}
\alias{as,matrix,genlight-method}
\alias{as,data.frame,genlight-method}
\alias{as,list,genlight-method}
\alias{coerce,matrix,genlight-method}
\alias{coerce,data.frame,genlight-method}
\alias{coerce,list,genlight-method}
\alias{as.list.genlight}
\alias{cbind.genlight}
\alias{rbind.genlight}
% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
% \alias{,genlight-method}
%%%%
\title{Formal class "genlight"}
\description{
  The class \code{genlight} is a formal (S4) class for storing a genotypes
  of binary SNPs in a compact way, using a bit-level coding scheme.
  This storage is most efficient with haploid data, where the memory
  taken to represent data can be reduced more than 50 times. However,
  \code{genlight} can be used for any level of ploidy, and still remain an
  efficient storage mode.

  A \code{genlight} object can be constructed from vectors of integers
  giving the number of the second allele for each locus and each
  individual (see 'Objects of the class genlight' below).

  \code{genlight} stores multiple genotypes. Each genotype is stored
  as a \linkS4class{SNPbin} object.
}
\section{Objects from the class genlight}{
  \code{genlight} objects can be created by calls to \code{new("genlight",
    ...)}, where '...' can be the following arguments:
  \describe{
    \item{\code{gen}}{input genotypes, where each genotype is coded as a
      vector of numbers of the second allele. If a list, each slot of the
      list correspond to an individual; if a matrix or a data.frame, rows
      correspond to individuals and columns to SNPs. If individuals or
      loci are named in the input, these names will we stored in the
      produced object. All individuals are expected to have the same
      number of SNPs. Shorter genotypes are completed with NAs, issuing a
      warning.}
    \item{\code{ploidy}}{an optional vector of integers indicating the ploidy of the
      genotypes. Genotypes can therefore have different ploidy. If not
      provided, ploidy will be guessed from the data (as the
      maximum number of second alleles in each individual).}
    \item{\code{ind.names}}{an optional vector of characters giving the labels
      of the genotypes.}
    \item{\code{loc.names}}{an optional vector of characters giving the labels
      of the SNPs.}
    \item{\code{loc.all}}{an optional vector of characters indicating
      the alleles of each SNP; for each SNP, alleles must be coded by two
      letters separated by '/', e.g. 'a/t' is valid, but 'a  t' or 'a |t'
      are not.}
    \item{\code{chromosome}}{an optional factor indicating the
      chromosome to which each SNP belongs.}
    \item{\code{position}}{an optional vector of integers indicating the
      position of the SNPs.}
    \item{\code{other}}{an optional list storing miscellaneous information.}
  }
}
\section{Slots}{
  The following slots are the content of instances of the class
  \code{genlight}; note that in most cases, it is better to retrieve
  information via accessors (see below), rather than by accessing the
  slots manually.
  \describe{
    \item{\code{gen}:}{a list of genotypes stored as  \linkS4class{SNPbin} objects.}
    \item{\code{n.loc}:}{an integer indicating the number of SNPs of the
      genotype.}
    \item{\code{ind.names}:}{a vector of characters indicating the names of
      genotypes.}
    \item{\code{loc.names}:}{a vector of characters indicating the names of
      SNPs.}
    \item{\code{loc.all}:}{a vector of characters indicating the alleles
      of each SNP.}
    \item{\code{chromosome}:}{an optional factor indicating the
  chromosome to which each SNP belongs.}
    \item{\code{position}:}{an optional vector of integers indicating the
      position of the SNPs.}
    \item{\code{ploidy}:}{a vector of integers indicating the ploidy of each individual.}
    \item{\code{pop}:}{a factor indicating the population of each individual.}
    \item{\code{other}:}{a list containing other miscellaneous information.}
  }
}
\section{Methods}{
  Here is a list of methods available for \code{genlight} objects. Most of
  these methods are accessors, that is, functions which are used to
  retrieve the content of the object. Specific manpages can exist for
  accessors with more than one argument. These are indicated by a '*'
  symbol next to the method's name. This list also contains methods
  for conversion from \code{genlight} to other classes.
  \describe{
    \item{[}{\code{signature(x = "genlight")}: usual method to subset
      objects in R. Is to be applied as if the object was a matrix where
      genotypes were rows and SNPs were columns. Indexing can be done via
      vectors of signed integers or of logicals. See details for extra
      supported arguments.}
    \item{show}{\code{signature(x = "genlight")}: printing of the
      object.}
    
    \item{$}{\code{signature(x = "genlight")}: similar to the @ operator;
      used to access the content of slots of the object.}
    \item{$<-}{\code{signature(x = "genlight")}: similar to the @ operator;
      used to replace the content of slots of the object.}
    
    \item{nInd}{\code{signature(x = "genlight")}: returns the number of
      individuals in the object.}
    \item{nLoc}{\code{signature(x = "genlight")}: returns the number of
      SNPs in the object.}
    
    \item{names}{\code{signature(x = "genlight")}: returns the names of
      the slots of the object.}
    
    \item{indNames}{\code{signature(x = "genlight")}: returns the names of
      the individuals, if provided when the object was contructed.}
    \item{indNames<-}{\code{signature(x = "genlight")}: sets the names of
      the individuals using a character vector of length
      \code{nInd(x)}.}
    
    \item{locNames}{\code{signature(x = "genlight")}: returns the names of
      the loci, if provided when the object was contructed.}
    \item{locNames<-}{\code{signature(x = "genlight")}: sets the names of
      the SNPs using a character vector of length \code{nLoc(x)}.}
    
    \item{ploidy}{\code{signature(x = "genlight")}: returns the ploidy of
      the genotypes.}
    \item{ploidy<-}{\code{signature(x = "genlight")}: sets the ploidy of
      the individuals using a vector of integers of size \code{nInd(x)};
      if a single value is provided, the same ploidy is assumed for all
      individuals.}
    
    \item{NA.posi}{\code{signature(x = "genlight")}: returns the indices
      of missing values (NAs) as a list with one vector of integer for each individual.}

    \item{alleles}{\code{signature(x = "genlight")}: returns the names
      of the alleles of each SNPs, if provided when the object was
      contructed.}
    \item{alleles<-}{\code{signature(x = "genlight")}: sets the names
      of the alleles of each SNPs using a character vector of length
      \code{nLoc(x)}; for each SNP, two alleles must be provided,
      separated by a "/", e.g. 'a/t', 'c/a', etc.}
    
    \item{chromosome}{\code{signature(x = "genlight")}: returns a factor
      indicating the chromosome of each SNPs, or NULL if the information
      is missing.}
    \item{chromosome<-}{\code{signature(x = "genlight")}: sets the
      chromosome to which SNPs belong using a factor of length
      \code{nLoc(x)}.}
    
    \item{chr}{\code{signature(x = "genlight")}: shortcut for
      \code{chromosome}.}
    \item{chr<-}{\code{signature(x = "genlight")}: shortcut for
      \code{chromosome<-}.}
    
    \item{position}{\code{signature(x = "genlight")}: returns an integer
      vector indicating the position of each SNPs, or NULL if the
      information is missing.}
    \item{position<-}{\code{signature(x = "genlight")}: sets the
      positions of the SNPs using an integer vector of length
      \code{nLoc(x)}.}
    
    \item{pop}{\code{signature(x = "genlight")}: returns a factor
      indicating the population of each individual, if provided when the
      object was contructed.}
    \item{pop<-}{\code{signature(x = "genlight")}: sets the population
      of each individual using a factor of length \code{nInd(x)}.}
    \item{other}{\code{signature(x = "genlight")}: returns the content of
      the slot \code{@other}.}
    \item{other<-}{\code{signature(x = "genlight")}: sets the content of
      the slot \code{@other}.}
    
    \item{as.matrix}{\code{signature(x = "genlight")}: converts a
      \code{genlight} object into a matrix of integers, with individuals
      in rows and SNPs in columns. The S4 method 'as' can be used as
      well (e.g. as(x, "matrix")).}
    \item{as.data.frame}{\code{signature(x = "genlight")}: same as \code{as.matrix}.}
    \item{as.list}{\code{signature(x = "genlight")}: converts a
      \code{genlight} object into a list of genotypes coded as vector of
      integers (numbers of second allele). The S4 method 'as' can be
      used as well (e.g. as(x, "list")).}
    
    \item{cbind}{\code{signature(x = "genlight")}: merges several
      \linkS4class{genlight} objects by column, i.e. regroups data of
      identical individuals genotyped for different SNPs.}
    \item{rbind}{\code{signature(x = "genlight")}: merges several
      \linkS4class{genlight} objects by row, i.e. regroups data of
      different individuals genotyped for the same SNPs.}
  }
}
\details{
  === On the subsetting using \code{[} ===
  
  The function \code{[} accepts the following extra arguments:
  \describe{
    \item{treatOther}{a logical stating whether elements of the
      \code{@other} slot should be treated as well (TRUE), or not
      (FALSE). If treated, elements of the list are examined for a
      possible match of length (vectors, lists) or number of rows
      (matrices, data frames) with the number of individuals. Those who
      match are subsetted accordingly. Others are left as is, issuing a
      warning unless the argument \code{quiet} is set to TRUE.}
    \item{quiet}{a logical indicating whether warnings should be issued
      when trying to subset components of the \code{@other} slot which
      do not match the number of individuals (TRUE), or not (FALSE,
      default). }
    \item{\dots}{further arguments passed to the genlight constructor.}
    }
}
\author{Thibaut Jombart (\email{t.jombart@imperial.ac.uk})}
\seealso{
 Related class:\cr
  -  \code{\linkS4class{SNPbin}}, for storing individual genotypes of
  binary SNPs\cr
  
  -  \code{\linkS4class{genind}}, for storing other types of genetic markers. \cr
}
\examples{
\dontrun{
## TOY EXAMPLE ##
## create and convert data
dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
x <- new("genlight", dat)
x

## examine the content of the object
names(x)
x@gen
x@gen[[1]]@snp # bit-level coding for first individual

## conversions
as.list(x)
as.matrix(x)

## round trips - must return TRUE
identical(x, new("genlight", as.list(x))) # list
identical(x, new("genlight", as.matrix(x))) # matrix
identical(x, new("genlight", as.data.frame(x))) # data.frame

## test subsetting
x[c(1,3)] # keep individuals 1 and 3
as.list(x[c(1,3)])
x[c(1,3), 1:2] # keep individuals 1 and 3, loci 1 and 2
as.list(x[c(1,3), 1:2])
x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)] # same, using logicals
as.list(x[c(TRUE,FALSE), c(TRUE,TRUE,FALSE,FALSE)])


## REAL-SIZE EXAMPLE ##
## 50 genotypes of 1,000,000 SNPs
dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e6, prob=c(.5, .49, .01), replace=TRUE))
names(dat) <- paste("indiv", 1:length(dat))
print(object.size(dat), unit="aut") # size of the original data

x <- new("genlight", dat) # conversion
x
print(object.size(x), unit="au") # size of the genlight object
object.size(dat)/object.size(x) # conversion efficiency



#### cbind, rbind ####
a <- new("genlight", list(toto=rep(1,10), tata=rep(c(0,1), each=5), titi=c(NA, rep(1,9)) ))

ara <- rbind(a,a)
ara
as.matrix(ara)

aca <- cbind(a,a)
aca
as.matrix(aca)


#### subsetting @other ####
x <- new("genlight", list(a=1,b=0,c=1), other=list(1:3, letters,data.frame(2:4)))
x
other(x)
x[2:3]
other(x[2:3])
other(x[2:3, treatOther=FALSE])


#### seppop ####
pop(x) # no population info
pop(x) <- c("pop1","pop1", "pop2") # set population memberships
pop(x)
seppop(x)
}
}
\keyword{classes}
