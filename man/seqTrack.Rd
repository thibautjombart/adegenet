\name{seqTrack}
\alias{seqTrack}
\alias{seqTrack.default}
\alias{seqTrack.matrix}
\alias{plotSeqTrack}
\alias{get.likelihood}
\alias{get.likelihood.seqTrack}
\alias{seqTrack-class}
\alias{as.igraph.seqTrack}
\alias{plot.seqTrack}
%\alias{as,seqTrack,graphNEL-method}
%\alias{coerce,seqTrack,graphNEL-method}
\title{SeqTrack algorithm for reconstructing genealogies}
\description{
  The SeqTrack algorithm [1] aims at reconstructing genealogies of sampled
  haplotypes or genotypes for which a collection date is
  available. Contrary to phylogenetic methods which aims at reconstructing
  hypothetical ancestors for observed sequences, SeqTrack
  considers that ancestors and descendents are sampled together, and
  therefore infers ancestry relationships among the sampled sequences.

  This approach proved more efficient than phylogenetic approaches for
  reconstructing transmission trees in densely sampled disease outbreaks
  [1]. This implementation defines a generic function \code{seqTrack}
  with methods for specific object classes.

}
\usage{
seqTrack(\ldots)

\method{seqTrack}{matrix}(x, x.names, x.dates, best = c("min", "max"),
    prox.mat = NULL, mu = NULL, haplo.length = NULL, \dots)

\method{as.igraph}{seqTrack}(x, col.pal=redpal, \dots)

\method{plot}{seqTrack}(x, y=NULL, col.pal=redpal, \dots)

plotSeqTrack(x, xy, use.arrows=TRUE, annot=TRUE, labels=NULL, col=NULL,
                         bg="grey", add=FALSE, quiet=FALSE,
                         date.range=NULL, jitter.arrows=0, plot=TRUE, \dots)

get.likelihood(\ldots)

\method{get.likelihood}{seqTrack}(x, mu, haplo.length, \ldots)
}
\arguments{
  \item{x}{for seqTrack, a matrix giving weights to pairs of ancestries
    such that x[i,j] is the weight of 'i ancestor of j'. For
    plotSeqTrack and get.likelihood. seqTrack, a \code{seqTrack}
    object.}
  \item{x.names}{a character vector giving the labels of the haplotypes/genotypes}
  \item{x.dates}{a vector of collection dates for the sampled
    haplotypes/genotypes. Dates must have the POSIXct format. See
    \code{details} or \code{?as.POSIXct} for more information.}
  \item{best}{a character string matching 'min' or 'max', indicating
    whether genealogies should minimize or maximize the sum of weights
    of ancestries.}
  \item{prox.mat}{an optional matrix of proximities between
    haplotypes/genotypes used to resolve ties in the choice of
    ancestors, by picking up the 'closest' ancestor amongst possible
    ancestors, in the sense of \code{prox.mat}. \code{prox.mat[i,j]}
    must indicate a proximity for the relationship 'i ancestor to
    j'. For instance, if \code{prox.mat} contains spatial proximities,
    then \code{prox.mat[i,j]} gives a measure of how easy it is to
    migrate from location 'i' to 'j'.
  }
  \item{mu}{(optional) a mutation rate, per site and per day. When 'x'
    contains numbers of mutations, used to resolve ties using a maximum
    likelihood approach (requires \code{haplo.length} to be provided).}
 \item{haplo.length}{(optional) the length of analysed sequences in
    number of nucleotides. When 'x' contains numbers of mutations, used
    to resolve ties using a maximum likelihood approach (requires
    \code{mu} to be provided).}
  \item{y}{unused argument, for compatibility with 'plot'.}
  \item{col.pal}{a color palette to be used to represent weights using
   colors on the edges of the graph. See \code{?num2col}. Note that the
   palette is inversed by default.}
  \item{xy}{spatial coordinates of the sampled haplotypes/genotypes.}
  \item{use.arrows}{a logical indicating whether arrows should be used to
    represented ancestries (pointing from ancestor to descendent, TRUE),
    or whether segments shall be used (FALSE).}
  \item{annot}{a logical indicating whether arrows or segments
    representing ancestries should be annotated (TRUE) or not (FALSE).}
  \item{labels}{a character vector containing annotations of the
    ancestries. If left empty, ancestries are annotated by the descendent.}
  \item{col}{a vector of colors to be used for plotting ancestries.}
  \item{bg}{a color to be used as background.}
  \item{add}{a logical stating whether the plot should be added to
    current figure (TRUE), or drawn as a new plot (FALSE, default).}
  \item{quiet}{a logical stating whether messages other than errors
    should be displayed (FALSE, default), or hidden (TRUE).}
   \item{date.range}{a vector of length two with POSIXct format indicating
     the time window for which ancestries should be displayed.}
   \item{jitter.arrows}{a positive number indicating the amount of noise
   to be added to coordinates of arrows; useful when several arrows
   overlap. See \code{\link{jitter}}.}
  \item{plot}{a logical stating whether a plot should be drawn (TRUE,
    default), or not (FALSE). In all cases, the
    function invisibly returns plotting information.}
  \item{\ldots}{further arguments to be passed to other methods}
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\references{
  Jombart T, Eggo R, Dodd P, Balloux F (2010) Reconstructing disease
  outbreaks from genetic data: a graph approach. Heredity. doi:
  10.1038/hdy.2010.78.
}
\value{
  === output of seqTrack ===\cr
  seqTrack function returns data.frame with the class \code{seqTrack},
  in which each row is an inferred ancestry described by the following columns:
  - id: indices identifying haplotypes/genotypes\cr
  - ances: index of the inferred ancestor\cr
  - weight: weight of the inferred ancestries\cr
  - date: date of the haplotype/genotype\cr
  - ances.date: date of the ancestor

  === output of plotSeqTrack ===\cr
  This graphical function invisibly returns the coordinates of the
  arrows/segments drawn and their colors, as a data.frame.
}
\details{
  === Maximum parsimony genealogies ===\cr
  Maximum parsimony genealogies can be obtained easily using this
  implementation of seqTrack. One has to provide in \code{x} a matrix of
  genetic distances. The most straightforward distance is the number of
  differing nucleotides. See \code{\link[ape]{dist.dna}} in the ape
  package for a wide range of genetic distances between aligned
  sequences. The argument \code{best} should be set to "min" (its
  default value), so that the identified genealogy minimizes the total
  number of mutations. If \code{x} contains number of mutations, then
  \code{mu} and \code{haplo.length} should also be provided for
  resolving ties in equally parsimonious ancestors using maximum
  likelihood.

  === Likelihood of observed genetic differentiation ===\cr
  The probability of oberving a given number of mutations between a
  sequence and its ancestor can be computed using
  \code{get.likelihood.seqTrack}. Note that this is only possible
  if \code{x} contained number of mutations.


  === Plotting/converting seqTrack objects to graphs ===\cr
  seqTrack objects are best plotted as graphs. From adegenet_1.3-5
  onwards, seqTrack objects can be converted to \code{igraph} objects (from the
  package \code{igraph}), which can in turn be plotted and manipulated
  using classical graph tools. The plot method does this operation
  automatically, using colors to represent edge weights, and using
  time-ordering of the data from top (ancient) to bottom (recent).
}
\seealso{
  \code{\link[ape]{dist.dna}} in the ape package to compute pairwise
  genetic distances in aligned sequences.
}
\examples{
\dontrun{
if(require(ape && require(igraph))){
## ANALYSIS OF SIMULATED DATA ##
## SIMULATE A GENEALOGY
dat <- haploGen(seq.l=1e4, repro=function(){sample(1:4,1)}, gen.time=1, t.max=3)
plot(dat, main="Simulated data")

## SEQTRACK ANALYSIS
res <- seqTrack(dat, mu=0.0001, haplo.length=1e4) 
plot(res, main="seqTrack reconstruction")

## PROPORTION OF CORRECT RECONSTRUCTION
mean(dat$ances==res$ances,na.rm=TRUE)


## ANALYSIS OF PANDEMIC A/H1N1 INFLUENZA DATA ##
## note:
## this is for reproduction purpose only
## seqTrack is best kept for the analysis
## of densely sampled outbreaks, which
## is not the case of this dataset.
## 
dat <- read.csv(system.file("files/pdH1N1-data.csv",package="adegenet"))
ha <-  read.dna(system.file("files/pdH1N1-HA.fasta",package="adegenet"), format="fa")
na <- read.dna(system.file("files/pdH1N1-NA.fasta",package="adegenet"), format="fa")


## COMPUTE NUCLEOTIDIC DISTANCES
nbNucl <- ncol(as.matrix(ha)) + ncol(as.matrix(na))
D <- dist.dna(ha,model="raw")*ncol(as.matrix(ha)) +
dist.dna(na,model="raw")*ncol(as.matrix(na))
D <- round(as.matrix(D))


## MATRIX OF SPATIAL CONNECTIVITY
## (to promote local transmissions)
xy <- cbind(dat$lon, dat$lat)
temp <- as.matrix(dist(xy))
M <- 1* (temp < 1e-10)


## SEQTRACK ANALYSIS
dat$date <- as.POSIXct(dat$date)
res <- seqTrack(D, rownames(dat), dat$date, prox.mat=M, mu=.00502/365, haplo.le=nbNucl)


## COMPUTE GENETIC LIKELIHOOD
p <- get.likelihood(res, mu=.00502/365, haplo.length=nbNucl)
# (these could be shown as colors when plotting results)
# (but mutations will be used instead)


## EXAMINE RESULTS
head(res)
tail(res)
range(res$weight, na.rm=TRUE)
barplot(table(res$weight)/sum(!is.na(res$weight)), ylab="Frequency",
xlab="Mutations between inferred ancestor and descendent", col="orange")


## DISPLAY SPATIO-TEMPORAL DYNAMICS 
if(require(maps)){
myDates <- as.integer(difftime(dat$date, as.POSIXct("2009-01-21"), unit="day"))
myMonth <- as.POSIXct(
c("2009-02-01", "2009-03-01","2009-04-01","2009-05-01","2009-06-01","2009-07-01"))
x.month <-  as.integer(difftime(myMonth, as.POSIXct("2009-01-21"), unit="day"))


## FIRST STAGE:
## SPREAD TO THE USA AND CANADA
curRange <- as.POSIXct(c("2009-03-29","2009-04-25"))
par(bg="deepskyblue")
map("world", fill=TRUE, col="grey")
opal <- palette()
palette(rev(heat.colors(10)))
plotSeqTrack(res, round(xy), add=TRUE,annot=FALSE,lwd=2,date.range=curRange,
col=res$weight+1)
title(paste(curRange, collapse=" to "))
legend("bottom", lty=1, leg=0:8, title="number of mutations", col=1:9,
lwd=2, horiz=TRUE)


## SECOND STAGE:
## SPREAD WITHIN AMERICA, FIRST SEEDING OUTSIDE AMERICA
curRange <- as.POSIXct(c("2009-04-30","2009-05-07"))
par(bg="deepskyblue")
map("world", fill=TRUE, col="grey")
opal <- palette()
palette(rev(heat.colors(10)))
plotSeqTrack(res, round(xy), add=TRUE,annot=FALSE,lwd=2,
date.range=curRange, col=res$weight+1)
title(paste(curRange, collapse=" to "))
legend("bottom", lty=1, leg=0:8, title="number of mutations",
col=1:9,lwd=2, horiz=TRUE)


## THIRD STAGE:
## PANDEMIC
curRange <- as.POSIXct(c("2009-05-15","2009-05-25"))
par(bg="deepskyblue")
map("world", fill=TRUE, col="grey")
opal <- palette()
palette(rev(heat.colors(10)))
plotSeqTrack(res, round(xy), add=TRUE,annot=FALSE,lwd=2, date.range=curRange,
col=res$weight+1)
title(paste(curRange, collapse=" to "))
legend("bottom", lty=1, leg=0:8, title="number of mutations",
col=1:9,lwd=2, horiz=TRUE)

}
}
}
}