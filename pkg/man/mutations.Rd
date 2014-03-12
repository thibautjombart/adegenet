\encoding{UTF-8}
\name{findMutations}
\alias{findMutations}
\alias{findMutations.DNAbin}
\alias{graphMutations}
\alias{graphMutations.DNAbin}
\title{Identify mutations between DNA sequences}
\description{
  The function \code{findMutations} identifies mutations (position and
  nature) of pairs of aligned DNA sequences. The function
  \code{graphMutations} does the same thing but plotting mutations on a
  directed graph.\cr

  Both functions are generics, but the only methods implemented in
  adegenet so far is for \code{\link[ape]{DNAbin}} objects.
}
\usage{
findMutations(\dots)

\method{findMutations}{DNAbin}(x, from=NULL, to=NULL, allcomb=TRUE, \dots)

graphMutations(\dots)

\method{graphMutations}{DNAbin}(x, from=NULL, to=NULL, allcomb=TRUE, plot=TRUE,
               curved.edges=TRUE, \dots)

}
\arguments{
  \item{x}{a \code{DNAbin} object containing aligned sequences, as a matrix.}
  \item{from}{a vector indicating the DNA sequences from which mutations
    should be found. If \code{NULL}, all sequences are considered (i.e.,
    \code{1:nrow(x)}).}
  \item{to}{a vector indicating the DNA sequences to which mutations
    should be found. If \code{NULL}, all sequences are considered (i.e.,
    \code{1:nrow(x)}).}
  \item{allcomb}{a logical indicating whether all combinations of
    sequences (from and to) should be considered (TRUE, default), or not
    (FALSE).}
  \item{plot}{a logical indicating whether the graph should be plotted.}
  \item{curved.edges}{a logical indicating whether the edges of the
    graph should be curved.}
  \item{\dots}{further arguments to be passed to other methods. Used in
    \code{graphMutations} where it is passed to the plot method for
    \code{igraph} objects.}
}
\value{
  For \code{findMutations}, a named list indicating the mutations from
  one sequence to another. For each comparison, a three-column matrix is
  provided, corresponding to the nucleotides in first and second
  sequence, and a summary of the mutation provided as:
  [position]:[nucleotide in first sequence]->[nucleotide in second
  sequence].

  For \code{graphMutations}, a graph with the class \code{igraph}.
}
\seealso{
  The \code{\link{fasta2DNAbin}} to read fasta alignments with minimum
  RAM use.
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}.
 }
\examples{
\dontrun{
data(woodmouse)

## mutations between first 3 sequences
findMutations(woodmouse[1:3,])

## mutations from the first to sequences 2 and 3
findMutations(woodmouse[1:3,], from=1)

## same, graphical display
g <- graphMutations(woodmouse[1:3,], from=1)

## some manual checks
as.character(woodmouse)[1:3,35]
as.character(woodmouse)[1:3,36]
as.character(woodmouse)[1:3,106]

}
}
