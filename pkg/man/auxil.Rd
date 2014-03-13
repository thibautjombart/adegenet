\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{checkType}
\alias{.rmspaces}
\alias{.genlab}
\alias{.readExt}
\alias{corner}
\alias{num2col}
\alias{fac2col}
\alias{any2col}
\alias{transp}
\alias{bluepal}
\alias{redpal}
\alias{greenpal}
\alias{greypal}
\alias{flame}
\alias{azur}
\alias{seasun}
\alias{lightseasun}
\alias{deepseasun}
\alias{spectral}
\alias{wasp}
\alias{funky}

\title{ Auxiliary functions for adegenet}
\description{
  adegenet implements a number of auxiliary procedures that might be of
  interest for users. These include graphical tools to translate
  variables (numeric or factors) onto a color scale, adding transparency
  to existing colors, pre-defined color palettes, extra functions to
  access documentation, and low-level treatment of character vectors.

  These functions are mostly auxiliary procedures used internally in
  adegenet.\cr

  These items include:
  \itemize{
    \item \code{num2col}: translates a numeric vector into colors. 
    \item \code{fac2col}: translates a factor into colors. 
    \item \code{any2col}: translates a vector of type numeric, character
    or factor into colors. 
    \item \code{transp}: adds transparency to a vector of colors. Note that
    transparent colors are not supported on some graphical devices.
    \item \code{corner}: adds text to a corner of a figure. 
    \item \code{checkType}: checks the type of markers being used in a
    function and issues an error if appropriate.
    \item \code{.rmspaces}: remove peripheric spaces in a character string. 
    \item \code{.genlab}: generate labels in a correct alphanumeric ordering. 
    \item \code{.readExt}: read the extension of a given file. 
  }

  Color palettes include:
  \itemize{
    \item \code{bluepal}: white -> dark blue
    \item \code{redpal}: white -> dark red
    \item \code{greenpal}: white -> dark green
    \item \code{greypal}: white -> dark grey
    \item \code{flame}: gold -> red
    \item \code{azur}: gold -> blue
    \item \code{seasun}: blue -> gold -> red
    \item \code{lightseasun}: blue -> gold -> red (light variant)
    \item \code{deepseasun}: blue -> gold -> red (deep variant)
    \item \code{spectral}: red -> yellow -> blue (RColorBrewer variant)
    \item \code{wasp}: gold -> brown -> black
    \item \code{funky}: many colors
  }
}
\seealso{
  The R package RColorBrewer, proposing a nice selection of color palettes.
}
\usage{
.genlab(base, n)
corner(text, posi="topleft",  inset=0.1, \dots)
num2col(x, col.pal=heat.colors, reverse=FALSE,
        x.min=min(x,na.rm=TRUE), x.max=max(x,na.rm=TRUE),
        na.col="transparent")
fac2col(x, col.pal=funky, na.col="transparent", seed=NULL)
any2col(x, col.pal=seasun, na.col="transparent")
transp(col, alpha=.5)
}
\arguments{
  \item{base}{a character string forming the base of the labels}
  \item{n}{the number of labels to generate}
  \item{text}{a character string to be added to the plot}
  \item{posi}{a character matching any combinations of "top/bottom" and
    "left/right".}
  \item{inset}{a vector of two numeric values (recycled if needed)
    indicating the inset, as a fraction of the plotting region.}
  \item{\dots}{further arguments to be passed to \code{\link{text}}}
  \item{x}{a numeric vector (for \code{num2col}) or a vector converted
    to a factor (for \code{fac2col}).}
  \item{col.pal}{a function generating colors according to a given palette.}
  \item{reverse}{a logical stating whether the palette should be
    inverted (TRUE), or not (FALSE, default).}
  \item{x.min}{the minimal value from which to start the color scale}
  \item{x.max}{the maximal value from which to start the color scale}
  \item{na.col}{the color to be used for missing values (NAs)}
  \item{seed}{a seed for R's random number generated, used to fix the
    random permutation of colors in the palette used; if NULL, no
    randomization is used and the colors are taken from the palette
    according to the ordering of the levels.}
  \item{col}{a vector of colors}
  \item{alpha}{a numeric value between 0 and 1 representing the alpha
    coefficient; 0: total transparency; 1: no transparency.}
}
\value{
  For \code{.genlab}, a character vector of size "n".
  \code{num2col} and \code{fac2col} return a vector of
  colors. \code{any2col} returns a list with the following components:
  \code{$col} (a vector of colors), \code{$leg.col} (colors for the
  legend), and \code{$leg.txt} (text for the legend).
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{

.genlab("Locus-",11)

## transparent colors using "transp"
plot(rnorm(1000), rnorm(1000), col=transp("blue",.3), pch=20, cex=4)


## numeric values to color using num2col
plot(1:100, col=num2col(1:100), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=bluepal), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=flame), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=wasp), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=azur,rev=TRUE), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=spectral), pch=20, cex=4)

## factor as colors using fac2col
dat <- cbind(c(rnorm(50,8), rnorm(100), rnorm(150,3),
rnorm(50,10)),c(rnorm(50,1),rnorm(100),rnorm(150,3), rnorm(50,5)))
fac <- rep(letters[1:4], c(50,100,150,50))
plot(dat, col=fac2col(fac), pch=19, cex=4)
plot(dat, col=transp(fac2col(fac)), pch=19, cex=4)
plot(dat, col=transp(fac2col(fac,seed=2)), pch=19, cex=4)

## use of any2col
x <- factor(1:10)
col.info <- any2col(x, col.pal=funky)
plot(x, col=col.info$col, main="Use of any2col on a factor")
legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")

x <- 100:1
col.info <- any2col(x, col.pal=wasp)
barplot(x, col=col.info$col, main="Use of any2col on a numeric")
legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")

}
\keyword{manip}