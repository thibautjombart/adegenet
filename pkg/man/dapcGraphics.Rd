\encoding{UTF-8}
\name{dapc graphics}
\alias{scatter.dapc}
\alias{assignplot}
\alias{compoplot}
\title{Graphics for Discriminant Analysis of Principal Components (DAPC)}
\description{
  These functions provide graphic outputs for Discriminant Analysis of
  Principal Components (DAPC, Jombart et al. 2010). See \code{?dapc} for
  details about this method. DAPC graphics are detailed in the DAPC tutorial
  accessible using \code{vignette("adegenet-dapc")}.

  These functions all require an object of class \code{dapc}
  (the ".dapc" can be ommitted when calling the functions):\cr
  - \code{scatter.dapc}: produces scatterplots of principal components (or
    'discriminant functions'), with a screeplot of eigenvalues as inset.\cr
  - \code{assignplot}: plot showing the probabilities of assignment of
  individuals to the different clusters.\cr
 - \code{compoplot}: barplot showing the probabilities of assignment of
  individuals to the different clusters.\cr

}
\usage{
\method{scatter}{dapc}(x, xax=1, yax=2, grp=x$grp, col=seasun(length(levels(grp))),
                       pch=20, bg="white", solid=.7, scree.da=TRUE,
                       scree.pca=FALSE, posi.da="bottomright",
                       posi.pca="bottomleft", bg.inset="white", ratio.da=.25,
                       ratio.pca=.25, inset.da=0.02, inset.pca=0.02,
                       inset.solid=.5, onedim.filled=TRUE, mstree=FALSE, lwd=1,
                       lty=1, segcol="black", legend=FALSE, posi.leg="topright",
                       cleg=1, txt.leg=levels(grp), cstar = 1, cellipse = 1.5,
                       axesell = FALSE, label = levels(grp), clabel = 1, xlim =
                       NULL, ylim = NULL, grid = FALSE, addaxes = TRUE, origin =
                       c(0,0), include.origin = TRUE, sub = "", csub = 1, possub =
                       "bottomleft", cgrid = 1, pixmap = NULL, contour = NULL, area
                       = NULL, \ldots)

assignplot(x, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

compoplot(x, only.grp=NULL, subset=NULL, new.pred=NULL, col=NULL, lab=NULL,
          legend=TRUE, txt.leg=NULL, ncol=4, posi=NULL, cleg=.8, bg=transp("white"), ...)
}
\arguments{
  \item{x}{a \code{dapc} object.}
  \item{xax,yax}{\code{integers} specifying which principal components of DAPC
    should be shown in x and y axes. }
  \item{grp}{a factor defining group membership for the individuals. The scatterplot
    is optimal only for the default group, i.e. the one used in the
    DAPC analysis.}
  \item{col}{a suitable color to be used for groups. The specified vector
    should match the number of groups, not the number of individuals.}
  \item{pch}{a \code{numeric} indicating the type of point to be used to indicate
    the prior group of individuals (see \code{\link{points}} documentation for
    more details); one value is expected for each group; recycled if necessary.}
  \item{bg}{the color used for the background of the scatterplot.}
  \item{solid}{a value between 0 and 1 indicating the alpha level for the colors of
    the plot; 0=full transparency, 1=solid colours.}
  \item{scree.da}{a logical indicating whether a screeplot of Discriminant Analysis
    eigenvalues should be displayed in inset (TRUE) or not (FALSE).}
  \item{scree.pca}{a logical indicating whether a screeplot of Principal Component Analysis
    eigenvalues should be displayed in inset (TRUE) or not (FALSE); retained axes are
    displayed in black.}
  \item{posi.da}{the position of the inset of DA eigenvalues; can match any
    combination of "top/bottom" and "left/right".}
  \item{posi.pca}{the position of the inset of PCA eigenvalues; can match any
    combination of "top/bottom" and "left/right".}
  \item{bg.inset}{the color to be used as background for the inset plots.}
  \item{ratio.da}{the size of the inset of DA eigenvalues as a proportion of the
    current plotting region.}
  \item{ratio.pca}{the size of the inset of PCA eigenvalues as a proportion of the
    current plotting region.}
   \item{inset.da}{a vector with two numeric values (recycled if needed) indicating
     the inset to be used for the screeplot of DA eigenvalues as a proportion of the
     current plotting region; see \code{?add.scatter} for more details.}
   \item{inset.pca}{a vector with two numeric values (recycled if needed) indicating
     the inset to be used for the screeplot of PCA eigenvalues as a proportion of the
     current plotting region; see \code{?add.scatter} for more details.}
   \item{inset.solid}{a value between 0 and 1 indicating the alpha level for the colors of
    the inset plots; 0=full transparency, 1=solid colours.}
  \item{onedim.filled}{a logical indicating whether curves should be filled when
    plotting a single discriminant function (TRUE), or not (FALSE).}
 \item{mstree}{a logical indicating whether a minimum spanning tree linking the
   groups and based on the squared distances between the groups inside the entire
   space should added to the plot (TRUE), or not (FALSE).}
  \item{lwd,lty,segcol}{the line width, line type, and segment colour to be used for
    the minimum spanning tree.}
  \item{legend}{a logical indicating whether a legend for group colours should added
    to the plot (TRUE), or not (FALSE).}
  \item{posi.leg}{the position of the legend for group colours; can match any
    combination of "top/bottom" and "left/right", or a set of x/y coordinates stored
    as a list (\code{locator} can be used).}
  \item{cleg}{a size factor used for the legend.}
  \item{cstar,cellipse,axesell,label,clabel,xlim,ylim,grid,addaxes,origin,include.origin,sub,csub,possub,cgrid,pixmap,contour,area}{arguments
    passed to \code{\link[ade4]{s.class}}; see \code{?s.class} for more informations}
  \item{only.grp}{a \code{character} vector indicating which groups should be
    displayed. Values should match values of \code{x$grp}. If \code{NULL}, all
    results are displayed}
  \item{subset}{\code{integer} or \code{logical} vector indicating which
    individuals should be displayed. If \code{NULL}, all
    results are displayed}
  \item{new.pred}{an optional list, as returned by the \code{predict} method for
    \code{dapc} objects; if provided, the individuals with unknown groups are added
    at the bottom of the plot. To visualize these individuals only, specify
    \code{only.grp="unknown"}.}
  \item{cex.lab}{a \code{numeric} indicating the size of labels.}
  \item{lab}{a vector of characters (recycled if necessary) of labels
    for the individuals; if left to NULL, the row names of \code{x$tab}
    are used.}
  \item{txt.leg}{a character vector indicating the text to be used in
    the legend; if not provided, group names stored in \code{x$grp} are
    used.}
  \item{ncol}{an integer indicating the number of columns of the legend,
    defaulting to 4.}
  \item{posi}{a characther string indicating the position of the legend; can match
    any combination of "top/bottom" and "left/right". See \code{?legend}.}
  \item{\ldots}{further arguments to be passed to other functions. For
    \code{scatter}, arguments passed to \code{points}; for \code{compoplot},
    arguments passed to \code{barplot}.}
}
\details{
  See the documentation of \code{\link{dapc}} for more information about the method.
}
\value{
  All functions return the matched call.\cr
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  - \code{\link{dapc}}: implements the DAPC.

  - \code{\link{find.clusters}}: to identify clusters without prior.

  - \code{\link{dapcIllus}}: a set of simulated data illustrating the DAPC

  - \code{\link{eHGDP}}, \code{\link{H3N2}}: empirical datasets illustrating
  DAPC
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
data(H3N2)
dapc1 <- dapc(H3N2, pop=H3N2$other$epid, n.pca=30,n.da=6)

## defautl plot ##
scatter(dapc1)

## showing different scatter options ##
## remove internal segments and ellipses, different pch, add MStree
scatter(dapc1, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2, posi.da="topleft")

## only ellipse, custom labels, use insets
scatter(dapc1, cell=2, pch="", cstar=0, posi.pca="topleft", posi.da="topleft", scree.pca=TRUE,
inset.pca=c(.01,.3), lab=paste("year\n",2001:2006), axesel=FALSE, col=terrain.colors(10))

## without ellipses, use legend for groups
scatter(dapc1, cell=0, cstar=0, scree.da=FALSE, clab=0, cex=3,
solid=.4, bg="white", leg=TRUE, posi.leg="topleft")

## only one axis
scatter(dapc1,1,1,scree.da=FALSE, legend=TRUE, solid=.4,bg="white")



## example using genlight objects ##
## simulate data
x <- glSim(50,4e3-50, 50, ploidy=2)
x
plot(x)

## perform DAPC
dapc2 <- dapc(x, n.pca=10, n.da=1)
dapc2

## plot results
scatter(dapc2, scree.da=FALSE, leg=TRUE, txt.leg=paste("group",
c('A','B')), col=c("red","blue"))

## SNP contributions
loadingplot(dapc2$var.contr)
loadingplot(tail(dapc2$var.contr, 100), main="Loading plot - last 100 SNPs")



## assignplot / compoplot ##
assignplot(dapc1, only.grp=2006)

data(microbov)
dapc3 <- dapc(microbov, n.pca=20, n.da=15)
compoplot(dapc3, lab="")
}
}
\keyword{multivariate}