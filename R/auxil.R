###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################

#' Auxiliary functions for adegenet
#'
#' adegenet implements a number of auxiliary procedures that might be of
#' interest for users. These include graphical tools to translate variables
#' (numeric or factors) onto a color scale, adding transparency to existing
#' colors, pre-defined color palettes, extra functions to access documentation,
#' and low-level treatment of character vectors.
#'
#' These functions are mostly auxiliary procedures used internally in
#' adegenet.\cr
#'
#' These items include:
#' \itemize{
#' \item \code{num2col}: translates a numeric vector into colors.
#'
#' \item \code{fac2col}: translates a factor into colors.
#'
#' \item \code{any2col}: translates a vector of type numeric, character or factor into colors.
#'
#' \item \code{transp}: adds transparency to a vector of colors. Note that transparent colors are not supported on some graphical devices.
#'
#' \item \code{corner}: adds text to a corner of a figure.
#'
#' \item \code{checkType}: checks the type of markers being used in a function and issues an error if appropriate.
#'
#' \item \code{.rmspaces}: remove peripheric spaces in a character string.
#'
#' \item \code{.genlab}: generate labels in a correct alphanumeric ordering.
#'
#' \item \code{.readExt}: read the extension of a given file.
#' }
#'
#' Color palettes include:
#' \itemize{
#'
#' \item \code{bluepal}: white -> dark blue
#'
#' \item \code{redpal}: white -> dark red
#'
#' \item \code{greenpal}: white -> dark green
#'
#' \item \code{greypal}: white -> dark grey
#'
#' \item \code{flame}: gold -> red
#'
#' \item \code{azur}: gold -> blue
#'
#' \item \code{seasun}: blue -> gold -> red
#'
#' \item \code{lightseasun}: blue -> gold -> red (light variant)
#'
#' \item \code{deepseasun}: blue -> gold -> red (deep variant)
#'
#' \item \code{spectral}: red -> yellow -> blue (RColorBrewer variant)
#'
#' \item \code{wasp}: gold -> brown -> black \item \code{funky}: many colors }
#'
#' @name Auxiliary functions
#' @aliases checkType .rmspaces .genlab .readExt corner num2col fac2col any2col
#' transp bluepal redpal greenpal greypal flame azur seasun lightseasun
#' deepseasun spectral wasp funky
#' @docType methods
#' @param base a character string forming the base of the labels
#' @param n the number of labels to generate
#' @param text a character string to be added to the plot
#' @param posi a character matching any combinations of "top/bottom" and
#' "left/right".
#' @param inset a vector of two numeric values (recycled if needed) indicating
#' the inset, as a fraction of the plotting region.
#' @param \dots further arguments to be passed to \code{\link{text}}
#' @param x a numeric vector (for \code{num2col}) or a vector converted to a
#' factor (for \code{fac2col}).
#' @param col.pal a function generating colors according to a given palette.
#' @param reverse a logical stating whether the palette should be inverted
#' (TRUE), or not (FALSE, default).
#' @param x.min the minimal value from which to start the color scale
#' @param x.max the maximal value from which to start the color scale
#' @param na.col the color to be used for missing values (NAs)
#' @param seed a seed for R's random number generated, used to fix the random
#' permutation of colors in the palette used; if NULL, no randomization is used
#' and the colors are taken from the palette according to the ordering of the
#' levels.
#' @param col a vector of colors
#' @param alpha a numeric value between 0 and 1 representing the alpha
#' coefficient; 0: total transparency; 1: no transparency.
#' @return For \code{.genlab}, a character vector of size "n".  \code{num2col}
#' and \code{fac2col} return a vector of colors. \code{any2col} returns a list
#' with the following components: \code{$col} (a vector of colors),
#' \code{$leg.col} (colors for the legend), and \code{$leg.txt} (text for the
#' legend).
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso The R package RColorBrewer, proposing a nice selection of color
#' palettes.
#' @keywords manip
#' @examples
#'
#'
#' .genlab("Locus-",11)
#'
#' ## transparent colors using "transp"
#' plot(rnorm(1000), rnorm(1000), col=transp("blue",.3), pch=20, cex=4)
#'
#'
#' ## numeric values to color using num2col
#' plot(1:100, col=num2col(1:100), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=bluepal), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=flame), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=wasp), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=azur,rev=TRUE), pch=20, cex=4)
#' plot(1:100, col=num2col(1:100, col.pal=spectral), pch=20, cex=4)
#'
#' ## factor as colors using fac2col
#' dat <- cbind(c(rnorm(50,8), rnorm(100), rnorm(150,3),
#' rnorm(50,10)),c(rnorm(50,1),rnorm(100),rnorm(150,3), rnorm(50,5)))
#' fac <- rep(letters[1:4], c(50,100,150,50))
#' plot(dat, col=fac2col(fac), pch=19, cex=4)
#' plot(dat, col=transp(fac2col(fac)), pch=19, cex=4)
#' plot(dat, col=transp(fac2col(fac,seed=2)), pch=19, cex=4)
#'
#' ## use of any2col
#' x <- factor(1:10)
#' col.info <- any2col(x, col.pal=funky)
#' plot(x, col=col.info$col, main="Use of any2col on a factor")
#' legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")
#'
#' x <- 100:1
#' col.info <- any2col(x, col.pal=wasp)
#' barplot(x, col=col.info$col, main="Use of any2col on a numeric")
#' legend("bottomleft", fill=col.info$leg.col, legend=col.info$leg.txt, bg="white")
#'
#'
#'
#######################
# Function rmspaces
#######################
# removes spaces and tab at the begining and the end of each element of charvec
.rmspaces <- function(charvec){
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
    return(charvec)
}





###################
# Function readExt
###################
.readExt <- function(char){
    temp <- as.character(char)
    temp <- unlist(strsplit(char,"[.]"))
    res <- temp[length(temp)]
    return(res)
}





###################
# Function .genlab
###################
# recursive function to have labels of constant length
# base = a character string
# n = number of labels
.genlab <- function(base, n) {
  f1 <- function(cha,n){
    if(nchar(cha)<n){
      cha <- paste("0",cha,sep="")
      return(f1(cha,n))
    } else {return(cha)}
  }
  w <- as.character(1:n)
  max0 <- max(nchar(w))
  w <- sapply(w, function(cha) f1(cha,max0))
  return(paste(base,w,sep=""))
}






#' Website and tutorials for adegenet
#'
#' These functions allow to access documentation for adegenet available
#' online.\cr
#'
#' These functions include: \itemize{ \item \code{adegenetWeb}: opens the
#' adegenet website in a web navigator.  \item \code{adegenetTutorial}: opens a
#' specified tutorial.  }
#'
#'
#' @name Website and tutorials
#' @aliases adegenetWeb adegenetTutorial adegenetTutorials
#' @docType methods
#' @param which a character string indicating the type of tutorial to open
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @examples
#'
#' \dontrun{
#' ## this opens the adegenet website
#' adegenetWeb()
#' adegenetTutorial("dapc")
#' }
#'

adegenetWeb <- function(){
    cat("Opening url \"http://adegenet.r-forge.r-project.org/\" ...\n")
    browseURL("http://adegenet.r-forge.r-project.org/")
}


adegenetTutorial <- function(which=c("basics","spca","dapc","genomics")){
    which <- match.arg(which)
    if(which=="basics"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf"
        cat("\n")
        cat("  >> Opening the general introduction to adegenet.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="spca"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf"
        cat("\n")
        cat("  >> Opening the sPCA tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="dapc"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf"
        cat("\n")
        cat("  >> Opening the DAPC tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="genomics"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf"
        cat("\n")
        cat("  >> Opening the genomics tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }

    browseURL(url)

    return(invisible(NULL))
}





############
# checkType
############
##
## WARNING: this does not work with S4 methods
##
checkType <- function(x){
    if(is.character(x)){
        markType <- x
    } else {
        markType <- x@type
    }

    if(markType=="codom") return() # always ok for codominant markers

    currCall <- as.character(sys.call(sys.parent()))[1]
    currFunction <- sub("[[:space:]]*[(].*","",currCall)
    if(currFunction==".local"){
        warning("Current call not found - stopping check (please report this warning).")
        return()
    }

    ## names of functions which are ok for dominant markers
    PAOk <- c("genind","genpop","genind2genpop","summary","df2genind", "genind2df",
                 "truenames","seppop","na.replace","nLoc","scaleGen","spca","selpop")

    PAWarn <- c("df2genind")

    ## function exists but is experimental
    if(currFunction %in% PAWarn){
        msg <- paste(currFunction,"is implemented but experimental presence/absence markers")
        warning(msg)
        return()
    }

    ## function not implemented
    if(! currFunction %in% PAOk){
        msgError <- paste(currFunction,"is not implemented for presence/absence markers")
        stop(msgError)
    } else return() # else, ok.
} # end checkType






##########
## transp
##########
## AUXIL FUNCTION TO USE TRANSPARENT COLORS
transp <- function(col, alpha=.5){
    res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
}



##########
## corner
##########
## AUXIL FUNCTION TO ADD LETTER TO A PLOT
corner <- function(text, posi="topleft",  inset=0.1, ...){
    oxpd <- par("xpd")
    on.exit(par(xpd=oxpd))
    par(xpd=TRUE)
    myUsr <- par("usr")
    xrange <- myUsr[1:2]
    yrange <- myUsr[3:4]
    x.size <- abs(diff(xrange))
    y.size <- abs(diff(yrange))
    inset <- rep(inset, length=2)
    x.inset <- inset[1]
    y.inset <- inset[2]

    if(length(grep("top", posi))==1){
        y <- yrange[2] - y.size*y.inset
    } else {
        y <- yrange[1] + y.size*y.inset
    }

    if(length(grep("right", posi))==1){
        x <- xrange[2] - x.size*x.inset
    } else {
        x <- xrange[1] + x.size*x.inset
    }

    text(x, y, lab=text, ...)
}





###########
## num2col
###########
## translate numeric values into colors of a palette
num2col <- function(x, col.pal=heat.colors, reverse=FALSE,
                    x.min=min(x,na.rm=TRUE), x.max=max(x,na.rm=TRUE),
                    na.col="transparent"){
    ## if(any(is.na(x))) warning("NAs detected in x")
    x[x < x.min] <- x.min
    x[x > x.max] <- x.max
    x <- x-x.min # min=0
    x.max <- x.max-x.min # update x.max
    x <- x/x.max # max=1
    x <- round(x*100)
    x[x<=0] <- 1
    if(!reverse) {
        pal <- col.pal(100)
    } else {
        pal <- rev(col.pal(100))
    }

    res <- pal[x]
    res[is.na(res)] <- na.col
    return(res)
}





###########
## fac2col
###########
## translate a factor into colors of a palette
## colors are randomized based on the provided seed
fac2col <- function(x, col.pal=funky, na.col="transparent", seed=NULL){
    ## get factors and levels
    x <- factor(x)
    lev <- levels(x)
    nlev <- length(lev)

    ## get colors corresponding to levels
    if(!is.null(seed)){
        set.seed(seed)
        newseed <- round(runif(1,1,1e9))
        on.exit(set.seed(newseed))
        col <- sample(col.pal(nlev))
    } else {
        col <- col.pal(nlev)
    }

    ## get output colors
    res <- rep(na.col, length(x))
    res[!is.na(x)] <- col[as.integer(x[!is.na(x)])]

    ## return
    return(res)
}


###########
## any2col
###########
any2col <- function(x, col.pal=seasun, na.col="transparent"){
    ## handle numeric data
    if(is.numeric(x)){
        col <- num2col(x, col.pal=col.pal, na.col=na.col)
        leg.col <- num2col(pretty(x), x.min=min(x, na.rm=TRUE),
                           x.max=max(x, na.rm=TRUE), col.pal=col.pal,
                           na.col=na.col)
        leg.txt <- pretty(x)
    } else{ ## handle factor
        x <- factor(x)
        col <- fac2col(x, col.pal=col.pal, na.col=na.col)
        leg.col <- col.pal(length(levels(x)))
        leg.txt <- levels(x)
    }

    return(list(col=col, leg.col=leg.col, leg.txt=leg.txt))
} # end any2col



## pre-defined palettes ##
## mono color
bluepal <- colorRampPalette(c("#F7FBFF","#DEEBF7","#C6DBEF",
                              "#9ECAE1","#6BAED6","#4292C6",
                              "#2171B5","#08519C","#08306B"))
redpal <- colorRampPalette(c("#FFF5F0","#FEE0D2","#FCBBA1",
                             "#FC9272","#FB6A4A","#EF3B2C",
                             "#CB181D","#A50F15","#67000D"))
greenpal <- colorRampPalette(c("#F7FCF5","#E5F5E0","#C7E9C0",
                               "#A1D99B","#74C476","#41AB5D",
                               "#238B45","#006D2C","#00441B"))
greypal <- colorRampPalette(c("#FFFFFF","#F0F0F0","#D9D9D9",
                              "#BDBDBD","#969696","#737373",
                              "#525252","#252525","#000000"))

## bi-color
flame <- colorRampPalette(c("gold","red3"))
azur <- colorRampPalette(c("gold","royalblue"))

## tri-color
seasun <- colorRampPalette(c("blue","gold","red"))
lightseasun <- colorRampPalette(c("deepskyblue2","gold","red1"))
deepseasun <- colorRampPalette(c("blue2","gold","red2"))
wasp <-  colorRampPalette(c("yellow2","brown", "black"))
spectral <- colorRampPalette(c("#D53E4F","#F46D43","#FDAE61",
                               "#FEE08B","#FFFFBF","#E6F598",
                               "#ABDDA4","#66C2A5","#3288BD"))

## psychedelic
funky <- colorRampPalette(c("#A6CEE3","#1F78B4","#B2DF8A",
                            "#33A02C","#FB9A99","#E31A1C",
                            "#FDBF6F","#FF7F00","#CAB2D6",
                            "#6A3D9A","#FFFF99","#B15928"))
