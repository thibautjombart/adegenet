###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################


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





#'
#' Functions to access online resources for adegenet
#'
#' These functions simply open websites or documents
#' available online providing resources for adegenet.
#'
#' \itemize{
#' \item adegenetWeb opens adegenet's website
#' \item adegenetTutorial opens adegenet tutorials
#' \item adegenetIssues opens the issue page on github;
#' this is used to report a bug or post a feature request.
#' }
#'
#' @param which a character string indicating which tutorial to open (see details)
#'
#' @details
#'
#' Available tutorials are:
#' \itemize{
#' \item 'basics': general introduction to adegenet; covers basic data structures, import/export, handling, and a number of population genetics methods
#' \item 'spca': spatial genetic structures using the spatial Principal Component Analysis
#' \item 'dapc': population structure using the Discriminant Analysis of Principal Components
#' \item 'genomics': handling large genome-wide SNP data using adegenet
#' \item 'strata': introduction to hierarchical population structure in adegenet
#' \item 'snapclust': introduction to fast maximum-likelihood genetic clustering using snapclust
#' }
#'
#' @export
#'
#' @rdname web
#'
adegenetWeb <- function(){
    cat("Opening url \"http://adegenet.r-forge.r-project.org/\" ...\n")
    browseURL("http://adegenet.r-forge.r-project.org/")
}




#' @rdname web
#' @export
adegenetTutorial <- function(which = c("basics","spca","dapc","genomics","strata","snapclust")){

    ## which <- match.arg(which)
    which <- which[1]
    choices <- c("basics","spca","dapc","genomics","strata","snapclust")
    if (!which %in% choices) {
        stop("Unknown tutorial")
    }
    if(which=="basics"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-basics.pdf"
        cat("\n")
        cat("  >> Opening the general introduction to adegenet.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="spca"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-spca.pdf"
        cat("\n")
        cat("  >> Opening the sPCA tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="dapc"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-dapc.pdf"
        cat("\n")
        cat("  >> Opening the DAPC tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="genomics"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-genomics.pdf"
        cat("\n")
        cat("  >> Opening the genomics tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
   if(which=="strata"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-strata.pdf"
        cat("\n")
        cat("  >> Opening the strata tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }
    if(which=="snapclust"){
        url <- "https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-snapclust.pdf"
        cat("\n")
        cat("  >> Opening the snapclust tutorial.\n")
        cat("  >> Seeking url: ",url,"\n ", sep="")
        cat("\n")
    }

    browseURL(url)

    return(invisible(NULL))
}




#' @rdname web
#' @export
adegenetIssues <- function(){
    cat("Opening url \"https://github.com/thibautjombart/adegenet/issues\" ...\n")
    browseURL("https://github.com/thibautjombart/adegenet/issues")
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



.listenToTheBrood <- function(){
    char <- c("?","??","?!?!?")
    for(i in 1:3){
        cat("\nGrind",  char[i], " (y/N): ")
        x <- readLines(con = getOption('adegenet.testcon'), n=1)
        if(x!="y") {
            cat("\n =( \n")
            return(invisible())
        }
    }

    cat("\nGRIIIINNND !!!\n")
    cat("\n\\m/ ^o^ \\m/\n\n")
    url.vec <- c("http://thebrooduk.bandcamp.com/album/swallowed-by-the-earth-single",
                 "http://thebrooduk.bandcamp.com/album/the-brood")
    Sys.sleep(2)
    browseURL(sample(url.vec,1))
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



## viridis
virid <- colorRampPalette(c("#440154FF", "#482173FF", "#433E85FF", "#38598CFF",
                            "#2D708EFF", "#25858EFF", "#1E9B8AFF", "#2BB07FFF",
                            "#51C56AFF", "#85D54AFF", "#C2DF23FF", "#FDE725FF"))


## reorder colors for hybrids
hybridpal <- function(col.pal = virid) {
    function(n) {
        if (n < 3) {
            return(col.pal(n))
        }
        col.pal(n)[c(1, n, 2:(n-1))]
    }
}
