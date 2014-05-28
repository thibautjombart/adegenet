## AVOID:
## airDistPlot.dist: no visible binding for global variable ‘groups’
## pairDistPlot.dist: no visible binding for global variable ‘distance’
if(getRversion() >= "2.15.1")  utils::globalVariables(c("groups","distance"))



##############
## GENERICS ##
##############
pairDistPlot <-  function (x, ...) UseMethod("pairDistPlot")
pairDist <-  function (x, ...) UseMethod("pairDistPlot")





#############
## DEFAULT ##
#############
pairDistPlot.default <- function(x, ...){
    stop(paste("No method for objects of class",class(x)))
} # end pairDistPlot.default


pairDist.default <- function(x, grp, within=FALSE, sep="-", ...){
    temp <- pairDistPlot(x=x, grp=grp, within=within, sep=sep,
                         data=TRUE, violin=FALSE, boxplot=FALSE, jitter=FALSE)
    return(temp$data)
}


##########
## DIST ##
##########
##
## this is the basic method
##
pairDistPlot.dist <- function(x, grp, within=FALSE, sep="-", data=TRUE, violin=TRUE, boxplot=TRUE,
                              jitter=TRUE, ...){
    ## CHECKS ##
    if(attr(x, "Size")!=length(grp)) stop("inconsistent length for grp")
    grp <- factor(grp)
    K <- length(levels(K))
    N <- length(grp)


    ## GET DATA FOR OUTPUT AND PLOTTING ##
    ## get groups of pairwise comparisons ##
    allCombs <- combn(N, 2)
    d <- as.vector(x)

    ## remove within if needed ##
    if(!within){
        toKeep <- grp[allCombs[1,]] != grp[allCombs[2,]]
        allCombs <- allCombs[,toKeep,drop=FALSE]
        d <- d[toKeep]
    }

    ## get group-group ##
    d.grp <- paste(grp[allCombs[1,]], grp[allCombs[2,]], sep=sep)


    ## BUILD OUTPUT ##
    out <- list()

    ## data ##
    fig.dat <- data.frame(distance=d, groups=d.grp)
    if(data){
        out$data <- fig.dat
    }

    ## plots ##
    base <- ggplot(data=fig.dat)

    ## violinplot
    if(violin){
        out$violin <- base + geom_violin(aes(x=groups, y=distance, fill=groups), alpha=.5) +
            coord_flip() + guides(fill=FALSE) + labs(x="",y="Pairwise distances")
    }

    ## boxplot
    if(boxplot){
        out$boxplot <- base + geom_boxplot(aes(x=groups, y=distance, fill=groups), alpha=.5) +
            coord_flip() + guides(fill=FALSE) + labs(x="",y="Pairwise distances")
    }

    ## jitter
    if(jitter){
        out$jitter <- base + geom_jitter(aes(x=groups, y=distance, colour=groups), alpha=.2) +
            coord_flip() + guides(colour=FALSE) + labs(x="",y="Pairwise distances")
    }


    return(out)

} # end pairDistPlot.dist







############
## MATRIX ##
############
pairDistPlot.matrix <- function(x, grp, within=FALSE, sep="-", data=TRUE, violin=TRUE, boxplot=TRUE,
                                jitter=TRUE, ...){
    ## CHECKS ##
    if(nrow(x) != ncol(x)) stop("x is not a square matrix")

    ## RETURN ##
    out <- pairDistPlot(as.dist(x), grp=grp, within=within, sep=sep,
                        data=data, violin=violin, boxplot=boxplot, jitter=jitter, ...)

    return(out)
} # end pairDistPlot.matrix






############
## GENIND ##
############
pairDistPlot.genind <- function(x, grp, within=FALSE, sep="-", data=TRUE, violin=TRUE, boxplot=TRUE,
                                jitter=TRUE, ...){
    ## CHECKS ##
    if(missing(grp)){
        if(!is.null(pop(x))) {
            grp <- pop(x)
        } else {
            stop("grp is missing with no population defined in x")
        }
    }


    ## RETURN ##
    D <- dist(x@tab)^2
    out <- pairDistPlot(D, grp=grp, within=within, sep=sep,
                        data=data, violin=violin, boxplot=boxplot, jitter=jitter, ...)

    return(out)
} # end pairDistPlot.matrix






############
## DNAbin ##
############
pairDistPlot.DNAbin <- function(x, grp, within=FALSE, sep="-", data=TRUE, violin=TRUE, boxplot=TRUE,
                                jitter=TRUE, ...){

    ## RETURN ##
    D <- dist.dna(x, ...)
    out <- pairDistPlot(D, grp=grp, within=within, sep=sep,
                        data=data, violin=violin, boxplot=boxplot, jitter=jitter, ...)

    return(out)
} # end pairDistPlot.matrix






