

#################
## findMutations
#################

## GENERIC
findMutations <- function(...){
    UseMethod("findMutations")
}



## METHOD FOR DNABIN
findMutations.DNAbin <- function(x, from=NULL, to=NULL, allcomb=TRUE, ...){
    ## CHECKS ##
    if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
    x <- as.matrix(x)

    ## function to pull out mutations from sequence a to b ##
    NUCL <- c('a','t','g','c')
    f1 <- function(a,b){
        seqa <- as.character(x[a,])
        seqb <- as.character(x[b,])
        temp <- which(seqa != seqb)
        ori <- seqa[temp]
        mut <- seqb[temp]
        names(ori) <- names(mut) <- temp
        toRemove <- !ori %in% NUCL | !mut %in% NUCL
        ori <- ori[!toRemove]
        mut <- mut[!toRemove]
        if(all(toRemove)) return(NULL)
        res <- data.frame(ori,mut)
        names(res) <- rownames(x)[c(a,b)]
        res$short <- paste(row.names(res),":",res[,1],"->",res[,2],sep="")
        return(res)
    }

    ## GET LIST OF PAIRS TO COMPARE ##
    ## handle from/to as character
    if(is.character(from)) from <- match(from, rownames(x))
    if(is.character(to)) to <- match(to, rownames(x))

    ## handle NULL
    if(is.null(from)) from <- 1:nrow(x)
    if(is.null(to)) to <- 1:nrow(x)

    ## get pairs
    if(allcomb){
        pairs <- expand.grid(to, from)[,2:1,drop=FALSE]
    } else {
        N <- max(length(from),length(to))
        from <- rep(from, length=N)
        to <- rep(to, length=N)
        pairs <- cbind(from, to)
    }

    ## remove unwanted comparisons
    pairs <- pairs[pairs[,1]!=pairs[,2],,drop=FALSE]

    ## GET NUMBER OF MUTATIONS ##
    out <- lapply(1:nrow(pairs), function(i) f1(pairs[i,1], pairs[i,2]))
    names(out) <- paste(rownames(x)[pairs[,1]], rownames(x)[pairs[,2]],sep="->")

    return(out)

} # end findMutations







##################
## graphMutations
##################

## GENERIC
graphMutations <- function(...){
    UseMethod("graphMutations")
}



## METHOD FOR DNABIN
graphMutations.DNAbin <- function(x, from=NULL, to=NULL, allcomb=TRUE, plot=TRUE, curved.edges=TRUE, ...){

    ## GET MUTATIONS ##
    x <- findMutations(x, from=from, to=to, allcomb=allcomb)

    ## GET GRAPH ##
    from <- gsub("->.*","",names(x))
    to <- gsub(".*->","",names(x))
    vnames <- sort(unique(c(from,to)))
    dat <- data.frame(from,to,stringsAsFactors=FALSE)
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(vnames, label=vnames))

    ## SET ANNOTATIONS FOR THE BRANCHES ##
    annot <- unlist(lapply(x, function(e) paste(e$short, collapse="\n")))
    E(out)$label <- annot
    E(out)$curved <- curved.edges

    ## PLOT / RETURN ##
    if(plot) plot(out, ...)

    return(out)
} # end graphMutations



