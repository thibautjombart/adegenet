#################
## findMutations
#################

## GENERIC


#' Identify mutations between DNA sequences
#' 
#' The function \code{findMutations} identifies mutations (position and nature)
#' of pairs of aligned DNA sequences. The function \code{graphMutations} does
#' the same thing but plotting mutations on a directed graph.\cr
#' 
#' Both functions are generics, but the only methods implemented in adegenet so
#' far is for \code{\link[ape]{DNAbin}} objects.
#' 
#' 
#' @aliases findMutations findMutations.DNAbin graphMutations
#' graphMutations.DNAbin
#' @param x a \code{DNAbin} object containing aligned sequences, as a matrix.
#' @param from a vector indicating the DNA sequences from which mutations
#' should be found. If \code{NULL}, all sequences are considered (i.e.,
#' \code{1:nrow(x)}).
#' @param to a vector indicating the DNA sequences to which mutations should be
#' found. If \code{NULL}, all sequences are considered (i.e.,
#' \code{1:nrow(x)}).
#' @param allcomb a logical indicating whether all combinations of sequences
#' (from and to) should be considered (TRUE, default), or not (FALSE).
#' @param plot a logical indicating whether the graph should be plotted.
#' @param curved.edges a logical indicating whether the edges of the graph
#' should be curved.
#' @param \dots further arguments to be passed to other methods. Used in
#' \code{graphMutations} where it is passed to the plot method for
#' \code{igraph} objects.
#' @return For \code{findMutations}, a named list indicating the mutations from
#' one sequence to another. For each comparison, a three-column matrix is
#' provided, corresponding to the nucleotides in first and second sequence, and
#' a summary of the mutation provided as: [position]:[nucleotide in first
#' sequence]->[nucleotide in second sequence].
#' 
#' For \code{graphMutations}, a graph with the class \code{igraph}.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}.
#' @seealso The \code{\link{fasta2DNAbin}} to read fasta alignments with
#' minimum RAM use.
#' @examples
#' 
#' \dontrun{
#' data(woodmouse)
#' 
#' ## mutations between first 3 sequences
#' findMutations(woodmouse[1:3,])
#' 
#' ## mutations from the first to sequences 2 and 3
#' findMutations(woodmouse[1:3,], from=1)
#' 
#' ## same, graphical display
#' g <- graphMutations(woodmouse[1:3,], from=1)
#' 
#' ## some manual checks
#' as.character(woodmouse)[1:3,35]
#' as.character(woodmouse)[1:3,36]
#' as.character(woodmouse)[1:3,106]
#' 
#' }
#' 
#' @export findMutations
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



