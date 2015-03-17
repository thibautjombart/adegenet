
#' Analyse the position of polymorphic sites
#'
#' These functions are used to describe the distribution of polymorphic sites
#' (SNPs) in an alignment.
#'
#' The function \code{snpposi.plot} plots the positions and density of SNPs in
#' the alignment.
#'
#' The function \code{snpposi.test} tests whether SNPs are randomly distributed
#' in the genome, the alternative hypothesis being that they are clustered.
#' This test is based on the distances of each SNP to the closest SNP. This
#' provides one measure of clustering for each SNP.  Different statistics can
#' be used to summarise these values (argument \code{stat}), but by default the
#' statistics used is the median.
#'
#' \code{snpposi.plot} and \code{snpposi.test} are generic functions with
#' methods for vectors of integers or numeric (indicating SNP position), and
#' for \code{\link[ape]{DNAbin}} objects.
#'
#'
#' @aliases snpposi.plot snpposi.plot.integer snpposi.plot.numeric
#' snpposi.plot.DNAbin snpposi.test snpposi.test.integer snpposi.test.numeric
#' snpposi.test.DNAbin
#' @param x a vector of integers or numerics containing SNP positions, or a set
#' of aligned sequences in a \code{DNAbin} object.
#' @param genome.size an integer indicating the length of genomes.
#' @param smooth a smoothing parameter for the density estimation; smaller
#' values will give more local peaks; values have to be positive but can be
#' less than 1.
#' @param col the color to be used for the plot; ignored if codon positions are
#' represented.
#' @param alpha the alpha level to be used for transparency (density curve).
#' @param codon a logical indicating if codon position should be indicated
#' (TRUE, default) or not.
#' @param start.at an integer indicating at which base of a codon the alignment
#' starts (defaults to 1); values other than 1, 2 and 3 will be ignored.
#' @param n.sim an integer indicating the number of randomizations to be used
#' in the Monte Carlo test.
#' @param stat a function used to summarize the measure of physical proximity
#' between SNPs; by default, the median is used.
#' @param \dots further arguments to be passed to the \code{integer} method.
#' @return A Monte Carlo test of the class \code{randtest}.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}.
#' @seealso The \code{\link{fasta2DNAbin}} to read fasta alignments with
#' minimum RAM use.
#' @examples
#'
#' if(require(ape)){
#' data(woodmouse)
#' snpposi.plot(woodmouse, codon=FALSE)
#' snpposi.plot(woodmouse)
#'
#' \dontrun{
#' snpposi.test(c(1,3,4,5), 100)
#' snpposi.test(woodmouse)
#' }
#' }
#'


## GENERIC
snpposi.test <- function(...){
    UseMethod("snpposi.test")
}





## METHOD FOR INTEGER - BASIC METHOD
snpposi.test.integer <- function(x, genome.size, n.sim=999, stat=median, ...){
    ## AUXILIARY FUNCTION ##
    ## computes the test statistics for one vector of SNP positions
    f1 <- function(snpposi, stat=median){
        temp <- as.matrix(dist(snpposi))
        diag(temp) <- 1e15
        out <- stat(apply(temp, 1, min))
        return(out)
    }

    ## GET OBSERVATION ##
    obs <- f1(x, stat=stat)

    ## GET SIMULATIONS ##
    n.snps <- length(x)
    sim <- sapply(1:n.sim, function(e) f1(sample(1:genome.size, n.snps, replace=FALSE), stat=stat))

    ## MAKE RANDTEST OUTPUT ##
    out <- as.randtest(obs=obs, sim=sim, alter="less")
    return(out)
} # end snpposi.test.integer





## METHOD FOR NUMERIC
snpposi.test.numeric <- function(x, ...){
    out <- snpposi.test(as.integer(x), ...)
    return(out)
}





## METHOD FOR DNABIN
snpposi.test.DNAbin <- function(x, ...){
    out <- snpposi.test(x=as.integer(seg.sites(x)),
                        genome.size=ncol(x), ...)
    return(out)
} # end snpposi.test.DNAbin











################
## snpposi.plot
################

## GENERIC
snpposi.plot <- function(...){
    UseMethod("snpposi.plot")
}




## METHOD FOR INTEGER - BASIC METHOD
snpposi.plot.integer <- function(x, genome.size, smooth=0.1, col="royalblue", alpha=.2,
                                 codon=TRUE, start.at=1, ...){
    ## IF WE REPRESENT DENSITY PER CODON POSITION ##
    if(codon){
        ## define base positions (1/2/3) ##
        fac <- rep(1:3, length=genome.size)
        if(start.at==2) fac <- c(2:3,fac)[1:genome.size]
        if(start.at==3) fac <- c(3,fac)[1:genome.size]
        fac <- factor(fac, levels=1:3)
        fac <- fac[x]

        ## make ggplot output ##
        out <- ggplot(data.frame(x=x, codon=fac), aes(x=x)) + xlim(0, genome.size)
        out <- out + geom_density(adjust=smooth, aes(fill=codon, colour=codon),alpha=I(alpha)) + geom_rug(aes(colour=codon),alpha=.7)
        out <- out + labs(x="Nucleotide position", title="Distribution of SNPs in the genome")
        out <- out + guides(fill=guide_legend(title="Codon position"), colour=guide_legend(title="Codon position"))
    } else {
        ## OTHERWISE, JUST ONE DENSITY ##
        ## make ggplot output ##
        out <- ggplot(data.frame(x=x), aes(x=x)) + xlim(0, genome.size)
        out <- out + geom_density(adjust=smooth, fill=transp(col,alpha=alpha), colour=col) + geom_rug(colour=col,alpha=.7)
        out <- out + labs(x="Nucleotide position", title="Distribution of SNPs in the genome")
    }


    ## return ##
    return(out)
} # end snpposi.plot.integer





## METHOD FOR NUMERIC
snpposi.plot.numeric <- function(x, ...){
    out <- snpposi.plot(as.integer(x), ...)
    return(out)
}




## METHOD FOR DNABIN
snpposi.plot.DNAbin <- function(x, ...){
    out <- snpposi.plot(x=as.integer(seg.sites(x)),
                        genome.size=ncol(x), ...)
    return(out)
} # end snpposi.plot.DNAbin

