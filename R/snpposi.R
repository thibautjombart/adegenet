################
## snpposi.test
################

## GENERIC
#' @export
snpposi.test <- function(...){
    UseMethod("snpposi.test")
}





## METHOD FOR INTEGER - BASIC METHOD
#' @method snpposi.test integer
#' @export
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
#' @method snpposi.test numeric
#' @export
snpposi.test.numeric <- function(x, ...){
    out <- snpposi.test(as.integer(x), ...)
    return(out)
}





## METHOD FOR DNABIN
#' @method snpposi.test DNAbin
#' @export
snpposi.test.DNAbin <- function(x, ...){
    out <- snpposi.test(x=as.integer(seg.sites(x)),
                        genome.size=ncol(x), ...)
    return(out)
} # end snpposi.test.DNAbin











################
## snpposi.plot
################

## GENERIC
#' @export
snpposi.plot <- function(...){
    UseMethod("snpposi.plot")
}




## METHOD FOR INTEGER - BASIC METHOD
#' @method snpposi.plot integer
#' @export
snpposi.plot.integer <- function(x, genome.size, smooth=0.1, col="royalblue", alpha=.2,
                                 codon=TRUE, start.at=1, ...){
    ## IF WE REPRESENT DENSITY PER CODON POSITION ##
    if(codon){
        ## define base positions (1/2/3) ##
        codon.posi <- ((2 + x) %% 3) + 1
        fac <- factor(codon.posi, levels=1:3)

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





#' @method snpposi.plot numeric
#' @export
snpposi.plot.numeric <- function(x, ...){
    out <- snpposi.plot(as.integer(x), ...)
    return(out)
}




#' @method snpposi.plot DNAbin
#' @export
snpposi.plot.DNAbin <- function(x, ...){
    out <- snpposi.plot(x=as.integer(seg.sites(x)),
                        genome.size=ncol(x), ...)
    return(out)
} # end snpposi.plot.DNAbin

