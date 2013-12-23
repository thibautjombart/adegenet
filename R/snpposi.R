################
## snpposi.test
################

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

