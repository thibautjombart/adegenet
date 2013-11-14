

simNeutralSNPs <- function(n.ind, n.all, pop.freq=1, LD=TRUE, ploidy=1,
                           block.minsize=10, block.maxsize=1000, sort.pop=FALSE){
    ## CHECKS ##
    ## force population frequencies to 1
    pop.freq <- pop.freq/sum(pop.freq)
    K <- length(pop.freq)


    ## GET POPULATION FACTOR ##
    pop <- paste("pop", sample(1:K, n.ind, replace=TRUE, prob=pop.freq),sep=".")
    if(sort.pop) pop <- sort(pop)


    ## AUXILIARY FUNCTIONS ##

    ## this function simulates a block of allele frequencies for K pop
    ## without linkage
    ## nbAll is a number of alleles to simulate
    ## 'pop' is a factor indicating populations
    simBlock.NOLD <- function(nbAll, pop){
        ## function to simulate one allele ##
        simOneAll <- function(){
            f <- runif(K)
            names(f) <- unique(pop)
            return(rbinom(n=n.ind, prob=f[pop], size=ploidy))
        }

        ## simulate all allele frequencies ##
        out <- replicate(nbAll, simOneAll())
        return(new("genlight", out))
    } # end simBlock.NOLD


    ## this function simulates a block of allele frequencies for K pop
    ## without linkage
    ## nbAll is a number of alleles to simulate
    ## 'pop' is a factor indicating populations
    simBlock.LD <- function(nbAll, pop){
        ## get master allele ##
        f.ori <- runif(K)
        names(f.ori) <- unique(pop)
        prob.ori <- f.ori[pop]

        ## function to simulate one allele ##
        simOneAll <- function(){
            return(rbinom(n=n.ind, prob=prob.ori, size=ploidy))
        }
        ## simulate all allele frequencies ##
        out <- replicate(nbAll, simOneAll())
        out <- matrix(out, ncol=nbAll)
        return(new("genlight", out))
    } # end simBlock.LD


    ## SIMULATE ALL DATA ##
    if(!LD){ # no LD
        out <- simBlock.NOLD(n.all, pop)
    } else { # with LD
        ## determine blocks ##
        block.sizes <- round(runif(1,min=block.minsize,max=block.maxsize))
        while(sum(block.sizes)<n.all){
            block.sizes <- c(block.sizes, round(runif(1,min=block.minsize,max=block.maxsize)))
        }
        block.sizes <- block.sizes[-length(block.sizes)]
        if(sum(block.sizes)<n.all) block.sizes <- c(block.sizes, n.all-sum(block.sizes))

        ## simulate all blocks ##
        temp <- lapply(block.sizes, simBlock.LD, pop)

        ## put blocks together ##
        out <- temp[[1]]
        if(length(temp)>2){
            for(i in 2:length(temp)){
                out <- cbind(out, temp[[i]])
            }
        }
    }

    return(out)
} # end simNeutralSNPs







## ##########
## ## glSim
## ##########
## glSim <- function(n.ind, n.snp.nonstruc, n.snp.struc=0, grp.size=round(n.ind/2), ploidy=1, alpha=0, block.size=NULL, LD=FALSE, ...){

##     ## BASIC CHECKS ##
##     if( any(c(n.ind, n.snp.nonstruc+n.snp.struc) <1)) stop("null numbers of individuals and/or SNPs requested")
##     ## alpha parameter
##     if(alpha>0.5){
##         alpha <- 0.5
##         warning("alpha cannot exceed 0.5 - changing the value to 0.5 (total forced asymmetry)")
##     }
##     if(alpha<0){
##         alpha <- 0
##         warning("alpha cannot be lower than 0 - changing the value to 0 (no forced asymmetry)")
##     }

##     ## handle block size
##     if(is.null(block.size)){
##         block.size <- n.snp.nonstruc + n.snp.struc
##     }

##     ## handle group sizes
##     grpA.size <- grp.size
##     if(grpA.size >= n.ind) stop("grpA.size is >= n.ind")
##     grpB.size <- n.ind - grpA.size


##     ## AUXIL FUNCTIONS ##
##     if(LD) { # LD - use mvrnorm
##         ## if(!require(MASS)) stop("MASS package is missing.")
##         QUANT <- qnorm(seq(0,1, le=ploidy+2), 0,1) # quantiles needed for continuous->discrete

##         f1 <- function(n,p){
##             Sig <- matrix(runif(p^2), p, p)
##             Sig <- t(Sig) %*% Sig/2 # get the covariance matrix
##             temp <- mvrnorm(n, rep(0,p), Sig) # continuous data
##             temp <- matrix(as.integer(cut(temp, breaks=QUANT))-1, nrow=n, ncol=p)
##             return(new("genlight", temp, ploidy=ploidy, ...))
##         }

##         if(n.snp.struc > 0){
##             if(alpha>1e-10) warning("assymetry parameter alpha ignored when LD used")
##             f2 <- function(n,p){
##                 temp <-rbind(f1(grpA.size,p), f1(grpB.size,p))
##                 ploidy(temp) <- ploidy
##                 return(temp)
##             }
##         }
##     } else { # no LD - use rbinom
##         ## draw p snp for i indiv and convert into a genlight - no structure
##         f1 <- function(n,p){
##             temp <- sapply(1:p, function(i) rbinom(n, ploidy, runif(1)))
##             if(n==1) {temp <- matrix(temp,nrow=1)}
##             return(new("genlight", temp, ploidy=ploidy, ...))
##         }

##         ## draw p snp for i indiv and convert into a genlight - differences between 2 groups
##         if(n.snp.struc > 0){
##             f2 <- function(n,p){
##                 probA <- runif(p, min=0, max=0.5-alpha)
##                 probB <- 1-probA
##                 tempA <- sapply(probA, function(i) rbinom(grpA.size, ploidy, i) )
##                 if(grpA.size==1) {tempA <- matrix(tempA,nrow=1)}
##                 tempB <- sapply(probB, function(i) rbinom(grpB.size, ploidy, i) )
##                 if(grpB.size==1) {tempB <- matrix(tempB,nrow=1)}
##                 return(new("genlight", rbind(tempA,tempB), ploidy=ploidy, ...))
##             }
##         }
##     }



##     ## NON-STRUCTURED DATA ##
##     ## generate data
##     if(n.snp.nonstruc <= block.size){ # no need to use blocks
##         res.ns <- f1(n.ind, n.snp.nonstruc)
##     } else { # proceed by blocks of SNPs
##         res.ns <- f1(n.ind, block.size)
##         snp.to.go <- n.snp.nonstruc - nLoc(res.ns) # number of SNPs left to simulate
##         while(snp.to.go > 0){
##             nb.snp.tmp <- min(snp.to.go, block.size)
##             res.ns <- cbind(res.ns, f1(n.ind, nb.snp.tmp))
##             snp.to.go <- n.snp.nonstruc - nLoc(res.ns) # number of SNPs left to simulate
##         }
##     } # end non-structured data


##     ## STRUCTURED DATA ##
##     if(n.snp.struc > 0){
##         if(n.snp.struc <= block.size){ # no need to use blocks
##             res.s <- f2(n.ind, n.snp.struc)
##         } else { # proceed by blocks of SNPs
##             res.s <- f2(n.ind, block.size)
##             snp.to.go <- n.snp.struc - nLoc(res.s) # number of SNPs left to simulate
##             while(snp.to.go > 0){
##                 nb.snp <- min(snp.to.go, block.size)
##                 res.s <- cbind(res.s, f2(n.ind, nb.snp))
##                 snp.to.go <- n.snp.struc - nLoc(res.s) # number of SNPs left to simulate
##             }

##         }
##     } # end structured data

##     ## RETURN RESULTS ##
##     res <- res.ns
##     if(n.snp.struc > 0){
##         res <- cbind(res,res.s)
##     }
##     if(!is.null(grp.size)){
##         pop(res) <- rep(c("A","B"), c(grpA.size, grpB.size))
##     }
##     indNames(res) <- paste("ind", 1:n.ind)

##     return(res)

## } # end glSim

