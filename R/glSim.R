

##########
## glSim
##########
glSim <- function(n.ind, n.snp.nonstruc, n.snp.struc = 0, grp.size = c(0.5, 0.5), k = NULL,
                    pop.freq = NULL, ploidy = 1, alpha = 0, parallel = FALSE,
                    LD = TRUE, block.minsize = 10, block.maxsize = 1000, theta = NULL,
                    sort.pop = FALSE, ...){
  
  
  
  ## BASIC CHECKS ##
  if( any(c(n.ind, n.snp.nonstruc+n.snp.struc) <1)) stop("null numbers of individuals and/or SNPs requested")
  ## alpha parameter
  if(alpha>0.5){
    alpha <- 0.5
    warning("alpha cannot exceed 0.5 - changing the value to 0.5 (total forced asymmetry)")
  }
  if(alpha<0){
    alpha <- 0
    warning("alpha cannot be lower than 0 - changing the value to 0 (no forced asymmetry)")
  }
  
  ## handle group sizes  
  if(length(grp.size)!=2) stop("grp.size should be a vector of length 2")
  grp.size <- grp.size/sum(grp.size)
  grpA.size <- grp.size[[1]]*n.ind
  if(grpA.size >= n.ind) stop("grpA.size is >= n.ind")
  grpB.size <- n.ind - grpA.size
  
  
  # handle pop.freq for k populations
  if(is.null(k) & is.null(pop.freq)){
    pop.freq <- 1
  }
  if(!is.null(k) & !is.null(pop.freq)){
    if(k != length(pop.freq)){
      warning("k != length(pop.freq), length(pop.freq) will be taken as k")
    }
  }
  if(!is.null(k) & is.null(pop.freq)){
    warning("pop.freq will be the result of randomly assorting individuals into k pops")
    pops <-c(1:k) 
    popBaseline <- rep(pops, c(10))
    popBaseline <- factor(sample(popBaseline, length(popBaseline)))
    popExtra <- factor(sample(pops, (n.ind - length(popBaseline)), replace=TRUE))
    pop <- c(popBaseline, popExtra)
    pop <- factor(pop)
    pop.freq <- as.vector(unlist(sapply(pops, function(e) sum(pop==e)))) 
  }
  
  
  n.all <- n.snp.nonstruc
  
  
  
  
  
  simNeutralSNPs <- function(n.ind, n.all, pop.freq=1, LD=TRUE, ploidy=1,
                             block.minsize=10, block.maxsize=1000, theta=NULL,
                             sort.pop=FALSE, parallel=parallel){
    ## CHECKS ##
    ## force population frequencies to 1
    pop.freq <- pop.freq/sum(pop.freq)
    K <- length(pop.freq)
    if(!is.null(theta) && theta<1e-14) theta <- NULL
    
    
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
      return(new("genlight", out, parallel=parallel))
    } # end simBlock.NOLD
    
    
    ## this function simulates a block of allele frequencies for K pop
    ## with linkage
    ## nbAll is a number of alleles to simulate
    ## 'pop' is a factor indicating populations
    
    if(is.null(theta)){
      ## FUNCTION WITHOUT THETA ##
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
        return(out)
      } # end simBlock.LD
    } else {
      ## FUNCTION WITH THETA ##
      ## function to alter frequencies ##
      tweak.freq <- function(f){
        f[f<.5] <- suppressWarnings(f[f<.5] + runif(length(f[f<.5]), 0, theta))
        f[f>.5] <- suppressWarnings(f[f>.5] - runif(length(f[f>.5]), 0, theta))
        f[f<0] <- 0
        f[f>1] <- 1
        return(f)
      }
      
      simBlock.LD <- function(nbAll, pop){
        ## get master allele ##
        f.ori <- runif(K)
        names(f.ori) <- unique(pop)
        prob.ori <- f.ori[pop]
        
        ## function to simulate one allele ##
        simOneAll <- function(){
          return(rbinom(n=n.ind, prob=tweak.freq(prob.ori), size=ploidy))
        }
        
        ## simulate all allele frequencies ##
        out <- replicate(nbAll, simOneAll())
        out <- matrix(out, ncol=nbAll)
        return(out)
      }
    } # end function with theta
    
    
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
      if(length(temp)>1){
        for(i in 2:length(temp)){
          out <- cbind(out, temp[[i]])
        }
      }
      out <- new("genlight", out, parallel=parallel)
    }
    out@other <- list(factor(pop))
    names(out@other) <- "ancestral.pops"
    outpop <- list(out, pop)
    return(outpop)
  } # end simNeutralSNPs
  
  
  ## carry out fn to get non-structural SNPs:
  res.ns <- simNeutralSNPs(n.ind, n.snp.nonstruc, pop.freq=pop.freq, ploidy=ploidy, 
                           LD=LD, block.minsize=block.minsize, block.maxsize=block.maxsize, theta=theta,
                           sort.pop=sort.pop, parallel=FALSE)
  
  
  
  ## GET STRUCTURAL SNPS:
  
  ## GET PHEN FACTOR ##
  phen <- paste("phen", sample(rep(c("A", "B"), c(grpA.size, grpB.size)), 
                               n.ind, replace=TRUE), sep=".")
  if(!sort.pop) phen <- sort(phen)           
  
  ## function simulating structural SNPs for phenotypic groups A and B 
  f2 <- function(){
    probA <- runif(1, min=0, max=0.5-alpha)  # generates probA for SNP p
    probB <- 1 - probA                       # generates probB for SNP p
    # get vector of probabilities by phenotype at SNP p               
    phenProbs <- c(probA, probB)
    names(phenProbs) <- unique(phen)
    # draw SNP p for all individuals according to the phenProbs for that SNP
    rbinom(n=rep(1,n.ind), prob=phenProbs[phen], size=ploidy)
  }
  
  if(n.snp.struc > 0){
    struct <- replicate(n.snp.struc, f2())
    struct <- matrix(struct, ncol=n.snp.struc)
    phen <- factor(phen)  
    struct <- new("genlight", struct, ploidy=ploidy, parallel=FALSE)
    struct@pop <- phen
  }  # end snp.struc
  
  
  pop <- res.ns[[2]]
  res.ns <- res.ns[[1]]
  if(n.snp.struc>0){
    res <- cbind(res.ns, struct)
    res@pop <- phen
  }
  else{
    res <- res.ns
  }
  res@other <- list(factor(pop))
  names(res@other) <- "ancestral.pops"
  return(res)
}   
# 
