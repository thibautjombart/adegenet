## ############
## ## haploPop
## ############
## ##
## ## Simulate only SNPs, allow reverse mutations.
## ##
## ## - haplo.length: length of simulated haplotypes
## ## - mu: substitution rate / nucleotide / year
## ## - n.steps: number of generations to simulate
## ##
## haploPop <- function(n.steps=20, ini.obj=NULL, ini.haplo=NULL, haplo.length=1e6, mu=1e-5, n.snp.ini=1,
##                      birth.func=function(){ sample(0:3, 1, prob=c(.05, .45, .35, .15))},
##                      max.pop.size=function(){1e4}, max.nb.pop=30, ini.pop.size=10, regen=FALSE,
##                      p.new.pop=function(){1e-4}, death.func=function(age){age>1},
##                      quiet=FALSE, allow.reverse=TRUE) {


##     ## SOME CHECKS
##     ## if(is.numeric(ini.pop.size)){
##     ##     ini.pop.size.val <- ini.pop.size
##     ##     ini.pop.size <- function(){ini.pop.size.val}
##     ## }

##     if(is.numeric(max.pop.size)){
##         max.pop.size.val <- max.pop.size
##         max.pop.size <- function(){max.pop.size.val}
##     }

##      if(is.numeric(p.new.pop)){
##         p.new.pop.val <- p.new.pop
##         p.new.pop <- function(){p.new.pop.val}
##     }

##     if(is.numeric(birth.func)){
##         birth.func.val <- birth.func[1]
##         birth.func <- function(){birth.func.val}
##     }

##     if(is.numeric(death.func)){
##         death.func.val <- death.func[1]
##         death.func <- function(age){age>death.func.val}
##     }


##     ## GLOBAL VARIABLES ##
##     SNP.POOL <- 1:haplo.length
##     vecS <- 1 # will be redefined later, but needed for evolveOnePop definition

##     ## AUXILIARY FUNCTIONS ##
##     if(allow.reverse){
##         createMutations <- function(N){ # L:genome length; N: pop size
##             nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
##             return( sample(SNP.POOL, size=nb.mutations, replace=TRUE) )
##         }
##     } else {
##         createMutations <- function(N){ # L:genome length; N: pop size
##             nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
##             res <- sample(SNP.POOL, size=nb.mutations, replace=TRUE)
##             SNP.POOL <<- setdiff(SNP.POOL, res)# update pool of SNPs
##             return(res)
##         }
##     }

##     ## clean reverse mutations
##     cleanRes <- function(vec){
##         temp <- table(vec)
##         return( as.integer(names(temp)[temp %% 2 != 0]) )
##     }


##     ## assign mutation to haplotypes
##     assignMutations <- function(myPop, mutations){ # mypop: list of `haplotypes'; mutations: vector of SNPs
##         if(length(mutations)==0 | length(myPop)==0) return(myPop)
##         id <- sample(1:length(myPop), size=length(mutations), replace=TRUE)
##         mutations <- split(mutations, id)

##         ## function to merge new mutations - handle reverse case
##         f1 <- function(a,b){
##             revMut <- intersect(a,b)
##             if(length(revMut)==0) return(c(a,b))
##             return(setdiff(c(a ,b), revMut))
##         }

##         ##myPop[as.integer(names(mutations))] <- mapply(c, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)
##          myPop[as.integer(names(mutations))] <- mapply(f1, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)

##         return(myPop)
##     } # end assignMutations


##     if(!regen){
##         ## VERSION FOR NO REGENERATION OF SUSCEPTIBLES
##         evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
##             ## strains get older
##             myAge <- myAge + 1
##             ## toKill <- death.func(myAge)
##             ## myPop[toKill] <- NULL
##             ## myAge <- myAge[!toKill]

##             ## generate new strains for new generation
##             sampSize <- round(min( length(myPop)*birth.func(), myS)) # number of strains for next step
##             if(sampSize<1){ # if no new strains
##                 ## old strains die
##                 toKill <- death.func(myAge)
##                 myPop[toKill] <- NULL
##                 myAge <- myAge[!toKill]
##                 return(list(pop=myPop, S=myS, age=myAge))
##             } # if there are new strains, do...
##             newGen <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains for new generations
##             newGen <- assignMutations(newGen, createMutations(sampSize)) # mutate strains
##             newAge <- rep(0, sampSize) # new ages for newborns

##             ## old strains die
##             toKill <- death.func(myAge)
##             myPop[toKill] <- NULL
##             myAge <- myAge[!toKill]

##             ## merge old and new generation
##             myPop <- c(myPop,newGen)
##             myAge <- c(myAge, newAge)

##             ## possibly create one or more new pop
##             if((length(listPop) < max.nb.pop) & (p.new.pop()>0)) { # total number of pop. limitation
##                 nbNewPop <- rbinom(1, length(myPop), prob=p.new.pop())
##             } else {
##                 nbNewPop <- 0
##             }
##             if(nbNewPop>0){
##                 ## newPop <- sample(listPop, size=nbNewPop, replace=TRUE) # wrong
##                 newPop <- lapply(sample(myPop, size=nbNewPop, replace=TRUE), as.list)
##                 listPop <<- c(listPop, newPop)
##                 vecS <<- c(vecS, replicate(nbNewPop, max.pop.size()) )
##                 listAges <<- c(listAges, replicate(nbNewPop, 0, simplify=FALSE) )
##             } # end new pop
##             return(list(pop=myPop, S=myS-sampSize, age=myAge))
##         } # end no regen version
##     } else { ## REGEN VERSION
##         evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
##             ## strains get older
##             myAge <- myAge + 1
##             ## toKill <- death.func(myAge)
##             ## myPop[toKill] <- NULL
##             ## myAge <- myAge[!toKill]
##             myS <- max.pop.size() ## DIFFERENCE between the two versions of the function

##             ## generate new strains for new generation
##             sampSize <- round(min( length(myPop)*birth.func(), myS)) # number of strains for next step
##             if(sampSize<1){ # if no sample
##                 ## old strains die
##                 toKill <- death.func(myAge)
##                 myPop[toKill] <- NULL
##                 myAge <- myAge[!toKill]
##                 return(list(pop=myPop, S=myS, age=myAge))
##             }
##             newGen <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains for new generations
##             newGen <- assignMutations(newGen, createMutations(sampSize)) # mutate strains
##             newAge <- rep(0, sampSize) # new ages for newborns

##             ## old strains die
##             toKill <- death.func(myAge)
##             myPop[toKill] <- NULL
##             myAge <- myAge[!toKill]

##             ## merge old and new generation
##             myPop <- c(myPop,newGen)
##             myAge <- c(myAge, newAge)

##             ## possibly create one or more new pop
##             if((length(listPop) < max.nb.pop) & (p.new.pop()>0)) { # total number of pop. limitation
##                 nbNewPop <- rbinom(1, length(myPop), prob=p.new.pop())
##             } else {
##                 nbNewPop <- 0
##             }
##             if(nbNewPop>0){
##                 ## newPop <- sample(listPop, size=nbNewPop, replace=TRUE) # wrong
##                 newPop <- lapply(sample(myPop, size=nbNewPop, replace=TRUE), as.list)
##                 listPop <<- c(listPop, newPop)
##                 vecS <<- c(vecS, replicate(nbNewPop, max.pop.size()) )
##                 listAges <<- c(listAges, replicate(nbNewPop, 0, simplify=FALSE) )
##             } # end new pop
##             return(list(pop=myPop, S=myS, age=myAge)) ## DIFFERENCE between the two versions of the function
##         } # end no regen version
##     } ## end evolveOnePop (both versions)



##     ## INITIATE SIMULATIONS ##
##     ## INITIALIZE FROM SCRATCH
##     vecS <- max.pop.size() # susceptibles

##     if(is.null(ini.obj)){
##         ##vecS <- max.pop.size() -  n.snp.ini # susceptibles
##         if(is.null(ini.haplo)) {
##             haplo.ini <- sample(SNP.POOL, n.snp.ini, replace=TRUE)
##         } else {
##             haplo.ini <- ini.haplo
##         }

##         ANCES <- haplo.ini
##         listPop <- list()
##         listPop[[1]] <- lapply(1:ini.pop.size, function(i) haplo.ini) # contains only one population of identical clones to start with
##         listAges <- list() # will contain vectors of ages of haplotypes (a time of appearance, age=0)
##         listAges[[1]] <- rep(0, ini.pop.size)
##     } else { ## INITIALIZE WITH PROVIDED OBJECT
##         if(!inherits(ini.obj, "haploPop")) stop("x is not a haploPop object")
##         ##vecS <- ini.obj$S
##         ANCES <- attr(ini.obj, "ances")
##         listPop <- ini.obj$pop
##         listAges <- ini.obj$ages
##     }

##     ## MAKE SIMULATIONS ##

##     ## evolve all populations
##     i <- 1L
##     if(!quiet){
##         cat("\nSimulating populations of haplotypes through time: \n")
##     }
##     ##while((sum(vecS)>0) & (i<(n.steps+1))){ # evolve all generations
##     while(i<(n.steps+1)){ # evolve all generations
##         i <- i + 1L # update iterator
##         if(!quiet){
##             catStep <- max(round(n.steps/100), 10)
##             cat(ifelse((i %% catStep)==0, paste(" ...", i), ""))
##         }


##         ## make populations evolve of one generation
##         ##idx <- which(vecS>0) # make sure that new pop won't evolve this time ! leads to not dying
##         idx <- 1:length(listPop)  # make sure that new pop won't evolve this time
##         if(length(idx)>0){
##             for(j in idx){
##                 temp <- evolveOnePop(listPop[[j]], vecS[j], listAges[[j]])
##                 listPop[[j]] <- temp$pop
##                 vecS[j] <- temp$S
##                 listAges[[j]] <- temp$age
##             }
##         }

##         ## ## purge non-susceptible pop
##         ## listPop <- listPop[vecS>0]
##         ## vecS <- vecS[vecS>0]

##         ## purge empty populations
##         toKeep <- sapply(listPop, length)>0
##         listPop <- listPop[toKeep]
##         vecS <- vecS[toKeep]
##         listAges <- listAges[toKeep]

##         ## stop if all pop go extinct
##         if(length(listPop)==0L){
##             if(!quiet) cat("\n All populations went extinct at time",i,"\n")
##             return(invisible(NULL))
##         }

##         ## FOR DEBUGGING
##         ## cat("\n=== ",i," ===")
##         ## cat("\nlistPop")
##         ## print(listPop)
##         ## cat("\nvecS")
##         ## print(vecS)
##         ## cat("\nlistAges")
##         ## print(listAges)
##         ## END DEBUGGING
##     } # end while

##     if(!quiet){
##         cat("\n... done! \n")
##     }

##     ## END OF SIMULATIONS ##


##     ## CLEAN RESULTS ##
##     ## handle reverse mutations
##     ## if(clean.haplo){
##     ##     if(!quiet){
##     ##         cat("\n... Cleaning haplotypes (handling reverse mutations)\n")
##     ##     }

##     ##     cleanRes <- function(vec){
##     ##         temp <- table(vec)
##     ##         return(sort(as.integer(names(temp)[temp %% 2 != 0])))
##     ##     }

##     ##     for(i in 1:length(listPop)){
##     ##         listPop[[i]] <- lapply(listPop[[i]], cleanRes)
##     ##     }

##     ##     if(!quiet){
##     ##         cat("\n... done! \n")
##     ##     }
##     ## }

##     ## RETURN RESULTS ##
##     res <- list(pop=listPop, ages=listAges, S=vecS)
##     class(res) <- "haploPop"
##     res$call <- match.call()
##     attr(res,"ances") <- ANCES # ancestral genotype
##     return(res)

## } # end haploPop






## ##################
## ## print.haploPop
## ##################
## print.haploPop <- function(x, ...){
##     myCall <- x$call

##     cat("\n== haploPop object ==\n")
##     cat("\nNumber of populations :", length(x$pop))

##     N <- sum(sapply(x$pop,length))
##     cat("\nNumber of haplotypes :", N)

##     N.mut <- length(unique(unlist(x$pop)))
##     cat("\nNumber of mutations :", N.mut)

##     N.empty <- sum(sapply(x$pop, function(e) length(e)==0))
##     cat("\nNumber of unmutated genotypes :", N.empty)

##     if( (length(x$pop) == length(x$ages)) & (length(x$pop) == length(x$S)) ){
##         cat("\nSlot lengths consistency:  OK\n")
##     } else {
##         cat("\nSlot lengths consistency: !! NOT OK !!\n")
##     }
## } # end print.haploPop






## ##################
## ## summary.haploPop
## ##################
## summary.haploPop <- function(object, ...){
##     x <- object$pop
##     myCall <- x$call
##     x$call <- NULL
##     res <- list()

##     ## cat("\t\n=======================================")
##     ## cat("\t\n= simulated populations of haplotypes =")
##     ## cat("\t\n=          (haploPop object)          =")
##     ## cat("\t\n=======================================\n")

##     cat("\nNumber of populations :", length(x))

##     cat("\nPopulation sizes :\n")
##     temp <- sapply(x,length)
##     names(temp) <- 1:length(temp)
##     print(temp)
##     res$pop.size <- temp

##     cat("\nNumber of SNPs per population :\n")
##     temp <- sapply(x,function(e) length(unique(unlist(e))))
##     names(temp) <- 1:length(temp)
##     print(temp)
##     res$n.snp <- temp

##     return(invisible(res))
## } # end print.haploPop






## ##################
## ## sample.haploPop
## ##################
## sample.haploPop <- function(x, n, n.pop=NULL, keep.pop=TRUE){
##     if(!inherits(x, "haploPop")) stop("x is not a haploPop object")
##     x$call <- NULL

##     if(!is.null(n.pop)){ # pre-treatment: reduce to n.pop populations with same size
##         ## kEEP ONLY SOME POP
##         popToKeep <- sample(which(sapply(x$pop, length) > n), n.pop, replace=FALSE) # keep n.pop large enough populations
##         if(length(popToKeep)==0L) stop("No population is big enough for this sampling.")
##         x$pop <- x$pop[popToKeep]
##         x$ages <- x$ages[popToKeep]
##         x$S <- x$S[popToKeep]

##         ## MAKE THEM THE SAME SIZE
##         popSizes <- sapply(x$pop, length)
##         for(i in 1:n.pop){
##             idx <- sample(1:popSizes[i], n, replace=FALSE)
##             x$pop[[i]] <- x$pop[[i]][idx]
##             x$ages[[i]] <- x$ages[[i]][idx]
##         }

##     } # end pop pre-treatment

##     if(keep.pop){
##         popSizes <- sapply(x$pop, length)
##         pop.id <- rep(1:length(x$pop), popSizes)
##     }

##     x$pop <- unlist(x$pop, recursive=FALSE)
##     x$ages <- unlist(x$ages, recursive=FALSE)

##     idx <- sample(1:length(x$pop), n, replace=FALSE)
##     res <- list(pop=list(), ages=list() )

##     if(keep.pop){
##         res$pop <- split(x$pop[idx], pop.id[idx])
##         res$ages <- split(x$ages[idx], pop.id[idx])
##     } else {
##         res$pop[[1]] <- x$pop[idx]
##         res$ages[[1]] <- x$ages[idx]
##     }

##     res$S <- rep(n, length(res$pop))

##     class(res) <- "haploPop"
##     attr(res, "ances") <- attr(x, "ances")
##     return(res)
## } # end sample.haploPop






## ###############
## ## dist.haploPop
## ###############
## dist.haploPop <- function(x, add.root=TRUE, res.type=c("dist","matrix")){
##     if(!inherits(x, "haploPop")) stop("x is not a haploPop object")
##     res.type <- match.arg(res.type)
##     ANCES <- attr(x,"ances")

##     x <- unlist(x$pop, recursive=FALSE)

##     ## handle root
##     if(add.root){ # add the root
##        x <- c(ANCES, x)
##     }

##     n <- length(x)

##     f1 <- function(a,b){
##         return(sum(!union(unlist(a),unlist(b)) %in% intersect(unlist(a),unlist(b))))
##     }

##     ## res <- outer(x, x, FUN=f1)
##     res <- matrix(0, ncol=n, nrow=n)
##     for(i in 1:(n-1)){
##         for(j in (i+1):n){
##             res[i,j] <- f1(x[[i]], x[[j]])
##         }
##     }

##     res <- res+t(res)

##     if(res.type=="dist"){
##         res <- as.dist(res)
##     }
##     return(res)
## } # end dist.haploPop






## ###############
## ## plot.haploPop
## ###############
## plot.haploPop <- function(x, y=NULL, type="unrooted", size.limit=300, show.pop=TRUE, col=NULL,
##                           transp=TRUE, tip.cex=2, method=c("nj", "bionj", "fastme.bal", "fastme.ols"), ...){
##     ## CHECKS ##
##     if(!require(ape)) stop("ape package is required")
##     if(!inherits(x, "haploPop")) stop("x is not a haploPop object")
##     method <- match.arg(method)

##     N <- sum(sapply(x$pop,length))

##     if(N > size.limit) {
##         stop("tree exceeds size limit")
##     }


##     ## PLOT TREE ##
##     f1 <- get(method)
##     if(method %in% c("nj","bionj")){
##         tre <- root(f1(dist.haploPop(x)),"1")
##     } else {
##         tre <- f1(dist.haploPop(x))
##      }
##     plot(tre, type=type, ...)
##     xy <- get("last_plot.phylo", envir = .PlotPhyloEnv)


##     ## SHOW POPULATIONS ##
##     if(!is.null(col)){
##         if(is.integer(col) | is.numeric(col)) {
##             col <- palette()[col]
##         }
##         if(transp){
##             transp <- function(col, alpha=.5){
##                 res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
##                 return(res)
##             }

##             col <- transp(col)
##         }
##           points(xy$xx[2:(N+1)], xy$yy[2:(N+1)], pch=20, col=col, cex=tip.cex)

##     } else if(show.pop){
##         nPop <- length(x$pop)
##         popSizes <- sapply(x$pop, length)
##         pop.id <- rep(1:length(x$pop), popSizes)
##         opal <- palette()
##         on.exit(palette(opal))
##         if(nPop>1){
##             pop.col <- rainbow(nPop)
##         } else {
##             pop.col <- c("red","red")
##         }

##         if(transp){
##             transp <- function(col, alpha=.5){
##                 res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
##                 return(res)
##             }

##             pop.col <- transp(pop.col)
##         }

##         palette(pop.col)
##         points(xy$xx[2:(N+1)], xy$yy[2:(N+1)], pch=20, col=pop.id, cex=tip.cex)
##     }


##     ## SHOW ROOT ##
##     points(xy$xx[1], xy$yy[1], pch=20, cex=3)

##     return(invisible(tre))
## } # end plot.haploPop
























## ##########################################################################
## ##########################################################################
## ##########################################################################
## ##########################################################################
## ##########################################################################
## ##########################################################################
## ##########################################################################























## ############
## ## haploPopDiv
## ############
## haploPopDiv <- function(n.steps=20, ini.obj=NULL, ini.haplo=NULL, haplo.length=1e6, mu=1e-5, n.snp.ini=1,
##                         birth.func=function(){ sample(0:3, 1, prob=c(.05, .45, .35, .15))},
##                         max.pop.size=function(){1e4}, max.nb.pop=30, ini.pop.size=10, regen=FALSE,
##                         p.new.pop=function(){1e-4}, death.func=function(age){age>1},
##                         quiet=FALSE, allow.reverse=TRUE,
##                         track=c("div", "distRoot", "freq","nbMut"), root.haplo=NULL, samp.size=50) {


##     ## SOME CHECKS
##     ## if(is.numeric(ini.pop.size)){
##     ##     ini.pop.size.val <- ini.pop.size
##     ##     ini.pop.size <- function(){ini.pop.size.val}
##     ## }

##     track <- match.arg(track)

##     if(is.numeric(max.pop.size)){
##         max.pop.size.val <- max.pop.size
##         max.pop.size <- function(){max.pop.size.val}
##     }

##      if(is.numeric(p.new.pop)){
##         p.new.pop.val <- p.new.pop
##         p.new.pop <- function(){p.new.pop.val}
##     }

##     if(is.numeric(birth.func)){
##         birth.func.val <- birth.func[1]
##         birth.func <- function(){birth.func.val}
##     }

##     if(is.numeric(death.func)){
##         death.func.val <- death.func[1]
##         death.func <- function(age){age>death.func.val}
##     }


##     ## GLOBAL VARIABLES ##
##     SNP.POOL <- 1:haplo.length
##     vecS <- 1 # will be redefined later, but needed for evolveOnePop definition

##     ## AUXILIARY FUNCTIONS ##
##     if(allow.reverse){
##         createMutations <- function(N){ # L:genome length; N: pop size
##             nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
##             return( sample(SNP.POOL, size=nb.mutations, replace=TRUE) )
##         }
##     } else {
##         createMutations <- function(N){ # L:genome length; N: pop size
##             nb.mutations <- sum(rbinom(N, size=haplo.length, prob=mu))
##             res <- sample(SNP.POOL, size=nb.mutations, replace=TRUE)
##             SNP.POOL <<- setdiff(SNP.POOL, res)# update pool of SNPs
##             return(res)
##         }
##     }


##     ## assign mutation to haplotypes
##     assignMutations <- function(myPop, mutations){ # mypop: list of `haplotypes'; mutations: vector of SNPs
##         if(length(mutations)==0 | length(myPop)==0) return(myPop)
##         id <- sample(1:length(myPop), size=length(mutations), replace=TRUE)
##         mutations <- split(mutations, id)

##         ## function to merge new mutations - handle reverse case
##         f1 <- function(a,b){
##             revMut <- intersect(a,b)
##             if(length(revMut)==0) return(c(a,b))
##             return(setdiff(c(a ,b), revMut))
##         }

##         ##myPop[as.integer(names(mutations))] <- mapply(c, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)
##          myPop[as.integer(names(mutations))] <- mapply(f1, myPop[as.integer(names(mutations))], mutations, SIMPLIFY=FALSE)

##         return(myPop)
##     } # end assignMutations


##     if(!regen){
##         ## VERSION FOR NO REGENERATION OF SUSCEPTIBLES
##         evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
##             ## strains get older
##             myAge <- myAge + 1
##             ## toKill <- death.func(myAge)
##             ## myPop[toKill] <- NULL
##             ## myAge <- myAge[!toKill]

##             ## generate new strains for new generation
##             sampSize <- round(min( length(myPop)*birth.func(), myS)) # number of strains for next step
##             if(sampSize<1){ # if no sample
##                 ## old strains die
##                 toKill <- death.func(myAge)
##                 myPop[toKill] <- NULL
##                 myAge <- myAge[!toKill]
##                 return(list(pop=myPop, S=myS, age=myAge))
##             }
##             newGen <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains for new generations
##             newGen <- assignMutations(newGen, createMutations(sampSize)) # mutate strains
##             newAge <- rep(0, sampSize) # new ages for newborns

##             ## old strains die
##             toKill <- death.func(myAge)
##             myPop[toKill] <- NULL
##             myAge <- myAge[!toKill]

##             ## merge old and new generation
##             myPop <- c(myPop,newGen)
##             myAge <- c(myAge, newAge)

##             ## possibly create one or more new pop
##             if((length(listPop) < max.nb.pop) & (p.new.pop()>0)) { # total number of pop. limitation
##                 nbNewPop <- rbinom(1, length(myPop), prob=p.new.pop())
##             } else {
##                 nbNewPop <- 0
##             }
##             if(nbNewPop>0){
##                 ## newPop <- sample(listPop, size=nbNewPop, replace=TRUE) # wrong
##                 newPop <- lapply(sample(myPop, size=nbNewPop, replace=TRUE), as.list)
##                 listPop <<- c(listPop, newPop)
##                 vecS <<- c(vecS, replicate(nbNewPop, max.pop.size()) )
##                 listAges <<- c(listAges, replicate(nbNewPop, 0, simplify=FALSE) )
##             } # end new pop
##             return(list(pop=myPop, S=myS-sampSize, age=myAge))
##         } # end no regen version
##     } else { ## REGEN VERSION
##         evolveOnePop <- function(myPop, myS, myAge){ # myPop: pop to evolve; myS: nb of susceptible in the pop; myAge: vector of ages
##             ## strains get older
##             myAge <- myAge + 1
##             ## toKill <- death.func(myAge)
##             ## myPop[toKill] <- NULL
##             ## myAge <- myAge[!toKill]
##             myS <- max.pop.size() ## DIFFERENCE between the two versions of the function

##             ## generate new strains for new generation
##             sampSize <- round(min( length(myPop)*birth.func(), myS)) # number of strains for next step
##             if(sampSize<1){ # if no sample
##                 ## old strains die
##                 toKill <- death.func(myAge)
##                 myPop[toKill] <- NULL
##                 myAge <- myAge[!toKill]
##                 return(list(pop=myPop, S=myS, age=myAge))
##             }
##             newGen <- myPop[sample(1:length(myPop), sampSize, replace=TRUE)] # sample strains for new generations
##             newGen <- assignMutations(newGen, createMutations(sampSize)) # mutate strains
##             newAge <- rep(0, sampSize) # new ages for newborns

##             ## old strains die
##             toKill <- death.func(myAge)
##             myPop[toKill] <- NULL
##             myAge <- myAge[!toKill]

##             ## merge old and new generation
##             myPop <- c(myPop,newGen)
##             myAge <- c(myAge, newAge)

##             ## possibly create one or more new pop
##             if((length(listPop) < max.nb.pop) & (p.new.pop()>0)) { # total number of pop. limitation
##                 nbNewPop <- rbinom(1, length(myPop), prob=p.new.pop())
##             } else {
##                 nbNewPop <- 0
##             }
##             if(nbNewPop>0){
##                 ## newPop <- sample(listPop, size=nbNewPop, replace=TRUE) # wrong
##                 newPop <- lapply(sample(myPop, size=nbNewPop, replace=TRUE), as.list)
##                 listPop <<- c(listPop, newPop)
##                 vecS <<- c(vecS, replicate(nbNewPop, max.pop.size()) )
##                 listAges <<- c(listAges, replicate(nbNewPop, 0, simplify=FALSE) )
##             } # end new pop
##             return(list(pop=myPop, S=myS, age=myAge)) ## DIFFERENCE between the two versions of the function
##         } # end no regen version
##     } ## end evolveOnePop (both versions)


##     ## INITIATE SIMULATIONS ##
##     ## INITIALIZE FROM SCRATCH
##     vecS <- max.pop.size() # susceptibles

##     if(is.null(ini.obj)){
##         ## vecS <- max.pop.size() -  n.snp.ini # susceptibles

##         if(is.null(ini.haplo)) {
##             haplo.ini <- sample(SNP.POOL, n.snp.ini, replace=TRUE)
##         } else {
##             haplo.ini <- ini.haplo
##         }
##         ANCES <- haplo.ini
##         listPop <- list()
##         listPop[[1]] <- lapply(1:ini.pop.size, function(i) haplo.ini) # contains only one population of identical clones to start with
##         listAges <- list() # will contain vectors of ages of haplotypes (a time of appearance, age=0)
##         listAges[[1]] <- rep(0, ini.pop.size)
##     } else { ## INITIALIZE WITH PROVIDED OBJECT
##         if(!inherits(ini.obj, "haploPop")) stop("x is not a haploPopDiv object")
##         ## vecS <- ini.obj$S
##         ANCES <- attr(ini.obj, "ances")
##         listPop <- ini.obj$pop
##         listAges <- ini.obj$ages
##     }

##     ## function getting pairwise distances
##     if(track=="div"){
##         fRes <- function(list.pop){
##             list.pop <- list(pop=list.pop) # kludge needed for dist.haploPop
##             class(list.pop) <- "haploPop" # kludge needed for dist.haploPop
##             N <- sum(sapply(list.pop$pop, length))
##             if(N<2) return(0)
##             if(N > samp.size){
##                 return(dist.haploPop(sample.haploPop(list.pop, samp.size, keep.pop=FALSE), add.root=FALSE)) # do not include the root in distances.
##             } else {
##                 return(dist.haploPop(list.pop, add.root=FALSE))
##             }
##         } # end fRes
##     }

##     ## function getting distances to the root
##     if(track=="distRoot"){
##         if(is.null(root.haplo)) {
##             root.haplo <- ANCES
##         }
##         fRes <- function(list.pop){
##             list.pop <- list(pop=list.pop) # kludge needed for sample.haploPop
##             class(list.pop) <- "haploPop" # kludge needed for sample.haploPop
##             N <- sum(sapply(list.pop$pop, length))
##             if(N<1) return(0)
##             if(N > samp.size){
##                 list.pop <- sample.haploPop(list.pop, samp.size, keep.pop=FALSE)
##             }

##             res <- sapply(unlist(list.pop$pop, recursive=FALSE), function(e) sum(!e %in% root.haplo))
##             return(res)
##         } # end fRes
##     }

##     ## function getting allele absolute frequencies
##     if(track=="freq"){
##         fRes <- function(list.pop){
##             res <- table(unlist(list.pop))
##             return(res)
##         } # end fRes
##     }

##     ## function getting allele absolute frequencies
##     if(track=="nbMut"){
##         fRes <- function(list.pop){
##             list.pop <- list(pop=list.pop) # kludge needed for sample.haploPop
##             class(list.pop) <- "haploPop" # kludge needed for sample.haploPop
##             N <- sum(sapply(list.pop$pop, length))
##             if(N<1) return(0)
##             if(N > samp.size){
##                 list.pop <- sample.haploPop(list.pop, samp.size, keep.pop=FALSE)
##             }

##             return( length(unique(unlist(list.pop))) )
##         } # end fRes

##     }



##     res <- list(div=list(), popSize=integer())
##     res$div[[1]] <- fRes(listPop)
##     res$popSize[1] <- sum(sapply(listPop, length))


##     ## MAKE SIMULATIONS ##

##     ## evolve all populations
##     i <- 1L
##     if(!quiet){
##         cat("\nSimulating populations of haplotypes through time: \n")
##     }
##     ##while((sum(vecS)>0) & (i<(n.steps+1))){ # evolve all generations
##     while(i<(n.steps+1)){ # evolve all generations
##         i <- i + 1L # update iterator
##         if(!quiet){
##             catStep <- max(round(n.steps/100), 10)
##             cat(ifelse((i %% catStep)==0, paste(" ...", i), ""))
##         }


##         ## make populations evolve of one generation
##         ##idx <- which(vecS>0) # make sure that new pop won't evolve this time ! leads to not dying
##         idx <- 1:length(listPop)  # make sure that new pop won't evolve this time
##         if(length(idx)>0){
##             for(j in idx){
##                 temp <- evolveOnePop(listPop[[j]], vecS[j], listAges[[j]])
##                 listPop[[j]] <- temp$pop
##                 vecS[j] <- temp$S
##                 listAges[[j]] <- temp$age
##             }

##         }

##         ## ## purge non-susceptible pop
##         ## listPop <- listPop[vecS>0]
##         ## vecS <- vecS[vecS>0]

##         ## purge empty populations
##         toKeep <- sapply(listPop, length)>0
##         listPop <- listPop[toKeep]
##         vecS <- vecS[toKeep]
##         listAges <- listAges[toKeep]

##         ## stop if all pop go extinct
##         if(length(listPop)==0L){
##             if(!quiet) cat("\n All populations went extinct at time",i,"\n")
##             return(res)
##         }

##         res$div[[i]] <- fRes(listPop)
##         res$popSize[i] <- sum(sapply(listPop, length))

##         ## FOR DEBUGGING
##         ## cat("\n=== ",i," ===")
##         ## cat("\nlistPop")
##         ## print(listPop)
##         ## cat("\nvecS")
##         ## print(vecS)
##         ## cat("\nlistAges")
##         ## print(listAges)
##         ## END DEBUGGING
##     } # end while

##     if(!quiet){
##         cat("\n... done! \n")
##     }

##     ## END OF SIMULATIONS ##

##     ## STORE HAPLOPOP OBJECT
##      obj <- list(pop=listPop, ages=listAges, S=vecS)
##     class(obj) <- "haploPop"
##     obj$call <- match.call()
##     attr(obj,"ances") <- ANCES # ancestral genotype

##     if(!quiet) cat("\nStored haploPop object in 'last.haploPop'\n")
##     assign("last.haploPop", obj, envir= .GlobalEnv)


##     ## RETURN RES
##     return(res)

## } # end haploPopDiv


