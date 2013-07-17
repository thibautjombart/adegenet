
############
## xvalDapc
############

xvalDapc <- function (x, ...) UseMethod("xvalDapc")

xvalDapc.data.frame <- function(x, grp, n.pca.max, n.da=NULL, training.set = 0.9,
                                result=c("groupMean","overall"),
                                center=TRUE, scale=FALSE, n.pca=NULL, n.rep=10, ...){

    ## CHECKS ##
    grp <- factor(grp)
    n.pca <- n.pca[n.pca>0]
    result <- match.arg(result)
    if(is.null(n.da)) {
        n.da <- length(levels(grp))-1
    }

    ##added:
    if(missing(training.set)){
      training.set <- 0.9}
    else{
      training.set <- training.set}

    ## GET TRAINING SET SIZE ##
    N <- nrow(x)
    
    ##deleted: N.training <- round(N*training.set)
    ##added: 
    groups <- levels(grp)
    if(all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
      N.training <- round(N*training.set)}
    else{
      groups1 <- (levels(grp))[(as.vector(which.min((lapply(groups, function(e) sum(as.vector(unclass(grp==e))))))))]
      popmin <- length(which(grp%in%groups1))
      training.set2 <- ((popmin - 1)/popmin)
      N.training <- round(N*training.set2)}
    

    ## GET FULL PCA ##
    if(missing(n.pca.max)) n.pca.max <- min(dim(x))
    pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
    n.pca.max <- min(n.pca.max,pcaX$rank,N.training-1)

    ## DETERMINE N.PCA IF NEEDED ##
    if(is.null(n.pca)){
        n.pca <- round(pretty(1:n.pca.max,10))
    }
    n.pca <- n.pca[n.pca>0 & n.pca<(N.training-1)]

    ## FUNCTION GETTING THE % OF ACCURATE PREDICTION FOR ONE NUMBER OF PCA PCs ##
    ## n.pca is a number of retained PCA PCs
    ##deleted: VOID.GRP <- FALSE # will be TRUE if empty group happened
    get.prop.pred <- function(n.pca){
        f1 <- function(){
            ##added:
                  if(all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
                    toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), size=(round(training.set*sum(as.vector(unclass(grp==e))))))))}
                  else{
                    toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), size=(round(training.set2*sum(as.vector(unclass(grp==e))))))))}
            ##deleted: toKeep <- sample(1:N, N.training)
            ##deleted: if(!(all(table(grp[toKeep])>0) & all(table(grp[-toKeep])>0))) VOID.GRP <<- TRUE
            temp.pca <- pcaX
            temp.pca$li <- temp.pca$li[toKeep,,drop=FALSE]
            temp.dapc <- suppressWarnings(dapc(x[toKeep,,drop=FALSE], grp[toKeep], n.pca=n.pca, n.da=n.da, dudi=temp.pca))
            temp.pred <- predict.dapc(temp.dapc, newdata=x[-toKeep,,drop=FALSE])
            if(result=="overall"){
                out <- mean(temp.pred$assign==grp[-toKeep])
            }
            if(result=="groupMean"){
                out <- mean(tapply(temp.pred$assign==grp[-toKeep], grp[-toKeep], mean), na.rm=TRUE)
            }
            return(out)
        }
        return(replicate(n.rep, f1()))
    }


    ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
    res.all <- unlist(lapply(n.pca, get.prop.pred))
    ##deleted: if(VOID.GRP) warning("At least one group was absent from the training / validating sets.\nTry using smaller training sets.")
    res <- data.frame(n.pca=rep(n.pca, each=n.rep), success=res.all)
    return(res)
} # end xvalDapc.data.frame


xvalDapc.matrix <- xvalDapc.data.frame






## #############
## ## discriVal
## #############

## discriVal <- function (x, ...) UseMethod("discriVal")

## discriVal.data.frame <- function(x, grp, n.pca.max, n.da=NULL, center=TRUE, scale=FALSE, n.pca=NULL, ...){

##     ## CHECKS ##
##     grp <- factor(grp)
##     n.pca <- n.pca[n.pca>0]
##     if(is.null(n.da)) {
##         n.da <- length(levels(grp))-1
##     }

##     ## GET FULL PCA ##
##     if(missing(n.pca.max)) n.pca.max <- min(dim(x))
##     pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
##     n.pca.max <- min(n.pca.max,pcaX$rank)

##     ## DETERMINE N.PCA IF NEEDED ##
##     if(is.null(n.pca)){
##         n.pca <- round(pretty(1:n.pca.max,10))
##     }
##     n.pca <- n.pca[n.pca>0 & n.pca<n.pca.max]

##     ## FUNCTION GETTING THE TOTAL DISCRIMINATION (SUM OF EIGENVALUES) FOR ONE GIVEN NB OF PCA PCs ##
##     ## n.pca is a number of retained PCA PCs
##     get.totdiscr <- function(n.pca){
##             temp.dapc <- suppressWarnings(dapc(x, grp, n.pca=n.pca, n.da=n.da, dudi=pcaX))
##             return(sum(temp.dapc$eig))
##     }


##     ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
##     res.all <- sapply(n.pca, get.totdiscr)
##     res <- data.frame(n.pca=n.pca, success=res.all)
##     return(res)
## } # end discriVal.data.frame


## discriVal.matrix <- discriVal.data.frame




## There's a bunch of problems down there, commenting it for nowé
## xval.dapc <- function(object, n.pca, n.da, training.set = 90, ...){
##   training.set = training.set/100
##   kept.id <- unlist(tapply(1:nInd(object), pop(object), function(e) {pop.size = length(e); pop.size.train = round(pop.size * training.set); sample(e, pop.size.train, replace=FALSE)})) # this can't work: nInd/pop not defined for DAPC objects
##   training <- object[kept.id]
##   validating <- object[-kept.id]
##   post = vector(mode = 'list', length = n.pca)
##   asgn = vector(mode = 'list', length = n.pca)
##   ind = vector(mode = 'list', length = n.pca)
##   mtch = vector(mode = 'list', length = n.pca)
##   for(i in 1:n.pca){
##     dapc.base = dapc(training, n.pca = i, n.da = 15) # Why 15??
##     dapc.p = predict.dapc(dapc.base, newdata = validating)
##     match.prp = mean(as.character(dapc.p$assign)==as.character(pop(validating)))
##     post[[i]] = dapc.p$posterior
##     asgn[[i]] = dapc.p$assign
##     ind[[i]] = dapc.p$ind.score
##     mtch[[i]] = match.prp
##   }
##   res = list(assign = asgn, posterior = post, ind.score = ind, match.prp = mtch)
##   return(res)
## } # end of xval.dapc

## xval.genind  <- function(object, n.pca, n.da, training.set = 90, ...){
##   res = xval.dapc(object = object, n.pca = n.pca, n.da = n.da, training.set = training.set)
##   return(res)
## }
## ###############
## ## randtest.dapc
## ###############
## ##randtest.dapc <- function(x, nperm = 999, ...){

## ##} # end randtest.dapc




######## TESTS IN R #######

## TEST PREDICT.DAPC ##
## data(sim2pop)
## temp <- seppop(sim2pop)
## temp <- lapply(temp, function(e) hybridize(e,e,n=30)) # force equal pop sizes
## hyb <- hybridize(temp[[1]], temp[[2]], n=30)
## newdat <- repool(temp[[1]], temp[[2]], hyb)
## pop(newdat) <- rep(c("pop A", "popB", "hyb AB"), c(30,30,30))


## ##dapc1 <- dapc(newdat[1:61],n.pca=10,n.da=1)
## dapc1 <- dapc(newdat[1:60],n.pca=2,n.da=1)
## scatter(dapc1)
## hyb.pred <- predict(dapc1, newdat[61:90])

## scatter(dapc1)
## points(hyb.pred$ind.scores, rep(.1, 30))

## assignplot(dapc1, new.pred=hyb.pred)
## title("30 indiv popA, 30 indiv pop B, 30 hybrids")
