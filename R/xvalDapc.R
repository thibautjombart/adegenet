
##############
## xvalDapc ##
##############

# Return randomly sampled indices from a group.
# @param e group name
# @param vector of group assignments per sample
# @param training.set fraction of samples to be kept for validation
.group_sampler <- function(e, grp, training.set){
  group_e   <- grp == e 
  N_group_e <- sum(group_e)
  if (N_group_e < 2){
    # If the size of the group is less than two, then leave the whole thing in 
    # the training set
    return(which(group_e))
  } else {
    samp_group <- which(group_e)
    samp_size <- round(training.set * N_group_e)
    return(sample(samp_group, size = samp_size))    
  }
}


# Function to subsample the data. This is to be used as the ran.gen function
# in the boot function. DATA is a data frame or matrix containing the samples,
# GRP is the group identities of the samples, PCA is the result of dudi.pca
# on the full data set, KEEP is the subset of samples based on the training
# set (mle). Note that this function only has two inputs, dat, and mle.
.boot_group_sampler <- function(dat = list(DATA = NULL, GRP = NULL, PCA = NULL, 
                                           KEEP = NULL), 
                                mle = NULL){
  to_keep    <- unlist(lapply(levels(dat$GRP), .group_sampler, dat$GRP, mle))
  dat$PCA$li <- dat$PCA$li[to_keep, , drop = FALSE]
  dat$KEEP   <- to_keep
  return(dat)
}


# Function to pass to the "statistic" parameter of boot. This will subset the
# data, calculate the DAPC, give the predictions and return the results.
.boot_dapc_pred <- function(x, n.pca = n.pca, n.da = n.da, result = "overall"){
  if (length(x$KEEP) == nrow(x$DATA)){
    out <- 1
  } else {
    new_dat   <- x$DATA[-x$KEEP, ,drop = FALSE]
    train_dat <- x$DATA[x$KEEP, ,drop = FALSE]
    
    new_grp   <- x$GRP[-x$KEEP]
    train_grp <- x$GRP[x$KEEP]      
    
    temp.dapc <- suppressWarnings(dapc(train_dat, train_grp, dudi = x$PCA, 
                                       n.pca = n.pca, n.da = n.da))
    temp.pred <- predict.dapc(temp.dapc, newdata = new_dat)
    if (result=="overall"){
      out <- mean(temp.pred$assign == new_grp)
    }
    if (result=="groupMean"){
      out <- mean(tapply(temp.pred$assign == new_grp, new_grp, mean), na.rm = TRUE)
    }
  }
  return(out)
}


# Function that will actually do the bootstrapping. It will return a numeric
# vector with the successes ratios. Note that the ellipses are used to pass
# parameters to boot. When implemented in the xvalDapc function, this will allow
# the user to implement this in parallel. 
# @param n.pca number of pcs
# @param x data frame/matrix with samples in rows
# @param n.da number of das
# @param groups the names of each of the populations
# @param grp factor of group assignments per sample
# @param training.set fraction of samples used for training
# @param training.set2 NULL or largest possible fraction that can be obtained
# @param pcaX principal componenets
# @param result user's choice of result type
# @param reps the number of replicates per number of retained PCs
# @param ... methods to be passed on to boot such as parallel and ncores
.get.prop.pred <- function(n.pca, x, n.da, groups, grp, training.set, 
                           training.set2, pcaX, result = "overall", reps = 100,
                           ...){
  
  groups_ge_ten <- vapply(groups, function(e) sum(grp == e) >= 10, logical(1))
  bootlist      <- list(DATA = x, GRP = grp, PCA = pcaX, KEEP = 1:nrow(x))
  if (all(groups_ge_ten) == TRUE){
    out <- boot::boot(bootlist, .boot_dapc_pred, sim = "parametric", R = reps,
                      ran.gen = .boot_group_sampler, mle = training.set, 
                      n.pca = n.pca, n.da = n.da, result = result, ...)$t
  } else {
    out <- boot::boot(bootlist, .boot_dapc_pred, sim = "parametric", R = reps,
                      ran.gen = .boot_group_sampler, mle = training.set2, 
                      n.pca = n.pca, n.da = n.da, result = result, ...)$t      
  }
  return(as.vector(out))
}




xvalDapc <- function(x, grp, n.pca.max = 300, n.da = NULL, training.set = 0.9, 
                     result = c("groupMean", "overall"), center = TRUE, scale = FALSE, 
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE, ...){
  
  ## CHECKS ##
  grp <- factor(grp)
  n.pca <- n.pca[n.pca > 0]
  if(!is.null(n.da)){
    n.da <- n.da
  }else{
    n.da <- length(levels(grp))-1}
  #   if(missing(n.da)){
  #   n.da <- length(levels(grp))-1}
  #   if(is.null(n.da)){
  #     n.da <- length(levels(grp))-1} # need to fix this to make interactive n.da selection an option! 
  #   else{
  #     n.da <- n.da}
  if(missing(training.set)){
    training.set <- 0.9
  }else{
    training.set <- training.set}
  if(missing(n.rep)){
    n.rep <- 30
  }else{
    n.rep <- n.rep}
  
  if(missing(result)){
    result <- "groupMean"
  }else{ 
    if(length(result) > 1) result <- result[1]
    result <- result}
  
  
  ## GET TRAINING SET SIZE ##
  N <- nrow(x)
  groups <- levels(grp)
  if (all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
    N.training <- round(N*training.set)
  } else {
    groups1 <- (levels(grp))[(as.vector(which.min((lapply(groups, function(e) sum(as.vector(unclass(grp==e))))))))]
    popmin <- length(which(grp%in%groups1))
    if(popmin==1){
      ## exclude smallest group; proceed with second smallest group:
      counter <- 0
      while(popmin==1){
        groups.temp <- factor(as.vector(grp[-which(grp %in% groups1)]), exclude=groups1)
        groups1.ori <- groups1
        groups1 <- levels(groups.temp)[(as.vector(which.min((lapply(levels(groups.temp), 
                                                              function(e) sum(as.vector(unclass(groups.temp==e))))))))]
        popmin <- length(which(groups.temp%in%groups1))
        groups1 <- c(groups1.ori, groups1)
        counter <- sum(counter, 1)           
      }
      if(counter==1){        
        msg <- "1 group has only 1 member so it cannot be represented in both training and validation sets."
      }else{
        msg <- paste(counter, "groups have only 1 member: these groups cannot be represented in both training and validation sets.") 
                      
      }      
      warning(msg)
    }
    training.set2 <- ((popmin - 1)/popmin)
    N.training <- round(N*training.set2)   
  }
  
  
  ## GET FULL PCA ##
  if(missing(n.pca.max)) n.pca.max <- min(dim(x))
  pcaX <- dudi.pca(x, nf=n.pca.max, scannf=FALSE, center=center, scale=scale)
  n.pca.max <- min(n.pca.max, pcaX$rank, N.training-1) # re-defines n.pca.max (so user's input may not be the value used...)
  
  ## DETERMINE N.PCA IF NEEDED ##
  if(n.pca.max < 10){
    runs <- n.pca.max
  }else{
    runs <- 10}
  
  if(is.null(n.pca)){
    n.pca <- round(pretty(1:n.pca.max, runs))
  }

  n.pca <- n.pca[n.pca>0 & n.pca<(N.training-1) & n.pca<n.pca.max]

  ## FUNCTION GETTING THE % OF ACCURATE PREDICTION FOR ONE NUMBER OF PCA PCs ##
  ## n.pca is a number of retained PCA PCs
  # get.prop.pred <- function(n.pca){
  #   f1 <- function(){
  #     if (all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
  #       toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), 
  #                                                          size=(round(training.set*sum(as.vector(unclass(grp==e))))))))
  #     } else {
  #       ## check if any groups have only 1 member
  #       grp.n.1 <- levels(grp)[which(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))==1)]                
  #       if(length(grp.n.1)!=0){
  #         ## for groups with 1 member, keep that member in the training set
  #         toKeep <- which(grp %in% grp.n.1) 
  #         toKeep <- c(toKeep, unlist(lapply(groups[-which(groups %in% grp.n.1)], function(e) sample(which(grp==e), 
  #                                                     size=(round(training.set2*sum(as.vector(unclass(grp==e)))))))))
  #       }else{
  #         ## if no group has only 1 member, proceed normally
  #         toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), 
  #                                                            size=(round(training.set2*sum(as.vector(unclass(grp==e))))))))
  #       }
  #       }
  #     temp.pca <- pcaX
  #     temp.pca$li <- temp.pca$li[toKeep,,drop=FALSE]
  #     temp.dapc <- suppressWarnings(dapc(x[toKeep,,drop=FALSE], grp[toKeep], 
  #                                        n.pca=n.pca, n.da=n.da, dudi=temp.pca))
  #     temp.pred <- predict.dapc(temp.dapc, newdata=x[-toKeep,,drop=FALSE])
  #     if(result=="overall"){
  #       out <- mean(temp.pred$assign==grp[-toKeep])
  #     }
  #     if(result=="groupMean"){
  #       out <- mean(tapply(temp.pred$assign==grp[-toKeep], grp[-toKeep], mean), na.rm=TRUE)
  #     }
  #     return(out)
  #   }
  #   return(replicate(n.rep, f1()))
  # } # end get.prop.pred
  
  
  ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
  res.all <- unlist(lapply(n.pca, .get.prop.pred, x, n.da, groups, grp,
                           training.set, training.set2, pcaX, result, 
                           n.rep, ...))
  xval <- data.frame(n.pca=rep(n.pca, each=n.rep), success=res.all)    
  
  
  n.pcaF <- as.factor(xval$n.pca)
  successV <- as.vector(xval$success)
  pca.success <- tapply(successV, n.pcaF, mean)
  n.opt <- which.max(tapply(successV, n.pcaF, mean))
  
  
  ###### MSE-BASED OPTIMAL n.pca SELECTION:
  temp <- seq(from=1, to=length(xval$n.pca), by=n.rep)
  orary <-c(temp+(n.rep-1))
  index <-c(1:length(temp))
  lins <-sapply(index, function(e) seq(from=temp[e], to=orary[e]))
  lin <-c(1:ncol(lins))
  col <-successV
  cait<-sapply(lin, function(e) ((col[lins[,e]])-1)^2)
  FTW <-sapply(lin, function(e) sum(cait[,e])/n.rep)
  RMSE <- sqrt(FTW)
  names(RMSE) <- xval$n.pca[temp]
  best.n.pca <- names(which.min(RMSE))
  
  # DAPC
  n.pca <- as.integer(best.n.pca)
  n.da <- nlevels(grp)-1
  dapc1 <- suppressWarnings(dapc(x, grp, n.pca=n.pca, n.da=n.da))
  
  # PLOT CROSS-VALIDATION RESULTS
  snps <- x
  phen <- grp
  random <- replicate(300, mean(tapply(sample(phen)==phen, phen, mean)))
  q.phen <- quantile(random, c(0.025,0.5,0.975))
  
  if(xval.plot==TRUE){
    smoothScatter(xval$n.pca, successV, nrpoints=Inf, pch=20, col=transp("black"),
                  ylim=c(0,1), xlab="Number of PCA axes retained",
                  ylab="Proportion of successful outcome prediction", 
                  main="DAPC Cross-Validation")
    abline(h=q.phen, lty=c(2,1,2))
  }
  
  
  # RESULTS
  xvalResults <- list(xval, q.phen, pca.success, (names(n.opt)), RMSE, best.n.pca, dapc1)
  names(xvalResults)[[1]] <- "Cross-Validation Results"
  names(xvalResults)[[2]] <- "Median and Confidence Interval for Random Chance"
  names(xvalResults)[[3]] <- "Mean Successful Assignment by Number of PCs of PCA"
  names(xvalResults)[[4]] <- "Number of PCs Achieving Highest Mean Success"
  names(xvalResults)[[5]] <- "Root Mean Squared Error by Number of PCs of PCA"
  names(xvalResults)[[6]] <- "Number of PCs Achieving Lowest MSE"
  names(xvalResults)[[7]] <- "DAPC"
  
  
  return(xvalResults)
  
} # end xvalDapc.data.frame


xvalDapc.data.frame <- xvalDapc
xvalDapc.matrix <- xvalDapc.data.frame

