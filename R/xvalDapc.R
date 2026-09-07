

#' @export
#' 
xvalDapc <- function (x, ...) UseMethod("xvalDapc")

##############
## xvalDapc ##
##############

# Return randomly sampled indices from a group.
# @param e group name
# @param vector of group assignments per sample
# @param training.set fraction of samples to be kept for training
.group_sampler <- function(e, grp, training.set){
  group_e   <- grp == e 
  N_group_e <- sum(group_e)
  if (N_group_e < 2){
    # If the size of the group is less than two, then leave the whole thing in 
    # the training set
    return(which(group_e))
  } else {
    samp_group <- which(group_e)
    samp_size <- max(1L, min(N_group_e - 1L,
                            round(training.set * N_group_e)))
    return(sample(samp_group, size = samp_size))    
  }
}


# Function to subsample the data. This is to be used as the ran.gen function
# in the boot function. DATA is a data frame or matrix containing the samples,
# GRP is the group identities of the samples, KEEP is the training subset.
# CENTER and SCALE specify the transformation fitted within each training
# set (mle). Note that this function only has two inputs, dat, and mle.
.boot_group_sampler <- function(dat = list(DATA = NULL, GRP = NULL,
                                           KEEP = NULL), 
                                mle = NULL){
  to_keep    <- unlist(lapply(levels(dat$GRP), .group_sampler, dat$GRP, mle))
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
    # Fit all preprocessing to the training individuals only.
    if (all(vapply(as.data.frame(train_dat), function(z) {
      length(unique(z)) < 2L
    }, logical(1)))) return(NA_real_)
    train_pca <- dudi.pca(train_dat, nf = n.pca, scannf = FALSE,
                          center = x$CENTER, scale = x$SCALE)
    # dapc() otherwise silently reduces n.pca to the available rank.
    # LDA's pooled within-group covariance has rank at most N - K.
    lda_rank_max <- nrow(train_dat) - nlevels(droplevels(train_grp))
    if (min(train_pca$rank, sum(train_pca$eig > 1e-14),
            lda_rank_max) < n.pca) return(NA_real_)
    temp.dapc <- suppressWarnings(dapc(train_dat, train_grp, dudi = train_pca,
                                       n.pca = n.pca, n.da = n.da))
    temp.pred <- predict.dapc(temp.dapc, newdata = new_dat)
    if (identical(result, "overall")){
      out <- mean(temp.pred$assign == new_grp)
    }
    if (identical(result, "groupMean")){
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
# @param center,scale training-set preprocessing options
# @param result user's choice of result type
# @param reps the number of replicates per number of retained PCs
# @param ... methods to be passed on to boot such as parallel and ncores
####################
## .get.prop.pred ##
####################
.get.prop.pred <- function(n.pca, x, n.da, groups, grp, training.set, 
                           center = TRUE, scale = FALSE,
                           result = "overall", reps = 100,
                           ...){
  
  bootlist <- list(DATA = x, GRP = grp, KEEP = seq_len(nrow(x)),
                   CENTER = center, SCALE = scale)

  out <- boot::boot(bootlist, .boot_dapc_pred, sim = "parametric", R = reps,
                      ran.gen = .boot_group_sampler, mle = training.set, 
                      n.pca = n.pca, n.da = n.da, result = result, ...)$t      
  #   }
  return(as.vector(out))
} # end .get.prop.pred




##############
#' @method xvalDapc default
#' @export
##############
xvalDapc.default <- function(x, grp, n.pca.max = 300, n.da = NULL, training.set = 0.9, 
                     result = c("groupMean", "overall"), center = TRUE, scale = FALSE, 
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE, ...){
  
  ## CHECKS ##
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || any(dim(x) == 0L))
    stop("x must be a non-empty numeric matrix or data frame.")
  if (anyNA(x))
    stop("xvalDapc requires complete data; missing values are not imputed. ",
         "Imputation before cross-validation can use validation information; ",
         "use a workflow that fits imputation within each training split.")
  if (any(!is.finite(x))) stop("x must contain finite values.")
  grp <- factor(grp)
  if (length(grp) != nrow(x) || anyNA(grp) || nlevels(grp) < 2L)
    stop("grp must identify at least two groups, with one non-missing label per row.")
  if (length(training.set) != 1L || !is.finite(training.set) ||
      training.set <= 0 || training.set >= 1)
    stop("training.set must be between 0 and 1, exclusively.")
  if (length(n.rep) != 1L || !is.finite(n.rep) || n.rep < 1 || n.rep != floor(n.rep))
    stop("n.rep must be a positive integer.")
  result <- match.arg(result)
  n.pca <- n.pca[n.pca > 0]
  if(!is.null(n.da)){
    n.da <- n.da
  }else{
    n.da <- length(levels(grp))-1}
  #   if(missing(n.da)){
  #   n.da <- length(levels(grp))-1}
  #   if(is.null(n.da)){
  #     n.da <- length(levels(grp))-1} # want to fix this to make interactive n.da selection an option! 
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
  ## identify the sizes of groups
  group.n <- as.vector(table(grp))
  if (all(group.n == 1L)) stop("No individuals are available for validation.")
  ## identify the smallest group size
  popmin <- min(group.n)
  
  ## check if any groups are of length 1:
  if(popmin == 1){
    singles <- which(group.n == 1)
    counter <- length(singles)
    
    ## get msg to print
    if(counter == 1){
      msg <- "1 group has only 1 member so it cannot be represented in both training and validation sets."
    }else{
      msg <- paste(counter, "groups have only 1 member: these groups cannot be represented in both training and validation sets.") 
    }
    warning(msg)
    
    ## exclude groups of length 1
    popmin <- min(group.n[-which(group.n==1)])        
  } # end if popmin ==1
  
  ## get training.set2 
  ## (ie. the max proportion we can use as training.set | smallest group)
  ## to be used as argument to .get.prop.pred
  training.set2 <- (popmin - 1)/popmin
  ## update training.set if needs reduction to accommodate small groups
  if(training.set2 < training.set) training.set <- training.set2
  
  ## get N.training | training.set 
  N.training <- sum(ifelse(group.n == 1L, 1L,
                          pmax(1L, pmin(group.n - 1L,
                                       round(group.n * training.set)))))
   
  
  
  ## Bound the candidate grid without fitting a full-data PCA. ##
  if(missing(n.pca.max)) n.pca.max <- min(dim(x))
  if (!length(n.pca.max) || any(!is.finite(n.pca.max)) || any(n.pca.max < 1))
    stop("n.pca.max must contain positive finite values.")
  lda.rank.max <- N.training - nlevels(grp)
  n.pca.max <- floor(min(max(n.pca.max), ncol(x), lda.rank.max))
  if (n.pca.max < 1L) stop("Too few training individuals for PCA.")
    
  ## DETERMINE N.PCA IF NEEDED ##
  if(n.pca.max < 10){
    runs <- n.pca.max
  }else{
    runs <- 10}
  
  if(is.null(n.pca)){
    n.pca <- round(pretty(1:n.pca.max, runs))
  }

  n.pca <- sort(unique(n.pca[is.finite(n.pca) & n.pca > 0 &
                            n.pca == floor(n.pca) & n.pca <= n.pca.max]))
  if (!length(n.pca)) stop("No feasible positive integer n.pca values.")
  
  
  ## GET %SUCCESSFUL OF ACCURATE PREDICTION FOR ALL VALUES ##
  res.all <- unlist(lapply(n.pca, .get.prop.pred, x, n.da, groups, grp,
                           training.set, center, scale, result,
                           n.rep, ...))
  xval <- data.frame(n.pca=rep(n.pca, each=n.rep), success=res.all)    
  valid <- tapply(is.finite(xval$success), xval$n.pca, all)
  if (!any(valid))
    stop("No candidate PC count is supported in every training split; reduce n.pca.")
  if (any(!valid))
    warning("PC counts exceeding the PCA or within-group LDA rank of a ",
            "training split were excluded ",
            "from selection: ", paste(names(valid)[!valid], collapse = ", "))
  
  
  n.pcaF <- as.factor(xval$n.pca)
  successV <- as.vector(xval$success)
  pca.success <- tapply(successV, n.pcaF, mean)
  pca.success[!valid] <- NA_real_
  n.opt <- which.max(tapply(successV, n.pcaF, mean))
  
  
  ###### MSE-BASED OPTIMAL n.pca SELECTION:
  RMSE <- sqrt(tapply((successV - 1)^2, n.pcaF, mean))
  RMSE[!valid] <- NA_real_
  ## if more than one n.pc give highest success, choose the largest one
  best.n.pca <- names(which(RMSE == min(RMSE, na.rm = TRUE)))
  if(length(best.n.pca) > 1) best.n.pca <- best.n.pca[length(best.n.pca)]
  
  # DAPC
  n.pca <- as.integer(best.n.pca)
  dapc1 <- suppressWarnings(dapc(x, grp, n.pca=n.pca, n.da=n.da,
                                center=center, scale=scale))
  
  # PLOT CROSS-VALIDATION RESULTS
  snps <- x
  phen <- grp
  random <- replicate(300, mean(tapply(sample(phen)==phen, phen, mean)))
  q.phen <- quantile(random, c(0.025,0.5,0.975))
  
  if(xval.plot==TRUE){
    plot_ok <- is.finite(successV)
    if (length(unique(xval$n.pca[plot_ok])) > 1L &&
        length(unique(successV[plot_ok])) > 1L) {
      smoothScatter(xval$n.pca[plot_ok], successV[plot_ok], nrpoints=Inf, pch=20, col=transp("black"),
                  ylim=c(0,1), xlab="Number of PCA axes retained",
                  ylab="Proportion of successful outcome prediction", 
                  main="DAPC Cross-Validation")
    } else {
      plot(xval$n.pca[plot_ok], successV[plot_ok], pch=20,
           ylim=c(0,1), xlab="Number of PCA axes retained",
           ylab="Proportion of successful outcome prediction",
           main="DAPC Cross-Validation")
    }
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


#' @method xvalDapc data.frame
#' @export
xvalDapc.data.frame <- xvalDapc.default
#' @method xvalDapc matrix
#' @export
xvalDapc.matrix <- xvalDapc.data.frame
#' @method xvalDapc genlight
#' @export
xvalDapc.genlight <- function(x, ...){
  xvalDapc.matrix(as.matrix(x), ...)
}
#' @method xvalDapc genind
#' @export
xvalDapc.genind <- function(x, ...){
  xvalDapc.matrix(tab(x), ...)
}
