
##############
## xval: temp xvalDapc ##
##############


xvalDapc <- function(x, grp, n.pca.max = 300, n.da = NULL, training.set = 0.9, 
                     result = c("groupMean", "overall"), center = TRUE, scale = FALSE, 
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE, ...){
  
  ## CHECKS ##
  grp <- factor(grp)
  n.pca <- n.pca[n.pca > 0]
  if(missing(n.da)){
  n.da <- length(levels(grp))-1}
  else{
    n.da <- n.da}
  if(missing(training.set)){
    training.set <- 0.9}
  else{
    training.set <- training.set}
  if(missing(n.rep)){
    n.rep <- 30}
  else{
    n.rep <- n.rep}
  
  
  ## GET TRAINING SET SIZE ##
  N <- nrow(x)
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
  n.pca.max <- min(n.pca.max, pcaX$rank, N.training-1)
  
  ## DETERMINE N.PCA IF NEEDED ##
  if(n.pca.max < 10){
    runs <- n.pca.max}
  else{
    runs <- 10}
  
  if(is.null(n.pca)){
    n.pca <- round(pretty(1:n.pca.max, runs))
  }
  n.pca <- n.pca[n.pca>0 & n.pca<(N.training-1)]
  
  ## FUNCTION GETTING THE % OF ACCURATE PREDICTION FOR ONE NUMBER OF PCA PCs ##
  ## n.pca is a number of retained PCA PCs
  get.prop.pred <- function(n.pca){
    f1 <- function(){
      if(all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
        toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), 
                                                           size=(round(training.set*sum(as.vector(unclass(grp==e))))))))}
      else{
        toKeep <- unlist(lapply(groups, function(e) sample(which(grp==e), 
                                                           size=(round(training.set2*sum(as.vector(unclass(grp==e))))))))}
      temp.pca <- pcaX
      temp.pca$li <- temp.pca$li[toKeep,,drop=FALSE]
      temp.dapc <- suppressWarnings(dapc(x[toKeep,,drop=FALSE], grp[toKeep], 
                                         n.pca=n.pca, n.da=n.da, dudi=temp.pca))
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
    print(abline(h=q.phen, lty=c(2,1,2)))
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



