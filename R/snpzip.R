
##############
### snpzip ###
##############


snpzip <- function(snps, y, plot=TRUE, xval.plot=FALSE, loading.plot=FALSE,
                   method=c("complete","single","average","centroid",
                            "mcquitty","median","ward"), ...) {
  
  ## dapc input prompts only SNP selection function
  if(class(y)=="dapc"){ 
    
    dapc1 <- y
    phen <- 0
    
    
    if(xval.plot==TRUE){
      warning("cross-validation not performed when x is a dapc object; xval.plot will not be shown")
      xval.plot=FALSE
    }
  }
  
  ## snps, phen input prompts cross-validation, DAPC, and SNP selection functions
  else{
    
    if(missing(y)){
      stop("phen argument needed")
    }
    phen <- y
    
  }
  
  
  if(class(y)!="dapc"){
    
    ################################################
    ######## Stratified Cross-Validation ###########
    ################################################
    
    xvalDapc <- function(x, grp, n.pca.max = 200, n.da = NULL, 
                         training.set = 0.9, result = "groupMean", 
                         center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 30, ...){
      
      ## CHECKS ##
      grp <- factor(grp)
      n.pca <- n.pca[n.pca>0]
      n.da <- length(levels(grp))-1
      if(missing(training.set)){
        training.set <- 0.9}
      else{
        training.set <- training.set}
      if(missing(n.rep)){
        n.rep<-30}
      else{
        n.rep<-n.rep}
      
      
      ## GET TRAINING SET SIZE ##
      N <- nrow(x)
      groups <- levels(grp)
      if(all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
        N.training <- round(N*training.set)}
      else{
        groups1 <- (levels(grp))[(as.vector(which.min((lapply(groups, 
                                                              function(e) sum(as.vector(unclass(grp==e))))))))]
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
        runs<- 10}
      
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
      
      
      ################################################
      #####  MSE Calculation and n.pca Selection #####
      ################################################
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
      
      
      
      ################################################
      #################  DAPC  #######################
      ################################################
      n.pcaF <- as.factor(xval$n.pca)
      successV <- as.vector(xval$success)
      pca.success <- tapply(successV, n.pcaF, mean)
      n.opt <- which.max(tapply(successV, n.pcaF, mean))
      n.pca <- as.integer(best.n.pca)
      n.da <- nlevels(grp)-1
      dapc1 <- dapc(x, grp, n.pca=n.pca, n.da=n.da)
      
      answerme <- list(n.da, n.pca, dapc1, xval, successV, RMSE, best.n.pca, pca.success, n.opt)
      return(answerme)
    }
    
    
    x <- snps
    grp <- phen
    
    XVAL <- xvalDapc(x, grp, n.da=NULL, training.set=0.9, result="groupMean",
                     center=TRUE, scale=FALSE, n.pca=NULL, ...)
    n.da <- XVAL[[1]]
    n.pca <- XVAL[[2]]
    dapc1 <- XVAL[[3]]
    xval <- XVAL[[4]]
    successV <- XVAL[[5]]
    RMSE <- XVAL[[6]]
    best.n.pca <- XVAL[[7]]
    pca.success <- XVAL[[8]]
    n.opt <- XVAL[[9]]
    
    
    ################################################
    ######## Show Cross-Validation Results #########
    ################################################
    snps <- x
    phen <- grp
    
    if(xval.plot==TRUE){
      par(ask=TRUE)
      random <- replicate(300, mean(tapply(sample(phen)==phen, phen, mean)))
      q.phen <- quantile(random, c(0.025,0.5,0.975))
      smoothScatter(xval$n.pca, successV, nrpoints=Inf, pch=20, col=transp("black"),
                    ylim=c(0,1), xlab="Number of PCA axes retained",
                    ylab="Proportion of successful outcome prediction", 
                    main="DAPC Cross-Validation")
      abline(h=q.phen,lty=c(2,1,2))
      
      xvalResults <- list(xval, q.phen, pca.success, (names(n.opt)), RMSE, best.n.pca)
      names(xvalResults)[[1]] <- "Cross-Validation Results"
      names(xvalResults)[[2]] <- "Median and Confidence Interval for Random Chance"
      names(xvalResults)[[3]] <- "Mean Successful Assignment by Number of PCs of PCA"
      names(xvalResults)[[4]] <- "Number of PCs Achieving Highest Mean Success"
      names(xvalResults)[[5]] <- "Root Mean Squared Error by Number of PCs of PCA"
      names(xvalResults)[[6]] <- "Number of PCs Achieving Lowest MSE"
      print(xvalResults)
    }
    
  } # end of snps, phen section
  
  
  ################################################
  #############  Plot DAPC Results   #############
  ################################################
  
  if(plot==TRUE){
    myCol <- colorRampPalette(c("blue", "gold", "red"))
    scatter(dapc1, bg="white", scree.da=FALSE, scree.pca=TRUE, posi.pca="topright",
            col=myCol((dapc1$n.da)+1), legend=TRUE, posi.leg="topleft")
    title("DAPC")}
  
  ################################################
  ###### Select Cluster of Structural SNPS #######
  ################################################
  
  if(missing(method)){
    method <- "ward"
  }
  else{ 
    method <- method}
  
  selector <- function(dapc1, dimension){
    
    z <- dapc1$var.contr[,dimension]
    xTotal <- dapc1$var.contr[,dimension]
    toto <- which(xTotal%in%tail(sort(xTotal), 2000))
    z <- sapply(toto, function(e) xTotal[e])
    
    D <- dist(z)
    clust <- hclust(D,method)
    pop <- factor(cutree(clust,k=2,h=NULL))
    m <- which.max(tapply(z,pop,mean))
    maximus <- which(pop==m)
    maximus <- as.vector(unlist(sapply(maximus, function(e) toto[e])))
    popvect <- as.vector(unclass(pop))
    n.snp.selected <- sum(popvect==m)
    sel.snps <- snps[,maximus]
    
    selection <- c((ncol(snps)-ncol(snps[,-maximus])), ncol(snps[,-maximus]))
    
    resultat <- list(selection, maximus, dimnames(sel.snps)[[2]], dapc1$var.contr[maximus, dimension])
    
    names(resultat)[[1]] <- "Number of selected vs. unselected alleles"
    names(resultat)[[2]] <- "List of selected alleles"
    names(resultat)[[3]] <- "Names of selected alleles"
    names(resultat)[[4]] <- "Contributions of selected alleles to discriminant axis"
    gc()
    return(resultat)
  }

    if(dapc1$n.da==1){
      features <- selector(dapc1, dimension=1)      
    }
    else{
      dimensions <- c(1:dapc1$n.da)
      features <- lapply(dimensions, function(e)
                        selector(dapc1, dimension=e))     
    }
    
    ################################################
    ####  Calculate success in discrimination   ####
    ################################################
    
    # overall
    grp <- dapc1$grp
    ass <- dapc1$assign
    dapc.success.overall <- length(which(ass==grp)) / length(grp)
    
    # by group        
    
    # give grp and ass numbered factor levels
    GRP <- factor(grp, levels=levels(grp), labels=c(0:(nlevels(grp)-1)))
    ASS <- factor(ass, levels=levels(ass), labels=c(0:(nlevels(ass)-1)))
    
    # make each level of those factors into one element of a list
    index <- c(0:(nlevels(grp) - 1))
    ASSIGN <- sapply(index, function(e) which(ASS==e))
    GROUP <- sapply(index, function(e) which(GRP==e))
    
    index2 <- c(1:(nlevels(grp)))
    dapc.success.byGroup <- sum(sapply(index2, function(e) 
      (length(which(ASSIGN[[e]]%in%GROUP[[e]]))) / 
                                         length(GROUP[[e]]))) / length(index2)
    
    dapc.success <- c(dapc.success.overall, dapc.success.byGroup)
    
    ################################################
    #### Loading Plot Delineating SNP Clusters  ####
    ################################################
    
    if(loading.plot==TRUE){
      if(dapc1$n.da==1){
      par(ask=TRUE)
      maximus <- features[[2]]
      decimus <- abs(dapc1$var.contr[maximus][(which.min(dapc1$var.contr[maximus]))])-0.000001
      meridius <- loadingplot(dapc1$var.contr[,1], threshold=c(decimus))
      }
      else{
        par(ask=TRUE)
        maximus <- features[[1]][[2]] # not sure this is the right way to deal with the lapply
                                      # plus need to specify that this is only the loadingplot for dimension 1 ... 
        decimus <- abs(dapc1$var.contr[maximus][(which.min(dapc1$var.contr[maximus]))])-0.000001
        meridius <- loadingplot(dapc1$var.contr[,1], threshold=c(decimus)) 
      }
    }
    
    ################################################
    ########## Return snpzip Results  ##############
    ################################################
    
    if(class(y)=="dapc"){
      
      answer <- list(dapc1$n.pca, features)
      names(answer)[[1]] <- "Number of PCs of PCA retained"
      names(answer)[[2]] <- "FS"
      return(answer) 
    }
    
    else{ 
      answer <- list(best.n.pca, features, dapc.success, dapc1)
      names(answer)[[1]] <- "Number of PCs of PCA retained"
      names(answer)[[2]] <- "FS"
      names(answer)[[3]] <- "Discrimination success overall & by group"
      names(answer)[[4]] <- "DAPC"
      return(answer)
    }
    
    par(ask=FALSE)
  } # end snpzip
  
  
  
  
  
  