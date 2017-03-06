#' Monte Carlo test prior to spca
#'
#'
#' @export
#'
#' @rdname spca.randtest
#'
#' @author original code by Valeria Montano with tweaks by Thibaut Jombart
#'
spca.randtest <- function(...) {
  UseMethod("spca.randtest")
}







spca.randtest.matrix <-function(obj, obs, nperm=100, p=0.05){

  if(!is.genind(obj) & !is.genpop(obj)){
    stop("obj must be a genind or genpop object")
  }

  if(!inherits(obs, "spca")){
    stop("obs must be an spca object")
  }


  rand <-function(obj){

    obj@tab<-obj@tab[sample(1:nrow(obj@tab)),]
    X <- scaleGen(obj, center = TRUE, scale = FALSE, truenames = TRUE)
    pcaX <- dudi.pca(X, center = FALSE, scale = FALSE, scannf = FALSE)
    sp_sim <- multispati(dudi = pcaX, listw = obs$lw, scannf = FALSE, nfposi = 1, nfnega = 1)
    return(sp_sim$eig)

  }

  sims<-sapply(1:nperm, function(e) rand(obj))

  if (is.list(sims)) sims<-do.call(cbind,sims)
  psims<-apply(sims,2,function(e) sum(e[which(e>0)]))
  nsims<-apply(sims,2,function(e) sum(e[which(e<0)]))

  test<-function(rval, eig){

    if (eig>0){pval<-as.randtest(rval, eig, alter="greater")
    }else{
      pval<-as.randtest(rval, eig, alter="less")}
    return(pval$pvalue)

  }

  psims<-as.randtest(psims, sum(obs$eig[which(obs$eig>0)]), alter="greater")$pvalue
  nsims<-as.randtest(nsims, sum(obs$eig[which(obs$eig<0)]), alter="less")$pvalue

  pvalues<-list("p-value_global_pattern"=c(), "p_value_local_pattern"=c())
  pvalues[[1]]<-psims
  pvalues[[2]]<-nsims

  if(psims <= p || nsims <= p){

    peig<-matrix(0, nrow=length(obs$eig), ncol=2)
    peig[,1]<-obs$eig
    peig[,2]<-sapply(obs$eig, function(e) test(sims[which(obs$eig==e),], e))
    pvalues[[3]]<-peig

  }else{pvalues[[3]]<-NULL}

  names(pvalues[1:2])<-
    if(length(pvalues)>2) colnames(pvalues[[3]])<-c("eigenvalue", "p-value")
  return(pvalues)
}

