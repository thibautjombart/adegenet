
#' Cross-validation for Discriminant Analysis of Principal Components (DAPC)
#'
#' The function \code{xvalDapc} performs stratified cross-validation of DAPC
#' using varying numbers of PCs (and keeping the number of discriminant
#' functions fixed); \code{xvalDapc} is a generic with methods for
#' \code{data.frame} and \code{matrix}.\cr
#'
#' The Discriminant Analysis of Principal Components (DAPC) relies on dimension
#' reduction of the data using PCA followed by a linear discriminant analysis.
#' How many PCA axes to retain is often a non-trivial question. Cross
#' validation provides an objective way to decide how many axes to retain:
#' different numbers are tried and the quality of the corresponding DAPC is
#' assessed by cross-validation: DAPC is performed on a training set, typically
#' made of 90\% of the observations (comprising 90\% of the observations in
#' each subpopulation) , and then used to predict the groups of the 10\% of
#' remaining observations.  The current method uses the average prediction
#' success per group (result="groupMean"), or the overall prediction success
#' (result="overall"). The number of PCs associated with the lowest Mean
#' Squared Error is then retained in the DAPC.
#'
#' @aliases xvalDapc xvalDapc.data.frame xvalDapc.matrix
#' @param x \code{a data.frame} or a \code{matrix} used as input of DAPC.
#' @param grp a \code{factor} indicating the group membership of individuals.
#' @param n.pca.max maximum number of PCA components to retain.
#' @param n.da an \code{integer} indicating the number of axes retained in the
#' Discriminant Analysis step. If \code{NULL}, n.da defaults to 1 less than the
#' number of groups.
#' @param training.set the proportion of data (individuals) to be used for the
#' training set; defaults to 0.9 if all groups have >= 10 members; otherwise,
#' training.set scales automatically to the largest proportion that still
#' ensures all groups will be present in both training and validation sets.
#' @param result a character string; "groupMean" for group-wise assignment
#' sucess, or "overall" for an overall mean assignment success; see details.
#' @param center a \code{logical} indicating whether variables should be
#' centred to mean 0 (TRUE, default) or not (FALSE). Always TRUE for
#' \linkS4class{genind} objects.
#' @param scale a \code{logical} indicating whether variables should be scaled
#' (TRUE) or not (FALSE, default). Scaling consists in dividing variables by
#' their (estimated) standard deviation to account for trivial differences in
#' variances.
#' @param n.pca an \code{integer} vector indicating the number of different
#' number of PCA axes to be retained for the cross validation; if \code{NULL},
#' this will be dertermined automatically.
#' @param n.rep the number of replicates to be carried out at each level of PC
#' retention; defaults to 30.
#' @param xval.plot a logical indicating whether a plot of the cross-validation
#' results should be generated.
#' @param \dots further arguments to be passed to other methods.
#' @return A \code{list} containing seven items, and a \code{plot} of the
#' results.  The first is a \code{data.frame} with two columns, the first
#' giving the number of PCs of PCA retained in the corresponding DAPC, and the
#' second giving the proportion of successful group assignment for each
#' replicate.  The second item gives the mean and confidence interval for
#' random chance.  The third gives the mean successful assignment at each level
#' of PC retention.  The fourth indicates which number of PCs is associated
#' with the highest mean success.  The fifth gives the Root Mean Squared Error
#' at each level of PC retention.  The sixth indicates which number of PCs is
#' associated with the lowest MSE.  The seventh item contains the DAPC carried
#' out with the optimal number of PCs, determined with reference to MSE.
#'
#' If \code{xval.plot=TRUE} a scatterplot of the results of cross-validation
#' will be displayed.
#' @author Caitlin Collins \email{caitlin.collins12@@imperial.ac.uk}, Thibaut
#' Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link{dapc}}
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant
#' analysis of principal components: a new method for the analysis of
#' genetically structured populations. BMC Genetics11:94.
#' doi:10.1186/1471-2156-11-94
#' @keywords multivariate
#' @examples
#'
#' \dontrun{
#' ## CROSS-VALIDATION ##
#' data(sim2pop)
#' xval <- xvalDapc(sim2pop@@tab, pop(sim2pop), n.pca.max=100, n.rep=3)
#' xval
#'
#' }
#'

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
    result <- result}


  ## GET TRAINING SET SIZE ##
  N <- nrow(x)
  groups <- levels(grp)
  if(all(lapply(groups, function(e) sum(as.vector(unclass(grp==e))))>=10)==TRUE){
    N.training <- round(N*training.set)
  }else{
    groups1 <- (levels(grp))[(as.vector(which.min((lapply(groups, function(e) sum(as.vector(unclass(grp==e))))))))]
    popmin <- length(which(grp%in%groups1))
    training.set2 <- ((popmin - 1)/popmin)
    N.training <- round(N*training.set2)}


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
  } # end get.prop.pred


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

