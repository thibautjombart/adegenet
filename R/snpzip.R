############
## snpzip ##
############

snpzip <- function(snps,phen,plot=TRUE,pca.plot=FALSE,method=c("complete","single","average","centroid","mcquitty","median","ward")) {

## DETERMINE N.DA & N.PCA ##

n.da <- nlevels(phen)-1
xval <- xvalDapc(snps,phen,70,n.da=n.da,training.set=0.8,result=c("groupMean"),n.pca=NULL,n.rep=30)
  ## xvalDapc will be changed to sxvalDapc & training.set to be removed from argument list)
n.pcaF <- as.factor(xval$n.pca)
successV <- as.vector(xval$success)
pca.success <- tapply(successV,n.pcaF,mean)
n.opt <- which.max(tapply(successV,n.pcaF,mean))
n.pca <- n.opt

## SHOW CROSS VALIDATION RESULTS ##

if(pca.plot==TRUE){
print(pca.success)
print(names(n.opt))
random <- replicate(300,mean(tapply(sample(phen)==phen,phen,mean)))
q.phen <- quantile(random,c(0.025,0.5,0.975))
print(q.phen)
smoothScatter(xval$n.pca,successV,nrpoint=Inf,pch=20,col=transp("black"),xlab="Number of PCA axes retained",ylab="Proportion of successful outcome prediction")
abline(h=q.phen,lty=c(2,1,2))}

## GET DAPC ##

dapc1 <- dapc(snps,phen,n.pca=n.pca,n.da=n.da)

## PLOT DAPC ##

if(plot==TRUE){
scatter(dapc1,bg="white",scree.da=FALSE,scree.pca=TRUE,posi.pca="topright",
col=c("royalblue","red"),legend=TRUE,posi.leg="topleft")
title("DAPC")}

## GET CLUSTER OF STRUCTURAL SNPS ##

x <- dapc1$var.contr
D <- dist(x)
clust <- hclust(D,method)
pop <- factor(cutree(clust,k=2,h=NULL))
m <- which.max(tapply(x,pop,mean))

## RETURN SNPZIP RESULTS ##

answer <- list(names(n.opt),table(pop),which(pop==m),dapc1$var.contr[pop==m])
names(answer)[[1]] <- "number of PCs of PCA retained"
names(answer)[[2]] <- "number of alleles in groups defined by contribution to DAPC"
names(answer)[[3]] <- "list of structuring alleles"
names(answer)[[4]] <- "contributions of structuring alleles to DAPC"
return(answer)
} # end snpzip









