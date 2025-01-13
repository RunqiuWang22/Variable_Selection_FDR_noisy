####read in imputed data
library(MASS)
library(glmnet)
library(Matrix)
library(hdme)
library(knockoff)
library(mice)
library(randomForest)
setwd("./")
source("estimation.r")


###different scale of measurement error:
for (scale in c(0.1,0.6)){
  ###different type of measurement
  for (type in c("W","X")){
    
impute_method="Min"
impute_Y="NA" 
####different estimation method: 
estimate_method="Lasso"
knockoff_method="Second" 
p=60
n=1000
effect=1
proMis=0.15
missing="MAR"
m=5
batch=1
rep=200


FDP=matrix(NA,nrow=rep,ncol=4)
power=matrix(NA,nrow=rep,ncol=4)
Z=matrix(NA,nrow=m*rep,ncol=p)
Ztilde=matrix(NA,nrow=m*rep,ncol=p)
Xmean=matrix(0,nrow=p,ncol=rep)
Xcov=matrix(NA,nrow=p,ncol=rep*p)
select=list()


for (i in 1:rep){
  try({
    set.seed(1111+batch*rep+i)
    datal<-datagen(n=n,p=p,rho=0.5,effect=effect,interc=-1) ##data without measurement error
    beta <- datal$beta
    data <- datal$data
    means<-colMeans(as.matrix(data[,-1])) ###calculate the mean
    cov <- cov(as.matrix(data[,-1]))  ###calculate the variance
    Tbeta <- ifelse(beta!=0,1,0) ###true beta
    S<-sum(beta!=0) ### number of nonzero element
    
    sigmae<-0.3^abs(outer(1:p,1:p,"-"))
    
    ### generating missing and measurement error
    datam <-datamissW(datal=datal, missing=missing, inputseed=1111+batch*50+i,proMis=proMis,sigmae=sigmae,scale=scale,type=type)
    
    ####perform imputation
    if (impute_method %in% c("Single", "Min")){
      mydata2 <- impdata(data=datam,method=impute_method,m=m)
    }
  
    if (impute_method %in% c("Single", "Min")){
      if (estimate_method %in% c("Lasso","Lasso1","Lasso_Order","RF","GMUL","GDS")) {
        select[[i]] <- myest(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso","CocoLasso_Order")) {
        select[[i]] <- myest(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso2","CocoLasso_Order2")) {
        select[[i]] <- myest(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=as.matrix(bdiag(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
    }
    
    if (impute_method %in% c("default","cart","rf")){
      
      if (estimate_method %in% c("Lasso","Lasso1","Lasso_Order","RF","GMUL","GDS")) {
        select[[i]] <- myest_mi(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso","CocoLasso_Order")) {
        select[[i]] <- myest_mi(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso2","CocoLasso_Order2")) {
        select[[i]] <- myest_mi(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=as.matrix(bdiag(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
    
    
    Z[(m*i-m+1):(m*i),] <- select[[i]]$Zmat
    Ztilde[(m*i-m+1):(m*i),] <- select[[i]]$Ztildemat
    Xmean[,i] <- apply(select[[i]]$Xknock,2,mean)
    Xcov[,(((p*(i-1)+1)):(p*i))] <- cov(select[[i]]$Xknock)
    }
    
    for (k in 1:4) { 
      if (length(select[[i]]$myselect[[k]])==1){
        if (is.na(select[[i]]$myselect[[k]])==T){ ### no index return
          FDP[i,k] <- 0
          power[i,k] <- 0}
        if (is.na(select[[i]]$myselect[[k]])==F) {
          FDP[i,k] <- length(setdiff(select[[i]]$myselect[[k]],which(Tbeta>0)))/max(length(select[[i]]$myselect[[k]]),1)
          power[i,k] <- length(intersect(select[[i]]$myselect[[k]],which(Tbeta>0)))/S}
      }
      if (length(select[[i]]$myselect[[k]])>1) {
        FDP[i,k] <- length(setdiff(select[[i]]$myselect[[k]],which(Tbeta>0)))/max(length(select[[i]]$myselect[[k]]),1)
        power[i,k] <- length(intersect(select[[i]]$myselect[[k]],which(Tbeta>0)))/S
      }
    }
  })
}

#save the results (need modify filename to incorporate all parameters above)
File1 <- sprintf("result/FDP_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
File2 <- sprintf("result/power_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
write.csv(FDP,  File1 )
write.csv(power, File2)

  }
}
