library(MASS)
library(glmnet)
library(Matrix)
library(hdme)
library(MASS)
library(glmnet)
library(knockoff)
library(randomForest)
library(Rcpp)
setwd("./")
source("estimation.r")
source("data_gen.r")
####different estimation method:  
for (estimate_method in c("Lasso","Lasso_Order","RF","GMUS","GDS","CocoLasso")){
####different knockoff method: x
 for (knockoff_method in c("Second")){
  ###different dimensional of p:
  for (p in c(60,120)){
   ###different sample size of n:
   for (n in c(1000)){
    ###different effect size
    for (effect in c(0.5,1.5)){
      ###different scale of measurement error:
      for (scale in c(0.6,1)){
batch=1
rep=200



FDP=matrix(NA,nrow=rep,ncol=4)
power=matrix(NA,nrow=rep,ncol=4)
Z=matrix(NA,nrow=p,ncol=rep)
Ztilde=matrix(NA,nrow=p,ncol=rep)
Xmean=matrix(NA,nrow=p,ncol=rep)
Xcov=matrix(NA,nrow=p,ncol=rep*p)
select=list()


for (i in 1:rep){
  try({
set.seed(1111+batch*rep+i)
ori_data<-datagen(n=n,p=p,rho=0.5,effect=effect,interc=-1) ##data without measurement error
means=colMeans(ori_data$data[,-1]) ###calculate the mean
cov = cov(ori_data$data[,-1])  ###calculate the variance
Tbeta = ifelse(ori_data$beta!=0,1,0) ###true beta set
S=sum(ori_data$beta!=0) ### number of nonzero element

sigmae<-0.3^abs(outer(1:p,1:p,"-"))
data<-datamea(data=ori_data$data,sigmae=sigmae,scale=scale) ##data with measurement error


    if (estimate_method %in% c("Lasso","Lasso_Order","RF","GMUS","GDS")) {
      select[[i]] <- myest(data=data,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Max",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
    }
    
    if (estimate_method %in% c("CocoLasso")) {
      select[[i]] <- myest(data=data,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Max",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
    }
    
    
 Z[,i]= select[[i]]$Z
 Ztilde[,i]=select[[i]]$Ztilde
 Xmean[,i] <- apply(select[[i]]$Xknock,2,mean)
 Xcov[,(((p*(i-1)+1)):(p*i))] <- cov(select[[i]]$Xknock)
    
    for (k in 1:4) { 
      FDP[i,k] <- length(setdiff(select[[i]]$myselect[[k]],which(Tbeta>0)))/max(length(select[[i]]$myselect[[k]]),1)
      power[i,k] <- length(intersect(select[[i]]$myselect[[k]],which(Tbeta>0)))/S
    }
  })
  
}

#save the results (need modify filename to incorporate all parameters above)
File1 <- sprintf("result/FDP_breast_match_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
File2 <- sprintf("result/power_breast_match_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
write.csv(FDP,  File1 )
write.csv(power, File2)


      }
    }
   }
  }
 }
}






