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
source("data_gen.r")
####different imputation method: three method, default, tree and random forest
for (impute_method in c("default","cart","rf")){
  ####different imputation method: with or without Y
  for (impute_Y in c("yes","no")){
    ####different estimation method: 
    for (estimate_method in c("Lasso","Lasso_Order","RF","GMUS","GDS","CocoLasso")){
      ####different knockoff method: 
      for (knockoff_method in c("Second")){
        ###different dimensional of p:
        for (p in c(60)){
          ###different sample size of n:
          for (n in c(1000)){
            ###different effect size
            for (effect in c(1)){
              ####different missing proportion
              for (proMis in c(0.15)){
                ####different betak setting: 1:all, 2:previous ones for monotone missing
                for (missing in c("MAR")){
                  ###different scale of measurement error:
                  for (scale in c(0.1,0.6)){
                    ###different type of measurement
                    for (type in c("W")){
                      ###imputation time for multiple imputation
                      for (m in c(5)){

batch=1
rep=200



FDP=matrix(NA,nrow=rep,ncol=4)
power=matrix(NA,nrow=rep,ncol=4)
Z1=Z2=matrix(NA,nrow=m*rep,ncol=p)
Ztilde1=Ztilde2=matrix(NA,nrow=m*rep,ncol=p)
Xmean=matrix(0,nrow=p,ncol=rep)
Xcov=matrix(NA,nrow=p,ncol=rep*p)
select1=select2=list()

for (i in 1:rep){
  try({
    set.seed(1111+batch*rep+i)
    datal<-datagen(n=n,p=p,rho=0.5,effect=effect,interc=-1) ##data without measurement error
    datal1 <- datal$datal1
    datal2 <- datal$datal2
  
    beta1 <- datal1$beta1
    data1 <- datal1$data1
    
    beta2 <- datal2$beta2
    data2 <- datal2$data2
    
    means1<-colMeans(as.matrix(data1[,-1])) ###calculate the mean
    cov1 <- cov(as.matrix(data1[,-1]))  ###calculate the variance
    
    means2<-colMeans(as.matrix(data2[,-2])) ###calculate the mean
    cov2 <- cov(as.matrix(data2[,-2]))  ###calculate the variance
    
    Tbeta=rep(0,p)
    Tbeta[which(beta1!=0)[which(beta1!=0)==which(beta2!=0)]]=1
    S<-sum(Tbeta!=0) ### number of nonzero element
    
    sigmae<-0.3^abs(outer(1:p,1:p,"-"))
    
    ### generating missing and measurement error
    datam1 <-datamissW(datal=datal1, missing=missing, inputseed=1111+batch*50+i,proMis=proMis,sigmae=sigmae,scale=scale,type=type)
    datam2 <-datamissW(datal=datal2, missing=missing, inputseed=1111+batch*50+i,proMis=proMis,sigmae=sigmae,scale=scale,type=type)
    
    ####perform imputation
    if (impute_method %in% c("Single", "Min")){
      mydata1 <- impdata(data=datam1,method=impute_method,m=m)
      mydata2 <- impdata(data=datam2,method=impute_method,m=m)
    }
    
    if (impute_Y=="yes"){
      if (impute_method=="default"){
        mydata1 <- impdata(data=datam1,method='Multi',m=m)	
        mydata2 <- impdata(data=datam2,method='Multi',m=m)	
      }
      else{
        mydata1 <- impdata(data=datam1,method=impute_method,m=m)
        mydata2 <- impdata(data=datam2,method=impute_method,m=m)
      }
    }
    if (impute_Y=="no"){
      if (impute_method=="default"){
        mydata3 <- impdata(data=datam1[,-1],method='Multi',m=m)	
        mydata4 <- impdata(data=datam2[,-1],method='Multi',m=m)	
      }
      else{
        mydata3 <- impdata(data=datam1[,-1],method=impute_method,m=m)
        mydata4 <- impdata(data=datam2[,-1],method=impute_method,m=m)
      }
      mydata3$Y <- rep(datam1$Y,5)
      mydata4$Y <- rep(datam2$Y,5)
      #order the columns
      col_order1 <- c(".imp", "Y", names(datam1)[-1])
      mydata1 <- mydata3[,col_order1]
      col_order2 <- c(".imp", "Y", names(datam2)[-1])
      mydata2 <- mydata4[,col_order2]
    }
    
    if (impute_method %in% c("Single", "Min")){
      if (estimate_method %in% c("Lasso","Lasso1","Lasso_Order","RF","GMUS","GDS")) {
        select1[[i]] <- myest(data=mydata1,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
        select2[[i]] <- myest(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso")) {
        select1[[i]] <- myest(data=mydata1,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
        select2[[i]] <- myest(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
    }
    
    if (impute_method %in% c("default","cart","rf")){
      
      if (estimate_method %in% c("Lasso","Lasso1","Lasso_Order","RF","GMUS","GDS")) {
        select1[[i]] <- myest_mi(data=mydata1,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
        select2[[i]] <- myest_mi(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      
      if (estimate_method %in% c("CocoLasso")) {
        select1[[i]] <- myest_mi(data=mydata1,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
        select2[[i]] <- myest_mi(data=mydata2,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Diff",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
      }
      

      
      
      Z1[(m*i-m+1):(m*i),] <- select1[[i]]$Zmat
      Ztilde1[(m*i-m+1):(m*i),] <- select1[[i]]$Ztildemat
      Z2[(m*i-m+1):(m*i),] <- select2[[i]]$Zmat
      Ztilde2[(m*i-m+1):(m*i),] <- select2[[i]]$Ztildemat
    }
    
  })
}



#save Z and Ztilde
File1 <- sprintf("result/Z1_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
File2 <- sprintf("result/Ztilde1_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
write.csv(Z1,File1)
write.csv(Ztilde1,File2)

File3 <- sprintf("result/Z2_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
File4 <- sprintf("result/Ztilde2_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
write.csv(Z2,File3)
write.csv(Ztilde2,File4)

                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}



