args=commandArgs(trailingOnly=TRUE)
iii=as.numeric(args[1])

####different imputation method: three method, default, tree and random forest
#impute_method_list=c("Single","Min", "default","cart","rf") #5
impute_method_list=c("default","cart","rf") #3
impute_method=impute_method_list[iii%%length(impute_method_list)+1]
iii=iii%/%length(impute_method_list)

####different imputation method: with or without Y
impute_Y_list=c("yes","no") #2
impute_Y=impute_Y_list[iii%%length(impute_Y_list)+1]
iii=iii%/%length(impute_Y_list)

####different estimation method: 
estimate_method_list=c("Lasso","Lasso1","Lasso_Order","RF","GMUL","GDS","CocoLasso", "CocoLasso_Order","CocoLasso2", "CocoLasso_Order2")#10
estimate_method=estimate_method_list[iii%%length(estimate_method_list)+1]
iii=iii%/%length(estimate_method_list)

####different knockoff method: 
knockoff_method_list=c("Second") #1
knockoff_method=knockoff_method_list[iii%%length(knockoff_method_list)+1]
iii=iii%/%length(knockoff_method_list)

###different dimensional of p:
p_list=c(60, 120) #2
p=p_list[iii%%length(p_list)+1]
iii=iii%/%length(p_list) 

###different sample size of n:
n_list=c(1000) #1
n=n_list[iii%%length(n_list)+1]
iii=iii%/%length(n_list) 

###different effect size
#effect_list=c(0.5,1) #1
effect_list=c(1)
effect=effect_list[iii%%length(effect_list)+1]
iii=iii%/%length(effect_list) 

####different missing proportion
#proMis_list=c(0.05,0.15) #2
proMis_list=c(0.15)
proMis=proMis_list[iii%%length(proMis_list)+1]
iii=iii%/%length(proMis_list)

####different betak setting: 1:all, 2:previous ones for monotone missing
missing_list=c("MAR") #1
missing=missing_list[iii%%length(missing_list)+1]
iii=iii%/%length(missing_list)

###different scale of measurement error:
scale_list=c(0.1,0.6) #2
scale=scale_list[iii%%length(scale_list)+1]
iii=iii%/%length(scale_list) 

###different type of measurement
type_list=c("X","W") #2
type=type_list[iii%%length(type_list)+1]
iii=iii%/%length(type_list) 

m_list=c(5)
m=m_list[iii%%length(m_list)+1]
iii=iii%/%length(m_list)

batch=iii+1
rep=50 #200 in total: remember to update!!!!!!

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
    
    if (impute_Y=="yes"){
      if (impute_method=="default"){
        mydata2 <- impdata(data=datam,method='Multi',m=m)	
      }
      else{
        mydata2 <- impdata(data=datam,method=impute_method,m=m)
      }
    }
    if (impute_Y=="no"){
      if (impute_method=="default"){
        mydata3 <- impdata(data=datam[,-1],method='Multi',m=m)	
      }
      else{
        mydata3 <- impdata(data=datam[,-1],method=impute_method,m=m)
      }
      mydata3$Y <- rep(datam$Y,5)
      #order the columns
      col_order <- c(".imp", "Y", names(datam)[-1])
      mydata2 <- mydata3[,col_order]
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

#save Z and Ztilde
File3 <- sprintf("result/Z_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
File4 <- sprintf("result/Ztilde_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
write.csv(Z,File3)
write.csv(Ztilde,File4)

#save xtidle mean and variance
File5 <- sprintf("result/Xmean_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
File6 <- sprintf("result/Xcov_%s_%s_%s_%s_%d_%d_%f_%f_%s_%f_%s_%d_%d.csv",impute_method,impute_Y,estimate_method,knockoff_method, p, n, effect, proMis,missing,scale,type,m,batch)
write.csv(Xmean,File5)
write.csv(Xcov,File6)



