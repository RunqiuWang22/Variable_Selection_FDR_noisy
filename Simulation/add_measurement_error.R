args=commandArgs(trailingOnly=TRUE)
iii=as.numeric(args[1])

####different estimation method:  
estimate_method_list=c("Lasso","Lasso1","Lasso_Order","RF","GMUL","GDS","CocoLasso", "CocoLasso_Order","CocoLasso2", "CocoLasso_Order2")#10
estimate_method=estimate_method_list[iii%%length(estimate_method_list)+1]
iii=iii%/%length(estimate_method_list)

####different knockoff method: x
#knockoff_method_list=c("Gaussian","Second") #2
knockoff_method_list=c("Second")
knockoff_method=knockoff_method_list[iii%%length(knockoff_method_list)+1]
iii=iii%/%length(knockoff_method_list)

###different dimensional of p:
#p_list=c(60, 120, 180, 240, 300) #5
p_list=c(60,120)
p=p_list[iii%%length(p_list)+1]
iii=iii%/%length(p_list) 

###different sample size of n:
#n_list=c(400, 1500,1000) #3
n_list=c(1000)
n=n_list[iii%%length(n_list)+1]
iii=iii%/%length(n_list) 

###different effect size
effect_list=c(0.5, 1,1.5) #3
effect=effect_list[iii%%length(effect_list)+1]
iii=iii%/%length(effect_list) 

###different scale of measurement error:
scale_list=c(0.01, 0.09, 0.25) #5
scale=scale_list[iii%%length(scale_list)+1]
iii=iii%/%length(scale_list) 


batch=iii+1
rep=50 ###200 in total: remember to update!!!!!! change to 20

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


    if (estimate_method %in% c("Lasso","Lasso1","Lasso_Order","RF","GMUL","GDS")) {
      select[[i]] <- myest(data=data,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Max",family="binomial",Sigma=NULL,delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
    }
    
    if (estimate_method %in% c("CocoLasso","CocoLasso_Order")) {
      select[[i]] <- myest(data=data,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Max",family="binomial",Sigma=cbind(rbind(scale*sigmae,scale*sigmae),rbind(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
    }
    
    if (estimate_method %in% c("CocoLasso2","CocoLasso_Order2")) {
      select[[i]] <- myest(data=data,q=c(0.1,0.2), method=knockoff_method,stat=estimate_method,flip="Max",family="binomial",Sigma=as.matrix(bdiag(scale*sigmae,scale*sigmae)),delta=NULL, offset=c(0,1),MEAN=means,COV=cov)
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

#save Z and Ztilde
File3 <- sprintf("result/Z_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
File4 <- sprintf("result/Ztilde_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
write.csv(Z,File3)
write.csv(Ztilde,File4)

#save xtidle mean and variance
File5 <- sprintf("result/Xmean_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
File6 <- sprintf("result/Xcov_%d_%d_%s_%s_%f_%f_%d.csv",p,n,estimate_method,knockoff_method,effect,scale,batch)
write.csv(Xmean,File5)
write.csv(Xcov,File6)







