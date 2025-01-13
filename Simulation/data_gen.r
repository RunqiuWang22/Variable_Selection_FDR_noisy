datagen<-function(n,p,rho,effect,interc){
  ###n sample size
  ###p number of predictors
  ##s proportion of non-zero element
  ###rho correlation parameter for design matrix X
  ###effect: effect scale multiplier
  Sigma<-rho^abs(outer(1:p,1:p,"-"))
  X<-mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
  s<-p*0.25 ### keep the proportion of sparsity at 25%
  b0=rep(c(3,1.5,0,0,2,0,0),s/3)
  beta=c(b0,rep(0,p-length(b0)))*effect
  beta<-beta*(-1)^(runif(length(beta))>0.5)
  logit <- function(x) {(1+exp(-x))^(-1)}
  Y<-rbinom(n, size = 1, prob = logit(X %*% beta+interc))
  data<-cbind(Y,X)
  data<-data.frame(data)
  names(data)=c("Y",paste("X",1:p,sep=""))
  return(list(data=data,beta=beta))
}

###generate data with measurement errors
datamea<- function(data,sigmae,scale){
  ###data without measurement error
  ###rhoe correlation parameter for design matrix e
  ###scale the scale for correlation matrix e
  Y<-data[,1]
  X<-data[,-1]
  n<-dim(X)[1]
  p<-dim(X)[2]
  W<-X+mvrnorm(n=n,mu=rep(0,p),Sigma=sigmae*scale)
  newdata<-cbind(Y,W)
  newdata
} 

###generate some missing data in X
datamiss = function(datal, missing,inputseed){
  set.seed(1234)
  ###MCAR: missing is completely random and unrelated to any observed or unobserved variables. 
  ###MAR: missing depends on observed responses, but are unrelated to the specific missing values
  beta = datal$beta
  data = datal$data
  X=as.matrix(data[,-1])
  Y=as.matrix(data[1])
  ###randomly choose missing column for true signals
  colindex_t = sample(which(beta!=0), size=S*(4/30),replace = F)
  ###randomly choose missing column for non-true signals 
  colindex_c = sample(which(beta==0), size=(p-S)*(4/30),replace = F)
  colindex = c(colindex_t, colindex_c)
  #generate R--missing indicator:1 means missing
  
  a<-list()
  R<-matrix(NA,nrow=n,ncol=p)
  logit <- function(x) (1+exp(-x))^(-1)
  prob<-matrix(NA,nrow=n,ncol=length(colindex))
  Xnew<-X
  ####create missing indicators
  if(missing%in%c("MCAR")){
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,0.2)
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("MAR")){
    set.seed(1234)
    for (k in 1:length(colindex)) {
      prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%runif(n=p, -2, 2))-5))  ### add a negative number
    }
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,prob[,k])
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("Monotone")){
    set.seed(1234)
    for (k in 1:length(colindex)) {
      prob[,k] <- as.vector(logit((cbind(1,X[,1:(colindex[k]-1)])%*%runif(n=length(1:(colindex[k]-1))+1, -2, 2))-5))  ### add a negative number
    }  
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,prob[,k])
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  mydata1 <- data.frame(Y, Xnew)
  return(mydata1)
}

###generate some missing data in X or W
datamissW = function(datal, missing,inputseed,proMis,sigmae,scale,type){
  ###datal: orginal data list without missing
  ###missing: missing type
  ###inputseed: each seed are different for different dataset
  ###proMis: missing proportion
  ###sigame: covariance matrix of sigmae for measurement error
  ###scale: scale for covariance matrix of sigmae for measurement error
  ###type: if type="X", using X matrix generating missing, if type="W", using W matri generating missing 
  set.seed(1234)
  ###MCAR: missing is completely random and unrelated to any observed or unobserved variables. 
  ###MAR: missing depends on observed responses, but are unrelated to the specific missing values
  beta = datal$beta
  S<-sum(beta!=0)
  data = datal$data ### original dataset
  dataW=datamea(data=data,sigmae=sigmae,scale=scale) ### dataset with measurement error
  if (type=="X"){
    X=as.matrix(data[,-1])
  }
  if (type=="W"){
    X=as.matrix(dataW[,-1])
  }
  Y=as.matrix(data[1])
  ###randomly choose missing column for true signals
  colindex_t = sample(which(beta!=0), size=S*(4/30),replace = F)
  ###randomly choose missing column for non-true signals 
  colindex_c = sample(which(beta==0), size=(p-S)*(4/30),replace = F)
  colindex = c(colindex_t, colindex_c)[order(c(colindex_t, colindex_c))]
  #generate R--missing indicator:1 means missing
  
  intercept=rep(0,length(colindex))
  a<-list()
  R<-matrix(NA,nrow=n,ncol=p)
  logit <- function(x) (1+exp(-x))^(-1)
  prob<-matrix(NA,nrow=n,ncol=length(colindex))
  Xnew<-X
  ####create missing indicators
  if(missing%in%c("MCAR")){
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,proMis)
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("MAR")){
    for (k in 1:length(colindex)) {
      set.seed(1234+k)
      prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%rnorm(n=p))-intercept[k]))  
    }
    aprob=apply(prob,2,mean)
    ### add a number to make sure the missing proportion is around a certain point
    for (k in 1:length(colindex)) {
      while (aprob[k]>proMis){
        set.seed(1234+k)
        intercept[k]=intercept[k]+0.0005
        prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%rnorm(n=p))-intercept[k]))  
        aprob[k]=mean(prob[,k])
      }
    }
    
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,prob[,k])
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("Monotone")){
    set.seed(1234)
    for (k in 1:length(colindex)) {
      set.seed(1234+k)
      prob[,k] <- as.vector(logit((cbind(1,X[,1:(colindex[k]-1)])%*%rnorm(n=length(1:(colindex[k]-1))+1))-intercept[k]))  ### add a negative number
    }  
    aprob=apply(prob,2,mean)
    ### add a number to make sure the missing proportion is around a certain point
    for (k in 1:length(colindex)) {
      while (aprob[k]>proMis){
        set.seed(1234+k)
        intercept[k]=intercept[k]+0.0005
        prob[,k] <- as.vector(logit((cbind(1,X[,1:(colindex[k]-1)])%*%rnorm(n=length(1:(colindex[k]-1))+1))-intercept[k])) 
        aprob[k]=mean(prob[,k])
      }
    }
    
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      a[[k]] <- rbinom(n,1,prob[,k])
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  mydata1 <- data.frame(Y, Xnew)
  
  if (type=="X"){
    mydata2=mydata1+datamea(data=mydata1,sigmae=sigmae,scale=scale) ###add measurement error
    mydata2$Y=mydata2$Y/2
  }
  if (type=="W"){
    mydata2=mydata1
  }
  
  
  
  return(mydata2)
}

####different way to handle missing data, data is in form of Y, X with no missing in Y
impdata<-function(data,method="Single",rr=log(2),m){
  ###Method: Single: single imputation by mean; Min: single imputation by min value minus rr; Multi: Multiple imputation; Ind: Create indicator for missing.
  if (method=="Single"){
    for (j in 2:ncol(data)){
      ss=which(is.na(data[,j]))
      if (length(ss)>0){
        data[ss,j]=mean(data[,j],na.rm=TRUE)
      }
    }
  }
  if (method=="Min"){
    for (j in 2:ncol(data)){
      ss=which(is.na(data[,j]))
      if (length(ss)>0){
        data[ss,j]=min(data[,j],na.rm=TRUE)-rr
      }
    }		
  }
  
  
  if (method%in%c("Multi","cart","rf")){
    if (sum(is.na(data))>0) {
      if (method=="Multi") {imp=mice(data,m=m)}
      else{imp=mice(data, m=m,method=method)}
      data=complete(imp,action="long")
      data=data[,-2]
    }
    
    else if (sum(is.na(data))==0){
      if (m==1) {data=data.frame(.imp=1,data)}
      if (m>1) {
        imp=data
        for (count in 1:(m-1)){
          imp=rbind(imp,data) #do we need to random sampling the data???
        }
        data=data.frame(.imp=rep(1:m,each=nrow(data)), imp)
      }
    }
  }
  
  
  if (method=="Ind"){
    newdata=NULL
    vname=NULL
    for (j in 2:ncol(data)){
      ss=which(is.na(data[,j]))
      if (length(ss)>0){
        data[ss,j]=mean(data[,j],na.rm=TRUE)
        vname=c(vname,paste(names(data)[j],"_mis",sep=""))
        newdata=cbind(newdata,as.numeric(is.na(data[,j])))
      }
    }
    vname=c(names(data),vname)
    data=data.frame(cbind(data,newdata))
    names(data)=vname				
  }	
  return(data)
}


generateX2 = function(Sigma_e,W,rhox){
  n = dim(W)[1]
  p = dim(W)[2]
  #Step 1: Sample a lot of WW from the empirical distribution of W and sample a lot of EE from the estimated distribution of measurement error independently, 
  NN=10000
  WW=EE=XX=matrix(NA,nrow=NN,ncol=p)
  for (j in 1:p){
    WW[,j] = W[sample(1:n,NN,replace = T),j]
  }
  
  EE = mvrnorm(n=NN,mu=rep(0,p),Sigma=Sigma_e)
  
  for (j in 1:p){
    XX[,j] = WW[,j] - EE[,j]
  }
  
  #Step 2:Sample n individual Z_i from multivariate normal distribution N(0,AR(\rho))
  Z<-mvrnorm(n=n,mu=rep(0,p),Sigma=rhox^abs(outer(1:p,1:p,"-")))
  #Step 3: For each i, j, compute X_{ij}= F_j^{-1}(\Phi(Z_{ij}))
  
  Xtilde=matrix(NA,nrow=n,ncol=p)
  for (i in 1:n) {
    for (j in 1:p) {
      Xtilde[i, j] <- quantile(XX[, j], pnorm(Z[i, j]))
    }
  }
  #F_XX[[j]] = ecdf(XX[,j])
  W_new = Xtilde +  mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma_e)
  return(Xtilde)
}

