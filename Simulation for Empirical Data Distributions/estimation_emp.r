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





datagen<-function(X,rho,effect,interc){
  ###n sample size
  ###p number of predictors
  ##s proportion of non-zero element
  ###rho correlation parameter for design matrix X
  ###effect: effect scale multiplier
  n=dim(X)[1]
  p=dim(X)[2]
  s<-p*0.3 ### keep the proportion of sparsity at 30% ###update
  b0=rep(c(3,1.5,0,0,2,0,0),s/3)
  beta=c(b0,1,rep(0,p-length(b0)-1))*effect
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

generate_intercept = function(datal, missing,proMis,sigmae,scale,type){
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
  colindex_t = sample(which(beta!=0), size=S*(1/10),replace = F)
  ###randomly choose missing column for non-true signals 
  colindex_c = sample(which(beta==0), size=(p-S)*(1/10),replace = F)
  colindex = c(colindex_t, colindex_c)[order(c(colindex_t, colindex_c))]
  
  ### generating missing
  intercept=rep(0,length(colindex))
  a<-list()
  R<-matrix(NA,nrow=n,ncol=p)
  logit <- function(x) (1+exp(-x))^(-1)
  prob<-matrix(NA,nrow=n,ncol=length(colindex))
  Xnew<-X
  if(missing%in%c("MCAR")){intercept=rep(0,length(colindex))}
  
  if(missing%in%c("MAR")){
    
    for (k in 1:length(colindex)) {
      set.seed(1234+k)
      prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%rep(1,p))-intercept[k]))  
    }
    aprob=apply(prob,2,mean)
    ### add a number to make sure the missing proportion is around a certain point
    for (k in 1:length(colindex)) {
      while (abs(aprob[k] - proMis) >0.01){
        set.seed(1234+k)
        intercept[k]=intercept[k] + ifelse(aprob[k] > 0.15, 0.05, -0.05)
        prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%rep(1,p))-intercept[k]))  
        aprob[k]=mean(prob[,k])
      }
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
      while (abs(aprob[k] - proMis) >0.01){
        set.seed(1234+k)
        intercept[k]=intercept[k] + ifelse(aprob[k] > 0.15, 0.05, -0.05)
        prob[,k] <- as.vector(logit((cbind(1,X[,1:(colindex[k]-1)])%*%rnorm(n=length(1:(colindex[k]-1))+1))-intercept[k])) 
        aprob[k]=mean(prob[,k])
      }
    }
  }
  
  return(intercept)
}

###generate some missing data in X or W
datamissW = function(datal, missing,inputseed,proMis,sigmae,scale,type,intercept){
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
  colindex_t = sample(which(beta!=0), size=S*(1/10),replace = F)
  ###randomly choose missing column for non-true signals 
  colindex_c = sample(which(beta==0), size=(p-S)*(1/10),replace = F)
  colindex = c(colindex_t, colindex_c)[order(c(colindex_t, colindex_c))]
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
      a[[k]] <- rbinom(n,1,proMis)
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("MAR")){
    
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      prob[,k] <- as.vector(logit((cbind(1,X[,-colindex[k]])%*%rep(1,p))-intercept[k]))  
      a[[k]] <- rbinom(n,1,prob[,k])
      R[,colindex[k]] <- a[[k]]
      Xnew[,colindex[k]] <- ifelse(R[,colindex[k]]==0, Xnew[,colindex[k]] ,NA)
    }
  }
  
  if(missing%in%c("Monotone")){
    
    set.seed(inputseed)
    for (k in 1:length(colindex)) {
      prob[,k] <- as.vector(logit((cbind(1,X[,1:(colindex[k]-1)])%*%rnorm(n=length(1:(colindex[k]-1))+1))-intercept[k])) 
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


####Knockoff, data are with 1 column of Y followed by p columns of X
myselectall=list() #return value
offset=c() #offset=1 or 0
q=c() #q=0.1 or 0.2
myest<-function(data,q,method,stat,flip,family,Sigma,delta,offset,MEAN=NULL,COV=NULL){
  ###Method: Fixed: fixed knockoff; Second: Second order model-X; Deep: deep knockoff
  ###flip: Max: signed max; Diff: difference
  ###offset: 0: Knockoff; 1: Knockoff+
  ###GMUL only support binary and poisson
  Y=data[,1]
  X=data[,-1]
  ###center columns
  X=scale(X,center=TRUE,scale=FALSE)
  p=ncol(X)
  ####knockoff creation
  if (method=="Fixed"){
    X<-X%*%diag(1/sqrt(colSums(X^2)))
    XKnock<-create.fixed(X)$Xk
  }
  if (method=="Second"){
    XKnock<-create.second_order(X)
  }	
  
  if (method=="Gaussian"){
    XKnock<-create.gaussian(X,mu=MEAN,Sigma=COV)
  }
  
  if (method=="Deep"){
    ###TBA:need get Python code work under R
  }
  fdata=data.frame(cbind(Y,X,XKnock))
  names(fdata)=c("Y",paste("X",1:p,sep=""),paste("Xk",1:p,sep=""))
  if (stat=="Lasso"){
    cvfit<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
    fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit$lambda.min)) 
    Z=abs(coef(fit)[1+(1:p)])
    Ztilde=abs(coef(fit)[1+p+(1:p)])
  }
  
  if (stat=="Lasso1"){
    cvfit<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) #updated
    fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit$lambda.min)/5) #updated lambda
    Z=abs(coef(fit)[1+(1:p)])
    Ztilde=abs(coef(fit)[1+p+(1:p)])
  }
  
  if (stat=="Lasso_Order"){
    fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y), nlambda = max(1000,5*p), family=family,standardize=FALSE,intercept=TRUE)
    Z=Ztilde=rep(NA,p)
    for (i in 1:p){
      Z[i]=fit$lambda[min(which(coef(fit)[1+i,]!=0))]
      Ztilde[i]=fit$lambda[min(which(coef(fit)[1+p+i,]!=0))]
      Z[which(is.na(Z))]=0
      Ztilde[which(is.na(Ztilde))]=0
    }
  }
  
  if (stat=="RF"){
    
    #fit=randomForest(as.factor(Y)~.,data=fdata,importance = TRUE)
    #tmp=importance(fit)
    fit=ranger(as.factor(Y)~.,data=fdata,importance = 'impurity')
    tmp = fit$variable.importance
    Z=tmp[1:p]
    Ztilde=tmp[p+(1:p)]
  }
  
  if (stat %in% c("CocoLasso","CocoLasso2", "CocoLasso3", "CocoLasso4")){
    cvfit_Lasso<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
    sum=0
    index=0
    while (sum==0){
    fit_Lasso<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit_Lasso$lambda[index+1])) #updated lambda
    sum=sum(abs(coef(fit_Lasso)[1+(1:p)]))
    index=index+1
    }
    if (sum<10^(-5)){ ## the lasso coefficient is 0, use the least square as the initial
      XX=cbind(X,XKnock)
      beta=solve(t(XX)%*%XX-Sigma)%*%t(XX)%*%Y
      sum=sum(abs(beta))
      radiix=seq(sum/1000,100*sum,length=1000)
      fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix)
      fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix[which.max(apply(fit1$betaCorr!=0,2,sum))])
      Z=abs(fit$betaCorr[1:p,1])
      Ztilde=abs(fit$betaCorr[p+(1:p),])
    } 
    
    if (sum>=10^(-5)) {
      #radiix=seq(log(sum/10),log(10*sum),length=10)
      eradiix=sum/10
      fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix)
      fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
      estimate=rbind(eradiix,fit1$betaCorr) ## The first row the radius, each column from 2:134 is estimated beta corresponded to each radius.
      count=0
      
      while (sum(abs(fit1$betaCorr)>10^(-3))<30 & count<1000) {
        count=count+1
        add_eradiix=eradiix*exp(0.5*count)
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
        #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
        add = rbind(add_eradiix,fit1$betaCorr)
        estimate=cbind(estimate,add)
        #write.csv(estimate,"estimate.csv")
      }
      
      ###use the one and the previous one to a check
      countx=count
      radiix=seq(log(eradiix*exp(0.5*(countx-1))),log(eradiix*exp(0.5*countx)),length=1000)
      estimate2=matrix(NA,ncol=length(radiix),nrow=dim(X)[2]*2+1)
      for (ii in 1:length(radiix)){
        fit2=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=exp(radiix[ii]))
        estimate2[,ii]=rbind(exp(radiix[ii]),fit2$betaCorr)
      }
      #write.csv(estimate2,"estimate2.csv")
      
      estimate=estimate2
      
        betaxx=estimate[-1,which(apply(estimate[-1,]!=0,2,sum)>=30)[1]]
        Z=abs(betaxx[1:p])
        Ztilde=abs(betaxx[p+(1:p)])
    }
  }
  
  if (stat %in% c("CocoLasso_Order", "CocoLasso_Order2","CocoLasso_Order3","CocoLasso_Order4")){
    cvfit_Lasso<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
    sum=0
    index=0
    while (sum==0){
      fit_Lasso<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit_Lasso$lambda[index+1])) #updated lambda
      sum=sum(abs(coef(fit_Lasso)[1+(1:p)]))
      index=index+1
    }
    
    if (sum<10^(-5)){ ## the lasso coefficient is 0, use the least square as the initial
    XX=cbind(X,XKnock)
    beta=solve(t(XX)%*%XX-Sigma)%*%t(XX)%*%Y
    sum=sum(abs(beta))
    radiix=seq(sum/1000,100*sum,length=1000)
    fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix)
    fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix[which.max(apply(fit1$betaCorr!=0,2,sum))])
    Z=abs(fit$betaCorr[1:p,1])
    Ztilde=abs(fit$betaCorr[p+(1:p),])
    } 
    
    if (sum>=10^(-5)) {
      #radiix=seq(log(sum/10),log(10*sum),length=10)
      eradiix=sum/10
      fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix)
      fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
      estimate=rbind(eradiix,fit1$betaCorr) ## The first row the radius, each column from 2:134 is estimated beta corresponded to each radius.
      
    
    ### find the upper bound
      countu=0
    while (sum(abs(fit1$betaCorr)>10^(-3))<30 & countu<1000 ) {
      countu=countu+1
      add_eradiix=eradiix*exp(0.5*countu)
      fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
      #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
      add = rbind(add_eradiix,fit1$betaCorr)
      estimate=cbind(estimate,add)
      #write.csv(estimate,"estimate.csv")
    }
      
    ### find the lower bound
      countl=0
      while (sum(abs(fit1$betaCorr)>10^(-3))<=1 & countl<countu ) {
        countl=countl+1
        add_eradiix=eradiix*exp(0.5*countl)
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
        #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
        add = rbind(add_eradiix,fit1$betaCorr)
        estimate=cbind(estimate,add)
        #write.csv(estimate,"estimate.csv")
      }
      
      ###use the one and lower bound to a check

      radiix=seq(log(eradiix*exp(0.5*(countl))),log(eradiix*exp(0.5*countu)),length=1000)
      estimate2=matrix(NA,ncol=length(radiix),nrow=dim(X)[2]*2+1)
      for (ii in 1:length(radiix)){
      fit2=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=exp(radiix[ii]))
      estimate2[,ii]=rbind(exp(radiix[ii]),fit2$betaCorr)
      }
      #write.csv(estimate2,"estimate2.csv")
      
      estimate=estimate2
      ### order cocoLasso
      ss=rep(NA,2*p) ###last no
      for (j in 2:(2*p+1)){
        aaa=which(estimate[j,]!=0)
        if (length(aaa)>=1){
          ss[j-1]=1/estimate[1,which(estimate[j,]!=0)[1]]
        }else{ss[j-1]=0}
      }
      Z=abs(ss[1:p])
      Ztilde=abs(ss[p+(1:p)])
    }
  }
  
  #by second derivative
  if (stat=="GMUS"){
    ### select the best delta
    fit1=gmus(W=cbind(X,XKnock),y=c(Y),family=family)
    param=coef(fit1)
    param=param[order(param$delta),]
    deriv <- function(x, y) diff(y) / diff(x)
    second.deriva=deriv(param$delta[-1], deriv(param$delta, param$nonzeros))
    
    fit=gmus(W=cbind(X,XKnock),y=c(Y),family=family,delta=param$delta[which.max(abs(second.deriva))]) 
    Z=abs(fit$beta[1:p,1])
    Ztilde=abs(fit$beta[p+(1:p),1])
  }
  if (stat=="GMUL"){
    ### select the best delta
    fit1=gmus(W=cbind(X,XKnock),y=c(Y),family=family)
    param=coef(fit1)
    param=param[order(param$delta),]
    deriv <- function(x, y) diff(y) / diff(x)
    second.deriva=deriv(param$delta[-1], deriv(param$delta, param$nonzeros))
    
    fit=gmus(W=cbind(X,XKnock),y=c(Y),family=family,delta=param$delta[which.max(abs(second.deriva))]) 
    Z=abs(fit$beta[1:p,1])
    Ztilde=abs(fit$beta[p+(1:p),1])
  }
  if (stat=="GDS_cv"){
    cvfit=cv_gds(X=cbind(X,XKnock),y=c(Y),family=family,no_lambda=max(1000,5*p))
    fit=gds(X=cbind(X,XKnock),y=c(Y),family=family, lambda = cvfit$lambda_min)
    Z=abs(fit$beta[1:p,1])
    Ztilde=abs(fit$beta[p+(1:p),1])
  }
  if (stat=="GDS"){
    fit=gds(X=cbind(X,XKnock),y=c(Y),family=family)
    Z=abs(fit$beta[1:p,1])
    Ztilde=abs(fit$beta[p+(1:p),1])
  }
  if (flip=="Max"){
    W=pmax(Z,Ztilde)*(-1)^as.numeric(Ztilde>=Z)
  }
  if (flip=="Diff"){
    W=Z-Ztilde
  }
  
  count=0
  for (offset0 in offset){
    for (q0 in q) {
      count=count+1
      mythred=knockoff.threshold(W,fdr=q0,offset=offset0)
      myselectall[[count]]=which(W>=mythred)
    }
  }
  return (list(myselect=myselectall,Z=Z,Ztilde=Ztilde,Xknock=XKnock))
}




#####with multiple imputation
####Knockoff, data are with 1 column of imp, 1 column of Y followed by p columns of X
myselectall=list() #return value
offset=c() #offset=1 or 0
q=c() #q=0.1 or 0.2
myest_mi<-function(data,q,method,stat,flip,family,Sigma,delta,offset,MEAN=NULL,COV=NULL){
  ###Mehtod: Fixed: fixed knockoff; Second: Second order model-X; Deep: deep knockoff
  ###flip: Max: signed max; Diff: difference
  ###offset: 0: Knockoff; 1: Knockoff+
  ###GMUL only support binary and poisson
  BB=max(data$.imp)
  p=ncol(data)-2
  Zmat=Ztildemat=matrix(data=NA,nrow=BB,ncol=p)
  Xknockmat=matrix(NA,nrow=dim(data)[1],ncol=p)
  for (bb in 1:BB){
    newdata=data[which(data$.imp==bb),-1] ### change data=data[which(data$.imp==bb),-1] to newdata: meet error
    Y=newdata[,1] 
    X=as.matrix(newdata[,-1])
    ###center columns
    X=scale(X,center=TRUE,scale=FALSE)
    ####knockoff creation
    if (method=="Fixed"){
      X<-X%*%diag(1/sqrt(colSums(X^2)))
      XKnock<-create.fixed(X)$Xk
    }
    if (method=="Second"){
      XKnock<-create.second_order(X)
    }	
    
    if (method=="Gaussian"){
      XKnock<-create.gaussian(X,mu=MEAN,Sigma=COV)
    }
    
    if (method=="Deep"){
      ###TBA:need get Pthon code work under R
    }
    fdata=data.frame(cbind(Y,X,XKnock))
    names(fdata)=c("Y",paste("X",1:p,sep=""),paste("Xk",1:p,sep=""))
    if (stat=="Lasso"){
      cvfit<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
      fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit$lambda.min)) 
      Z=abs(coef(fit)[1+(1:p)])
      Ztilde=abs(coef(fit)[1+p+(1:p)])
    }
    
    if (stat=="Lasso1"){
      cvfit<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) #updated
      fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit$lambda.min)/5) #updated lambda
      Z=abs(coef(fit)[1+(1:p)])
      Ztilde=abs(coef(fit)[1+p+(1:p)])
    }
    
    if (stat=="Lasso_Order"){
      fit<-glmnet(as.matrix(cbind(X,XKnock)),c(Y), nlambda = max(1000,5*p), family=family,standardize=FALSE,intercept=TRUE)
      Z=Ztilde=rep(NA,p)
      for (i in 1:p){
        Z[i]=fit$lambda[min(which(coef(fit)[1+i,]!=0))]
        Ztilde[i]=fit$lambda[min(which(coef(fit)[1+p+i,]!=0))]
        Z[which(is.na(Z))]=0
        Ztilde[which(is.na(Ztilde))]=0
      }
    }
    if (stat=="RF"){
      #fit=randomForest(as.factor(Y)~.,data=fdata)
      #tmp=importance(fit)
      fit=ranger(as.factor(Y)~.,data=fdata,importance = 'impurity')
      tmp = fit$variable.importance
      Z=tmp[1:p]
      Ztilde=tmp[p+(1:p)]
    }
    
    if (stat %in% c("CocoLasso","CocoLasso2", "CocoLasso3", "CocoLasso4")){
      cvfit_Lasso<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
      sum=0
      index=0
      while (sum==0){
        fit_Lasso<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit_Lasso$lambda[index+1])) #updated lambda
        sum=sum(abs(coef(fit_Lasso)[1+(1:p)]))
        index=index+1
      }
      
      if (sum<10^(-5)){ ## the lasso coefficient is 0, use the least square as the initial
        XX=cbind(X,XKnock)
        beta=solve(t(XX)%*%XX-Sigma)%*%t(XX)%*%Y
        sum=sum(abs(beta))
        radiix=seq(sum/1000,100*sum,length=1000)
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix)
        fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix[which.max(apply(fit1$betaCorr!=0,2,sum))])
        Z=abs(fit$betaCorr[1:p,1])
        Ztilde=abs(fit$betaCorr[p+(1:p),])
      } 
      
      if (sum>=10^(-5)) {
        #radiix=seq(log(sum/10),log(10*sum),length=10)
        eradiix=sum/10
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix)
        fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
        estimate=rbind(eradiix,fit1$betaCorr) ## The first row the radius, each column from 2:134 is estimated beta corresponded to each radius.
        count=0
        
        while (sum(abs(fit1$betaCorr)>10^(-3))<30 & count<1000) {
          count=count+1
          add_eradiix=eradiix*exp(0.5*count)
          fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
          #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
          add = rbind(add_eradiix,fit1$betaCorr)
          estimate=cbind(estimate,add)
          #write.csv(estimate,"estimate.csv")
        }
        
        ###use the one and the previous one to a check
        countx=count
        radiix=seq(log(eradiix*exp(0.5*(countx-1))),log(eradiix*exp(0.5*countx)),length=1000)
        estimate2=matrix(NA,ncol=length(radiix),nrow=dim(X)[2]*2+1)
        for (ii in 1:length(radiix)){
          fit2=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=exp(radiix[ii]))
          estimate2[,ii]=rbind(exp(radiix[ii]),fit2$betaCorr)
        }
        #write.csv(estimate2,"estimate2.csv")
        
        estimate=estimate2
        
        betaxx=estimate[-1,which(apply(estimate[-1,]!=0,2,sum)>=30)[1]]
        Z=abs(betaxx[1:p])
        Ztilde=abs(betaxx[p+(1:p)])
      }
    }
    
    if (stat %in% c("CocoLasso_Order", "CocoLasso_Order2","CocoLasso_Order3","CocoLasso_Order4")){
      cvfit_Lasso<-cv.glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE) 
      sum=0
      index=0
      while (sum==0){
        fit_Lasso<-glmnet(as.matrix(cbind(X,XKnock)),c(Y),family=family,standardize=FALSE,intercept=TRUE,lambda=(cvfit_Lasso$lambda[index+1])) #updated lambda
        sum=sum(abs(coef(fit_Lasso)[1+(1:p)]))
        index=index+1
      }
      
      ##try
      if (sum<10^(-5)){ ## the lasso coefficient is 0, use the least square as the initial
        XX=cbind(X,XKnock)
        beta=solve(t(XX)%*%XX-Sigma)%*%t(XX)%*%Y
        sum=sum(abs(beta))
        radiix=seq(sum/1000,100*sum,length=1000)
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix)
        fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=radiix[which.max(apply(fit1$betaCorr!=0,2,sum))])
        Z=abs(fit$betaCorr[1:p,1])
        Ztilde=abs(fit$betaCorr[p+(1:p),])
      } 
      
      if (sum>=10^(-5)) {
        #radiix=seq(log(sum/10),log(10*sum),length=10)
        eradiix=sum/10
        fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix)
        fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
        estimate=rbind(eradiix,fit1$betaCorr) ## The first row the radius, each column from 2:134 is estimated beta corresponded to each radius.
        
        
        ### find the upper bound
        countu=0
        while (sum(abs(fit1$betaCorr)>10^(-3))<30 & countu<1000 ) {
          countu=countu+1
          add_eradiix=eradiix*exp(0.5*countu)
          fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
          #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
          add = rbind(add_eradiix,fit1$betaCorr)
          estimate=cbind(estimate,add)
          #write.csv(estimate,"estimate.csv")
        }
        
        ### find the lower bound
        countl=0
        while (sum(abs(fit1$betaCorr)>10^(-3))<=1 & countl<countu ) {
          countl=countl+1
          add_eradiix=eradiix*exp(0.5*countl)
          fit1=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix)
          #fit=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=add_eradiix[which.max(apply(abs(fit1$betaCorr)>10^(-3),2,sum))])
          add = rbind(add_eradiix,fit1$betaCorr)
          estimate=cbind(estimate,add)
          #write.csv(estimate,"estimate.csv")
        }
        
        ###use the one and lower bound to a check
        
        radiix=seq(log(eradiix*exp(0.5*(countl))),log(eradiix*exp(0.5*countu)),length=1000)
        estimate2=matrix(NA,ncol=length(radiix),nrow=dim(X)[2]*2+1)
        for (ii in 1:length(radiix)){
          fit2=corrected_lasso(W=cbind(X,XKnock),y=c(Y),sigmaUU=Sigma,family=family,radii=exp(radiix[ii]))
          estimate2[,ii]=rbind(exp(radiix[ii]),fit2$betaCorr)
        }
        #write.csv(estimate2,"estimate2.csv")
        
        estimate=estimate2
        ### order cocoLasso
        ss=rep(NA,2*p) ###last no
        for (j in 2:(2*p+1)){
          aaa=which(estimate[j,]!=0)
          if (length(aaa)>=1){
            ss[j-1]=1/estimate[1,which(estimate[j,]!=0)[1]]
          }else{ss[j-1]=0}
        }
        Z=abs(ss[1:p])
        Ztilde=abs(ss[p+(1:p)])
      }
    }
    
    
    #by second derivative
    if (stat=="GMUS"){
      ### select the best delta
      fit1=gmus(W=cbind(X,XKnock),y=c(Y),family=family)
      param=coef(fit1)
      param=param[order(param$delta),]
      deriv <- function(x, y) diff(y) / diff(x)
      second.deriva=deriv(param$delta[-1], deriv(param$delta, param$nonzeros))
      
      fit=gmus(W=cbind(X,XKnock),y=c(Y),family=family,delta=param$delta[which.max(abs(second.deriva))]) 
      Z=abs(fit$beta[1:p,1])
      Ztilde=abs(fit$beta[p+(1:p),1])
    }
    if (stat=="GMUL"){
      ### select the best delta
      fit1=gmus(W=cbind(X,XKnock),y=c(Y),family=family)
      param=coef(fit1)
      param=param[order(param$delta),]
      deriv <- function(x, y) diff(y) / diff(x)
      second.deriva=deriv(param$delta[-1], deriv(param$delta, param$nonzeros))
      
      fit=gmus(W=cbind(X,XKnock),y=c(Y),family=family,delta=param$delta[which.max(abs(second.deriva))]) 
      Z=abs(fit$beta[1:p,1])
      Ztilde=abs(fit$beta[p+(1:p),1])
    }
    if (stat=="GDS_cv"){
      cvfit=cv_gds(X=cbind(X,XKnock),y=c(Y),family=family,no_lambda=max(1000,5*p))
      fit=gds(X=cbind(X,XKnock),y=c(Y),family=family, lambda = cvfit$lambda_min)
      Z=abs(fit$beta[1:p,1])
      Ztilde=abs(fit$beta[p+(1:p),1])
    }
    if (stat=="GDS"){
      fit=gds(X=cbind(X,XKnock),y=c(Y),family=family)
      Z=abs(fit$beta[1:p,1])
      Ztilde=abs(fit$beta[p+(1:p),1])
    }
    
    Zmat[bb,]=Z
    Ztildemat[bb,]=Ztilde
    Xknockmat[(((bb-1)*dim(XKnock)[1]+1)):(bb*dim(XKnock)[1]),]=XKnock
  }
  
  PP = ifelse(Ztildemat>Zmat,1,0)
  PP1 = ifelse(Ztildemat==Zmat,1/2,0)
  PP3 = PP+PP1*2
  p_value = apply(PP3,2,function(X){(sum(X)+1)/(length(X)+1)})
  
  if (flip=="Diff"){
    W=apply(abs(Zmat - Ztildemat),2,max)
  }
  
  seq_thres=function(W,p_value,fdr,offset){
    
  re_p_value=p_value[order(-W)]
  ratio.calculate = function(k) {(offset + sum(re_p_value[1:k] > 0.5))/max(1,sum(re_p_value[1:k] <= 0.5))}
  ratio=rep(NA,p)
  for (k in 1:p){
   ratio[k]=ratio.calculate(k)
  }
  ok = which(ratio <= fdr)
  if (length(ok)>0) {
    KK = max(ok)
    set=order(-W)[intersect(1:KK,which(re_p_value<=0.5))]
  }
  else {set=NULL}
  return(set)
  }
  
  count=0
  for (offset0 in offset){
    for (q0 in q) {
      count=count+1
      temp_result=seq_thres(W,p_value,fdr=q0,offset=offset0)
      if (length(temp_result)>0){
      myselectall[[count]]=temp_result
      }
      if (length(temp_result)==0){
        myselectall[[count]]=NA
      }
    }
  }
  return (list(myselect=myselectall,Zmat=Zmat,Ztildemat=Ztildemat,Xknock=Xknockmat))
}


