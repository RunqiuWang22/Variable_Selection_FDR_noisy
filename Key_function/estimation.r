####Knockoff, data are with 1 column of Y followed by p columns of X
myselectall=list() #return value
offset=c() #offset=1 or 0
q=c() #q=0.1 or 0.2
myest<-function(data,q,method,stat,flip,family,Sigma,delta,offset,MEAN=NULL,COV=NULL){
  ###Method: Fixed: fixed knockoff; Second: Second order model-X; Deep: deep knockoff
  ###flip: Max: signed max; Diff: difference
  ###offset: 0: Knockoff; 1: Knockoff+
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
    fit=randomForest(as.factor(Y)~.,data=fdata)
    tmp=importance(fit)
    Z=tmp[1:p,1]
    Ztilde=tmp[p+(1:p),1]
  }
  
  ###calculate the prediction value of Random forest
  if (stat=="RF_pred"){
    set.seed(1111)
    sample_index <- sample(seq_len(nrow(fdata)), size = 0.8 * nrow(fdata))
    training <- fdata[sample_index, ]
    testing <- fdata[-sample_index, ]
    fit=randomForest(as.factor(Y)~.,data=training)
    tmp=importance(fit)
    Z=tmp[1:p,1]
    Ztilde=tmp[p+(1:p),1]
    Y_hat = as.numeric(predict(fit,newdata=testing,type = "prob")[,2])
    YY = c(Y_hat,testing$Y)
  }
  
  if (stat=="CocoLasso"){
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
      fit=randomForest(as.factor(Y)~.,data=fdata)
      tmp=importance(fit)
      Z=tmp[1:p,1]
      Ztilde=tmp[p+(1:p),1]
    }
    
    ###calculate the prediction value of Random forest
    if (stat=="RF_pred"){
      set.seed(1111)
      sample_index <- sample(seq_len(nrow(fdata)), size = 0.8 * nrow(fdata))
      training <- fdata[sample_index, ]
      testing <- fdata[-sample_index, ]
      fit=randomForest(as.factor(Y)~.,data=training)
      tmp=importance(fit)
      Z=tmp[1:p,1]
      Ztilde=tmp[p+(1:p),1]
      Y_hat = as.numeric(predict(fit,newdata=testing,type = "prob")[,2])
      YY = c(Y_hat,testing$Y)
    }
    
    if (stat=="CocoLasso"){
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
