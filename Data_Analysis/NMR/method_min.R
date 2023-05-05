setwd("./")
source("estimation.r")

final_brease_p=read.csv("final_breast_match_nmr.csv",header=T)
final_Colore_p=read.csv("final_Colore_match_nmr.csv",header=T)

final_brease_p_AC=read.csv("final_breast_nmr.csv",header=T)
final_Colore_p_AC=read.csv("final_Colore_nmr.csv",header=T)

Sigma_e=read.csv("nmr_Sigma_e.csv",header=T)

final_brease_p=final_brease_p[,-1]
final_Colore_p=final_Colore_p[,-1]
final_brease_p_AC=final_brease_p_AC[,-1]
final_Colore_p_AC=final_Colore_p_AC[,-1]
Sigma_e=as.matrix(Sigma_e[,-1])

Sigma_ee=cbind(rbind(Sigma_e,Sigma_e),rbind(Sigma_e,Sigma_e))
p=dim(Sigma_e)[1]

###Changed to 20 parallel jobs each with 5 runs
rep=5 ### need to change in the server
args=commandArgs(trailingOnly=TRUE)
num=as.numeric(args[1])
iniseed=num*rep


###Different selection method
real=function(data,method,Sigmax,deltax,name){
  p=dim(data)[2]-1
  select=list()
  beta_select1=matrix(NA,nrow=p,ncol=rep)
  rownames(beta_select1)=names(data)[-1]
  beta_select2=matrix(NA,nrow=p,ncol=rep)
  rownames(beta_select2)=names(data)[-1]
  beta_select3=matrix(NA,nrow=p,ncol=rep)
  rownames(beta_select3)=names(data)[-1]
  beta_select4=matrix(NA,nrow=p,ncol=rep)
  rownames(beta_select4)=names(data)[-1]
  Z=matrix(NA,nrow=p,ncol=rep)
  Ztilde=matrix(NA,nrow=p,ncol=rep)
  seed=c()
  
  for (k in 1:rep){
    try({
    seed[k]=iniseed+k
    set.seed(seed[k])
    select[[k]]=myest(data=data,q=c(0.1,0.2),method="Second",stat=method,flip="Max",family="binomial",Sigma=Sigmax,delta=deltax,offset=c(0,1))
    
    Z[,k]= select[[k]]$Z
    Ztilde[,k]=select[[k]]$Ztilde
    rownames(Z)=names(data)[-1]
    rownames(Ztilde)=names(data)[-1]
    
    beta_select1[select[[k]]$myselect[[1]],k]=1#offset=0, q=0.1
    beta_select1[setdiff(1:p,select[[k]]$myselect[[1]]),k]=0
    beta_select2[select[[k]]$myselect[[2]],k]=1#offset=0, q=0.2
    beta_select2[setdiff(1:p,select[[k]]$myselect[[2]]),k]=0
    beta_select3[select[[k]]$myselect[[3]],k]=1#offset=1, q=0.1
    beta_select3[setdiff(1:p,select[[k]]$myselect[[3]]),k]=0
    beta_select4[select[[k]]$myselect[[4]],k]=1#offset=1, q=0.2
    beta_select4[setdiff(1:p, select[[k]]$myselect[[4]]),k]=0
  })
  }
  
  sumselect1=apply(beta_select1,1,sum)
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  
  sumselect2=apply(beta_select2,1,sum)
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  
  sumselect3=apply(beta_select3,1,sum)
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  
  sumselect4=apply(beta_select4,1,sum)
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  
  
  File1=sprintf("result/Z_%s_%s_%d.csv",method,name,iniseed)
  File2=sprintf("result/Ztilde_%s_%s_%d.csv",method,name,iniseed)
  File3=sprintf("result/beta1_%s_%s_%d.csv",method,name,iniseed)
  File4=sprintf("result/beta2_%s_%s_%d.csv",method,name,iniseed)
  File5=sprintf("result/beta3_%s_%s_%d.csv",method,name,iniseed)
  File6=sprintf("result/beta4_%s_%s_%d.csv",method,name,iniseed)
  File7=sprintf("result/sumbeta1_%s_%s_%d.csv",method,name,iniseed)
  File8=sprintf("result/sumbeta2_%s_%s_%d.csv",method,name,iniseed)
  File9=sprintf("result/sumbeta3_%s_%s_%d.csv",method,name,iniseed)
  File10=sprintf("result/sumbeta4_%s_%s_%d.csv",method,name,iniseed)
  
  write.csv(Z,File1)
  write.csv(Ztilde,File2)
  write.csv(beta_select1,File3)
  write.csv(beta_select2,File4)
  write.csv(beta_select3,File5)
  write.csv(beta_select4,File6)
  write.csv(sumselect1,File7)
  write.csv(sumselect2,File8)
  write.csv(sumselect3,File9)
  write.csv(sumselect4,File10)
  
}

### brease cancer 
#real(data=final_brease_p,method="Lasso",Sigmax=NULL,deltax=NULL,name='pb')
#real(data=final_brease_p,method="Lasso1",Sigmax=NULL,deltax=NULL,name='pb')
#real(data=final_brease_p,method="Lasso_Order",Sigmax=NULL,deltax=NULL,name='pb')
#real(data=final_brease_p,method="RF",Sigmax=NULL,deltax=NULL,name='pb')
#real(data=final_brease_p,method="GMUL",Sigmax=NULL,deltax=NULL,name='pb')
#real(data=final_brease_p,method="GDS",Sigmax=NULL,deltax=NULL,name='pb')
real(data=final_brease_p,method="CocoLasso",Sigmax=Sigma_ee,deltax=NULL,name='pb')
real(data=final_brease_p,method="CocoLasso2",Sigmax=as.matrix(bdiag(Sigma_e,Sigma_e)),deltax=NULL,name='pb')

### colore cancer
#real(data=final_Colore_p,method="Lasso",Sigmax=NULL,deltax=NULL,name='pc')
#real(data=final_Colore_p,method="Lasso1",Sigmax=NULL,deltax=NULL,name='pc')
#real(data=final_Colore_p,method="Lasso_Order",Sigmax=NULL,deltax=NULL,name='pc')
#real(data=final_Colore_p,method="RF",Sigmax=NULL,deltax=NULL,name='pc')
#real(data=final_Colore_p,method="GMUL",Sigmax=NULL,deltax=NULL,name='pc')
#real(data=final_Colore_p,method="GDS",Sigmax=NULL,deltax=NULL,name='pc')
real(data=final_Colore_p,method="CocoLasso",Sigmax=Sigma_ee,deltax=NULL,name='pc')
real(data=final_Colore_p,method="CocoLasso2",Sigmax=as.matrix(bdiag(Sigma_e,Sigma_e)),deltax=NULL,name='pc')

### brease cancer ALL control
#real(data=final_brease_p_AC,method="Lasso",Sigmax=NULL,deltax=NULL,name='pbAC')
#real(data=final_brease_p_AC,method="Lasso1",Sigmax=NULL,deltax=NULL,name='pbAC')
#real(data=final_brease_p_AC,method="Lasso_Order",Sigmax=NULL,deltax=NULL,name='pbAC')
#real(data=final_brease_p_AC,method="RF",Sigmax=NULL,deltax=NULL,name='pbAC')
#real(data=final_brease_p_AC,method="GMUL",Sigmax=NULL,deltax=NULL,name='pbAC')
#real(data=final_brease_p_AC,method="GDS",Sigmax=NULL,deltax=NULL,name='pbAC')
real(data=final_brease_p_AC,method="CocoLasso",Sigmax=Sigma_ee,deltax=NULL,name='pbAC')
real(data=final_brease_p_AC,method="CocoLasso2",Sigmax=as.matrix(bdiag(Sigma_e,Sigma_e)),deltax=NULL,name='pbAC')
### colore cancer ALL control
#real(data=final_Colore_p_AC,method="Lasso",Sigmax=NULL,deltax=NULL,name='pcAC')
#real(data=final_Colore_p_AC,method="Lasso1",Sigmax=NULL,deltax=NULL,name='pcAC')
#real(data=final_Colore_p_AC,method="Lasso_Order",Sigmax=NULL,deltax=NULL,name='pcAC')
#real(data=final_Colore_p_AC,method="RF",Sigmax=NULL,deltax=NULL,name='pcAC')
#real(data=final_Colore_p_AC,method="GMUL",Sigmax=NULL,deltax=NULL,name='pcAC')
#real(data=final_Colore_p_AC,method="GDS",Sigmax=NULL,deltax=NULL,name='pcAC')
real(data=final_Colore_p_AC,method="CocoLasso",Sigmax=Sigma_ee,deltax=NULL,name='pcAC')
real(data=final_Colore_p_AC,method="CocoLasso2",Sigmax=as.matrix(bdiag(Sigma_e,Sigma_e)),deltax=NULL,name='pcAC')
