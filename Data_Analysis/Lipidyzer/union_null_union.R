###read the rowname for multi
setwd("/Users/runqiuwang/Downloads/Dr.Dai's project/2022Spring_project/real/lipid/Multi")
read=read.csv("final_breast_GCMS.csv",header=T)
zname = names(read)[-c(1,2,3)]
### combine Z and Ztilde
### combine Z and Ztilde
com_multi = function(ZZ,stat,name){
  
  File1 = sprintf("%s_%s_%s_%d.csv",ZZ,stat,name,0)
  X = read.csv(File1, header=T)
  X = X[,-1]
  p=dim(X)[2]
  addmean=matrix(NA,nrow=5,ncol=p)
  
  Z=NULL
  for (i in 1:20){
    File2 = sprintf("%s_%s_%s_%d.csv",ZZ,stat,name,5*i-5)
    if (file.exists(File2)==T){
      add = read.csv(File2, header=T)      
      add = add[,-1]
      for (k in 1:5){
        addmean[k,] = as.vector(apply(add[((5*k-4):(5*k)),],2,mean))
      }
      Z = rbind(Z,addmean)
    }
  }  
  
  ### add names to the Z
  TZ = t(Z)
  rownames(TZ)=zname
  TZ
}

###Multi
###analytic sample (include y)
setwd("/Users/runqiuwang/Downloads/Dr.Dai's project/2022Spring_project/real/GCMS/Multi/result")
union1 = function(stat){
  Z_Lasso_b = com_multi(ZZ="Z",stat=stat,name="pbAC")
  Z_Lasso_c = com_multi(ZZ="Z",stat=stat,name="pcAC")
  Ztilde_Lasso_b = com_multi(ZZ="Ztilde",stat=stat,name="pbAC")
  Ztilde_Lasso_c = com_multi(ZZ="Ztilde",stat=stat,name="pcAC")
  
  p=dim(Z_Lasso_b)[1]
  mythred1=mythred2=mythred3=mythred4=list()
  myselect1=myselect2=myselect3=myselect4=list()
  Z=Ztilde=W=matrix(NA,nrow=p, ncol=100)
  beta_select1=beta_select2=beta_select3=beta_select4=matrix(0,nrow=p,ncol=100)
  rownames(beta_select1)=rownames(Z_Lasso_b)
  rownames(beta_select2)=rownames(Z_Lasso_b)
  rownames(beta_select3)=rownames(Z_Lasso_b)
  rownames(beta_select4)=rownames(Z_Lasso_b)
  
  for (i in 1:100){
    Z[,i]=Z_Lasso_b[,i]*Z_Lasso_c[,i]+Ztilde_Lasso_b[,i]*Ztilde_Lasso_c[,i]
    Ztilde[,i]=Z_Lasso_b[,i]*Ztilde_Lasso_c[,i]+Z_Lasso_c[,i]*Ztilde_Lasso_b[,i]
    W[,i]=Z[,i]-Ztilde[,i]
    
    mythred1[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=0)
    myselect1[[i]]=which(W[,i]>=mythred1[[i]])
    beta_select1[myselect1[[i]],i]=1
    
    mythred2[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=1)
    myselect2[[i]]=which(W[,i]>=mythred2[[i]])
    beta_select2[myselect2[[i]],i]=1
    
    mythred3[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=0)
    myselect3[[i]]=which(W[,i]>=mythred3[[i]])
    beta_select3[myselect3[[i]],i]=1
    
    mythred4[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=1)
    myselect4[[i]]=which(W[,i]>=mythred4[[i]])
    beta_select4[myselect4[[i]],i]=1
  }
  
  sumselect1=apply(beta_select1,1,sum)
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  paste1 = paste(names(sumselect1)," (", sumselect1, "%)",sep="")
  sumselect1=cbind(sumselect1,paste1)
  
  sumselect2=apply(beta_select2,1,sum)
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  paste2 = paste(names(sumselect2)," (", sumselect2, "%)",sep="")
  sumselect2=cbind(sumselect2,paste2)
  
  sumselect3=apply(beta_select3,1,sum)
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  paste3 = paste(names(sumselect3)," (", sumselect3, "%)",sep="")
  sumselect3=cbind(sumselect3,paste3)
  
  sumselect4=apply(beta_select4,1,sum)
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  paste4 = paste(names(sumselect4)," (", sumselect4, "%)",sep="")
  sumselect4=cbind(sumselect4,paste4)
  
  File7=sprintf("analytic/sumbeta1_%s.csv",stat)
  File8=sprintf("analytic/sumbeta2_%s.csv",stat)
  File9=sprintf("analytic/sumbeta3_%s.csv",stat)
  File10=sprintf("analytic/sumbeta4_%s.csv",stat)
  
  write.csv(sumselect1,File7)
  write.csv(sumselect2,File8)
  write.csv(sumselect3,File9)
  write.csv(sumselect4,File10)
}


union1(stat="Lasso")
union1(stat="Lasso1")
union1(stat="Lasso_Order")
union1(stat="RF")
union1(stat="GMUL")
union1(stat="GDS")
union1(stat="CocoLasso")
union1(stat="CocoLasso2")

###all sample (include y)
union2 = function(stat){
  Z_Lasso_b = com_multi(ZZ="Z",stat=stat,name="pbAllAC")
  Z_Lasso_c = com_multi(ZZ="Z",stat=stat,name="pcAllAC")
  Ztilde_Lasso_b = com_multi(ZZ="Ztilde",stat=stat,name="pbAllAC")
  Ztilde_Lasso_c = com_multi(ZZ="Ztilde",stat=stat,name="pcAllAC")
  
  p=dim(Z_Lasso_b)[1]
  mythred1=mythred2=mythred3=mythred4=list()
  myselect1=myselect2=myselect3=myselect4=list()
  Z=Ztilde=W=matrix(NA,nrow=p, ncol=100)
  beta_select1=beta_select2=beta_select3=beta_select4=matrix(0,nrow=p,ncol=100)
  rownames(beta_select1)=rownames(Z_Lasso_b)
  rownames(beta_select2)=rownames(Z_Lasso_b)
  rownames(beta_select3)=rownames(Z_Lasso_b)
  rownames(beta_select4)=rownames(Z_Lasso_b)
  
  for (i in 1:100){
    Z[,i]=Z_Lasso_b[,i]*Z_Lasso_c[,i]+Ztilde_Lasso_b[,i]*Ztilde_Lasso_c[,i]
    Ztilde[,i]=Z_Lasso_b[,i]*Ztilde_Lasso_c[,i]+Z_Lasso_c[,i]*Ztilde_Lasso_b[,i]
    W[,i]=Z[,i]-Ztilde[,i]
    
    mythred1[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=0)
    myselect1[[i]]=which(W[,i]>=mythred1[[i]])
    beta_select1[myselect1[[i]],i]=1
    
    mythred2[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=1)
    myselect2[[i]]=which(W[,i]>=mythred2[[i]])
    beta_select2[myselect2[[i]],i]=1
    
    mythred3[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=0)
    myselect3[[i]]=which(W[,i]>=mythred3[[i]])
    beta_select3[myselect3[[i]],i]=1
    
    mythred4[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=1)
    myselect4[[i]]=which(W[,i]>=mythred4[[i]])
    beta_select4[myselect4[[i]],i]=1
  }
  
  sumselect1=apply(beta_select1,1,sum)
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  paste1 = paste(names(sumselect1)," (", sumselect1, "%)",sep="")
  sumselect1=cbind(sumselect1,paste1)
  
  sumselect2=apply(beta_select2,1,sum)
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  paste2 = paste(names(sumselect2)," (", sumselect2, "%)",sep="")
  sumselect2=cbind(sumselect2,paste2)
  
  sumselect3=apply(beta_select3,1,sum)
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  paste3 = paste(names(sumselect3)," (", sumselect3, "%)",sep="")
  sumselect3=cbind(sumselect3,paste3)
  
  sumselect4=apply(beta_select4,1,sum)
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  paste4 = paste(names(sumselect4)," (", sumselect4, "%)",sep="")
  sumselect4=cbind(sumselect4,paste4)
  
  File7=sprintf("all/sumbeta1_%s.csv",stat)
  File8=sprintf("all/sumbeta2_%s.csv",stat)
  File9=sprintf("all/sumbeta3_%s.csv",stat)
  File10=sprintf("all/sumbeta4_%s.csv",stat)
  
  write.csv(sumselect1,File7)
  write.csv(sumselect2,File8)
  write.csv(sumselect3,File9)
  write.csv(sumselect4,File10)
}


union2(stat="Lasso")
union2(stat="Lasso1")
union2(stat="Lasso_Order")
union2(stat="RF")
union2(stat="GMUL")
union2(stat="GDS")
union2(stat="CocoLasso")
union2(stat="CocoLasso2")

###analytic (exclude Y)
union3 = function(stat){
  Z_Lasso_b = com_multi(ZZ="Z",stat=stat,name="pbAnalyticEYAC")
  Z_Lasso_c = com_multi(ZZ="Z",stat=stat,name="pcAnalyticEYAC")
  Ztilde_Lasso_b = com_multi(ZZ="Ztilde",stat=stat,name="pbAnalyticEYAC")
  Ztilde_Lasso_c = com_multi(ZZ="Ztilde",stat=stat,name="pcAnalyticEYAC")
  
  p=dim(Z_Lasso_b)[1]
  mythred1=mythred2=mythred3=mythred4=list()
  myselect1=myselect2=myselect3=myselect4=list()
  Z=Ztilde=W=matrix(NA,nrow=p, ncol=100)
  beta_select1=beta_select2=beta_select3=beta_select4=matrix(0,nrow=p,ncol=100)
  rownames(beta_select1)=rownames(Z_Lasso_b)
  rownames(beta_select2)=rownames(Z_Lasso_b)
  rownames(beta_select3)=rownames(Z_Lasso_b)
  rownames(beta_select4)=rownames(Z_Lasso_b)
  
  for (i in 1:100){
    Z[,i]=Z_Lasso_b[,i]*Z_Lasso_c[,i]+Ztilde_Lasso_b[,i]*Ztilde_Lasso_c[,i]
    Ztilde[,i]=Z_Lasso_b[,i]*Ztilde_Lasso_c[,i]+Z_Lasso_c[,i]*Ztilde_Lasso_b[,i]
    W[,i]=Z[,i]-Ztilde[,i]
    
    mythred1[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=0)
    myselect1[[i]]=which(W[,i]>=mythred1[[i]])
    beta_select1[myselect1[[i]],i]=1
    
    mythred2[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=1)
    myselect2[[i]]=which(W[,i]>=mythred2[[i]])
    beta_select2[myselect2[[i]],i]=1
    
    mythred3[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=0)
    myselect3[[i]]=which(W[,i]>=mythred3[[i]])
    beta_select3[myselect3[[i]],i]=1
    
    mythred4[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=1)
    myselect4[[i]]=which(W[,i]>=mythred4[[i]])
    beta_select4[myselect4[[i]],i]=1
  }
  
  sumselect1=apply(beta_select1,1,sum)
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  paste1 = paste(names(sumselect1)," (", sumselect1, "%)",sep="")
  sumselect1=cbind(sumselect1,paste1)
  
  sumselect2=apply(beta_select2,1,sum)
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  paste2 = paste(names(sumselect2)," (", sumselect2, "%)",sep="")
  sumselect2=cbind(sumselect2,paste2)
  
  sumselect3=apply(beta_select3,1,sum)
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  paste3 = paste(names(sumselect3)," (", sumselect3, "%)",sep="")
  sumselect3=cbind(sumselect3,paste3)
  
  sumselect4=apply(beta_select4,1,sum)
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  paste4 = paste(names(sumselect4)," (", sumselect4, "%)",sep="")
  sumselect4=cbind(sumselect4,paste4)
  
  File7=sprintf("analyticEY/sumbeta1_%s.csv",stat)
  File8=sprintf("analyticEY/sumbeta2_%s.csv",stat)
  File9=sprintf("analyticEY/sumbeta3_%s.csv",stat)
  File10=sprintf("analyticEY/sumbeta4_%s.csv",stat)
  
  write.csv(sumselect1,File7)
  write.csv(sumselect2,File8)
  write.csv(sumselect3,File9)
  write.csv(sumselect4,File10)
}


union3(stat="Lasso")
union3(stat="Lasso1")
union3(stat="Lasso_Order")
union3(stat="RF")
union3(stat="GMUL")
union3(stat="GDS")
union3(stat="CocoLasso")
union3(stat="CocoLasso2")

###all (exclude Y)
union4 = function(stat){
  Z_Lasso_b = com_multi(ZZ="Z",stat=stat,name="pbAllEYAC")
  Z_Lasso_c = com_multi(ZZ="Z",stat=stat,name="pcAllEYAC")
  Ztilde_Lasso_b = com_multi(ZZ="Ztilde",stat=stat,name="pbAllEYAC")
  Ztilde_Lasso_c = com_multi(ZZ="Ztilde",stat=stat,name="pcAllEYAC")
  
  p=dim(Z_Lasso_b)[1]
  rep=dim(Z_Lasso_b)[2]
  mythred1=mythred2=mythred3=mythred4=list()
  myselect1=myselect2=myselect3=myselect4=list()
  Z=Ztilde=W=matrix(NA,nrow=p, ncol=rep)
  beta_select1=beta_select2=beta_select3=beta_select4=matrix(0,nrow=p,ncol=rep)
  rownames(beta_select1)=rownames(Z_Lasso_b)
  rownames(beta_select2)=rownames(Z_Lasso_b)
  rownames(beta_select3)=rownames(Z_Lasso_b)
  rownames(beta_select4)=rownames(Z_Lasso_b)
  
  for (i in 1:rep){
    Z[,i]=Z_Lasso_b[,i]*Z_Lasso_c[,i]+Ztilde_Lasso_b[,i]*Ztilde_Lasso_c[,i]
    Ztilde[,i]=Z_Lasso_b[,i]*Ztilde_Lasso_c[,i]+Z_Lasso_c[,i]*Ztilde_Lasso_b[,i]
    W[,i]=Z[,i]-Ztilde[,i]
    
    mythred1[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=0)
    myselect1[[i]]=which(W[,i]>=mythred1[[i]])
    beta_select1[myselect1[[i]],i]=1
    
    mythred2[[i]]=knockoff.threshold(W[,i],fdr=0.1,offset=1)
    myselect2[[i]]=which(W[,i]>=mythred2[[i]])
    beta_select2[myselect2[[i]],i]=1
    
    mythred3[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=0)
    myselect3[[i]]=which(W[,i]>=mythred3[[i]])
    beta_select3[myselect3[[i]],i]=1
    
    mythred4[[i]]=knockoff.threshold(W[,i],fdr=0.2,offset=1)
    myselect4[[i]]=which(W[,i]>=mythred4[[i]])
    beta_select4[myselect4[[i]],i]=1
  }
  
  sumselect1=apply(beta_select1,1,sum)
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  paste1 = paste(names(sumselect1)," (", sumselect1, "%)",sep="")
  sumselect1=cbind(sumselect1,paste1)
  
  sumselect2=apply(beta_select2,1,sum)
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  paste2 = paste(names(sumselect2)," (", sumselect2, "%)",sep="")
  sumselect2=cbind(sumselect2,paste2)
  
  sumselect3=apply(beta_select3,1,sum)
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  paste3 = paste(names(sumselect3)," (", sumselect3, "%)",sep="")
  sumselect3=cbind(sumselect3,paste3)
  
  sumselect4=apply(beta_select4,1,sum)
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  paste4 = paste(names(sumselect4)," (", sumselect4, "%)",sep="")
  sumselect4=cbind(sumselect4,paste4)
  
  File7=sprintf("allEY/sumbeta1_%s.csv",stat)
  File8=sprintf("allEY/sumbeta2_%s.csv",stat)
  File9=sprintf("allEY/sumbeta3_%s.csv",stat)
  File10=sprintf("allEY/sumbeta4_%s.csv",stat)
  
  write.csv(sumselect1,File7)
  write.csv(sumselect2,File8)
  write.csv(sumselect3,File9)
  write.csv(sumselect4,File10)
}


union4(stat="Lasso")
union4(stat="Lasso1")
union4(stat="Lasso_Order")
union4(stat="RF")
union4(stat="GMUL")
union4(stat="GDS")
union4(stat="CocoLasso")
union4(stat="CocoLasso2")


