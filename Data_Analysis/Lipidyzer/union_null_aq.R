setwd("~/Downloads/Result_COCO1/result_newtry2")
###read the rowname for multi
read=read.csv("~/Downloads/Result_COCO1/result_newtry2/final_breast_aq.csv",header=T)
zname = names(read)[-c(1,2)]
p=length(zname)
### combine Z and Ztilde
com_multi = function(ZZ,stat,cancer,name){
  
  addmean=matrix(NA,nrow=5,ncol=p)
  
  Z=NULL
  for (i in 1:20){
    File2 = sprintf("%s_%s_%s_Analytic_%s_%d.csv",ZZ,stat,cancer,name,i)
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

union1 = function(stat){
  Z_Lasso_b = com_multi(ZZ="Z",stat=stat,cancer='breast',name="aq")
  Z_Lasso_c = com_multi(ZZ="Z",stat=stat,cancer='Colore',name="aq")
  Ztilde_Lasso_b = com_multi(ZZ="Ztilde",stat=stat,cancer='breast',name="aq")
  Ztilde_Lasso_c = com_multi(ZZ="Ztilde",stat=stat,cancer='Colore',name="aq")
  
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
  
  File7=sprintf("analytic_aq/sumbeta1_%s.csv",stat)
  File8=sprintf("analytic_aq/sumbeta2_%s.csv",stat)
  File9=sprintf("analytic_aq/sumbeta3_%s.csv",stat)
  File10=sprintf("analytic_aq/sumbeta4_%s.csv",stat)
  
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
union1(stat="CocoLasso2")

