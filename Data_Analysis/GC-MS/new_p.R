setwd("/Volumes/Runqiu Wang/Dr.Dai's project/2022Spring_project/real/GCMS/Multi/")

### combine Z and Ztilde
whole = function(ZZ,stat,cancer,name){
  ###read the rowname for multi ### for composition
  read=read.csv("final_breast_GCMS.csv",header=T)
  zname = names(read)[-c(1,2,3)]
  p=length(zname)
  
  Z=NULL
  for (m in 1:20){
    i=(m-1)*5
    File2 = sprintf("result3/%s_%s_%s_%d.csv",ZZ,stat,cancer,i)
    if (file.exists(File2)==T){
      add = read.csv(File2, header=T)      
      add = add[,-1]
      Z = rbind(Z,add)
    }
  }
  
  ### add names to the Z
  colnames(Z)=zname
  Z
}
#xx=whole(ZZ="Z",stat="Lasso_Order",cancer="pbAC",name="p")
seq_thres=function(fdr,offset,pvalue,W){
  p=length(pvalue)
  #ratio.calculate = function(k) {(offset + sum(pvalue[1:k] >= 0.5))/max(1,sum(re_p_value[1:k] <= 0.5))}
  #ratio.calculate1 = function(k) {(offset + sum(pvalue[1:k] >= 0.5))/max(1,sum(re_p_value1[1:k] <= 0.5))}
  #ratio.calculate2 = function(k) {(offset + sum(re_p_value1[1:k] > 0.5))/max(1,sum(re_p_value1[1:k] <= 0.5))}
  ratio.calculate3 = function(k) {(offset + sum(pvalue[1:k] > 0.5))/max(1,sum(pvalue[1:k] <= 0.5))}
  
  ratio=rep(NA,p)
  set=c()
  for (k in 1:p){
    ratio[k]=ratio.calculate3(k)
  }
  ok = which(ratio <= fdr)
  if (length(ok)>0) {
    KK = max(ok)
    set=order(-W)[intersect(1:KK,which(pvalue<=0.5))]
  }
  else {set=NULL} 
  return(set)
}


###analytic sample (include y)

new <- function(stat,name,cancer){
  ZmatW = whole(ZZ="Z",stat=stat,cancer=cancer,name=name)
  ZtildematW = whole(ZZ="Ztilde",stat=stat,cancer=cancer,name=name)
  p=dim(ZmatW)[2]
  rep=dim(ZmatW)[1]
  beta_select1=beta_select2=beta_select3=beta_select4=matrix(0,nrow=p,ncol=rep/5)
  rownames(beta_select1)=rownames(beta_select2)=rownames(beta_select3)=rownames(beta_select4)=colnames(ZmatW)
  
  q=c(0.1,0.2)
  offset=c(0,1)
  m=5
  for (i in 1:(rep/m)){
    Zmat=ZmatW[(m*i-m+1):(m*i),]
    Ztildemat=ZtildematW[(m*i-m+1):(m*i),]
    #### calculate for the result
    PP = ifelse(Ztildemat>Zmat,1,0)
    PP1 = ifelse(Ztildemat==Zmat,1/2,0)
    PP2 = PP+PP1
    PP3 = PP+PP1*2
    
    #p_value = apply(PP,2,function(X){(sum(X)+1)/(length(X)+1)})
    # p_value1 = apply(PP2,2,function(X){(sum(X)+1)/(length(X)+1)})
    p_value2 = apply(PP3,2,function(X){(sum(X)+1)/(length(X)+1)})
    
    W=apply(abs(Zmat - Ztildemat),2,max)
    
    
    #re_p_value=p_value[order(-W)]
    # re_p_value1=p_value1[order(-W)]
    re_p_value2=p_value2[order(-W)]
    
    myselectall=list()
    count=0
    for (offset0 in offset){
      for (q0 in q) {
        count=count+1
        temp_result=seq_thres(fdr=q0,offset=offset0,pvalue=re_p_value2,W=W)
        if (length(temp_result)>0)
        {myselectall[[count]]=temp_result}
        if (length(temp_result)==0)
        {myselectall[[count]]=NA}  ### if didn't select, set as NA
      }
    }
    
    beta_select1[myselectall[[1]],i]=1
    beta_select2[myselectall[[2]],i]=1
    beta_select3[myselectall[[3]],i]=1
    beta_select4[myselectall[[4]],i]=1
    
  }
  
  
  sumselect1=apply(beta_select1,1,sum)/(rep/m)*100
  sumselect1=sumselect1[order(sumselect1,decreasing = T)]
  paste1 = paste(names(sumselect1)," (", sumselect1, "%)",sep="")
  sumselect1=cbind(sumselect1,paste1)
  
  sumselect2=apply(beta_select2,1,sum)/(rep/m)*100
  sumselect2=sumselect2[order(sumselect2,decreasing = T)]
  paste2 = paste(names(sumselect2)," (", sumselect2, "%)",sep="")
  sumselect2=cbind(sumselect2,paste2)
  
  sumselect3=apply(beta_select3,1,sum)/(rep/m)*100
  sumselect3=sumselect3[order(sumselect3,decreasing = T)]
  paste3 = paste(names(sumselect3)," (", sumselect3, "%)",sep="")
  sumselect3=cbind(sumselect3,paste3)
  
  sumselect4=apply(beta_select4,1,sum)/(rep/m)*100
  sumselect4=sumselect4[order(sumselect4,decreasing = T)]
  paste4 = paste(names(sumselect4)," (", sumselect4, "%)",sep="")
  sumselect4=cbind(sumselect4,paste4)
  
  #File7=sprintf("sumbeta1_%s_%s_%s.csv",stat,name,cancer)
  #File8=sprintf("sumbeta2_%s_%s_%s.csv",stat,name,cancer)
  # File9=sprintf("sumbeta3_%s_%s_%s.csv",stat,name,cancer)
  # File10=sprintf("sumbeta4_%s_%s_%s.csv",stat,name,cancer)
  
  #write.csv(sumselect1,File7)
  #write.csv(sumselect2,File8)
  #write.csv(sumselect3,File9)
  # write.csv(sumselect4,File10)
  return(list(sumselect1=sumselect1, sumselect2=sumselect2, sumselect3=sumselect3, sumselect4=sumselect4))
}

#Lasso = new(stat="Lasso",name="p",cancer="pbAC")$sumselect1[,2]
#Lasso_Order = new(stat="Lasso_Order",name="p",cancer="pbAC")$sumselect1[,2]
#GMUL = new(stat="GMUL",name="p",cancer="pbAC")$sumselect1[,2]
#GDS = new(stat="GDS",name="p",cancer="pbAC")$sumselect1[,2]
#CocoLasso3 = new(stat="CocoLasso3",name="p",cancer="pbAC")$sumselect1[,2]
#RF = new(stat="RF",name="p",cancer="pbAC")$sumselect1[,2]
#result_breast_p = cbind(Lasso,Lasso_Order,GMUL,GDS,CocoLasso3,RF)
#colnames(result_breast_p) = c("Lasso","Lasso_Order","GMUL","GDS", "CocoLasso3","RF")
#write.csv(result_breast_p,"result3/result_breast_p.csv")

#Lasso = new(stat="Lasso",name="p",cancer="pcAC")$sumselect1[,2]
#Lasso_Order = new(stat="Lasso_Order",name="p",cancer="pcAC")$sumselect1[,2]
#GMUL = new(stat="GMUL",name="p",cancer="pcAC")$sumselect1[,2]
#GDS = new(stat="GDS",name="p",cancer="pcAC")$sumselect1[,2]
#CocoLasso3 = new(stat="CocoLasso3",name="p",cancer="pcAC")$sumselect1[,2]
#RF = new(stat="RF",name="p",cancer="pcAC")$sumselect1[,2]
#result_Colore_p = cbind(Lasso,Lasso_Order,GMUL,GDS,CocoLasso3,RF)
#colnames(result_Colore_p) = c("Lasso","Lasso_Order","GMUL","GDS", "CocoLasso3","RF")
#write.csv(result_Colore_p,"result3/result_Colore_p.csv")

Lasso = new(stat="Lasso",name="p",cancer="pb")$sumselect1[,2]
Lasso_Order = new(stat="Lasso_Order",name="p",cancer="pb")$sumselect1[,2]
GMUL = new(stat="GMUL",name="p",cancer="pb")$sumselect1[,2]
GDS = new(stat="GDS",name="p",cancer="pb")$sumselect1[,2]
CocoLasso3 = new(stat="CocoLasso3",name="p",cancer="pb")$sumselect1[,2]
RF = new(stat="RF",name="p",cancer="pb")$sumselect1[,2]
result_breast_p = cbind(Lasso,Lasso_Order,GMUL,GDS,CocoLasso3,RF)
colnames(result_breast_p) = c("Lasso","Lasso_Order","GMUL","GDS", "CocoLasso3","RF")
write.csv(result_breast_p,"result_breast_match.csv")

Lasso = new(stat="Lasso",name="p",cancer="pc")$sumselect1[,2]
Lasso_Order = new(stat="Lasso_Order",name="p",cancer="pc")$sumselect1[,2]
GMUL = new(stat="GMUL",name="p",cancer="pc")$sumselect1[,2]
GDS = new(stat="GDS",name="p",cancer="pc")$sumselect1[,2]
CocoLasso3 = new(stat="CocoLasso3",name="p",cancer="pc")$sumselect1[,2]
RF = new(stat="RF",name="p",cancer="pc")$sumselect1[,2]
result_Colore_p = cbind(Lasso,Lasso_Order,GMUL,GDS,CocoLasso3,RF)
colnames(result_Colore_p) = c("Lasso","Lasso_Order","GMUL","GDS", "CocoLasso3","RF")
write.csv(result_Colore_p,"result_Colore_match.csv")
