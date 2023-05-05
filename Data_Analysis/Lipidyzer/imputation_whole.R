args=commandArgs(trailingOnly=TRUE)
iii=as.numeric(args[1])

### different imputation method
impute_method_list=c("whole","wholeEY", "mean", "min")
method=impute_method_list[iii%%length(impute_method_list)+1]
iii=iii%/%length(impute_method_list)

### two type of dataset
data_list=c("p","aq")
data_name=data_list[iii%%length(data_list)+1]
iii=iii%/%length(data_list)

mm=iii+1 #mm should be 1 to 5

library(mice)
setwd("~/Downloads/Dr.Dai's project/2022Spring_project/real/lipid/MI/") #use m=1 for imputation,generate different seed for m=2,3,4,5

impdata=function(data){
  set.seed(mm+111)
  imp=mice(data,m=1)
  data=complete(imp,action="long")
  data=data[,-2]
  return(data)
}

final_p=read.csv("final_p.csv")
final_p=final_p[,-1]

final_aq=read.csv("final_aq.csv")
final_aq=final_aq[,-1]


if (data_name=='p') {
  if (method == "whole") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=final_p[which(final_p$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=final_p[which(final_p$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=final_p[which(final_p$casetype=="Breast ca"),]
    data4=final_p[which(final_p$match %in% bre$match),]
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=final_p[which(final_p$casetype=="Colorecal ca"),]
    data5=final_p[which(final_p$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    
    data2x=impdata(data=data2[,-c(dim(data2)[2],(dim(data2)[2]-1))]) #remove casetype, match
    data3x=impdata(data=data3[,-c(dim(data3)[2],(dim(data3)[2]-1))])
    data4x=impdata(data=data4[,-c(dim(data4)[2],(dim(data4)[2]-1))])
    data5x=impdata(data=data5[,-c(dim(data5)[2],(dim(data5)[2]-1))])
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
    
    data2$.imp=mm
    data3$.imp=mm
    data4$.imp=mm
    data5$.imp=mm
    
    File2=sprintf("data2_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data2,File2,row.names = F)
    
    File3=sprintf("data3_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data3,File3,row.names = F)
    
    File4=sprintf("data4_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data4,File4,row.names = F)
    
    File5=sprintf("data5_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data5,File5,row.names = F)
  }
  
  
  if (method == "AnalyticEY") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=final_p[which(final_p$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=final_p[which(final_p$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=final_p[which(final_p$casetype=="Breast ca"),]
    data4=final_p[which(final_p$match %in% bre$match),]
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=final_p[which(final_p$casetype=="Colorecal ca"),]
    data5=final_p[which(final_p$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    
    data2x=impdata(data=data2[,-c(1,dim(data2)[2],(dim(data2)[2]-1))]) #remove Y, casetype, match
    data3x=impdata(data=data3[,-c(1,dim(data3)[2],(dim(data3)[2]-1))])
    data4x=impdata(data=data4[,-c(1,dim(data4)[2],(dim(data4)[2]-1))])
    data5x=impdata(data=data5[,-c(1,dim(data5)[2],(dim(data5)[2]-1))])
    
    data2x$Y=data2$Y #give the original value of Y to imputation data
    data3x$Y=data3$Y
    data4x$Y=data4$Y
    data5x$Y=data5$Y
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
    
    data2$.imp=mm
    data3$.imp=mm
    data4$.imp=mm
    data5$.imp=mm
    
    File2=sprintf("data2_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data2,File2,row.names = F)
    
    File3=sprintf("data3_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data3,File3,row.names = F)
    
    File4=sprintf("data4_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data4,File4,row.names = F)
    
    File5=sprintf("data5_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data5,File5,row.names = F)
  }
  
  
  if (method =="All") {
    #use all data to do the imputation
    datax=impdata(data=final_p[,-c(dim(final_p)[2],(dim(final_p)[2]-1))]) #remove casetype, match
    datax$casetype=rep(final_p$casetype,1)
    datax$match=rep(final_p$match,1)
    
    datax$.imp=mm
    File=sprintf("datax_%d_%s_%s.csv",mm,data_name,method)
    write.csv(datax,File,row.names = F)
  }
  
  
  if (method =="AllEY") {
    #use all data to do the imputation
    datax=impdata(data=final_p[,-c(1,dim(final_p)[2],(dim(final_p)[2]-1))]) #remove Y, casetype, match
    datax$Y=rep(final_p$Y,1)
    
    datax$casetype=rep(final_p$casetype,1)
    datax$match=rep(final_p$match,1)
    datax$.imp=mm
    File=sprintf("datax_%d_%s_%s.csv",mm,data_name,method)
    write.csv(datax,File,row.names = F)
    
  }
}

if (data_name=='aq') {
  if (method == "Analytic") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=final_aq[which(final_aq$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=final_aq[which(final_aq$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=final_aq[which(final_aq$casetype=="Breast ca"),]
    data4=final_aq[which(final_aq$match %in% bre$match),]
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=final_aq[which(final_aq$casetype=="Colorecal ca"),]
    data5=final_aq[which(final_aq$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    
    data2x=impdata(data=data2[,-c(dim(data2)[2],(dim(data2)[2]-1))]) #remove casetype, match
    data3x=impdata(data=data3[,-c(dim(data3)[2],(dim(data3)[2]-1))])
    data4x=impdata(data=data4[,-c(dim(data4)[2],(dim(data4)[2]-1))])
    data5x=impdata(data=data5[,-c(dim(data5)[2],(dim(data5)[2]-1))])
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
    
    data2$.imp=mm
    data3$.imp=mm
    data4$.imp=mm
    data5$.imp=mm
    
    File2=sprintf("data2_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data2,File2,row.names = F)
    
    File3=sprintf("data3_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data3,File3,row.names = F)
    
    File4=sprintf("data4_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data4,File4,row.names = F)
    
    File5=sprintf("data5_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data5,File5,row.names = F)
  }
  
  
  if (method == "AnalyticEY") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=final_aq[which(final_aq$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=final_aq[which(final_aq$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=final_aq[which(final_aq$casetype=="Breast ca"),]
    data4=final_aq[which(final_aq$match %in% bre$match),]
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=final_aq[which(final_aq$casetype=="Colorecal ca"),]
    data5=final_aq[which(final_aq$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    
    data2x=impdata(data=data2[,-c(1,dim(data2)[2],(dim(data2)[2]-1))]) #remove Y, casetype, match
    data3x=impdata(data=data3[,-c(1,dim(data3)[2],(dim(data3)[2]-1))])
    data4x=impdata(data=data4[,-c(1,dim(data4)[2],(dim(data4)[2]-1))])
    data5x=impdata(data=data5[,-c(1,dim(data5)[2],(dim(data5)[2]-1))])
    
    data2x$Y=data2$Y #give the original value of Y to imputation data
    data3x$Y=data3$Y
    data4x$Y=data4$Y
    data5x$Y=data5$Y
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
    
    data2$.imp=mm
    data3$.imp=mm
    data4$.imp=mm
    data5$.imp=mm
    
    File2=sprintf("data2_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data2,File2,row.names = F)
    
    File3=sprintf("data3_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data3,File3,row.names = F)
    
    File4=sprintf("data4_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data4,File4,row.names = F)
    
    File5=sprintf("data5_%d_%s_%s.csv",mm,data_name,method)
    write.csv(data5,File5,row.names = F)
  }
  
  
  if (method =="All") {
    #use all data to do the imputation
    datax=impdata(data=final_aq[,-c(dim(final_aq)[2],(dim(final_aq)[2]-1))]) #remove casetype, match
    datax$casetype=rep(final_aq$casetype,1)
    datax$match=rep(final_aq$match,1)
    
    datax$.imp=mm
    File=sprintf("datax_%d_%s_%s.csv",mm,data_name,method)
    write.csv(datax,File,row.names = F)
  }
  
  
  if (method =="AllEY") {
    #use all data to do the imputation
    datax=impdata(data=final_aq[,-c(1,dim(final_aq)[2],(dim(final_aq)[2]-1))]) #remove Y, casetype, match
    datax$Y=rep(final_aq$Y,1)
    
    datax$casetype=rep(final_aq$casetype,1)
    datax$match=rep(final_aq$match,1)
    datax$.imp=mm
    File=sprintf("datax_%d_%s_%s.csv",mm,data_name,method)
    write.csv(datax,File,row.names = F)
  }
}


