setwd("/Users/runqiuwang/Downloads/Dr.Dai's project/2022Spring_project/real/")
source("Rcode/estimation.r")
final_p=read.csv("LCMS/Multi/final_p.csv")
final_p=final_p[,-1]

final_aq=read.csv("LCMS/Multi/final_aq.csv")
final_aq=final_aq[,-1]

MI = function(method,data,name) {
  if (method == "Analytic") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=data[which(data$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=data[which(data$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=data[which(data$casetype=="Breast ca"),]
    data4=data[which(data$match %in% bre$match),]
    data4=data4[order(data4$match),]
    #cc=tapply(data4$match,data4$match,length)
    data4=data4[which(data4$match!=67),] 
    #delete match=67 since it only have case group
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=data[which(data$casetype=="Colorecal ca"),]
    data5=data[which(data$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    data2x=impdata(data=data2[,-c(1,dim(data2)[2],(dim(data2)[2]-1))],method="Multi",rr=log(2)) #delete ID, casetype, match
    data3x=impdata(data=data3[,-c(1,dim(data3)[2],(dim(data3)[2]-1))],method="Multi",rr=log(2))
    data4x=impdata(data=data4[,-c(1,dim(data4)[2],(dim(data4)[2]-1))],method="Multi",rr=log(2))
    data5x=impdata(data=data5[,-c(1,dim(data5)[2],(dim(data5)[2]-1))],method="Multi",rr=log(2))
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
  }
  
  if (method =="All") {
    #use all data to do the imputation
    datax=impdata(data=data[,-c(dim(data)[2],(dim(data)[2]-1))],method="Multi",rr=log(2))
    datax$casetype=rep(data$casetype,5)
    datax$match=rep(data$match,5)
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=datax[which(datax$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=datax[which(datax$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=datax[which(datax$casetype=="Breast ca"),]
    data4=datax[which(datax$match %in% bre$match),]
    data4=data4[order(data4$match),]
    #cc=tapply(data4$match,data4$match,length)
    data4=data4[which(data4$match!=67),] 
    #delete match=67 since it only have case group
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=datax[which(datax$casetype=="Colorecal ca"),]
    data5=datax[which(datax$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    data2=data2[,-c(dim(data2)[2],(dim(data2)[2]-1))]
    data3=data3[,-c(dim(data3)[2],(dim(data3)[2]-1))]
    data4=data4[,-c(dim(data4)[2],(dim(data4)[2]-1))]
    data5=data5[,-c(dim(data5)[2],(dim(data5)[2]-1))]
  }
  
  if (method == "AnalyticEY") {
    #categorize into group then do the imputation
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=data[which(data$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=data[which(data$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=data[which(data$casetype=="Breast ca"),]
    data4=data[which(data$match %in% bre$match),]
    data4=data4[order(data4$match),]
    #cc=tapply(data4$match,data4$match,length)
    data4=data4[which(data4$match!=67),] 
    #delete match=67 since it only have case group
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=data[which(data$casetype=="Colorecal ca"),]
    data5=data[which(data$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    data2x=impdata(data=data2[,-c(1,dim(data2)[2],(dim(data2)[2]-1))],method="Multi",rr=log(2)) #exclude Y during the imputation
    data3x=impdata(data=data3[,-c(1,dim(data3)[2],(dim(data3)[2]-1))],method="Multi",rr=log(2))
    data4x=impdata(data=data4[,-c(1,dim(data4)[2],(dim(data4)[2]-1))],method="Multi",rr=log(2))
    data5x=impdata(data=data5[,-c(1,dim(data5)[2],(dim(data5)[2]-1))],method="Multi",rr=log(2))
    
    data2x$Y=rep(data2$Y,5)
    data3x$Y=rep(data3$Y,5)
    data4x$Y=rep(data4$Y,5)
    data5x$Y=rep(data5$Y,5)
    
    data2=data2x
    data3=data3x
    data4=data4x
    data5=data5x
  }
  
  
  if (method =="AllEY") {
    #use all data to do the imputation
    datax=impdata(data=data[,-c(1,dim(data)[2],(dim(data)[2]-1))],method="Multi",rr=log(2)) #exlude Y during the imputation
    datax$Y=rep(data$Y,5)
    datax$casetype=rep(data$casetype,5)
    datax$match=rep(data$match,5)
    
    
    ###separate different cancer: using common controls
    ### breast cancer vs control
    data2=datax[which(datax$casetype!="Colorecal ca"),]
    ### Colorecal cancer vs control
    data3=datax[which(datax$casetype!="Breast ca"),]
    
    ###separate different cancer with matched group
    ### breast cancer 
    bre=datax[which(datax$casetype=="Breast ca"),]
    data4=datax[which(datax$match %in% bre$match),]
    data4=data4[order(data4$match),]
    #cc=tapply(data4$match,data4$match,length)
    data4=data4[which(data4$match!=67),] 
    #delete match=67 since it only have case group
    data4=data4[order(data4$match),]
    
    ### colore cancer
    col=datax[which(datax$casetype=="Colorecal ca"),]
    data5=datax[which(datax$match %in% col$match),]
    #cc=tapply(data5$match,data5$match,length)
    data5=data5[order(data5$match),]
    
    data2=data2[,-c(dim(data2)[2],(dim(data2)[2]-1))]
    data3=data3[,-c(dim(data3)[2],(dim(data3)[2]-1))]
    data4=data4[,-c(dim(data4)[2],(dim(data4)[2]-1))]
    data5=data5[,-c(dim(data5)[2],(dim(data5)[2]-1))]
  }
  
  
  File2=sprintf("LCMS/Multi/data_breast_%s_%s.csv",method,name)
  write.csv(data2, File2)
  File3=sprintf("LCMS/Multi/data_Colore_%s_%s.csv",method,name)
  write.csv(data3, File3)
  File4=sprintf("LCMS/Multi/data_breast_match_%s_%s.csv",method,name)
  write.csv(data4, File4)
  File5=sprintf("LCMS/Multi/data_Colore_match_%s_%s.csv",method,name)
  write.csv(data5, File5)
}

#MI(method="All",data=final_p,name='p')
#MI(method="AnalyticEY",data=final_p,name='p')
#MI(method="AllEY",data=final_p,name='p')

#MI(method="All",data=final_aq,name='aq')
#MI(method="AnalyticEY",data=final_aq,name='aq')
#MI(method="AllEY",data=final_aq,name='aq')

whole = function(way,data,name) {
  data=data[,-c(dim(data)[2],(dim(data)[2]-1))]
  
  if (way=='mean') {
    data1=impdata(data=data,method="Single",rr=log(2))
  }
  
  if (way=='min') {
    data1=impdata(data=data,method="Min",rr=log(2))
  }
  
  if (way=='MI') {
    data1=impdata(data=data,method="Multi",rr=log(2))
  }
  
  if (way=='MIEY') {
    data1=impdata(data=data[,-1],method="Multi",rr=log(2))
    data1$Y=rep(data$Y,5)
  }
  
  File6=sprintf("LCMS/Multi/data_whole_%s_%s.csv",way,name)
  write.csv(data1, File6)
}
whole(way='mean',data=final_aq,name='aq')
whole(way='min',data=final_aq,name='aq')
whole(way='MI',data=final_aq,name='aq')
whole(way='MIEY',data=final_aq,name='aq')

whole(way='mean',data=final_p,name='p')
whole(way='min',data=final_p,name='p')
whole(way='MI',data=final_p,name='p')
whole(way='MIEY',data=final_p,name='p')



