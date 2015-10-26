#setwd("d:\\Users\\xzhou\\Desktop\\")
library(data.table)
data=fread("input_curve.csv",skip=4,header=F)
shell=fread("input_curve.csv",nrows=4,header=F)
g.start=1e-3
v.start=0.5


shell=shell[,seq(2,ncol(shell),by=2),with=F]
shell=data.table(t(shell))
setnames(shell,names(shell),c("Audience","Media","Frequency","cpp"))

g=rep(0,ncol(data))
v=rep(0,ncol(data))
k=rep(0,ncol(data))
max_reach=rep(0,ncol(data))
mape=rep(0,ncol(data))

for (i in seq(1,ncol(data),2)){
  print(i)
  d=data[[i+1]]
  id=data[[i]]
  dataset=data.frame(d=d[!is.na(d)],id=id[!is.na(id)])
  k.start <- max(dataset$d)
  control1 <- nls.control(maxiter= 10000, minFactor= 1e-30, warnOnly= FALSE,tol=1e-05)
  nl.reg <- try(nls(d ~ k * ((1-exp(-g * id))**v),data=dataset,start= list(k=k.start,g=g.start,v=v.start),
                    control= control1),silent=T)
  if (class(nl.reg)=="try-error") break else {
    k[i]=coef(nl.reg)[1]
    g[i]=coef(nl.reg)[2]
    v[i]=coef(nl.reg)[3]
    max_reach[i]=min(max(dataset$d),coef(nl.reg)[1])
    mape.temp=abs(resid(nl.reg))[-1]/dataset$d[-1]
    mape[i]=mean(mape.temp[mape.temp!=Inf])
  }
}

if (class(nl.reg) =="try-error") {
  curve.name=as.vector(as.matrix(shell[i,1:3,with=F]))
  curve.name1=paste(curve.name,collapse=", ")
  print(paste("Error: Curve: ",curve.name1," cannot be fitted due to unusual shape(s). Please contact Support for help.")) 
} else{
  shell$k=k[k!=0]
  shell$g=g[g!=0]
  shell$v=v[v!=0] 
  shell[,inf_point:=log(v)/g]
  shell$inf_point[shell$inf_point<0]=0
  shell$max_reach=max_reach[max_reach!=0]
  shell$mape=mape[mape!=0]
  write.csv(shell,"output.csv",row.names=F)
}
