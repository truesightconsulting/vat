#######################################################################################
# just one optm with true optm method
# Load in setup files
# does two opt's, one without cstr and one with
#######################################################################################
setwd("C:\\Users\\876036-mzhou\\Desktop\\opt_0_10_06_2015_06_54_04_05183\\")
start=Sys.time()
#######################################################################################
# OPTM w/o constraint
#######################################################################################
print("Optimization w/o constraint")
print("Optimization Initialization")
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(doSNOW)))
suppressMessages(suppressWarnings(library(nloptr)))
max.level=1e+10 # max constrain if missing
no.cluster=6 # no. of cores to use for large no. of curves

ex.curve=fread("modelinput_curves.csv")
ex.budget=fread("input_budget.csv")
ex.min=fread("input_constrain_min.csv")
ex.max=fread("input_constrain_max.csv")
ex.aud=fread("input_audience.csv")
ex.freq=fread("input_freq.csv")
ex.media=fread("input_media.csv")
ex.cpp=fread("modelinput_cpp.csv")
ex.shell=fread("modelinput_shell.csv")
ex.shell.rollup=fread("modelinput_shell_rollup.csv")
ex.dupe=fread("input_dupefactor.csv")

#######################################################################################
# OPTM initialization
#######################################################################################
budget=ex.budget$Budget
max.level=budget

# merge all the info with curve
curve=merge(ex.curve,ex.cpp[,c("aud_num","media_num","cpp"),with=F],
            by=c("aud_num","media_num"),all.x=T)
curve=merge(curve,ex.media[,c("media_num","flag_media"),with=F],
            by=c("media_num"),all.x=T)
curve=merge(curve,ex.freq[,c("freq_num","flag_freq"),with=F],
            by=c("freq_num"),all.x=T)
curve=merge(curve,ex.aud[,c("aud_num","flag_aud"),with=F],
            by=c("aud_num"),all.x=T)

# filter out non-selected curves 
flag=curve$flag_aud+curve$flag_media+curve$flag_freq
curve=curve[flag==3,]

# create some varaibles for optm
curve$r_grs=curve$r_net=curve$sp_next=curve$r_grs_next=rep(0,nrow(curve))
curve$g1=curve$g/curve$cpp
curve$inf_point1=curve$inf_point/curve$cpp

# filter out non-selected shell
index=ex.shell$shell_num %in% curve$shell_num
shell=ex.shell[index,]

# merge shell with min and max
ex.max$max_spend=rep(max.level,nrow(ex.max))
ex.min$min_spend=rep(0,nrow(ex.max))

shell=merge(shell,ex.max[,c("shell_num","max_spend"),with=F],
            by="shell_num",all.x=T)
shell=merge(shell,ex.min[,c("shell_num","min_spend"),with=F],
            by="shell_num",all.x=T)

# calc max spend from max reach
curve$max_reach[is.na(curve$max_reach)]=curve$k[is.na(curve$max_reach)]
max_reach_sp=log((1-(curve$max_reach/curve$k)^curve$v))/(-curve$g1)
max_reach_sp[max_reach_sp==Inf]=max.level
curve$max_reach_sp=max_reach_sp
shell1=merge(shell[,c("shell_num","max_spend","min_spend"),with=F],curve[,c("shell_num","max_reach_sp","inf_point1"),with=F],by="shell_num",all.x=T)
max_spend=pmin(shell1$max_spend,shell1$max_reach_sp)
shell2=data.table(shell_num=shell1$shell_num,max_spend=max_spend)
min_spend=pmax(shell1$min_spend,shell1$inf_point1)
shell2=data.table(shell2,min_spend=min_spend,inf_point1=shell1$inf_point1,min_spend_t=shell1$min_spend)
shell=merge(shell[,!c("max_spend","min_spend"),with=F],shell2,by="shell_num",all.x=T)
# Ignore the turn point spend for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
shell$min_spend=shell$min_spend_t

# merge curve with shell 
shell$sp_current=shell$min_spend
curve=merge(curve, shell[,c("shell_num","max_spend","min_spend","sp_current"),with=F],
            by="shell_num",all.x=T)


# missing curve check
if (nrow(curve)==0){
  stop("Error: There is no curve existing under selected dimension mix.")
}else{
  index=curve$k==0 | curve$g==0 | curve$v==0
  if (sum(index)!=0) {
    #curve[index,c("Audience","Frequency","Media"),with=F]
    print("Warning: There is curve(s) missing. The result below drops those missing curve(s).")
    curve=curve[!index]
  }
  
  # single curve simulation
  if (nrow(curve)==1){
    curve$sp_current=budget
    curve$g1=curve$g/curve$cpp
    curve$r_grs=curve$k*((1-exp(-curve$g1*curve$sp_current))^curve$v)
    col1=c("Media","Budget","Allocation","Gross reach","Net reach","Total 30s GRPs")
    curve$grp=curve$sp_current/curve$cpp
    output1=curve[,c("sp_current","r_grs","grp"),with=F]
    output1=output1[1,]
    output2=data.table(Media=c("Budget","Gross reach","Total 30s GRPs"),v1=as.vector(as.matrix(output1)))
    setnames(output2,"v1",curve[["Media"]])
    write.table(output2,"output_alloc_net_net.csv",sep=",")
  }else{
    # check some posibble errors
    if (sum(curve$sp_current)>budget){
      stop("Error: Total current spend is greater than total planed budget.")
    }else if (sum(curve$sp_current>curve$max_spend) !=0) {
      stop("Error: One or more media channel's spend is greater than its maximum constrain.")
    }else{
      #######################################################################################
      # OPTM initialization
      #######################################################################################
      print("Optimization")
      x0=budget*curve$k/sum(curve$k)
      fn=function(x){
        #x=curve$sp_current
        curve=curve[,r_grs:=k*((1-exp(-g1*x)))**v]
        result.random.shell=curve[,c("media_num","r_grs"),with=F]
        # compute random scenario
        for.ran=result.random.shell[order(media_num)]
        n=nrow(for.ran)
        no.media=1:n
        r.dupe.box=rep(0,n)
        r.dupe.box.media=matrix(0,nr=n,nc=n,dimnames=list(for.ran$media_num,1:n))
        r.dupe.box.media[,1]=for.ran$r_grs
        for (j in 2:n){
          #j=2
          combo=combn(no.media,j)
          r.combo=matrix(for.ran$r_grs[combo],nr=j)
          r.combo.prod=apply(r.combo,2,prod)
          r.combo.tran=data.table(t(rbind(combo,r.combo.prod)))
          r.combo.melt=data.table(melt(r.combo.tran,id="r.combo.prod"))
          r.combo.meida=r.combo.melt[,list(r.prod=sum(r.combo.prod)),by=c("value")]
          r.dupe.box.media[,j]=r.combo.meida$r.prod
          r.dupe.box[j]=sum(r.combo.prod)
        }
        -(sum(for.ran$r_grs)+sum((-1)^(1:n+1)*r.dupe.box))
      }
      hin<- function(x) {
        h <- numeric(1)
        h[1] <- budget -sum(x)
        h[2] <- sum(x)-budget
        return(h)
      }
      a=curve$min_spend
      b=curve$max_spend
      b[a==b]=b[a==b]+0.01
      x0=budget*curve$k/sum(curve$k)
      x0[x0>b]=b[x0>b]
      x0[x0<a]=a[x0<a]
      optm <- cobyla(x0, fn, hin = hin, lower=a, upper=b,
                     nl.info = F, control = list(maxeval = 10000,ftol_abs=1e-6,xtol_rel=1e-6))
      curve$sp_current=optm$par
      
      #######################################################################################
      # OPTM Rollup Random
      #######################################################################################
      result.random.shell=curve[,c("media_num","shell_num","Media","max_spend","min_spend","sp_current","r_grs"),
                                with=F]
      # compute random scenario
      for.ran=result.random.shell[order(media_num)]
      no.media=1:nrow(curve)
      
      r.dupe.box=rep(0,nrow(curve))
      r.dupe.box.media=matrix(0,nr=nrow(curve),
                              nc=nrow(curve),
                              dimnames=list(for.ran$media_num,
                                            1:nrow(curve)))
      r.dupe.box.media[,1]=for.ran$r_grs
      
      for (j in 2:nrow(curve)){
        #j=2
        combo=combn(no.media,j)
        r.combo=matrix(for.ran$r_grs[combo],nr=j)
        r.combo.prod=apply(r.combo,2,prod)
        r.combo.tran=data.table(t(rbind(combo,r.combo.prod)))
        r.combo.melt=data.table(melt(r.combo.tran,id="r.combo.prod"))
        r.combo.meida=r.combo.melt[,list(r.prod=sum(r.combo.prod)),by=c("value")]
        r.dupe.box.media[,j]=r.combo.meida$r.prod
        r.dupe.box[j]=sum(r.combo.prod)
      }
      r_net_ran_total=sum(for.ran$r_grs)+sum((-1)^(1:nrow(curve)+1)*r.dupe.box)
      for.ran$r_net=apply(sweep(r.dupe.box.media,MARGIN=2,(-1)^(1:nrow(curve)+1),"*"),
                          1,sum)
      for.ran$r_dupe=for.ran$r_grs-for.ran$r_net
      r_dupe_ran_total=r_net_ran_total-sum(for.ran$r_net) #zx
      
      #######################################################################################
      # OPTM Rollup Max
      #######################################################################################
      result.max.shell=curve[,c("media_num","shell_num","Media","max_spend","min_spend","sp_current","r_grs"),
                             with=F]
      
      # compute max scenario
      for.max=result.max.shell[order(-r_grs)]
      for.max$r_net=rep(0,nrow(for.max))
      #r_grs_lag=c(for.max$r_grs[2:length(for.max$r_grs)],0)
      for.max$r_net[1]=for.max$r_grs[1]-for.max$r_grs[2]
      r_net_max_total=max(for.max$r_grs)
      r_dupe_max_total=max(for.max$r_grs[for.max$r_grs!=max(for.max$r_grs)])
      for.max=for.max[order(media_num)]
      for.max$r_dupe=for.max$r_grs-for.max$r_net #zx
      
      # save unconstraint results
      r_dupe_ran_total1=r_dupe_ran_total
      r_net_ran_total1=r_net_ran_total
      r_net_max_total1=r_net_max_total
      r_dupe_max_total1=r_dupe_max_total
    }# error check
  }#single curve check
}# zero curve check

#######################################################################################
# OPTM w/ constraint
#######################################################################################
print("OPTM w/ constraint")
print("Optimization Initialization")
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(doSNOW)))
suppressMessages(suppressWarnings(library(nloptr)))
max.level=1e+10 # max constrain if missing
no.cluster=6 # no. of cores to use for large no. of curves

ex.curve=fread("modelinput_curves.csv")
ex.budget=fread("input_budget.csv")
ex.min=fread("input_constrain_min.csv")
ex.max=fread("input_constrain_max.csv")
ex.aud=fread("input_audience.csv")
ex.freq=fread("input_freq.csv")
ex.media=fread("input_media.csv")
ex.cpp=fread("modelinput_cpp.csv")
ex.shell=fread("modelinput_shell.csv")
ex.shell.rollup=fread("modelinput_shell_rollup.csv")
ex.dupe=fread("input_dupefactor.csv")

#######################################################################################
# OPTM initialization
#######################################################################################
budget=ex.budget$Budget

# merge all the info with curve
curve=merge(ex.curve,ex.cpp[,c("aud_num","media_num","cpp"),with=F],
            by=c("aud_num","media_num"),all.x=T)
curve=merge(curve,ex.media[,c("media_num","flag_media"),with=F],
            by=c("media_num"),all.x=T)
curve=merge(curve,ex.freq[,c("freq_num","flag_freq"),with=F],
            by=c("freq_num"),all.x=T)
curve=merge(curve,ex.aud[,c("aud_num","flag_aud"),with=F],
            by=c("aud_num"),all.x=T)

# filter out non-selected curves 
flag=curve$flag_aud+curve$flag_media+curve$flag_freq
curve=curve[flag==3,]

# create some varaibles for optm
curve$r_grs=curve$r_net=curve$sp_next=curve$r_grs_next=rep(0,nrow(curve))
curve$g1=curve$g/curve$cpp
curve$inf_point1=curve$inf_point/curve$cpp

# filter out non-selected shell
index=ex.shell$shell_num %in% curve$shell_num
shell=ex.shell[index,]

# merge shell with min and max
ex.max$max_spend[is.na(ex.max$max_spend)]=max.level
ex.min$min_spend[is.na(ex.min$min_spend)]=0

shell=merge(shell,ex.max[,c("shell_num","max_spend"),with=F],
            by="shell_num",all.x=T)
shell=merge(shell,ex.min[,c("shell_num","min_spend"),with=F],
            by="shell_num",all.x=T)

# calc max spend from max reach
curve$max_reach[is.na(curve$max_reach)]=curve$k[is.na(curve$max_reach)]
max_reach_sp=log((1-(curve$max_reach/curve$k)^curve$v))/(-curve$g1)
max_reach_sp[max_reach_sp==Inf]=max.level
curve$max_reach_sp=max_reach_sp
shell1=merge(shell[,c("shell_num","max_spend","min_spend"),with=F],curve[,c("shell_num","max_reach_sp","inf_point1"),with=F],by="shell_num",all.x=T)
max_spend=pmin(shell1$max_spend,shell1$max_reach_sp)
shell2=data.table(shell_num=shell1$shell_num,max_spend=max_spend)
min_spend=pmax(shell1$min_spend,shell1$inf_point1)
shell2=data.table(shell2,min_spend=min_spend,max_reach_sp=shell1$max_reach_sp,min_spend_t=shell1$min_spend)
shell=merge(shell[,!c("max_spend","min_spend"),with=F],shell2,by="shell_num",all.x=T)
# Ignore the turn point spend for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
shell$min_spend=shell$min_spend_t

# merge curve with shell 
shell$sp_current=shell$min_spend
curve=merge(curve, shell[,c("shell_num","max_spend","min_spend","sp_current"),with=F],
            by="shell_num",all.x=T)
if (sum(curve$max_reach_sp<curve$min_spend) !=0) {
  print("Warning: One or more media channel's minimum constrain is greater than maximum reach spend. The minimum spend will be replaced by the maximum reach spend")
}
curve$min_spend=pmin(curve$max_reach_sp,curve$min_spend)
curve$sp_current=curve$min_spend

# missing curve check
if (nrow(curve)==0){
  stop("Error: There is no curve existing under selected dimension mix.")
}else{
  index=curve$k==0 | curve$g==0 | curve$v==0
  if (sum(index)!=0) {
    #curve[index,c("Audience","Frequency","Media"),with=F]
    print("Warning: There is curve(s) missing. The result below drops those missing curve(s).")
    curve=curve[!index]
  }
  
  # single curve simulation
  if (nrow(curve)==1){
    curve$sp_current=budget
    curve$g1=curve$g/curve$cpp
    curve$r_grs=curve$k*((1-exp(-curve$g1*curve$sp_current))^curve$v)
    col1=c("Media","Budget","Allocation","Gross reach","Net reach","Total 30s GRPs")
    curve$grp=curve$sp_current/curve$cpp
    output1=curve[,c("sp_current","r_grs","grp"),with=F]
    output1=output1[1,]
    output2=data.table(Media=c("Budget","Gross reach","Total 30s GRPs"),v1=as.vector(as.matrix(output1)))
    setnames(output2,"v1",curve[["Media"]])
    write.table(output2,"output_alloc_net_net.csv",col.names = F,sep=",")
  }else{
    # check some posibble errors
    if (sum(curve$sp_current)>budget){
      stop("Error: Total current spend is greater than total planed budget.")
    }else if (sum(curve$sp_current>curve$max_spend) !=0) {
      stop("Error: One or more media channel's spend is greater than its maximum constrain.")
    }else{
      #######################################################################################
      # OPTM initialization
      #######################################################################################
      print("Optimization")
      fn=function(x){
        #x=curve$sp_current
        curve=curve[,r_grs:=k*((1-exp(-g1*x)))**v]
        result.random.shell=curve[,c("media_num","r_grs"),with=F]
        # compute random scenario
        for.ran=result.random.shell[order(media_num)]
        n=nrow(for.ran)
        no.media=1:n
        r.dupe.box=rep(0,n)
        r.dupe.box.media=matrix(0,nr=n,nc=n,dimnames=list(for.ran$media_num,1:n))
        r.dupe.box.media[,1]=for.ran$r_grs
        for (j in 2:n){
          #j=2
          combo=combn(no.media,j)
          r.combo=matrix(for.ran$r_grs[combo],nr=j)
          r.combo.prod=apply(r.combo,2,prod)
          r.combo.tran=data.table(t(rbind(combo,r.combo.prod)))
          r.combo.melt=data.table(melt(r.combo.tran,id="r.combo.prod"))
          r.combo.meida=r.combo.melt[,list(r.prod=sum(r.combo.prod)),by=c("value")]
          r.dupe.box.media[,j]=r.combo.meida$r.prod
          r.dupe.box[j]=sum(r.combo.prod)
        }
        -(sum(for.ran$r_grs)+sum((-1)^(1:n+1)*r.dupe.box))
      }
      hin<- function(x) {
        h <- numeric(1)
        h[1] <- budget -sum(x)
        h[2] <- sum(x)-budget
        return(h)
      }
      a=curve$min_spend
      b=curve$max_spend
      b[a==b]=b[a==b]+0.01
      x0=budget*curve$k/sum(curve$k)
      x0[x0>b]=b[x0>b]
      x0[x0<a]=a[x0<a]
      optm <- cobyla(x0, fn, hin = hin, lower=a, upper=b,
                     nl.info = F, control = list(maxeval = 10000,ftol_abs=1e-6,xtol_rel=1e-6))
      curve$sp_current=optm$par
      
      #######################################################################################
      # OPTM Rollup Random
      #######################################################################################
      result.random.shell=curve[,c("media_num","shell_num","Media","max_spend","min_spend","sp_current","r_grs"),
                                with=F]
      # compute random scenario
      for.ran=result.random.shell[order(media_num)]
      no.media=1:nrow(curve)
      
      r.dupe.box=rep(0,nrow(curve))
      r.dupe.box.media=matrix(0,nr=nrow(curve),
                              nc=nrow(curve),
                              dimnames=list(for.ran$media_num,
                                            1:nrow(curve)))
      r.dupe.box.media[,1]=for.ran$r_grs
      
      for (j in 2:nrow(curve)){
        #j=2
        combo=combn(no.media,j)
        r.combo=matrix(for.ran$r_grs[combo],nr=j)
        r.combo.prod=apply(r.combo,2,prod)
        r.combo.tran=data.table(t(rbind(combo,r.combo.prod)))
        r.combo.melt=data.table(melt(r.combo.tran,id="r.combo.prod"))
        r.combo.meida=r.combo.melt[,list(r.prod=sum(r.combo.prod)),by=c("value")]
        r.dupe.box.media[,j]=r.combo.meida$r.prod
        r.dupe.box[j]=sum(r.combo.prod)
      }
      r_net_ran_total=sum(for.ran$r_grs)+sum((-1)^(1:nrow(curve)+1)*r.dupe.box)
      for.ran$r_net=apply(sweep(r.dupe.box.media,MARGIN=2,(-1)^(1:nrow(curve)+1),"*"),
                          1,sum)
      for.ran$r_dupe=for.ran$r_grs-for.ran$r_net
      r_dupe_ran_total=r_net_ran_total-sum(for.ran$r_net) #zx
      
      #######################################################################################
      # OPTM Rollup Max
      #######################################################################################
      result.max.shell=curve[,c("media_num","shell_num","Media","max_spend","min_spend","sp_current","r_grs"),
                                with=F]
      
      # compute max scenario
      for.max=result.max.shell[order(-r_grs)]
      for.max$r_net=rep(0,nrow(for.max))
      #r_grs_lag=c(for.max$r_grs[2:length(for.max$r_grs)],0)
      for.max$r_net[1]=for.max$r_grs[1]-for.max$r_grs[2]
      r_net_max_total=max(for.max$r_grs)
      r_dupe_max_total=max(for.max$r_grs[for.max$r_grs!=max(for.max$r_grs)])
      for.max=for.max[order(media_num)]
      for.max$r_dupe=for.max$r_grs-for.max$r_net #zx
      
      # calc final results
      r_dupe_ran_total=max(r_dupe_ran_total,r_dupe_ran_total1)
      r_dupe_max_total=max(r_dupe_max_total,r_dupe_max_total1)
      r_net_max_total=min(r_net_max_total,r_net_max_total1)
      r_net_ran_total=min(r_net_ran_total,r_net_ran_total1)
      
      # compute the final dupe factor
      for (k in 2:ncol(ex.dupe)) set(ex.dupe, j=k, value=as.numeric(ex.dupe[[k]]))
      dupe=melt(ex.dupe,id.vars="Media")
      dupe=dupe[(!is.na(dupe$value)),]
      dupe=merge(dupe,ex.media,by="Media",all.x=T)
      dupe=dupe[dupe$flag_media==1,]
      setnames(dupe,"Media","Media1")
      setnames(dupe,"variable","Media")
      dupe=merge(dupe,ex.media,by="Media",all.x=T)
      dupe=dupe[dupe$flag_media.y==1,]
      dupe.temp=dupe[,c("Media1","media_num.y","value"),with=F]
      setnames(dupe.temp,"Media1","Media")
      setnames(dupe.temp,"media_num.y","media_num.x")
      dupe.temp1=rbind(dupe[,c("Media","media_num.x","value"),with=F],dupe.temp)
      setnames(dupe.temp1,"media_num.x","media_num")
      dupe.temp1=data.table(dupe.temp1)
      dupe.media=dupe.temp1[,list(dupe=mean(value)),by=c("media_num")]
      
      # rescale spend, net reach and dupe reach
      for.max1=copy(for.max)
      setnames(for.max1,"sp_current","sp_max")
      setnames(for.max1,"r_net","r_net_max")
      setnames(for.max1,"r_dupe","r_dupe_max")#zx
      
      for.ran1=for.ran[,!c("media_num","shell_num","max_spend","min_spend","r_grs"),with=F]
      setnames(for.ran1,"sp_current","sp_ran")
      setnames(for.ran1,"r_net","r_net_ran")
      setnames(for.ran1,"r_dupe","r_dupe_ran")#zx
      final=merge(for.max1,for.ran1,by="Media",all=T)
      final=merge(final,dupe.media,by="media_num",all.x=T)
      
      # re-calc actual budget
      budget=min(sum(final$sp_max),sum(final$sp_ran))
      
      # net reach
      if (r_net_max_total>=r_net_ran_total){
        min=final$r_net_ran
        max=final$r_net_max
      }else{
        max=final$r_net_ran
        min=final$r_net_max
      }
      
      r_net_adj=min+(1-final$dupe)*(max-min)
      final$r_net_adj=r_net_adj
      
      # spend 
      if (r_net_max_total>=r_net_ran_total){
        min=final$sp_ran
        max=final$sp_max
      }else{
        max=final$sp_ran
        min=final$sp_max
      }
      
      sp_adj=min+(1-final$dupe)*(max-min)
      final$sp_adj=sp_adj
      
      # rescale spend
      final=merge(final,curve[,c("media_num","k","g1","v","cpp"),with=F],by=c("media_num"),all.x=T)
      if(sum(final$sp_adj)!=budget){
        final$sp_adj=final$sp_adj+((budget-sum(final$sp_adj))*final$sp_adj/sum(final$sp_adj))
      }
      final$flag_min=final$sp_adj<final$min_spend
      final$flag_max=final$sp_adj>final$max_spend
      final$flag=final$flag_min+final$flag_max
      final$sp_adj[final$flag==1]=
        (final$min_spend*final$flag_min+final$flag_max*final$max_spend)[final$flag==1]
      if(sum(final$sp_adj)!=budget){
        final$sp_adj[final$flag!=1]=final$sp_adj[final$flag!=1]+
          ((budget-sum(final$sp_adj))*final$sp_adj[final$flag!=1]/sum(final$sp_adj[final$flag!=1]))
      }
      final$allo=final$sp_adj/budget
      final$r_grs=final$k*((1-exp(-final$g1*final$sp_adj))^final$v)
      final$grp=final$sp_adj/final$cpp
      
      #total net and dupe reach 
      r_net_adj_total=(1-mean(dupe.temp1$value))*(max(c(r_net_max_total,r_net_ran_total))-
                                                    min(c(r_net_max_total,r_net_ran_total)))+min(c(r_net_max_total,r_net_ran_total))
      
      
      r_dupe_adj_total=r_net_adj_total-sum(final$r_net_adj)
      
      
      final.total=data.table(r_net_max_total,r_dupe_max_total,r_net_ran_total,
                             r_dupe_ran_total,r_net_adj_total,r_dupe_adj_total)
      
      
      # for output dupe matrix
      index=ex.dupe$Media %in% ex.media$Media[ex.media$flag_media==1]
      dupe=ex.dupe[index]
      dupe=dupe[,c(T,index),with=F]
      dupe1=plyr::join(dupe,final[,c("Media","r_grs"),with=F],by=c("Media"),type="left")
      dupe.m1=matrix(rep(dupe1$r_grs,nrow(dupe)),nr=nrow(dupe),byrow=F)#zx
      dupe.m2=matrix(rep(dupe1$r_grs,nrow(dupe)),nr=nrow(dupe),byrow=T)#zx
      dupe.m3=data.matrix(dupe[,!"Media",with=F])
      dupe.min=dupe.m1*dupe.m2
      dupe.max=matrix(0,nr=nrow(dupe.m1),nc=ncol(dupe.m1))
      for (i in 1:nrow(dupe.m1)){
        index=which(dupe.m1[i,]<dupe.m2[i,])
        dupe.max[i,]=dupe.m2[i,]
        dupe.max[i,][index]=dupe.m1[i,][index]  
      }
      dupe.final=dupe.min+dupe.m3*(dupe.max-dupe.min)
      dupe.final=cbind(dupe$Media,dupe.final)
      colnames(dupe.final)[1]="Media"
      
      # for output
      for.output=data.frame(final[,c("Media","sp_adj","allo","r_grs","r_net_adj","grp"),with=F])
      # several adjustments for two channel scenario
      r_n=for.output$r_grs
      if(sum(r_n!=0)==2){
        index=dupe.final[,-1]!=0
        index[is.na(index)]=F
        r_dupe_adj_total=as.numeric(dupe.final[,-1][(!is.na(dupe.final[,-1]))&index])
        r_net_adj_total=sum(for.output$r_grs)-r_dupe_adj_total
        for.output$r_net_adj[r_n!=0]=for.output$r_grs[r_n!=0]-r_dupe_adj_total
      }
      
      x=c("Total",apply(for.output[,-1],2,sum))
      x[names(for.output)=="r_net_adj"]=r_net_adj_total
      for.output=rbind(for.output,c("All other duplicated reach",c("","","",r_dupe_adj_total,"")))
      for.output=rbind(for.output,x)
      names(for.output)=c("Media","Budget","Allocation","Gross reach","Net reach","Total 30s GRPs")
      for.output.t=data.frame(t(for.output))
      names(for.output.t)=for.output$Media
      for.output.t=data.frame(for.output.t[curve$Media],for.output.t[c("All other duplicated reach","Total")])
      
      write.table(for.output.t,"output_alloc_net_net.csv",col.names = F,sep=",")
      write.csv(dupe.final,"output_dupe_net_net.csv",row.names=F)
    }# error check
  }#single curve check
}# zero curve check
end=Sys.time()-start
print(paste("Modeling time: ",round(end[[1]],digit=2),attr(end,"units"),sep=""))
