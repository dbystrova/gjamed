rm(list=ls())
#setwd("~/Documents/GitHub/gjamed")
setwd("~/GitHub/gjamed")
########## Simulating the data
library(MASS)
library(repmis)
library(gjam)
library(rlist)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(NLRoot)
library(Rcpp)
library(plyr)
library(ggplot2)
library(ggsn)
library(parallel)detectCores()
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


generate_data<-function(Sp=50,nsamples=500,qval=20,Ktrue=4){
  S<-Sp
  n<- nsamples
  q<- qval
  env<-runif(-50,50,n=n)
  X<-cbind(1,poly(env,2)) #nxK
  idx<-sample(S)
  B_0<-seq(0,100,length.out=S)[idx]
  B_1<-seq(0,100,length.out=S)[idx]
  B_2<-seq(0,100,length.out=S)[idx]
  B<-cbind(B_0,B_1,B_2) #SxK
  L<-X%*%t(B) #nxS
  K_t<- Ktrue
  cat("True number of clusters : ",K_t,"\n")
  A<-matrix(NA,nrow=ceiling(K_t),ncol=q)
  for(i in 1:ceiling(K_t)){
    A[i,]<-mvrnorm(n = 1,rep(0,q), Sigma=3*diag(q)) #Nxq short and skinny
  }
  idx<-sample((1:ceiling(K_t)),S,replace=T) 
  Lambda<-A[idx,] #Sxr tall and skinny
  Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  Sigma_true<-Sigma
  Y<-vector()
   for(i in 1:n){
   Y<-rbind(Y,mvrnorm(1, mu=L[i,], Sigma=Sigma))
   }
  #Y<-mvrnorm(n = n, mu=rep(0,S), Sigma=Sigma)
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env1","env2")
  
  return(list(xdata=xdata, Y=Y,idx=idx,S_true=Sigma_true))
}


simulation_fun_oneDS<-function(data_set,Sp, Ntr, rval,nsamples=500, Ktrue,q=20, it=1000, burn=500,type="GJAM"){
  S<-Sp
  n<- nsamples
  r <- rval
  iterations<-it
  #env<-runif(-50,50,n=n)
  #X<-cbind(1,poly(env,2)) #nxK
  #idx<-sample(S)
  #B_0<-seq(0,100,length.out=S)[idx]
  #B_1<-seq(0,100,length.out=S)[idx]
  #B_2<-seq(0,100,length.out=S)[idx]
  #B<-cbind(B_0,B_1,B_2) #SxK
  #L<-X%*%t(B) #nxS
  
  K=sum(S/(S+(1:S)-1)) #104, his prior number of clusters when alpha=S
  if(type=="GJAM"){ cat("Prior expected number of clusters : ",K,"\n")}
  else{cat("Prior expected number of clusters : ",Ktrue,"\n")}
  K_t= Ktrue
  cat("True number of clusters : ",K_t,"\n")
  # A<-matrix(NA,nrow=ceiling(K_t),ncol=q)
  # # sig=matrix(runif(n=q*q),ncol=q)
  # for(i in 1:ceiling(K_t)){
  #   A[i,]<-mvrnorm(n = 1,rep(0,q), Sigma=3*diag(q)) #Nxq short and skinny
  # }
  # idx<-sample((1:ceiling(K_t)),S,replace=T) 
  # Lambda<-A[idx,] #Sxr tall and skinny
  # Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  # Sigma_true<-Sigma
  # Y<-mvrnorm(n = n, mu=rep(0,S), Sigma=Sigma)
  # 
  # xdata<-as.data.frame(X[,-1])
  # colnames(xdata)<-c("env1","env2")
  xdata<-data_set$xdata
  Y<-data_set$Y
  idx<- data_set$idx
  Sigma_true<- data_set$S_true
  formula<-as.formula(~env1+env2)
  if(type=="GJAM"){
    rl <- list(r = r, N =Ntr-1)
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-gjam(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
    alpha.DP<-S
    pk_chains<-NULL
  }
  
  
  if(type=="0"){
    #func<-function(x) {sum(x/(x+(1:S)-1))-K_t}
    #alpha.DP<-.bisec(func,0.01,100)
    alpha.DP<-S
    rl <- list(r = r, N =Ntr-1,alpha.DP=alpha.DP)
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-.gjam0(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
    pk_chains<- fit$chains$pk_g
  }
  if(type=="1"){
    func<-function(x) {sum(x/(x+(1:S)-1))-K_t}
    alpha.DP<-.bisec(func,0.01,100)
    shape=((alpha.DP)^2)/20
    rate=alpha.DP/20
    rl  <- list(r = r, N = Ntr-1, rate=rate,shape=shape)
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-.gjam_1(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
    pk_chains<- fit$chains$pk_g
  }
  if(type=="2"){
    func<-function(x) {sum(x/(x+(1:S)-1))-K_t}
    alpha.DP<-.bisec(func,0.01,100)
    shape=((alpha.DP)^2)/20
    rate=alpha.DP/20
    rl  <- list(r = r, N = Ntr-1, rate=rate,shape=shape,V=1)
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-.gjam_2(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
    pk_chains<- fit$chains$pk_g
  }
  if(type=="3"){
    eps=0.05
    # alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    # colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    # alp_sig$sigma=seq(0.05,0.5,length.out = 20)
    # #loop to run bisection on a grid for sigma
    # for(i in 1:20){
    #   ####corrected formula : added -1
    #   func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S) -1)/(x+c(1:S) -1))-1) - K_t}
    #   alp_sig[i,"alpha"]<-.bisec(func,0.0001,100)
    #   N_eps<-floor(.compute_tau_mean_large_dim(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps) + 2*.compute_tau_var_large_dim(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps))
    #   ifelse(N_eps<=150,alp_sig[i,"is_less_150"]<-T,alp_sig[i,"is_less_150"]<-F)
    #   N_eps
    # }
    # if(sum(alp_sig$is_less_150==T)==0) cat("!! no choice under N=150, need to recheck!!!")
    # k<-min(which(alp_sig$is_less_150==T)) #max sigma s.t. N<150
    # sigma_py<-alp_sig[k,"sigma"]
    # alpha.PY<-alp_sig[k,"alpha"]
    
    ##CHOOSING PARAMETERS ALPHA and SIGMA
    sigma_py<-0.25
    funcPY_root<-function(x) {(x/sigma_py)*(prod((x+sigma_py+c(1:S) -1)/(x+c(1:S) -1))-1) - K_t}
    alpha.PY<-.bisec(funcPY_root,0.0001,100)
    
    
    N_eps<-floor(.compute_tau_mean_large_dim(sigma_py,alpha.PY,eps) + 2*.compute_tau_var_large_dim(sigma_py,alpha.PY,eps))
    rl   <- list(r = r, N = N_eps, sigma_py=sigma_py, alpha=alpha.PY)
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-.gjam_3(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
    pk_chains<- fit$chains$pk_g
    alpha.DP<-alpha.PY
    Ntr<-N_eps+1
  }
  if(type=="4"){
    eps=0.1
    alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    alp_sig$sigma=seq(0.05,0.4,length.out = 20)
    #loop to run bisecetion on a grid for sigma
    for(i in 1:20){
      ####corrected added  -1
      func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S)-1)/(x+c(1:S) -1))-1) - K_t}
      alp_sig[i,"alpha"]<-.bisec(func,0.01,100)
      N_eps<-floor(.compute_tau_mean_large_dim(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps) + 2*.compute_tau_var_large_dim(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps))
      ifelse(N_eps<=150,alp_sig[i,"is_less_150"]<-T,alp_sig[i,"is_less_150"]<-F)
      N_eps
    }
    
    if(sum(alp_sig$is_less_150==T)==0) cat("!! no choice under N=150, need to recheck!!!")
    
    k<-max(which(alp_sig$is_less_150==T)) #max sigma s.t. N<150
    sigma_py<-alp_sig[k,"sigma"]
    alpha.PY<-alp_sig[k,"alpha"]
    #fixing hyperparameters
    ro.disc=1-2* sigma_py
    shape=((alpha.PY)^2)/10
    rate=alpha.PY/10
    # 95% quantile of alpha
    alpha.max=qgamma(.95, shape=shape, rate=rate)
    alpha.max_val<-5
    sigma_py_max<-0.5
    N_eps<-floor(.compute_tau_mean_large_dim(sigma_py_max,alpha.max_val,eps) + 2*.compute_tau_var_large_dim(sigma_py_max,alpha.max_val,eps))
    
    
    rl   <- list(r = r, N = N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
    ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
    fit<-.gjam_4(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.PY_g
    sigma.chains<-fit$chains$discount.PY_g
    pk_chains<- fit$chains$pk_g
    Ntr<-N_eps+1
    alpha.DP<-alpha.PY
    
  }
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  df<-as.data.frame(trace)
  df$iter<-1:it
  #plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
  p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
    labs(title=paste0("Trace plot for the number of groups K for S=",S," r=",r," true K=",K_t," type=",type), caption=paste0("Number of iterations: ",it," burnin: ",burn,"number of samples: ",nsamples))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = Ktrue,color = "red")
  #plot(p)
  #####Weights plot
  if(type%in%c("0","1","2","3","4")){
    pk<- apply(fit$chains$pk_g[-c(1:burn),],2,mean)
    last_pk<- round(pk[Ntr-1],3)
    df_weights <- data.frame(matrix(NA, nrow = Ntr-1, ncol =1))
    df_weights$pw<-pk 
    df_weights$tr<-1:(Ntr-1)
    pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
      geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
      geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Weights for the case: S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr, " pN=",last_pk), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
    # pl_weigths
  }
  #####Alpha plot
  if(type%in%c("1","2","4")){ 
    df_alpha <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
    df_alpha$alpha<- alpha.chains[-c(1:burn)]
    df_alpha$type<- "posterior"
    df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
    #df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
    alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
    df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)
    
    df_alpha_prior$type<- "prior"
    df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
    ###Compute mean
    mu <- ddply(df_alpha_all, "type", summarise, grp.mean=mean(alpha))
    mu$grp.mean[which(mu$type=="prior")]=alpha.DP
    p_alpha_1<- ggplot(df_alpha_all, aes(x=alpha, color=type)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
      geom_density()+labs(title=paste0("Distribution alpha: S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
    # p_alpha_1
    
    p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
      geom_density(color="red")+labs(title=paste0("Posterior distribution for alpha"), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples," S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
      scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
    # p_alpha_2
  }  
  #######Sigma plot
  if(type%in%c("4")){ 
    df_sigma <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
    df_sigma$sigma<- sigma.chains[-c(1:burn)]
    ###Compute mean
    mu <- ddply(df_sigma, type, summarise, grp.mean=mean(sigma))
    p_sigma<- ggplot(df_sigma, aes(x=sigma)) + geom_vline(data=mu, aes(xintercept=grp.mean),linetype="dashed")+
      geom_density()+labs(title=paste0("Distribution sigma: S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
    # p_sigma
    
    dfs<-as.data.frame(sigma.chains)
    dfs$iter<-1:it
    #plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
    p_trace_sigma<-ggplot(dfs, aes(y=sigma.chains, x=iter)) + geom_point() + 
      labs(title=paste0("Trace plot for the sigma for S=",S," r=",r," true K=",K_t," type=",type), caption=paste0("Number of iterations: ",it," burnin: ",burn,"number of samples: ",nsamples))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
    # p_trace_sigma
  } 
  
  
  ###TOPDF
  # pdf(paste0("plots/Plot-clusters_S",S,"r",r,"Ktr",K_t,"mod_type_",type,".pdf"))
  # plot(p)
  # if(type%in%c("1","2","4")){ 
  #   plot(p_alpha_1)
  #   plot(p_alpha_2)
  # }
  # if(type%in%c("4")){ 
  #   plot(p_sigma)
  #   plot(p_trace_sigma)
  # }
  # if(type%in%c("0","1","2","3","4")){ 
  #   plot(pl_weigths)
  # }
  #dev.off()
  # 
  # N_dim<-(it-burn)
  # Z<-array(dim=c(Sp,r,N_dim))
  # for(j in 1:N_dim){
  #   K<-fit$chains$kgibbs[j,]
  #   Z[,,j]  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)[K,]
  #   #sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
  # }
  # Lambda_mean<-apply(Z,c(1,2),mean)
  # err<-sum((Lambda_mean-Lambda)^2)/(S*q)
  # fit<-fit$rmspeAll
  #fit_er<-fit$rmspeAll
  N_dim<-(it-burn)
  sigma<-array(dim=c(Sp,Sp,N_dim))
  for(j in 1:N_dim){
    K<-fit$chains$kgibbs[j,]
    Z  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)
    sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
    
  }
  sigma_mean<-apply(sigma,c(1,2),mean)
  err<-sum((sigma_mean-Sigma_true)^2)/(Sp*Sp)
  rmspe<-fit$fit$rmspeAll
  
  ########Plot error graph
  M_m<- as.vector(sigma_mean)
  M_t<- as.vector(Sigma_true)
  df_vs<- as.data.frame(M_m)
  df_vs$M_t<- M_t
  plot_vs <- ggplot(df_vs)+
    aes(x = M_t, y = M_m)+geom_point()+geom_smooth(method = "lm") +  geom_abline(intercept = -min(M_t), slope = 1, color="red")+
    labs(title=paste0(" True vs estimated covariance parameters"), caption=paste0("Number of iterations: ",it," burnin: ",burn,"number of samples: ",nsamples,"S=",S," r=",r," true K=",K_t," type=",type))+ xlab("True")+ylab("Estimated")+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
  
  #plot(plot_vs)
  
  pdf(paste0("plots/Plot-clusters_S",S,"r",r,"Ktr",K_t,"mod_type_",type,".pdf"))
  plot(p)
  if(type%in%c("1","2","4")){ 
    plot(p_alpha_1)
    plot(p_alpha_2)
  }
  if(type%in%c("4")){ 
    plot(p_sigma)
    plot(p_trace_sigma)
  }
  if(type%in%c("0","1","2","3","4")){ 
    plot(pl_weigths)
  }
  plot(plot_vs)
  
  dev.off()
  #chain=fit$chains$kgibbs,
  return(list(trace=trace,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,pk_val=pk, pkN=pk_chains[,Ntr-1], 
              coeff_t=Sigma_true,coeff_f=sigma_mean,
              err=err,fit=rmspe))
}

####### Just one possible test case
S_vec<-100
q<-20
r_vec<-5
Ktr<-10
n_samples<-500
l<-1
it<-5000
burn<-200

#####################################################################################
#one dataset only
data_set<- generate_data(Sp=S_vec,nsamples=n_samples,qval=q,Ktrue=Ktr)
save(data_set, file = paste0("data/DS_S_",S_vec,"_q_",q,"_n_500_",Ktr,"l_",l,".Rda"))
list0<-list2<-list3<-list4<-list5<-NULL
list0<-list.append(list0,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K_",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec, q=20,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="0")))
list2<-list.append(list2,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=150,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="1")))
#names(list2)<-paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l)
list3<-list.append(list3,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="2")))
#names(list3)<-paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l)
list4<-list.append(list4,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="3")))
#names(list4)<-paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l)
list5<-list.append(list5,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="4")))
#names(list5)<-paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l)
# x_1<-x_2<-x_3<-x_4<-x_5<-1:25
#  for(i in 0:4){
#  for(j in 1:39) assign(paste0("x_",i+1),c(get(paste0("x_",i+1)), (5*(5*j+i)+1):(5*(5*j+i)+5)))
#  }
# 
# list<-list(list0,list2,list3,list4,list5)
# table<-data.frame()
# for(i in 1:length(list)){
#   str<-names(list)[[i]]
#   tmp<-data.frame("trace"=list[[i]]$trace[get(paste0("x_",i))],"S"=rep(S_vec,length(list[[i]]$trace[get(paste0("x_",i))])),
#                   "type"=rep(i-1,length(list[[i]]$trace[get(paste0("x_",i))])),
#                   "x"=get(paste0("x_",i)))
#   table<-rbind(table,tmp)
# }

table<-data.frame()
for(i in 1:length(list)){
  str<-names(list)[[i]]
  tmp<-data.frame("trace"=list[[i]]$trace,"S"=rep(S_vec,length(list[[i]]$trace)),
                  "type"=rep(i-1,length(list[[i]]$trace)),
                  "x"=1:length(list[[i]]$trace))
  table<-rbind(table,tmp)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

p1<-ggplot(table[which(table$type=="0"),], aes(x=table$x[which(table$type=="0")],y=table$trace[which(table$type=="0")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p1
p2<-ggplot(table[which(table$type=="1"),], aes(x=table$x[which(table$type=="1")],y=table$trace[which(table$type=="1")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p2
p3<-ggplot(table[which(table$type=="2"),], aes(x=table$x[which(table$type=="2")],y=table$trace[which(table$type=="2")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p3
p4<-ggplot(table[which(table$type=="3"),], aes(x=table$x[which(table$type=="3")],y=table$trace[which(table$type=="3")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p4

p<-ggplot(table, aes(x=x,y=trace,col=as.factor(type)))+geom_point()+
   scale_color_manual(name = c(""), values = cols, labels=c("Original model","DP with prior on alpha 1","DP with prior on alpha 2","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))+
   labs(title="Trace of the number of clusters for the different models. S=100")+xlab("iterations")+theme_bw()
pdf("plots/POSTER_plot_all.pdf")
p
dev.off()


list0$err
list3$err
list4$err
list5$err


# 
# table_K_10_S_100_r_5_mean<-table
# save(table_K_10_S_100_r_5_mean,file="table_K_10_S_100_r_5_mean.rda")


#####################################################################################
#5 datasets
l<-1
data_set<- generate_data(Sp=S_vec,nsamples=n_samples,qval=q,Ktrue=Ktr)
#save(data_set, file = paste0("data/DS_S_",S_vec,"_q_",q,"_n_500_",Ktr,"l_",l,".Rda"))
list0<-assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K_",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec, q=20,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="0"))
list2<-assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=150,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="1"))
list3<-assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="2"))
list4<-assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="3"))
list5<-assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="4"))



list0<-list2<-list3<-list4<-list5<-NULL
for(l in 2:5){
  data_set<- generate_data(Sp=S_vec,nsamples=n_samples,qval=q,Ktrue=Ktr)
  #save(data_set, file = paste0("data/DS_S_",S_vec,"_q_",q,"_n_500_",Ktr,"l_",l,".Rda"))
  list0<-list.append(list0,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K_",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec, q=20,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="0")))
  list2<-list.append(list2,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=150,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="1")))
  list3<-list.append(list3,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="2")))
  list4<-list.append(list4,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="3")))
  list5<-list.append(list5,assign(paste0("S_",S_vec,"_r_",r_vec,"_N_150_n_500_K",Ktr,"l_",l),simulation_fun_oneDS(data_set,Sp=S_vec, Ntr=S_vec,q=20, rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="4")))
}

# 
# tmp_list<-list(trace=list0$trace,
#                idx=list0$idx,K=list0$K,
#                alpha=list0$alpha,alpha.chains=list0$alpha.chains,pk_val=list0$pk_val, pkN=list0$pkN, 
#                coeff_t=list0$coeff_f,coeff_f=list0$coeff_f,
#                err=list0$err,fit=list0$fit)
# list0<-list(tmp_list,list0[[12]],list0[[13]],list0[[14]],list0[[15]])
# 
# tmp_list<-list(trace=list2$trace,
#                idx=list2$idx,K=list2$K,
#                alpha=list2$alpha,alpha.chains=list2$alpha.chains,pk_val=list2$pk_val, pkN=list2$pkN, 
#                coeff_t=list2$coeff_f,coeff_f=list2$coeff_f,
#                err=list2$err,fit=list2$fit)
# list2<-list(tmp_list,list2[[12]],list2[[13]],list2[[14]],list2[[15]])
# 
# tmp_list<-list(trace=list3$trace,
#                idx=list3$idx,K=list3$K,
#                alpha=list3$alpha,alpha.chains=list3$alpha.chains,pk_val=list3$pk_val, pkN=list3$pkN, 
#                coeff_t=list3$coeff_f,coeff_f=list3$coeff_f,
#                err=list3$err,fit=list3$fit)
# list3<-list(tmp_list,list3[[12]],list3[[13]],list3[[14]],list3[[15]])
# 
# tmp_list<-list(trace=list4$trace,
#                idx=list4$idx,K=list4$K,
#                alpha=list4$alpha,alpha.chains=list4$alpha.chains,pk_val=list4$pk_val, pkN=list4$pkN, 
#                coeff_t=list4$coeff_f,coeff_f=list4$coeff_f,
#                err=list4$err,fit=list4$fit)
# list4<-list(tmp_list,list4[[12]],list4[[13]],list4[[14]],list4[[15]])
# 
# tmp_list<-list(trace=list5$trace,
#                idx=list5$idx,K=list5$K,
#                alpha=list5$alpha,alpha.chains=list5$alpha.chains,pk_val=list5$pk_val, pkN=list5$pkN, 
#                coeff_t=list5$coeff_f,coeff_f=list5$coeff_f,
#                err=list5$err,fit=list5$fit)
# list5<-list(tmp_list,list5[[12]],list5[[13]],list5[[14]],list5[[15]])
# table<-data.frame()
# for(i in 1:length(list)){
#   str<-names(list)[[i]]
#   tmp<-data.frame("trace"=list[[i]]$trace[get(paste0("x_",i))],"S"=rep(S_vec,length(list[[i]]$trace[get(paste0("x_",i))])),
#                   "type"=rep(i-1,length(list[[i]]$trace[get(paste0("x_",i))])),
#                   "x"=get(paste0("x_",i)))
#   table<-rbind(table,tmp)
# }


table_0<-data.frame()
for(i in 1:length(list0)){
  table_0[i,"err"]<-list0[[i]]$err
  table_0[i,"K_mean"]<-mean(list0[[i]]$trace)
  table_0[i,"K_last"]<-mean(list0[[i]]$trace[length(list0[[i]]$trace)])
  table_0[i,"type"]<-"0"
  table_0[i,"S"]<-S_vec
  table_0[i,"Ktr"]<-Ktr
}

table_2<-data.frame()
for(i in 1:length(list2)){
  table_2[i,"err"]<-list2[[i]]$err
  table_2[i,"K_mean"]<-mean(list2[[i]]$trace)
  table_2[i,"K_last"]<-mean(list2[[i]]$trace[length(list2[[i]]$trace)])
  table_2[i,"type"]<-"2"
  table_2[i,"S"]<-S_vec
  table_2[i,"Ktr"]<-Ktr
}
table_3<-data.frame()
for(i in 1:length(list3)){
  table_3[i,"err"]<-list3[[i]]$err
  table_3[i,"K_mean"]<-mean(list3[[i]]$trace)
  table_3[i,"K_last"]<-mean(list3[[i]]$trace[length(list3[[i]]$trace)])
  table_3[i,"type"]<-"3"
  table_3[i,"S"]<-S_vec
  table_3[i,"Ktr"]<-Ktr
}
table_4<-data.frame()
for(i in 1:length(list4)){
  table_4[i,"err"]<-list4[[i]]$err
  table_4[i,"K_mean"]<-mean(list4[[i]]$trace)
  table_4[i,"K_last"]<-mean(list4[[i]]$trace[length(list4[[i]]$trace)])
  table_4[i,"type"]<-"4"
  table_4[i,"S"]<-S_vec
  table_4[i,"Ktr"]<-Ktr
}
table_5<-data.frame()
for(i in 1:length(list5)){
  table_5[i,"err"]<-list5[[i]]$err
  table_5[i,"K_mean"]<-mean(list5[[i]]$trace)
  table_5[i,"K_last"]<-mean(list5[[i]]$trace[length(list5[[i]]$trace)])
  table_5[i,"type"]<-"5"
  table_5[i,"S"]<-S_vec
  table_5[i,"Ktr"]<-Ktr
}
table<-rbind(table_0,table_2,table_3,table_4,table_5)

#plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)
p_err<-ggplot(table)+geom_boxplot(aes(x=as.factor(type),y= sqrt(as.numeric(err)), fill = as.factor(type)))+
  scale_fill_discrete(name = c("Models"), labels=c("Original model","DP with prior on alpha 1","DP with prior on alpha 2","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))+
  labs(title=paste0("Error in estimating the covariance matrix. S=",S_vec," r=",r_vec," Ktr=",Ktr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",n_samples))+
  ylab("Frobenius norm of the error")+xlab("")
pdf(paste0("plots/with_mean_SIGMA_ERROR_S_",S_vec," r=",r_vec," Ktr=",Ktr,".pdf"))
p_err
dev.off()

p_K<-ggplot(table)+geom_boxplot(aes(x=as.factor(type),y= as.numeric(K_mean), fill = as.factor(type)))+
  scale_fill_discrete(name = c("Models"), labels=c("Original model","DP with prior on alpha 1","DP with prior on alpha 2","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))+
  labs(title=paste0("Posterior mean of the retrieved number of clusters. S=",S_vec," r=",r_vec," Ktr=",Ktr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",n_samples))+
  ylab("Posterior mean of the retrieved number of clusters")+xlab("")+geom_hline(yintercept = Ktr,col="red")
pdf(paste0("plots/with_mean_POST_MEAN_CLUST_S_",S_vec," r=",r_vec," Ktr=",Ktr,".pdf"))
p_K
dev.off()

p_K_last<-ggplot(table)+geom_boxplot(aes(x=as.factor(type),y= as.numeric(K_last), fill = as.factor(type)))+
  scale_fill_discrete(name = c("Models"), labels=c("Original model","DP with prior on alpha 1","DP with prior on alpha 2","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))+
  labs(title=paste0("Last traceplot element of the retrieved number of clusters. S=",S_vec," r=",r_vec," Ktr=",Ktr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",n_samples))+
  ylab("Last traceplot element of the retrieved number of clusters")+xlab("")+geom_hline(yintercept = Ktr,col="red")
pdf(paste0("plots/with_mean_CLUST_LAST_S_",S_vec," r=",r_vec," Ktr=",Ktr,".pdf"))
p_K_last
dev.off()


