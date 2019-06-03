rm(list=ls())
setwd("~/Documents/GitHub/gjamed")
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
  Y<-mvrnorm(n = n, mu=rep(0,S), Sigma=Sigma)
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
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
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
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
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
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
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
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
    fit<-.gjam_2(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
    pk_chains<- fit$chains$pk_g
  }
  if(type=="3"){
    eps=0.05
    alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    alp_sig$sigma=seq(0.05,0.5,length.out = 20)
    #loop to run bisecetion on a grid for sigma
    for(i in 1:20){
      ####corrected formula : added -1
      func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S) -1)/(x+c(1:S) -1))-1) - K_t}
      alp_sig[i,"alpha"]<-.bisec(func,0.01,100)
      N_eps<-floor(.compute_tau_mean(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps) + 2*.compute_tau_var(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps))
      ifelse(N_eps<=150,alp_sig[i,"is_less_150"]<-T,alp_sig[i,"is_less_150"]<-F)
      N_eps
    }
    if(sum(alp_sig$is_less_150==T)==0) cat("!! no choice under N=150, need to recheck!!!")
    k<-max(which(alp_sig$is_less_150==T)) #max sigma s.t. N<150
    sigma_py<-alp_sig[i,"sigma"]
    alpha.PY<-alp_sig[i,"alpha"]
    
    N_eps<-floor(.compute_tau_mean(sigma_py,alpha.PY,eps) + 2*.compute_tau_var(sigma_py,alpha.PY,eps))
    rl   <- list(r = r, N = N_eps, sigma_py=sigma_py, alpha=alpha.PY)
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
    fit<-.gjam_3(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
    pk_chains<- fit$chains$pk_g
    alpha.DP<-alpha.PY
    Ntr<-N_eps+1
  }
  if(type=="4"){
    eps=0.05
    alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    alp_sig$sigma=seq(0.05,0.4,length.out = 20)
    #loop to run bisecetion on a grid for sigma
    for(i in 1:20){
      ####corrected added  -1
      func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S)-1)/(x+c(1:S) -1))-1) - K_t}
      alp_sig[i,"alpha"]<-.bisec(func,0.01,100)
      N_eps<-floor(.compute_tau_mean(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps) + 2*.compute_tau_var(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps))
      ifelse(N_eps<=150,alp_sig[i,"is_less_150"]<-T,alp_sig[i,"is_less_150"]<-F)
      N_eps
    }
    
    if(sum(alp_sig$is_less_150==T)==0) cat("!! no choice under N=150, need to recheck!!!")
    
    k<-max(which(alp_sig$is_less_150==T)) #max sigma s.t. N<150
    sigma_py<-alp_sig[i,"sigma"]
    alpha.PY<-alp_sig[i,"alpha"]
    #fixing hyperparameters
    ro.disc=1-2* sigma_py
    shape=((alpha.PY)^2)/10
    rate=alpha.PY/10
    # 95% quantile of alpha
    alpha.max=qgamma(.95, shape=shape, rate=rate)
    
    N_eps<-floor(.compute_tau_mean(sigma_py,alpha.max,eps) + 2*.compute_tau_var(sigma_py,alpha.max,eps))
    
    
    rl   <- list(r = r, N = N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
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

  return(list(trace=trace, chain=fit$chains$kgibbs,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,pk_chain=pk_chains, 
              coeff_t=Sigma_true,coeff_f=sigma,
              err=err,fit=rmspe))
}

####### Just one possible test case
#sim<-simulation_fun(Sp=50, Ntr=150, rval=3,nsamples=500, Ktrue=4,it=1000,burn=200,type=2)
# plot(as.vector(sim$coeff_t),as.vector(sim$coeff_f))
# x11()
# heatmap(sim$coeff_f)
# x11()
# heatmap(sim$coeff_t)
# plot(sim$trace)
# plot(sim$idx,sim$K)
#possible parameters to add in the loop:
# - type (T) of gjam to be fitted
# - n, the number of simulated normals
# - N, the truncation level
# - S, the number of species
# - K_t, the true number of clusters

###########Simulation for Continous data case :small S K=4###################################################


####Small S, N==S, n=500

list=list()
list2=list()
list3=list()
list4=list()
list5=list()
list0=list()
S_vec<-c(20,50,80,100)
r_vec<-c(5,10,15,20)
k<-1
it<-5000
burn<-2000
n_samples<-500
Ktr<-4
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    data_set<- generate_data(Sp=S_vec[i],nsamples=n_samples,qval=q,Ktrue=Ktrue)
    save(data_set, file = paste0("data/DS_S_",S_vec[i],"_q_",q,"_n_500_",Ktrue,"_.Rda"))
    
    list<-list.append(list,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i], q=20,rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="GJAM")))
    names(list)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list0<-list.append(list0,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i], q=20,rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="0")))
    names(list0)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list2<-list.append(list2,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=150,q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="1")))
    names(list2)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list3<-list.append(list3,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i],q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="2")))
    names(list3)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list4<-list.append(list4,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i],q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="3")))
    names(list4)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list5<-list.append(list5,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i],q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=4,it=it,burn=burn,type="4")))
    names(list5)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    k=k+1
  }
}


save(list, file = "ODSim_smallSK4_gjam.Rda")
save(list0, file = "ODSim_smallSK4_gjam0.Rda")
save(list2, file = "ODSim_smallSK4_type1.Rda")
save(list3, file = "ODSim_smallSK4_type2.Rda")
save(list4, file = "ODSim_smallSK4_type3.Rda")
save(list5, file = "ODSim_smallSK4_type4.Rda")








#########Testing K=8########################################################################################## 


####Small S, N==S, n=500

list=list()
list2=list()
list3=list()
list4=list()
list5=list()
list0=list()
S_vec<-c(20,50,80,100)
r_vec<-c(5,10,15,20)
k<-1
it<-5000
burn<-2000
n_samples<-500
q<-20
Ktr<-8
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    data_set<- generate_data(Sp=S_vec[i],nsamples=n_samples,qval=q,Ktrue=Ktr)
    save(data_set, file = paste0("data/DS_S_",S_vec[i],"_q_",q,"_n_500_",Ktr,"_.Rda"))
    
    list<-list.append(list,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_,150_n_500_Kt_",Ktr,"_TGJ"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i], q=20,rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="GJAM")))
    names(list)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list0<-list.append(list0,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_",Ktr,"_T0"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i], q=20,rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="0")))
    names(list0)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list2<-list.append(list2,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_",Ktr,"_T1"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=150,q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="1")))
    names(list2)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list3<-list.append(list3,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_",Ktr,"_T2"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=S_vec[i],q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="2")))
    names(list3)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list4<-list.append(list4,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_",Ktr,"_T3"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=150,q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="3")))
    names(list4)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    list5<-list.append(list5,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_",Ktr,"_T4"),simulation_fun_oneDS(data_set,Sp=S_vec[i], Ntr=150,q=20, rval=r_vec[j],nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="4")))
    names(list5)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    k=k+1
  }
}


save(list, file = "ODSim_smallSK8_gjam.Rda")
save(list0, file = "ODSim_smallSK8_gjam0.Rda")
save(list2, file = "ODSim_smallSK8_type1.Rda")
save(list3, file = "ODSim_smallSK8_type2.Rda")
save(list4, file = "ODSim_smallSK8_type3.Rda")
save(list5, file = "ODSim_smallSK8_type4.Rda")




#########Testing########################################################################################## 
#########Testing########################################################################################## 
#########Testing########################################################################################## 
#########Testing########################################################################################## 




listx=list()
S_vec<-c(20)
r_vec<-c(5)
k<-1
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    listx<-list.append(listx,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_1"),simulation_fun_paper(Sp=S_vec[i], Ntr=S_vec[i], q=20,rval=r_vec[j],nsamples=500, Ktrue=8,it=5000,burn=2000,type="1")))
    names(listx)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
    k=k+1
  }
}

pk<- apply(fit$chains$pk_g,2,mean)
#pk<-apply()
last_pk<- round(pk[Ntr-1],3)
#plot(1:(Ntr-1),pk,type='h', ylim=c(0,max(pk)),col=gray(.5),xlab=expression(italic(i)),ylab=expression(italic(p_k)) ,lwd=3)
#points(1:(Ntr-1),pk)
#mtext(paste0("Weights p_k for N=",Ntr," ,p_N=",last_pk),side=3)

df_weights <- data.frame(matrix(NA, nrow = Ntr-1, ncol =1))
df_weights$pw<-pk 
df_weights$tr<-1:(Ntr-1)
#df_weights$tr<- as.factor(df_weights$tr)
# p <- ggplot(df_weights, aes(x=tr,xend=tr,y=0,yend=pw)) +geom_segment()+ geom_point(size=0.1)
# p
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Weights for the case: S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr, " pN=",last_pk), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))

pl_weigths

df_alpha <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
df_alpha$alpha<- fit$chains$alpha.DP_g[-c(1:2000)]
df_alpha$type<- "posterior"
df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
df_alpha_prior$type<- "prior"
df_alpha<- rbind(df_alpha[-1,],df_alpha_prior[-1,])

mu <- ddply(df_alpha, "type", summarise, grp.mean=mean(alpha))


p_alpha<- ggplot(df_alpha, aes(x=alpha, color=type)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density()+labs(title=paste0("Distribution alpha: S=",S," ,r=",r," true gr K=",K_t," ,type=",type, " ,N=",Ntr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))

p_alpha
###########

sigma<- LtT2[[1]]$coeff_f
sigma_m<- apply(sigma,c(1,2),mean)
sigma_t<-LtT2[[1]]$coeff_t
M_m<- as.vector(sigma_m)
M_t<- as.vector(sigma_t)
plot(M_t,M_m)
df_vs<- as.data.frame(M_m)
df_vs$M_t<- M_t

plot_vs <- ggplot(df_vs)+
  aes(x = M_t, y = M_m)+geom_point()+geom_smooth(method = "lm") +  geom_abline(intercept = -min(M_t), slope = 1, color="red")
plot_vs

