#setwd("~/Bureau/PY_comparisons/gjamed-master")
rm(list=ls())
#setwd("~/Tesi/Code/modified_gjam/Gjam/")
setwd("~/Documents/GitHub/gjamed")
########## Simulating the data
library(MASS)
library(repmis)
library(gjam)
library(rlist)y
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(NLRoot)
library(Rcpp)
library(ggplot2)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


simulation_fun<-function(Sp, Ntr, rval,nsamples=500, Ktrue, it=1000, burn=500,type="GJAM"){
  S<-Sp
  n<- nsamples
  r <- rval
  iterations<-it
  env<-runif(-50,50,n=n)
  X<-cbind(1,poly(env,2)) #nxK
  idx<-sample(S)
  B_0<-seq(0,100,length.out=S)[idx]
  B_1<-seq(0,100,length.out=S)[idx]
  B_2<-seq(0,100,length.out=S)[idx]
  B<-cbind(B_0,B_1,B_2) #SxK
  L<-X%*%t(B) #nxS
  
  K=sum(S/(S+(1:S)-1)) #104, his prior number of clusters when alpha=S
  cat("Prior expected number of clusters : ",K,"\n")
  K_t= Ktrue
  cat("True number of clusters : ",K_t)
  A<-matrix(NA,nrow=ceiling(K_t),ncol=r)
  sig=matrix(runif(n=r*r),ncol=r)
  for(i in 1:ceiling(K_t)){
    A[i,]<-mvrnorm(n = 1, rep(20,r), Sigma=3*diag(r)) #Nxr short and skinny
  }
  idx<-sample((1:ceiling(K_t)),S,replace=T) 
  Lambda<-A[idx,] #Sxr tall and skinny
  Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  Y<-mvrnorm(n = 500, mu=rep(0,S), Sigma=Sigma)
  
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env1","env2")
  formula<-as.formula(~env1+env2)
  if(type=="GJAM"){
  rl <- list(r = r, N =Ntr-1)
  ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
  fit<-gjam(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
  alpha.chains<-NULL
  alpha.DP<-S
  }
  if(type=="0"){
    func<-function(x) {sum(x/(x+(1:S)-1))-K_t}
    alpha.DP<-.bisec(func,0.01,100)
    rl <- list(r = r, N =Ntr-1,alpha.DP=alpha.DP)
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
    fit<-.gjam0(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
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
  }
  if(type=="2"){
    func<-function(alpha) {sum(x/(x+(1:S)-1))-K_t}
    alpha.DP<-.bisec(func,0.01,100)
    shape=((alpha.DP)^2)/20
    rate=alpha.DP/20
    rl  <- list(r = r, N = Ntr-1, rate=rate,shape=shape,V=1)
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
    fit<-.gjam_2(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
  }
  if(type=="3"){
    eps=0.1
      
    alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    alp_sig$sigma=seq(0.05,0.5,length.out = 20)
    #loop to run bisecetion on a grid for sigma
    for(i in 1:20){
    func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S))/(x+c(1:S)))-1) - K_t}
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
  }
  if(type=="4"){
    
    eps=0.1
    
    alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
    colnames(alp_sig)<-c("alpha","sigma","is_less_150")
    alp_sig$sigma=seq(0.05,0.4,length.out = 20)
    #loop to run bisecetion on a grid for sigma
    for(i in 1:20){
      func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S))/(x+c(1:S)))-1) - K_t}
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
    alpha.chains<-fit$chains$alpha.DP_g
  }
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  df<-as.data.frame(trace)
  df$iter<-1:it
  #plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
  p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
    labs(title=paste0("Trace plot for the number of groups K for S=",S," r=",r," true K=",K_t), caption=paste0("Number of iterations ",it," burnin ",burn))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = Ktrue,color = "red")
  plot(p)
  pdf(paste0("Plot-clustersS",S,"r",r,"Ktr",K_t,".pdf"))
  plot(p)
  dev.off()
  
  
  N_dim<-(it-burn)
  Z<-array(dim=c(Sp,r,N_dim))
  for(j in 1:N_dim){
    K<-fit$chains$kgibbs[j,]
    Z[,,j]  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)[K,]
    #sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
  }
  Lambda_mean<-apply(Z,c(1,2),mean)
  err<-sum((Lambda_mean-Lambda)^2)/(S*r)
  fit<-fit$rmspeAll
  
  return(list(trace=trace, chain=fit$chains$kgibbs,
              coeff_t=Lambda,coeff_f=Lambda_mean,
              err=err,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,fit=fit))
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
list=list()
S_vec<-c(50,100,150,200)
r_vec<-c(3,5,10,15)
k<-1
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    list<-list.append(list,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T_0"),simulation_fun(Sp=S_vec[i], Ntr=150, rval=r_vec[j],nsamples=500, Ktrue=4,it=1000,burn=100,type="GJAM")))
    names(list)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_Kt_4_T0")
     k=k+1
    }
}

table<-data.frame()
for(i in 1:length(list)){
  str<-names(list)[[i]]
  tmp<-data.frame("trace"=list[[i]]$trace,"S"=rep(substr(str,regexpr("S",str)+2 , regexpr("_r",str)-1),length(list[[i]]$trace)),
                  "r"=rep(substr(str,regexpr("r",str)+2, regexpr("_N",str)-1),length(list[[i]]$trace)),
                  "N"=rep(substr(str,regexpr("N",str)+2 , regexpr("_n",str)-1),length(list[[i]]$trace)),
                  "n"=rep(substr(str,regexpr("n",str)+2 , regexpr("_Kt",str)-1),length(list[[i]]$trace)),
                  "Kt"=rep(substr(str,regexpr("Kt_",str)+3,nchar(str)) ,length(list[[i]]$trace)),
                  "It"=rep(length(list[[i]]$trace),length(list[[i]]$trace)),
                  "x"=1:length(list[[i]]$trace),
                  "err"=rep(list[[i]]$err,length(list[[i]]$trace)),
                  "fit"=rep(list[[i]]$fit,length(list[[i]]$trace)))
  table<-rbind(table,tmp)
}






table_r_3<-table[which(table$r==3),]
p_3<-ggplot(data=table_r_3,aes(x=table_r_3$x,y=table_r_3$trace,color=table_r_3$S))+geom_line()+
  labs(title="r=3, K_t=4")
table_r_5<-table[which(table$r==5),]
p_5<-ggplot(data=table_r_5,aes(x=table_r_5$x,y=table_r_5$trace,color=table_r_5$S))+geom_line()+
  labs(title="r=5, K_t=4")
table_r_10<-table[which(table$r==10),]
p_10<-ggplot(data=table_r_10,aes(x=table_r_10$x,y=table_r_10$trace,color=table_r_10$S))+geom_line()+
  labs(title="r=10, K_t=4")
table_r_15<-table[which(table$r==15),]
p_15<-ggplot(data=table_r_15,aes(x=table_r_10$x,y=table_r_15$trace,color=table_r_15$S))+geom_line()+
  labs(title="r=15, K_t=4")
grid.arrange(p_3,p_5,p_10,p_15)





