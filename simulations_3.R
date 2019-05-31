#test 2
#setwd("~/Bureau/PY_comparisons/gjamed-master")
rm(list=ls())
setwd("~/Tesi/Code/modified_gjam/Gjam/")


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
library(ggplot2)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


#See section 5.2 of Clark.
# We repeat this simulations by simulating the rows of PSi as standard multivariate normal \in R^S
# And then Sigma= (Psi*t(Psi))^(-1) which is equivalent to Sigma \sim IW(I,S).
# So there are no clusters. However, the strong reduction of the number of parameters leads to a much better
# estimation of the covariance matrix. We want to check if there's any advantage in giving a vague prior to alpha
# for DP and alpha and sigma for PY.
# The mean structure is defined by an intercept and a single predictor.

simulation_fun<-function(Sp, Ntr,nsamples=500, it=1000, burn=500, type="GJAM"){
  S<-Sp
  n<- nsamples
  iterations<-it
  X<-cbind(rep(1,n),scale(runif(0,100,n=n))) #intercept and environment, scaled. nxk
  idx<-sample(S)
  B_0<-scale(seq(0,100,length.out=S)[idx]) #coefficients, scaled (we will take Y=1 if V>0 and viceversa,
                                 #so we don't want the mean to be too far away from 0 or we will have all 0s or 1s)
  B_1<-scale(seq(0,100,length.out=S)[idx])
  B<-cbind(B_0,B_1) #SxK
  L<-X%*%t(B) #nxS

  Psi<- rmvnormRcpp(n=1,mu=rep(0,S),sigma=diag(S))
  for(i in 2:S) Psi<-rbind(Psi, rmvnormRcpp(n=1,mu=rep(0,S),sigma=diag(S)))
  
  Sigma_true=solveRcpp(Psi%*%t(Psi)) #Now sigma is an Inverse Wishart
  
  
  Y_cont<-L+rmvnormRcpp(n = n, mu=rep(0,S), sigma=Sigma_true)
  Y<- ifelse(Y_cont>0,1,0)
  
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env")
  formula<-as.formula(~env)
  if(type=="GJAM"){ 
    # I use gjam0 with alpha=S bc there's only one predictor, so problems with colMeans(chain) 
    rl <- list(r = r, N =Ntr-1,alpha.DP=S)
    ml<-list(ng=it,burnin=burn,typeNames='PA',reductList=rl)
    fit<-.gjam0(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-NULL
  }
  

  if(type=="1"){
    #Mean 1 variance 50, "vague" prior [my prof does not like a=b=0.001]
    shape=0.02
    rate=0.02
    rl  <- list(r = r, N = Ntr-1, rate=rate,shape=shape)
    ml<-list(ng=it,burnin=burn,typeNames='PA',reductList=rl)
    fit<-.gjam_1(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
  }
  if(type=="2"){
    #Mean 1 variance 50, "vague" prior [my prof does not like a=b=0.001]
    
    shape=0.02
    rate=0.02
    rl  <- list(r = r, N = Ntr-1, rate=rate,shape=shape,V=1)
    ml<-list(ng=it,burnin=burn,typeNames='PA',reductList=rl)
    fit<-.gjam_2(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
  }

  if(type=="4"){
    
    eps=0.1
    
    #fixing hyperparameters
    ro.disc=0.5 #equal prob to be a DP
    #Mean 1 variance 50, "vague" prior [my prof does not like a=b=0.001]
    shape=0.02
    rate=0.02
    # 95% quantile of alpha
    alpha.max=qgamma(.95, shape=shape, rate=rate)
    
    N_eps<-floor(.compute_tau_mean(sigma_py,alpha.max,eps) + 2*.compute_tau_var(sigma_py,alpha.max,eps))
    
    
    rl   <- list(r = r, N = N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
    
    ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
    fit<-.gjam_4(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
    alpha.chains<-fit$chains$alpha.DP_g
  }
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  
  N_dim<-(it-burn)
  sigma<-array(dim=c(Sp,Sp,N_dim))
  for(j in 1:N_dim){
    K<-fit$chains$kgibbs[j,]
    Z  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)
    sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
    
  }
  sigma_mean<-apply(sigma,c(1,2),mean)
  err<-sum((sigma_mean-Sigma_true)^2)/(S*S)
  rmspe<-fit$fit$rmspeAll
  
  #out-of-sample on species on 300 sites, 30% left of species then 10% left out species
  X_out<-cbind(rep(1,300),scale(runif(0,100,n=300))) #intercept and environment, scaled. nxk
  L_out<-X_out%*%t(B) #nxS
  Y_cont_out<-L_out+rmvnormRcpp(n = 300, mu=rep(0,S), sigma=Sigma_true)
  Y_out<- ifelse(Y_cont_out>0,1,0)
  xdata=as.data.frame(X_out[,-1])
  colnames(xdata)<-c("env")
  newdata_0.3   <- list(xdata=xdata,ydataCond = Y_out[,c(1:floor(S*0.7))])
  #newdata   <- list(ydataCond = Y[,c(1:floor(S*0.7))])
  pred_0.3        <- gjamPredict(output = fit, newdata = newdata_0.3)
  newdata_0.9   <- list(xdata=xdata,ydataCond = Y_out[,c(1:floor(S*0.9))])
  #newdata   <- list(ydataCond = Y[,c(1:floor(S*0.7))])
  pred_0.9        <- gjamPredict(output = fit, newdata = newdata_0.9)
  
  #tjur
  
  
  
  return(list(trace=trace, chain=fit$chains$kgibbs,
              coeff_t=Sigma_true,coeff_f=sigma,
              err=err,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,fit=rmspe,
              pred_0.9_err=,
              pred_0.3_err=
              ))
}

####### Just one possible test case
#sim<-simulation_fun(Sp=50, Ntr=150, rval=3,nsamples=500,it=1000,burn=200,type=2)
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
list=list()
S_vec<-
r_vec<-
k<-1
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    list<-list.append(list,assign(paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_T_0"),simulation_fun(Sp=S_vec[i], Ntr=150, rval=r_vec[j],nsamples=500,it=1000,burn=100,type="GJAM")))
    names(list)[[k]]<-paste0("S_",S_vec[i],"_r_",r_vec[j],"_N_150_n_500_T0")
    k=k+1
  }
}
