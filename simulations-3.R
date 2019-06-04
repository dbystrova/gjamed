setwd("~/Documents/GitHub/gjamed")
########## Simulating the data
rm(list=ls())
library(MASS)
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(NLRoot)
library(Rcpp)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")




simulation_fun<-function(Sp, Ntr, rval,nsamples=500, Ktrue, it=1000, burn=500){
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
    A[i,]<-mvrnorm(n = 1, rep(0,r), Sigma=3*diag(r)) #Nxr short and skinny
  }
  idx<-sample((1:ceiling(K_t)),S,replace=T) 
  Lambda<-A[idx,] #Sxr tall and skinny
  Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  Y<-mvrnorm(n = 500, mu=rep(0,S), Sigma=Sigma)
  
  rl <- list(r = r, N =Ntr-1)
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env1","env2")
  ml<-list(ng=it,burnin=burn,typeNames='CA',reductList=rl)
  formula<-as.formula(~env1+env2)
  fit<-gjam(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  df<-as.data.frame(trace)
  df$iter<-1:it
  #plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
  p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + geom_hline(yintercept=K_t,color="red")+
    labs(title=paste0("Trace plot for the number of groups K for S=",S,", r=",r,", true K=",K_t),caption=paste0("Number of iterations ",it,", burnin ",burn, ", Prior expected number of clusters: ",round(K,2)))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10),strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))

  plot(p)
  
  pdf(paste0("Plot-clustersS",S,"r",r,"Ktr",K_t,".pdf"))
  plot(p)
  dev.off()
  return(list(chains=fit$chains$kgibbs,trace=trace))
  
}

#######Simulation case1: S-500 , r=3, Ktrue=8
list=list()
S_vec<-c(20,50,100,150)
r_vec<-c(3,5,10,15)
table<-as.data.frame(matrix(NA,ncol=3))
colnames(table)<-c("trace","S","r")
k<-1
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    list<-list.append(list,assign(paste0("trace_S_",S_vec[i],"r",r_vec[j],"_N_150_n_500_Kt_4"),simulation_fun(Sp=S_vec[i], Ntr=S_vec[i], rval=r_vec[j],nsamples=500, Ktrue=4,it=2000,burn=500)))
    names(list)[[k]]<-paste0("trace_S_",S_vec[i],"r",r_vec[j],"_N_150_n_500_Kt_4")
    k=k+1
  }
}


save(list, file="simulationclustk4NsSs.Rdata")


#######Simulation case1: S-500 , r=3, Ktrue=8
list=list()
S_vec<-c(20,50,100,150)
r_vec<-c(10,15)
table<-as.data.frame(matrix(NA,ncol=3))
colnames(table)<-c("trace","S","r")
k<-1
for(i in 1:length(S_vec)){
  for(j in 1:length(r_vec)){
    list<-list.append(list,assign(paste0("trace_S_",S_vec[i],"r",r_vec[j],"_N_150_n_500_Kt_4"),simulation_fun(Sp=S_vec[i], Ntr=S_vec[i], rval=r_vec[j],nsamples=500, Ktrue=4,it=10000,burn=3000)))
    names(list)[[k]]<-paste0("trace_S_",S_vec[i],"r",r_vec[j],"_N_150_n_500_Kt_4")
    k=k+1
  }
}


save(list, file="simulationclustk4NsSsR1015Niter.Rdata")




##########################################################################################





