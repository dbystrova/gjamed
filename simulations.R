########## Simulating the data
rm(list=ls())
setwd("~/Tesi/Code/modified_gjam/Gjam/")
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


S<-150
n=500
r=3


env<-runif(-50,50,n=n)
X<-cbind(1,poly(env,2)) #nxK

idx<-sample(S)
B_0<-seq(0,100,length.out=S)[idx]
B_1<-seq(0,100,length.out=S)[idx]
B_2<-seq(0,100,length.out=S)[idx]
B<-cbind(B_0,B_1,B_2) #SxK

L<-X%*%t(B) #nxS

K=sum(S/(S+(1:S)-1)) #104, his prior number of clusters when alpha=S

K_t=

A<-matrix(NA,nrow=ceiling(K_t),ncol=r)
sig=matrix(runif(n=r*r),ncol=r)
for(i in 1:ceiling(K_t)){
  A[i,]<-mvrnorm(n = 1, rep(0,r), Sigma=3*diag(r)) #Nxr short and skinny
}

# A<-matrix(NA,nrow=,ncol=r)
# sig=matrix(runif(n=r*r),ncol=r)
# for(i in 1:100){
#   A[i,]<-mvrnorm(n = 1, rep(0,r), Sigma=3*diag(r)) #Nxr short and skinny
# }

#idx<-sample((1:100),S,replace=T) 
idx<-sample((1:ceiling(K_t)),S,replace=T) 
Lambda<-A[idx,] #Sxr tall and skinny
Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS


#Y<-L+mvrnorm(n = 500, mu=rep(0,S), Sigma=Sigma)
Y<-mvrnorm(n = 500, mu=rep(0,S), Sigma=Sigma)

################# Fit Gjam

rl <- list(r = 3, N = 150-1)
xdata<-as.data.frame(X[,-1])
colnames(xdata)<-c("env1","env2")
ml<-list(ng=1000,burnin=200,typeNames='CA',reductList=rl)
formula<-as.formula(~env1+env2)
fit<-gjam(formula,xdata,ydata=as.data.frame(Y),modelList = ml)

plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))


####

func<-function(x,K) {sum(x/(x+(1:S)-1))-K}

#just to check how the function look like
# alpha<-c(1:100)
# for(i in 1:100) kmean[i]<-func(alpha[i])
# plot(alpha,kmean)

K=ceiling(K_t)
func<-function(x) {sum(x/(x+(1:S)-1))-K}

#choose alpha to have the right K
alpha.DP<-.bisec(func,0.5,100)

rl0<-list(r = 3, N = 150-1,alpha.DP=alpha.DP)
ml0<-list(ng=1000,burnin=200,typeNames='CA',reductList=rl0)
fit0<-.gjam0(formula,xdata,ydata=as.data.frame(Y),modelList = ml0)

plot(apply(fit0$chains$kgibbs,1,function(x) length(unique(x))))
