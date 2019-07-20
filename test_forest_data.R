rm(list=ls())
setwd("~/GitHub/gjamed")
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(ggmcmc)
library(arm)
library(coda)
library(Rcpp)
library(ggplot2)
library(biomod2)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]
y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots
form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )
S<-ncol(treeYdata) #95
it<-5000
burn<-500
rl <- list(r = 8, N = S,alpha.DP=S)
rl1 <- list(r = 8, N = S, rate=0.1,shape=0.1)
rl2  <- list(r = 8, N = S,rate=0.1,shape=0.1,V=5) #here to modify N


N_eps<-floor(.compute_tau_mean(0.3,10,0.1) + 2*.compute_tau_var(0.3,10,0.1)) #alpha=10,sigma=0.3
rl3   <- list(r = 8, N = N_eps, sigma_py=0.3, alpha=2)


N_eps<-floor(.compute_tau_mean(0.5,10,0.1) + 2*.compute_tau_var(0.5,10,0.1))
rl4   <- list(r = 8, N = N_eps,rate=0.1,shape=0.1,V1=5,ro.disc=0.5) #here to modify N

#ml4   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl4) #change ml
ml4   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl4, holdoutIndex = 1:100)
#ml3   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl3) #change ml
ml3   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl3, holdoutIndex = 1:100)
#ml2   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl2) #change ml
ml2   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl2, holdoutIndex = 1:100) 
#ml1   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl1) #change ml
ml1   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl1, holdoutIndex = 1:100) 
#ml   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl) #change ml
ml   <- list(ng = it, burnin = burn, typeNames = 'DA', reductList = rl, holdoutIndex = 1:100) 

fit<-.gjam0(form, xdata = xdata, ydata = treeYdata, modelList = ml)
#save(fit,file="models_forest_data/fit.Rda")
save(fit,file="models_forest_data_OSS/fit.Rda")

fit1<-.gjam_1(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
#save(fit1,file="models_forest_data/fit1.Rda")
save(fit1,file="models_forest_data_OSS/fit1.Rda")

fit2<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
#save(fit2,file="models_forest_data/fit2.Rda")
save(fit2,file="models_forest_data_OSS/fit2.Rda")

fit3 <- .gjam_3(form,xdata,treeYdata,ml3)
#save(fit3,file="models_forest_data/fit3.Rda")
save(fit3,file="models_forest_data_OSS/fit3.Rda")

fit4<-.gjam_4(form, xdata = xdata, ydata = treeYdata, modelList = ml4)
#save(fit4,file="models_forest_data/fit4.Rda")
save(fit4,file="models_forest_data_OSS/fit4.Rda")


fit$fit$rmspeAll  #2.257202
fit1$fit$rmspeAll #2.205223
fit2$fit$rmspeAll #2.209276
fit3$fit$rmspeAll #2.203517
fit4$fit$rmspeAll #2.092136

fit$fit$DIC   #2510192
fit1$fit$DIC  #2496028
fit2$fit$DIC  #2505273
fit3$fit$DIC  #2496243
fit4$fit$DIC  #2460125


#check that alpha and sigma posteriors
#gjam1  
alpha<-mcmc(fit1$chains$alpha.DP_g)
plot(alpha)
acfplot(alpha)
cumuplot(alpha)

#gjam2
alpha<-mcmc(fit2$chains$alpha.DP_g)
plot(alpha)
acfplot(alpha)
cumuplot(alpha)

##gjam4
alpha<-mcmc(fit4$chains$alpha.PY_g)
alpha<-mcmc(fit4$chains$alpha.PY_g[seq(1,length(fit4$chains$alpha.PY_g),by=20)])
plot(alpha)
acfplot(alpha)
cumuplot(alpha)

discount<-mcmc(fit4$chains$discount.PY_g)
plot(discount)
acfplot(discount)
cumuplot(discount)

#check the convergence

#for sigma
gjam_mc<- mcmc(fit$chains$sgibbs)
gjam_mc1<- mcmc(fit1$chains$sgibbs) 
gjam_mc2<- mcmc(fit2$chains$sgibbs) 
gjam_mc3<- mcmc(fit3$chains$sgibbs) 
gjam_mc4<- mcmc(fit4$chains$sgibbs) 


par(mfrow=c(2,3))
hist(effectiveSize(gjam_mc), main="ess(sigma) gjam",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc1), main="ess(sigma) gjam1",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc2), main="ess(sigma) gjam2",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc3), main="ess(sigma) gjam3",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc4), main="ess(sigma) gjam4",lwd=2,col=gray(.6),breaks=100)

# for betas
beta_mcmc<-mcmc(fit$chains$bgibbs)
beta_mcmc1<-mcmc(fit1$chains$bgibbs)
beta_mcmc2<-mcmc(fit2$chains$bgibbs)
beta_mcmc3<-mcmc(fit3$chains$bgibbs)
beta_mcmc4<-mcmc(fit4$chains$bgibbs)

#nESS
par(mfrow=c(2,3))
hist(effectiveSize(beta_mcmc), main="ess(beta) gjam",lwd=2,col=gray(.6))
hist(effectiveSize(beta_mcmc1), main="ess(beta) gjam1",lwd=2,col=gray(.6))
hist(effectiveSize(beta_mcmc2), main="ess(beta) gjam2",lwd=2,col=gray(.6))
hist(effectiveSize(beta_mcmc3), main="ess(beta) gjam3",lwd=2,col=gray(.6))
hist(effectiveSize(beta_mcmc4), main="ess(beta) gjam4",lwd=2,col=gray(.6))



#check the traceplots of K
trace0<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
trace1<-apply(fit1$chains$kgibbs,1,function(x) length(unique(x)))
trace2<-apply(fit2$chains$kgibbs,1,function(x) length(unique(x)))
trace3<-apply(fit3$chains$kgibbs,1,function(x) length(unique(x)))
trace4<-apply(fit4$chains$kgibbs,1,function(x) length(unique(x)))

table<-data.frame()
table<-data.frame("trace"=c(trace0,
                            #trace1,
                            trace2,trace3,trace4),
                  "type"=c(rep("0",length(trace0)),
                           #rep("1",length(trace1)),
                           rep("2",length(trace2)),rep("3",length(trace3)),rep("4",length(trace4))),
                  "x"=rep(1:it,4))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(4)

#single traceplots - not useful
# p1<-ggplot(table[which(table$type=="0"),], aes(x=table$x[which(table$type=="0")],y=table$trace[which(table$type=="0")]))+geom_point()
# p1
# p2<-ggplot(table[which(table$type=="1"),], aes(x=table$x[which(table$type=="1")],y=table$trace[which(table$type=="1")]))+geom_point()
# p2
# p3<-ggplot(table[which(table$type=="2"),], aes(x=table$x[which(table$type=="2")],y=table$trace[which(table$type=="2")]))+geom_point()
# p3
# p4<-ggplot(table[which(table$type=="3"),], aes(x=table$x[which(table$type=="3")],y=table$trace[which(table$type=="3")]))+geom_point()
# p4

# traceplots altogether
p<-ggplot(table, aes(x=x,y=trace,col=as.factor(type)))+geom_point()+
  scale_color_manual(name = c(""), values = cols, labels=c("Original model",
                                                           #"DP with prior on alpha 1",
                                                           "DP with prior on alpha 2","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+theme_bw()
pdf("plot_forest_data/forest_data_trace_K.pdf")
p
dev.off()

#check the last weight
pk_chains0_last<- mcmc(fit$chains$pk_g[,ncol(fit$chains$pk_g)])
plot(pk_chains1_last)
pk_chains1_last<- mcmc(fit1$chains$pk_g[,ncol(fit1$chains$pk_g)])
plot(pk_chains1_last)
pk_chains2_last<- mcmc(fit2$chains$pk_g[,ncol(fit2$chains$pk_g)])
plot(pk_chains2_last)
pk_chains3_last<- mcmc(fit3$chains$pk_g[,ncol(fit3$chains$pk_g)])
plot(pk_chains3_last)
pk_chains4_last<- mcmc(fit4$chains$pk_g[,ncol(fit4$chains$pk_g)])
plot(pk_chains4_last)

#Out of sample prediction
sum((fit$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5987808
sum((fit1$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5961853
sum((fit2$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5979932
sum((fit3$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5977189
sum((fit4$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5983279

