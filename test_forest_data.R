rm(list=ls())
setwd("~/GitHub/gjamed")
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(ggplot2)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
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

rl <- list(r = 8, N = S)
rl1 <- list(r = 8, N = S, rate=0.1,shape=0.1)
rl2  <- list(r = 8, N = 20,rate=0.1,shape=0.1,V=1) #here to modify N


N_eps<-floor(.compute_tau_mean(0.3,10,0.1) + 2*.compute_tau_var(0.3,10,0.1)) #alpha=10,sigma=0.3
rl3   <- list(r = 8, N = N_eps, sigma_py=0.3, alpha=2)


N_eps<-floor(.compute_tau_mean(0.5,10,0.1) + 2*.compute_tau_var(0.5,10,0.1))
rl4   <- list(r = 8, N = N_eps,rate=0.1,shape=0.1,V1=1,ro.disc=0.5) #here to modify N

ml4   <- list(ng = 5000, burnin = 500, typeNames = 'DA', reductList = rl4) #change ml
ml3   <- list(ng = 5000, burnin = 500, typeNames = 'DA', reductList = rl3) #change ml
ml2   <- list(ng = 5000, burnin = 500, typeNames = 'DA', reductList = rl2) #change ml
ml1   <- list(ng = 5000, burnin = 500, typeNames = 'DA', reductList = rl1) #change ml
ml   <- list(ng = 5000, burnin = 500, typeNames = 'DA', reductList = rl) #change ml


fit<-gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit1<-.gjam_1(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
fit2<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
fit3 <- .gjam_3(form,xdata,treeYdata,ml3)
fit4<-.gjam_4(form, xdata = xdata, ydata = treeYdata, modelList = ml4)

fit$fit$rmspeAll
fit1$fit$rmspeAll
fit2$fit$rmspeAll
fit3$fit$rmspeAll
fit4$fit$rmspeAll

#check that alpha and sigma respected the truncation bounds
alpha<-mcmc(fit4$chains$alpha.PY_g)
plot(alpha)
acfplot(alpha)
cumuplot(alpha)

discount<-mcmc(fit4$chains$discount.PY_g)
plot(discount)
acfplot(discount)
cumuplot(discount)

#check the convergence
gjam_mc<- mcmc(fit$chains$sgibbs) #thinned sigma
s2s1<- mcmc(fit$chains$sgibbs[,2]) #s2s1 chain

#acf
#acfplot(gjam_mc)  ##autocor plot
ggs_autocorrelation(ggs(gjam_mc)) ###autocr

#traceplots
x11()
plot(gjam_mc) #traceplots

#nESS crazy low
xyplot(s2s1)
hist(effectiveSize(gjam_mc), main="ess(sigma)",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc_nd), main="ess(sigma) only non diagonal terms",lwd=2,col=gray(.6))
#cumuplot
cumuplot(s2s1)
x11()
cumuplot(gjam_mc)
dev.off()

beta_mcmc<-mcmc(mod_gjam_low$chains$bgibbs)
#acf
ggs_autocorrelation(ggs(beta_mcmc)) ###autocr
#nESS
hist(effectiveSize(beta_mcmc), main="ess(beta)",lwd=2,col=gray(.6))
#traceplots
x11()
plot(beta_mcmc)
dev.off()
#cumplots
x11()
cumuplot(beta_mcmc)
dev.off()



#check the traceplots of K
trace0<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
trace1<-apply(fit1$chains$kgibbs,1,function(x) length(unique(x)))
trace2<-apply(fit2$chains$kgibbs,1,function(x) length(unique(x)))
trace3<-apply(fit3$chains$kgibbs,1,function(x) length(unique(x)))
trace4<-apply(fit4$chains$kgibbs,1,function(x) length(unique(x)))

table<-data.frame()
table<-data.frame("trace"=c(trace0,trace1,trace2,trace3,trace4),
                  "type"=c(rep("0",length(trace0)),rep("1",length(trace1)),rep("2",length(trace2)),rep("3",length(trace3)),rep("4",length(trace4))))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(4)

p1<-ggplot(table[which(table$type=="0"),], aes(x=table$x[which(table$type=="0")],y=table$trace[which(table$type=="0")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p1
p2<-ggplot(table[which(table$type=="1"),], aes(x=table$x[which(table$type=="1")],y=table$trace[which(table$type=="1")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p2
p3<-ggplot(table[which(table$type=="2"),], aes(x=table$x[which(table$type=="2")],y=table$trace[which(table$type=="2")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p3
p4<-ggplot(table[which(table$type=="3"),], aes(x=table$x[which(table$type=="3")],y=table$trace[which(table$type=="3")]))+geom_point()+geom_hline(yintercept = Ktr,col="red")
p4

p<-ggplot(table, aes(x=x,y=trace,col=as.factor(type)))+geom_point()+
  scale_color_manual(name = c(""), values = cols, labels=c("Original model","DP with prior on alpha","PY with fixed alpha, sigma","PY with prior on alpha, sigma"))
p


#check the last weight
pk_chains0<- fit$chains$pk_g
pk_chains1- fit1$chains$pk_g
pk_chains2<- fit2$chains$pk_g
pk_chains3<- fit3$chains$pk_g
pk_chains4<- fit4$chains$pk_g
