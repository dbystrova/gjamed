rm(list=ls())
#setwd("~/Tesi/Code/modified_gjam/Gjam/")
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]

set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)



#train <- data[train_ind, ]
#test <- data[-train_ind, ]
## dim(train)  / 3712  131
## dim(test)  /  1591  131


formula <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots

rl <- list(r = 8, N = 20)
rl1 <- list(r = 8, N = 20,rate=10,shape=10)
rl2  <- list(r = 8, N = 20,rate=10,shape=10,V=1) #here to modify N
N_eps<-floor(.compute_tau_mean(0.3,2,0.1) + 2*.compute_tau_var(0.3,2,0.1))
rl3   <- list(r = 8, N = N_eps, sigma_py=0.3, alpha=2)
N_eps<-floor(.compute_tau_mean(0.5,10,0.1) + 2*.compute_tau_var(0.5,10,0.1))
rl4   <- list(r = 8, N = N_eps,rate=10,shape=10,V1=1,ro.disc=0.5) #here to modify N

ml4   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl4) #change ml
ml3   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl3) #change ml
ml2   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl2) #change ml
ml1   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl1) #change ml
ml   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl) #change ml

form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

fit<-gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit1<-.gjam_1(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
fit2<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
fit3 <- .gjam_3(form,xdata,treeYdata,ml3)
fit4<-.gjam_4(form, xdata = xdata, ydata = treeYdata, modelList = ml4)









fit$fit$rmspeAll
fit4$fit$rmspeAll

alpha<-fit$chains$alpha.PY_g[seq(from=200,to=length(fit$chains$alpha.PY_g),by=20)]
alpha<-mcmc(alpha)
plot(alpha)
acfplot(alpha)
cumuplot(alpha)

discount<-mcmc(fit$chains$discount.PY_g)
plot(discount)
acfplot(discount)
cumuplot(discount)

rl0  <- list(r = 8, N = 20,alpha.DP=) #here to modify N
rl1  <- list(r = 8, N = 20,rate=10,shape=10,V=0.01) #here to modify N
rl2  <- list(r = 8, N = 20,rate=10,shape=10,V=0.1) #here to modify N
rl3  <- list(r = 8, N = 20,rate=10,shape=10,V=1) #here to modify N
rl4  <- list(r = 8, N = 20,rate=10,shape=10,V=10) #here to modify N
rl5  <- list(r = 8, N = 20,rate=10,shape=10,V=100) #here to modify N

ml1   <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl0) #change ml
ml1   <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl1) #change ml
ml2  <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl2) #change ml
ml3   <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl3) #change ml
ml4   <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl4) #change ml
ml5   <- list(ng = 2000, burnin = 500, typeNames = 'DA', reductList = rl5) #change ml

fit0<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml0)
fit1<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
fit2<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
fit3<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml3)
fit4<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml4)
fit5<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml5)

alpha1<-fit1$chains$alpha.DP_g[seq(from=200,to=length(fit1$chains$alpha.DP_g),by=20)]
alpha1<-mcmc(alpha1)
plot(alpha1)

alpha2<-fit2$chains$alpha.DP_g[seq(from=200,to=length(fit2$chains$alpha.DP_g),by=20)]
alpha2<-mcmc(alpha2)
plot(alpha2)

alpha3<-fit3$chains$alpha.DP_g[seq(from=200,to=length(fit3$chains$alpha.DP_g),by=20)]
alpha3<-mcmc(alpha3)
plot(alpha3)

alpha4<-fit4$chains$alpha.DP_g[seq(from=200,to=length(fit4$chains$alpha.DP_g),by=20)]
alpha4<-mcmc(alpha4)
plot(alpha4)

alpha5<-fit5$chains$alpha.DP_g[seq(from=200,to=length(fit5$chains$alpha.DP_g),by=20)]
alpha5<-mcmc(alpha5)
plot(alpha5)





trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:1000
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 19,color = "red")
p


