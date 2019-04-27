
setwd("~/Tesi/Code/modified_gjam/Gjam")
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
source("R/simple_gjam_4.R")



#d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]


formula <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots

rl1 <- list(r = 8, N = 20,rate=10,shape=10)
rl2  <- list(r = 8, N = 20,rate=10,shape=10,V=5) #here to modify N
rl4   <- list(r = 8, N = 20,rate=10,shape=10,V1=5,V2=1,ro.disc=0.5) #here to modify N
N_eps<-floor(.compute_tau_mean(0.3,2,0.1) + 2*.compute_tau_var(0.3,2,0.1))
rl3   <- list(r = 8, N = N_eps, sigma_py=0.3, alpha=2)
ml   <- list(ng = 1000, burnin = 500, typeNames = 'DA', reductList = rl1) #change ml

form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

fit<-.gjam_4(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit<-.gjam_2(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit<-.gjam_1(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit <- .gjam_3(form,xdata,treeYdata,ml)


alpa<-fit$chains$alpha.PY_g[seq(from=1,to=length(fit$chains$alpha.PY_g),by=20)]
alpha<-mcmc(alpha)
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
thin(alpha)

discount<-mcmc(fit$chains$discount.PY_g)
plot(discount)
acfplot(discount)
cumuplot(discount)






.compute_tau_mean

