#rm(list=ls())
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(raster)
library(ggplot2)
library(rgdal)
library(biomod2)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


#setwd("~/Downloads/RFate-master/data_supplements/Bauges")


B_coords_xy<- load_object("DB.XY.RData")
#B__observations_XY<- load_object("DB.observations.xy.RData")

#load abundance data
PA_data<-load_object("DOM.mat.sites.species.PA.RData")
PA_data_df<- as.data.frame(PA_data)
PA_data_df$cite<- rownames(PA_data)

#raster stack for the 100.tif 
zone.name="ENV_VARIABLES"
zone.env.folder="EOBS_1970_2005"
zone.env.variables=c("bio_1_0","bio_12_0","bio_19_0","bio_8_0","slope")
env.files = list.files(path = paste0(zone.name, "/", zone.env.folder)
                       , pattern = paste0(paste0(zone.env.variables, ".img", collapse = "|"), "|", paste0(zone.env.variables, ".tif", collapse = "|")), full.names = TRUE)


##function from the package
getSDM_env = function(zone.name, zone.env.folder, zone.env.variables, maskSimul)
{
  env.files = list.files(path = paste0(zone.name, "/", zone.env.folder)
                         , pattern = paste0(paste0(zone.env.variables, ".img", collapse = "|")
                                            , "|"
                                            , paste0(zone.env.variables, ".tif", collapse = "|"))
                         , full.names = TRUE)
  zone.env.stk.CALIB = raster::stack(env.files)
  #maskSimul=raster("MASK_100m.tif")
  origin(maskSimul) = origin(zone.env.stk.CALIB)
  zone.env.stk.PROJ = stack(zone.env.stk.CALIB * maskSimul)
  names(zone.env.stk.PROJ) = names(zone.env.stk.CALIB)
  
  return(list(env.CALIB = zone.env.stk.CALIB
              , env.PROJ = zone.env.stk.PROJ))
}


new_B_env<-getSDM_env(zone.name, zone.env.folder, zone.env.variables, maskSimul=raster("MASK_100m.tif")) 
#B_ENV_proj<- load_object("Bauges.zone.env.stk.PROJ.RData")
B_env_raster<- new_B_env$env.PROJ

B_env<-as.data.frame(extract(B_env_raster, B_coords_xy))
B_env$cite<- rownames(B_coords_xy)


### merge environmental covariates and presence/abscence data by cite.
PA_env_df <- merge(B_env,PA_data_df,by="cite")
## delete cites with NA for environment
NAs_values<- is.na(PA_env_df$bio_1_0)&is.na(PA_env_df$bio_12_0)&is.na(PA_env_df$bio_19_0)&is.na(PA_env_df$bio_8_0)&is.na(PA_env_df$slope)
PA_env_df_1<- PA_env_df[!NAs_values,]
###  delete total 0 for environment
zeros_values<- (PA_env_df_1$bio_1_0==0)&(PA_env_df_1$bio_12_0==0)&(PA_env_df_1$bio_19_0==0)&(PA_env_df_1$bio_8_0==0)&(PA_env_df_1$slope==0)
PA_env_df_2<- PA_env_df_1[!zeros_values,]

## non missing data sets
#PA_non_miss_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,any),]
#PA_non_na_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,all),]

#### Keeping PA data with at least one 1/0
PA_env_df_3<- PA_env_df_2[which(!(rowSums(is.na(PA_env_df_2[,7:131]))==125)),]

PA_non_miss_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,any),]

####Separation test/train
#data<- PA_env_df_3
data<- PA_non_miss_env
smp_size <- floor(0.70 * nrow(data))
## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
## dim(train)  / 3712  131
## dim(test)  /  1591  131
########################################################################Group numbers
Species_names_groups<- read.csv("PFG_Bauges_Description_2017.csv", sep="\t")
# K=16 functional groups
#############################################################################Fitting the model 

it<-1000
burn<-500
holdout<- sample(seq_len(nrow(train)), size = 100)

###########GJAM standart model
y<- train[,7:131]
xdata<- train[,2:6]
formula <- as.formula( ~ bio_1_0+ bio_19_0+bio_8_0+ slope)
Ydata  <- gjamTrimY(y,50)$y             # at least 10 plots - re-group rare species
S<- ncol(Ydata)
rl <- list(r = 5, N = 50)
ml   <- list(ng = 100, burnin = 10, typeNames = 'PA', reductList = rl) #change ml
fit<-gjam(formula, xdata = xdata, ydata = Ydata, modelList = ml)
#save(fit,file="models_Bagues_data_OSS/fit.Rda")
#no Holdout
#save(fit,file="models_forest_data_OSS/fit.Rda")


####### Out of sample prediction  -DOESN't work : chol() error ??
y_test<- test[,7:131]
xdata_test<- test[,2:6]

new <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
p1  <- gjamPredict(output = fit, newdata = new)
plot(y_test, p1$sdList$yMu ,ylab = 'Predicted',cex=.1)
abline(0,1)

############################################
####Check the trace for number of groups. From the previous analysis we know that the number of functional groups is 16

trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:1000
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 16,color = "red")
p

####The Dirichlet process prior with conjugation############################################


##We don't need this model
rl1 <- list(r = 5, N = 50,rate=10,shape=10)
fit1<-.gjam_1(form, xdata = xdata, ydata = Ydata, modelList = ml1, holdoutIndex =holdout)





####The Dirichlet process prior : multinomial  + prior on alpha
##Computing the hyper-parameters for K=16
K=16
S=125
func<-function(x) {sum(x/(x+(1:S)-1))-16}
alpha.DP<-.bisec(func,0.01,100)
shape=((alpha.DP)^2)/20
rate=alpha.DP/20


rl2  <- list(r = 5, N = 50,rate=rate,shape=shape,V=1) #here to modify N
ml2   <- list(ng = 1000, burnin = 500, typeNames = 'PA', reductList = rl2) #change ml

fit2<-.gjam_2(formula, xdata = xdata, ydata = Ydata, modelList = ml2)
save(fit2,file="models_Bagues_data_OSS/fit2.Rda")


trace<-apply(fit2$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:1000
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 16,color = "red")
p

#################################################################################PY
K=16
eps<-0.1
sigma_py<-0.25
funcPY_root<-function(x) {(x/sigma_py)*(prod((x+sigma_py+c(1:S) -1)/(x+c(1:S) -1))-1) - K}
alpha.PY<-.bisec(funcPY_root,0.0001,100)
N_eps<-floor(.compute_tau_mean_large_dim(sigma_py,alpha.PY,eps) + 2*.compute_tau_var_large_dim(sigma_py,alpha.PY,eps))
rl3   <- list(r = 5, N = N_eps, sigma_py=sigma_py, alpha=alpha.PY)
ml3   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl3)


fit3 <- .gjam_3(formula,xdata,Ydata,ml3)



save(fit3,file="models_Bagues_data_OSS/fit3.Rda")





K<-16
eps=0.1
alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
colnames(alp_sig)<-c("alpha","sigma","is_less_150")
alp_sig$sigma=seq(0.05,0.4,length.out = 20)
#loop to run bisecetion on a grid for sigma
for(i in 1:20){
  ####corrected added  -1
  func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S)-1)/(x+c(1:S) -1))-1) - K}
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

rl4   <- list(r = 5, N =N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
ml4   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl4)

fit4<-.gjam_4(formula, xdata = xdata, ydata = Ydata, modelList = ml4)


save(fit4,file="models_Bagues_data_OSS/fit4.Rda")




########################################## Model comparison



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








#################################################################################Other models


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
#trace1<-apply(fit1$chains$kgibbs,1,function(x) length(unique(x)))
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
#pdf("plot_forest_data/forest_data_trace_K.pdf")
p
#dev.off()

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
#TJUR coefficient


sum((fit$prediction$ypredMu[1:100,]-Ydata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5987808
sum((fit1$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5961853
sum((fit2$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5979932
sum((fit3$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5977189
sum((fit4$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5983279



#################################################################################Other models - BIOMOD


## Not run:
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",package="biomod2"), row.names = 1)
# the presence/absences data for our species
myRespName <- 'VulpesVulpes'
myResp <- as.numeric(DataSpecies[,myRespName])
# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12) 
myExpl = stack( system.file( "external/bioclim/current/bio3.grd",package="biomod2"),
system.file( "external/bioclim/current/bio4.grd",package="biomod2"),
system.file( "external/bioclim/current/bio7.grd",package="biomod2"),
system.file( "external/bioclim/current/bio11.grd", package="biomod2"),
system.file( "external/bioclim/current/bio12.grd",package="biomod2"))
# 1. Formatting Data

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = myExpl,
                                    resp.xy = myRespXY,
                                    resp.name = myRespName)
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()
# 3. Doing Modelisation
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                     models = c('SRE'),
                                     models.options = myBiomodOption,
                                     NbRunEval=1,
                                     DataSplit=80,
                                     Prevalence=0.5,
                                     VarImport=0,
                                     models.eval.meth = c('TSS','ROC'),
                                     do.full.models=FALSE,
                                     modeling.id="test2")
# files have been created on hard drive
list.files(myRespName,all.files=TRUE,recursive=TRUE)
# remove properly the modeling objects and all the file saved on hard drive
RemoveProperly(myBiomodModelOut)
# check files had been removed
list.files(myRespName,all.files=TRUE,recursive=TRUE)
## End(Not run)

