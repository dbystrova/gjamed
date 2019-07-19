rm(list=ls())

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

## delete cites with NA/0  for environment

NAs_values<- is.na(PA_env_df$bio_1_0)&is.na(PA_env_df$bio_12_0)&is.na(PA_env_df$bio_19_0)&is.na(PA_env_df$bio_8_0)&is.na(PA_env_df$slope)

PA_env_df_1<- PA_env_df[!NAs_values,]
zeros_values<- (PA_env_df_1$bio_1_0==0)&(PA_env_df_1$bio_12_0==0)&(PA_env_df_1$bio_19_0==0)&(PA_env_df_1$bio_8_0==0)&(PA_env_df_1$slope==0)
PA_env_df_2<- PA_env_df_1[!zeros_values,]


PA_non_miss_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,any),]
PA_non_na_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,all),]
PA_env_df_3<- PA_env_df_2[which(!(rowSums(is.na(PA_env_df_2[,7:131]))==125)),]




smp_size <- floor(0.30 * nrow(PA_env_df_3))
## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(PA_env_df_3)), size = smp_size)

train <- PA_env_df_3[train_ind, ]
test <- PA_env_df_3[-train_ind, ]

########################################################################Group numbers
Species_names_groups<- read.csv("PFG_Bauges_Description_2017.csv", sep="\t")

#############################################################################Fitting the model 

y<- train[,7:131]
xdata<- train[,2:5]
formula <- as.formula( ~ bio_1_0 + bio_12_0 + I(bio_1_0^2) + I(bio_12_0^2) )

Ydata  <- gjamTrimY(y,10)$y             # at least 10 plots

rl <- list(r = 8, N = 50)


ml   <- list(ng = 1000, burnin = 500, typeNames = 'PA', reductList = rl) #change ml

fit<-gjam(formula, xdata = xdata, ydata = Ydata, modelList = ml)





trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:1000
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 19,color = "red")
p




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




