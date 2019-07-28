#rm(list=ls())
#setwd("~/Documents/GitHub/gjamed")
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
library(AUC)
#library(formattable)
library(gridExtra)
library(grid)

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
AB_data<-load_object("DOM.mat.sites.species.abund.RData")


### Colnames PA vs AB
S_PA<- colnames(PA_data)
S_AB<- colnames(AB_data)

setdiff(S_PA,S_AB) #14797"


PA_data_df<- as.data.frame(PA_data)
PA_data_df$cite<- rownames(PA_data)

AB_data_df<- as.data.frame(AB_data)
AB_data_df$cite<- rownames(AB_data)


PA_AB_common_plots<- intersect(AB_data_df$cite, PA_data_df$cite)



#L1<- apply(!is.na(PA_data_df[,2:126]),1, which)
#L2<- apply(!is.na(AB_data_df[,2:126]),1, which)


#spdf <- SpatialPoints(B_coords_xy,proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80
#                                                    +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))


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

#X<-  mask(zone.env.stk.CALIB, maskSimul)


rasmask = raster("MASK_100m.tif")
# File to environmental data rasteres
file.env = "ENV_VARIABLES/EOBS_1970_2005"
# 1) Get all the environmental variables in a single raster stack. So that you can extract them all afterwards.
env.bauges = stack(list.files(file.env,full.names = T))
# 1.a) set the projection same as the MASK
env.bauges = projectRaster(env.bauges, rasmask, res=100, method = "ngb")
# 1.b) Get only the values within MASK
env.bauges = mask(env.bauges, rasmask, maskvalue=0)
B_env_new<-as.data.frame(extract(env.bauges, B_coords_xy))
B_env_new$cite<- rownames(B_coords_xy)


#pdf("Bagues_10k_r5_bio12.pdf")
new_B_env<-getSDM_env(zone.name, zone.env.folder, zone.env.variables, maskSimul=raster("MASK_100m.tif")) 
B_ENV_proj<-load_object("Bauges.zone.env.stk.PROJ.RData")
B_env_raster<- new_B_env$env.PROJ
B_env<-as.data.frame(extract(B_env_raster, B_coords_xy))
B_env_true<-as.data.frame(extract(B_ENV_proj, B_coords_xy))
B_env$cite<- rownames(B_coords_xy)
#all(na.omit(B_env) == na.omit(B_env_true))

##Delete 0 values and NA's
NAs_values<- is.na(B_env$bio_1_0)&is.na(B_env$bio_12_0)&is.na(B_env$bio_19_0)&is.na(B_env$bio_8_0)&is.na(B_env$slope)
B_env_1<- B_env[!NAs_values,]
zeros_values<- (B_env_1$bio_1_0==0)&(B_env_1$bio_12_0==0)&(B_env_1$bio_19_0==0)&(B_env_1$bio_8_0==0)&(B_env_1$slope==0)
B_env_2<- B_env_1[!zeros_values,]

##Delete 0 values and NA's for B_env_new
B_env_new2<- B_env_new[,c("cite","bio_1_0","bio_12_0","bio_19_0","bio_8_0","slope")]
NAs_values<- is.na(B_env_new2$bio_1_0)&is.na(B_env_new2$bio_12_0)&is.na(B_env_new2$bio_19_0)&is.na(B_env_new2$bio_8_0)&is.na(B_env_new2$slope)
B_env_new3<- B_env_new2[!NAs_values,]
zeros_values<- (B_env_new3$bio_1_0==0)&(B_env_new3$bio_12_0==0)&(B_env_new3$bio_19_0==0)&(B_env_new3$bio_8_0==0)&(B_env_new3$slope==0)
B_env_new4<- B_env_new3[!zeros_values,]

####### Merge with presence_absence by cite
PA_env_df <- merge(B_env_2,PA_data_df,by="cite")
summary(PA_env_df)

AB_env_df <- merge(B_env_2,AB_data_df,by="cite")
summary(AB_env_df)


PA_env_df2 <- merge(B_env_new4,PA_data_df,by="cite")
summary(PA_env_df2)

AB_env_df2 <- merge(B_env_new4,AB_data_df,by="cite")
summary(AB_env_df2)
### merge environmental covariates and presence/abscence data by cite.
#PA_env_df <- merge(B_env,PA_data_df,by="cite")
## delete cites with NA for environment
#NAs_values<- is.na(PA_env_df$bio_1_0)&is.na(PA_env_df$bio_12_0)&is.na(PA_env_df$bio_19_0)&is.na(PA_env_df$bio_8_0)&is.na(PA_env_df$slope)
#PA_env_df_1<- PA_env_df[!NAs_values,]
###  delete total 0 for environment
#zeros_values<- (PA_env_df_1$bio_1_0==0)&(PA_env_df_1$bio_12_0==0)&(PA_env_df_1$bio_19_0==0)&(PA_env_df_1$bio_8_0==0)&(PA_env_df_1$slope==0)
#PA_env_df_2<- PA_env_df_1[!zeros_values,]

## non missing data sets
#PA_non_miss_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,any),]
#PA_non_na_env<- PA_env_df_2[!apply(is.na(PA_env_df_2[,7:131]),1,all),]

#### Keeping PA data with at least one 1/0
PA_env_df_3<- PA_env_df2[which(!(rowSums(is.na(PA_env_df2[,7:131]))==125)),]


#summary(PA_env_df_3)
# L1<- sapply(PA_env_df, function(x) sum(is.na(x)))
# formattable(summary(L1))
# L2<- sapply(AB_env_df, function(x) sum(is.na(x)))
# formattable(summary(L2))
###Convert from NA to 0
#PA_env_df_2_fil_o_y<- PA_env_df_3[,7:131]
#PA_env_df_2_fil_o_x<- PA_env_df_3[,1:6]
#PA_env_df_2_fil_o_y[is.na(PA_env_df_2_fil_o_y)] <- 0
#PA_env_df_2_0<- cbind(PA_env_df_2_fil_o_x,PA_env_df_2_fil_o_y)
#PA_env_df_2_0<-apply(PA_env_df_2[,7:131],1,function(x) x[is.na(x)] <- 0 )
#PA_env_df_2_00<- data.frame(PA_env_df_2[,1:6], apply(PA_env_df_2[,7:131],1, ) )
####Separation test/train
#data<- PA_env_df_3
#data<- PA_env_df_2_0
#smp_size <- floor(0.70 * nrow(data))
## set the seed to make your partition reproducible
#set.seed(123)
#train_ind <- sample(seq_len(nrow(data)), size = smp_size)
#train <- data[train_ind, ]
#test <- data[-train_ind, ]
## dim(train)  / 3712  131
## dim(test)  /  1591  131

env_data<- AB_env_df2[,1:6]
###duplicates
#duplicated_env<- env_data[,2:6] %>% duplicated()
#env_data_nodup<- env_data[!duplicated_env,]
###scaling
env_data_n<- as.data.frame(scale(env_data[,2:6]))
env_data_norm<- cbind(env_data[,1],env_data_n )
names(env_data_norm)<- c("cite", names(env_data_n))

#### Using abundance 
Ydat<- AB_env_df2[,c(1,7:130)]
Ydat_num<- Ydat[,2:125]
Ydat_num[Ydat_num> 0] <- 1
Ydat[,2:125]<- Ydat_num
AB_norm_PA<- merge(env_data_norm,Ydat,by="cite")
#summary(AB_norm_PA)

########################################################################Group numbers
Species_names_groups<- read.csv("PFG_Bauges_Description_2017.csv", sep="\t")

true_names<- as.data.frame(names(table(Species_names_groups$PFG)))
names(true_names)<- c("PFG")
true_names$K_n<- 1:16

Species_names_groups_num<- merge(Species_names_groups,true_names, by="PFG" )


# K=16 functional groups
#############################################################################Fitting the model 
data<- AB_norm_PA

set.seed(123)
smp_size <- floor(0.90 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]
#save(train, file = "sample_data_train_12.Rds")
#save(test, file = "sample_data_test_12.Rds")

#train<- load_object( "sample_data_train4.Rds")
#test<- load_object( "sample_data_test4.Rds")


it<-10000
burn<-5000
holdout<- sample(seq_len(nrow(train)), size = 100)


###################PCA
xdata_init<- data[,2:6]
#
##Normalization
#scaled.data <- as.data.frame(scale(xdata_init))
scaled.data<- xdata_init

###PCA 
xnew<- prcomp(scaled.data)
xdata_pca<- as.data.frame(xnew$x[,1:2])
eigs <- xnew$sdev^2
eigs[1] / sum(eigs) +eigs[2] / sum(eigs)


#### Matrix for Models Selection 



##########GJAM standart model
y<- train[,7:130]
#xdata<- xdata_pca[train_ind,]
xdata<- train[,2:6]

#xdata<- scaled.data
#xdata<- xdata_pca
#save(xdata, file = "sample_data_train_pca.Rds")

#formula <- as.formula( ~ PC1 +  PC2 + I(PC1^2)+  I(PC2^2))
#temp*deficit + I(temp^2) + I(deficit^2) 
formula <- as.formula( ~ bio_12_0+ I(bio_12_0^2) + slope + I(slope^2))
Ydata  <- gjamTrimY(y,10)$y             # at least 10 plots - re-group rare species
S<- ncol(Ydata)
rl <- list(r =5, N = S)
ml   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl,PREDICTX = F) #change ml
fit<-gjam(formula, xdata = xdata, ydata = Ydata, modelList = ml)


###Check the rank of X, should be number of columns
x <- model.matrix(formula, xdata)
qr(x)$rank

#save(fit,file="models_Bagues_data_OSS/fit_4.Rda")
#save(fit,file="models_Bagues_data_OSS/fit_12.Rda")
#no Holdout
#fit<- load_object("fit_4.Rda")
#fit$modelList$formula

#it2<- load_object("fit2_4.Rda")
#fit2$modelList$formula
#fit3<- load_object("fit3_4.Rda")
#fit3$modelList$ng
#fit4<- load_object("fit4_4.Rda")
#fit4$modelList$burnin
#fit<- load_object("models_Bagues_data_OSS/fit_12.Rda")
####### Out of sample prediction
Ykeep<- as.vector(colnames(Ydata))
y_test<- test[,c(Ykeep[1:(ncol(Ydata)-1)])]
#xdata_test<- test[,2:6]
xdata_test<- test[,2:6]
#save(xdata_test, file = "sample_data_test_pca.Rds")


new <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
p1  <- gjamPredict(output = fit, newdata = new)


AUC_GJAM<-vector()
for(i in 1:ncol(y_test)){ 
  label<- is.na(y_test[,i])
  predict<- p1$sdList$yMu[!label,i]
  test_value<- y_test[!label,i]
  if(sum(test_value)>0){
    AUC_GJAM<-c(AUC_GJAM,auc(roc(predict,factor(test_value))))
  }
}

mean(AUC_GJAM)



Tjur_GJAM<-vector()

for(k in 1:ncol(y_test)){
  label<- is.na(y_test[,k])
  predict<- p1$sdList$yMu[!label,k]
  test_value<- y_test[!label,k]
  indx <- test_value==1
  Tjur_GJAM <- c(Tjur_GJAM,(mean(predict[indx]) - mean(predict[!indx])))
}


mean(na.omit(Tjur_GJAM))



############################################
####Check the trace for number of groups. From the previous analysis we know that the number of functional groups is 16
### TOPLOT
trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:it
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 16,color = "red")
p

####The Dirichlet process prior with conjugation############################################
####The Dirichlet process prior : multinomial  + prior on alpha
##Computing the hyper-parameters for K=16
K=16
S=ncol(Ydata)
func<-function(x) {sum(x/(x+(1:S)-1))-16}
alpha.DP<-.bisec(func,0.01,100)
shape=((alpha.DP)^2)/20
rate=alpha.DP/20


rl2  <- list(r = 5, N = S,rate=rate,shape=shape,V=1) #here to modify N
ml2   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl2,PREDICTX = F) #change ml

fit2<-.gjam_2(formula, xdata = xdata, ydata = Ydata, modelList = ml2)
#save(fit2,file="models_Bagues_data_OSS/fit2_12.Rda")
#fit2<- load_object("models_Bagues_data_OSS/fit2_.Rda")


new <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
p2  <- gjamPredict(output = fit2, newdata = new)

AUC_GJAM2<-vector()
for(i in 1:ncol(y_test)){ 
  label<- is.na(y_test[,i])
  predict<- p2$sdList$yMu[!label,i]
  test_value<- y_test[!label,i]
  if(sum(test_value)>0){
    AUC_GJAM2<-c(AUC_GJAM2,auc(roc(predict,factor(test_value))))
  }
}

mean(AUC_GJAM2)

Tjur_GJAM2<-vector()

for(k in 1:ncol(y_test)){
  label<- is.na(y_test[,k])
  predict<- p2$sdList$yMu[!label,k]
  test_value<- y_test[!label,k]
  indx <- test_value==1
  Tjur_GJAM2 <- c(Tjur_GJAM2,(mean(predict[indx]) - mean(predict[!indx])))
}

mean(na.omit(Tjur_GJAM2))

### TO PLOT

trace<-apply(fit2$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:it
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
ml3   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl3,PREDICTX = F)


fit3 <- .gjam_3(formula,xdata,Ydata,ml3)



#save(fit3,file="models_Bagues_data_OSS/fit3_12.Rda")
#save(fit3,file="models_Bagues_data_OSS/fit3_2.Rda")

#fit3<- load_object("models_Bagues_data_OSS/fit3_12.Rda")

new <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
p3  <- gjamPredict(output = fit3, newdata = new)

AUC_PY1<-vector()
for(i in 1:ncol(y_test)){ 
  label<- is.na(y_test[,i])
  predict<- p3$sdList$yMu[!label,i]
  test_value<- y_test[!label,i]
  if(sum(test_value)>0){
    AUC_PY1<-c(AUC_PY1,auc(roc(predict,factor(test_value))))
  }
}

mean(AUC_PY1)



Tjur_PY1<-vector()

for(k in 1:ncol(y_test)){
  label<- is.na(y_test[,k])
  predict<- p3$sdList$yMu[!label,k]
  test_value<- y_test[!label,k]
  indx <- test_value==1
  Tjur_PY1 <- c(Tjur_PY1,(mean(predict[indx]) - mean(predict[!indx])))
}

mean(na.omit(Tjur_PY1))
# 
# k_mat_200<- as.data.frame(K_dmat)
# k_mat_200$n<- 200
# names(k_mat_200)<- c("alpha","sigma","PY","PY_A1","NG_A2","N")
# 
# 
# write.csv(k_mat_200, file = "K_200.csv")
# formattable(k_mat_200)
# 
############################################
####Check the trace for number of groups. From the previous analysis we know that the number of functional groups is 16
### TO PLOT
trace<-apply(fit3$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:it
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 16,color = "red")
p

####PY 4
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
ml4   <- list(ng = it, burnin = burn, typeNames = 'PA', reductList = rl4,PREDICTX = F)

fit4<-.gjam_4(formula, xdata = xdata, ydata = Ydata, modelList = ml4)
#fit4<- load_object("models_Bagues_data_OSS/fit4.Rda")


#save(fit4,file="models_Bagues_data_OSS/fit4_12.Rda")


new <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
p4  <- gjamPredict(output = fit4, newdata = new)

AUC_PY2<-vector()
for(i in 1:ncol(y_test)){ 
  label<- is.na(y_test[,i])
  predict<- p4$sdList$yMu[!label,i]
  test_value<- y_test[!label,i]
  if(sum(test_value)>0){
    AUC_PY2<-c(AUC_PY2,auc(roc(predict,factor(test_value))))
  }
}

mean(AUC_PY2)
######### TJUR
Tjur_PY2<-vector()

for(k in 1:ncol(y_test)){
  label<- is.na(y_test[,k])
  predict<- p1$sdList$yMu[!label,k]
  test_value<- y_test[!label,k]
  indx <- test_value==1
  Tjur_PY2 <- c(Tjur_PY2,(mean(predict[indx]) - mean(predict[!indx])))
}
mean(na.omit(Tjur_PY2))
############################################
####Check the trace for number of groups. From the previous analysis we know that the number of functional groups is 16
### TO PLOT

trace<-apply(fit4$chains$kgibbs,1,function(x) length(unique(x)))
df<-as.data.frame(trace)
df$iter<-1:it
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace, x=iter)) + geom_point() + 
  labs(title=paste0("Trace plot for the number of groups"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept = 16,color = "red")
p

########################################## Model comparison

#fit$fit$rmspeAll  #2.257202
#fit1$fit$rmspeAll #2.205223
#fit2$fit$rmspeAll #2.209276
#fit3$fit$rmspeAll #2.203517
#fit4$fit$rmspeAll #2.092136

fit$fit$DIC   #2510192
#fit1$fit$DIC  #2496028
fit2$fit$DIC  #2505273
fit3$fit$DIC  #2496243
fit4$fit$DIC  #2460125

AUC_data<- matrix(NA, nrow =length(AUC_GJAM), ncol =4)
AUC_data[,1]<- AUC_GJAM
AUC_data[,2]<- AUC_GJAM2
AUC_data[,3]<- AUC_PY1
AUC_data[,4]<- AUC_PY2
AUC_data_df<- as.data.frame(AUC_data)
names(AUC_data_df)<- c("GJAM","GJAM2","PY1","PY2")
AUC_fin<- melt(AUC_data_df)
### TO PLOT

p2<-ggplot(data=AUC_fin)+geom_boxplot(aes(y=as.numeric(value),x=as.factor(variable),fill=as.factor(variable)))+
  scale_y_continuous(name="AUC")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM2","PY1","PY2"))+theme_bw() + xlab("Models")  +theme_bw()
p2

AUC_fin_table<- as.data.frame(t(apply(AUC_data,2,mean)))
names(AUC_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
grid.newpage()
grid.table(AUC_fin_table)
grid.newpage()


Tjur_data<- matrix(NA, nrow =length(Tjur_GJAM), ncol =4)
Tjur_data[,1]<- Tjur_GJAM
Tjur_data[,2]<- Tjur_GJAM2
Tjur_data[,3]<- Tjur_PY1
Tjur_data[,4]<- Tjur_PY2
Tjur_data_df<- as.data.frame(Tjur_data)
names(Tjur_data_df)<- c("GJAM","GJAM2","PY1","PY2")
Tjur_fin<- melt(Tjur_data_df)


Tjur_fin_table<- as.data.frame(t(apply(na.omit(Tjur_data),2,mean)))
names(Tjur_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
grid.newpage()
grid.table(Tjur_fin_table)
grid.newpage()


p3<-ggplot(data=Tjur_fin)+geom_boxplot(aes(y=as.numeric(value),x=as.factor(variable),fill=as.factor(variable)))+
  scale_y_continuous(name="Tjur")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM2","PY1","PY2"))+ xlab("Models") +theme_bw()
p3

################################################################################ Partition

Species<- colnames(Ydata)[1: (ncol(Ydata)-1)]

True_clustering<- as.data.frame(Species)
names(True_clustering)<- c("CODE_CBNA")
True_clust<- merge(Species_names_groups_num,True_clustering, by="CODE_CBNA")
#True_clust<- True_clust[order(True_clust$CODE_CBNA),]


# 
# n <- ncol(Mat_new)
# M <- nrow(Mat_new)
# 
# if(sum(Mat_new %in% (1:n)) < n * M){
#   stop("All elements of cls must be integers in 1:nobs")
# }   


Mat<- fit4$chains$kgibbs[(burn+1):it,]
Mat_new<- matrix(NA, nrow=nrow(Mat),ncol=ncol(Mat))

for(i in 1:nrow(Mat) ){
  label<-unique(Mat[i,])
  change<- label[which(label>ncol(Mat))]
  all_n<- 1:ncol(Mat)
  admis<- all_n[which(!(all_n %in% label))]
  Mat_new[i,]<-Mat[i,]
  for(k in 1:length(change)){
    old_row<- Mat_new[i,]
    rep <-old_row==change[k]
    Mat_new[i,] <- replace(Mat_new[i,], rep, admis[k])
  }
}


tr_cl<-True_clust$K_n 
CM_DP1<- comp.psm(fit$chains$kgibbs[(burn+1):it,])
CM_DP2<- comp.psm(fit2$chains$kgibbs[(burn+1):it,])
CM_PY1<- comp.psm(fit3$chains$kgibbs[(burn+1):it,])
CM_PY2<- comp.psm(Mat_new)




mbind_DP1 <- minbinder(CM_DP1)
mbind_DP2 <- minbinder(CM_DP2)
mbind_PY1 <- minbinder(CM_PY1)
mbind_PY2 <- minbinder(CM_PY2)
ns<-ncol(fit$chains$kgibbs[(burn+1):it,])-1
# compare clusterings found by different methods with true grouping
Ar_D_DP1<- arandi(mbind_DP1$cl[1:ns], tr_cl)
Ar_D_DP2<-arandi(mbind_DP2$cl[1:ns], tr_cl)
Ar_D_PY1<-arandi(mbind_PY1$cl[1:ns], tr_cl)
Ar_D_PY2<-arandi(mbind_PY2$cl[1:ns], tr_cl)
Ar_D_fin_table<- as.data.frame(t(c(Ar_D_DP1,Ar_D_DP2,Ar_D_PY1,Ar_D_PY2)))
names(Ar_D_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
grid.newpage()
grid.table(Ar_D_fin_table)
grid.newpage()


vi.dist_DP1<- vi.dist(mbind_DP1$cl[1:ns], tr_cl)
vi.dist_DP2<- vi.dist(mbind_DP2$cl[1:ns], tr_cl)
vi.dist_PY1<- vi.dist(mbind_PY1$cl[1:ns], tr_cl)
vi.dist_PY2<- vi.dist(mbind_PY2$cl[1:ns], tr_cl)
VI_D_fin_table<- as.data.frame(t(c(vi.dist_DP1,vi.dist_DP2,vi.dist_PY1,vi.dist_PY2)))
names(VI_D_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
grid.newpage()

grid.table(VI_D_fin_table)
grid.newpage()





#################################################################################Other models



#gjam2
alpha<-mcmc(fit2$chains$alpha.DP_g)
alpha_ <- fit2$chains$alpha.DP_g[seq(burn, it, by=1)]
plot(alpha_, main="alpha DP multinomial")



#pdf("Bauges_data_trace_alpha_DP_by1.pdf")

trace_a<-alpha_
df<-as.data.frame(trace_a)
df$iter<-seq(burn, it, by=1)
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace_a, x=iter)) + geom_line() + 
  labs(title=bquote("Trace plot for " ~ alpha ~" for DP"))+
  theme_bw() +xlab("Value")+ ylab("Iteration")+ theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
p
#dev.off()



acfplot(alpha)
cumuplot(alpha)




##gjam2
alpha<-mcmc(fit4$chains$alpha.PY_g)
alpha_PY <- fit4$chains$alpha.PY_g[seq(burn, 15000, by=5)]

plot(alpha,main="alpha PY")
acfplot(alpha)
cumuplot(alpha)


pdf("Bauges_data_trace_alpha_PY_by5.pdf")


trace_a<-alpha_PY
df<-as.data.frame(trace_a)
df$iter<-seq(burn, 15000, by=5)
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace_a, x=iter)) + geom_line() + 
  labs(title=bquote("Trace plot for " ~ alpha ~" for PY"))+
  theme_bw() +xlab("Value")+ ylab("Iteration")+ theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
p
dev.off()


discount<-mcmc(fit4$chains$discount.PY_g)
sigma_PY <- fit4$chains$discount.PY_g[seq(burn, 15000, by=5)]

pdf("Bauges_data_trace_sigma_PY_by5.pdf")


trace_s<-sigma_PY
df<-as.data.frame(trace_s)
df$iter<-seq(burn, 15000, by=5)
#plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
p<-ggplot(df, aes(y=trace_s, x=iter)) + geom_line() + 
  labs(title=bquote("Trace plot for " ~ sigma ~" for PY"))+
  theme_bw() +xlab("Value")+ ylab("Iteration")+ ylim(c(0,1))+theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
p
dev.off()




plot(discount,main="sigma PY")
acfplot(discount)
cumuplot(discount)

#check the convergence

#for sigma
gjam_mc<- mcmc(fit$chains$sgibbs)
#gjam_mc1<- mcmc(fit1$chains$sgibbs) 
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
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+  theme_bw()+theme(plot.title=element_text(hjust =0.5),legend.position="top")+geom_hline(yintercept = 16,color = "red")
pdf("Bauges_data_trace_K.pdf")
p
dev.off()

#check the last weight
#pk_chains0_last<- mcmc(fit$chains$pk_g[,ncol(fit$chains$pk_g)])
#plot(pk_chains1_last)
#pk_chains1_last<- mcmc(fit1$chains$pk_g[,ncol(fit1$chains$pk_g)])
#plot(pk_chains1_last)
pk_chains2_last<- mcmc(fit2$chains$pk_g[,ncol(fit2$chains$pk_g)])
plot(pk_chains2_last)
pk_chains3_last<- mcmc(fit3$chains$pk_g[,ncol(fit3$chains$pk_g)])
plot(pk_chains3_last)
pk_chains4_last<- mcmc(fit4$chains$pk_g[,ncol(fit4$chains$pk_g)])
plot(pk_chains4_last)

#Out of sample prediction
#TJUR coefficient



pk_chains2_last<- mcmc(fit2$chains$pk_g[burn:it,ncol(fit2$chains$pk_g)])
plot(pk_chains2_last)
GJAM2_pk_last<- mean(pk_chains2_last)

pk_chains3_last<- mcmc(fit3$chains$pk_g[burn:it,ncol(fit3$chains$pk_g)])
plot(pk_chains3_last)
PY1_pk_last<- mean(pk_chains3_last)

pk_chains4_last<- mcmc(fit4$chains$pk_g[burn:it,ncol(fit4$chains$pk_g)])
plot(pk_chains4_last)
PY2_pk_last<- mean(pk_chains4_last)





###Not applicable here. 

#sum((fit$prediction$ypredMu[1:100,]-Ydata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5987808
#sum((fit1$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5961853
#sum((fit2$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5979932
#sum((fit3$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5977189
#sum((fit4$prediction$ypredMu[1:100,]-treeYdata[1:100,])^2)/sum(treeYdata[1:100,]^2) #0.5983279


# Final matrix
form<-c(formula)
Fin_all<-as.data.frame(matrix(NA,nrow=20,ncol=9))
names(Fin_all)<- c("Parameter","GJAM","GJAM2","PY1","PY2","r", "iterations", "burn","formula")
Fin_all$iterations<- it
Fin_all$burn<- burn
Fin_all$r<-5
Fin_all$formula<-as.character(form)
Fin_all[1,1]<- "DIC"
Fin_all[1,2:5]<- c(fit$fit$DIC,fit2$fit$DIC,fit3$fit$DIC,fit4$fit$DIC)
Fin_all[2,1]<- "mean AUC"
Fin_all[2,2:5]<- AUC_fin_table
Fin_all[3,1]<- "mean Tjur"
Fin_all[3,2:5]<- Tjur_fin_table
Fin_all[4,1]<- "mean p_N"
Fin_all[4,2:5]<- c(0,GJAM2_pk_last,PY1_pk_last,PY2_pk_last)

Fin_all[5,1]<- "VI dist"
Fin_all[5,2:5]<- VI_D_fin_table
Fin_all[6,1]<- "AR dist"
Fin_all[6,2:5]<- Ar_D_fin_table
Fin_all[,2:5]<- round(Fin_all[,2:5], 3)
grid.newpage()

grid.table(Fin_all[1:6,1:8])
grid.newpage()
grid.table(Fin_all[1:6,9])
grid.newpage()

###Sensitivity table
grid.newpage()

grid.table(fit$parameters$sensTable)
grid.newpage()

grid.table(fit2$parameters$sensTable)
grid.newpage()

grid.table(fit3$parameters$sensTable)
grid.newpage()

grid.table(fit4$parameters$sensTable)
grid.newpage()

grid.table(fit$inputs$designTable)
grid.newpage()

grid.table(fit2$inputs$designTable)
grid.newpage()

grid.table(fit3$inputs$designTable)
grid.newpage()

grid.table(fit4$inputs$designTable)
grid.newpage()



dev.off()
