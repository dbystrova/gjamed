rm(list=ls())
setwd("/Users/dariabystrova/Documents/GitHub/gjamed")
########## Simulating the data
library(MASS)
library(repmis)
library(gjam)
library(rlist)
library(truncnorm)
#library(coda)
library(RcppArmadillo)
library(arm)
library(NLRoot)
library(Rcpp)
library(plyr)
library(ggplot2)
#library(ggsn)
library(parallel)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions_mod.R")
source("R/simple_gjam_0.R")
source("R/simple_gjam_1.R")
source("R/simple_gjam_2.R")
source("R/simple_gjam_3.R")
source("R/simple_gjam_4.R")


generate_data<-function(Sp=50,nsamples=500,qval=20,Ktrue=4){
  S<-Sp
  n<- nsamples
  q<- qval
  env<-runif(-50,50,n=n)
  X<-cbind(1,poly(env,2)) #nxK
  idx<-sample(S)
  B_0<-seq(0,100,length.out=S)[idx]
  B_1<-seq(0,100,length.out=S)[idx]
  B_2<-seq(0,100,length.out=S)[idx]
  B<-cbind(B_0,B_1,B_2) #SxK
  L<-X%*%t(B) #nxS
  K_t<- Ktrue
  cat("True number of clusters : ",K_t,"\n")
  A<-matrix(NA,nrow=ceiling(K_t),ncol=q)
  for(i in 1:ceiling(K_t)){
    A[i,]<-mvrnorm(n = 1,rep(0,q), Sigma=3*diag(q)) #Nxq short and skinny
  }
  idx<-sample((1:ceiling(K_t)),S,replace=T) 
  Lambda<-A[idx,] #Sxr tall and skinny
  Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  Sigma_true<-Sigma
  Y<-mvrnorm(n = n, mu=rep(0,S), Sigma=Sigma)
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env1","env2")
  
  return(list(xdata=xdata, Y=Y,idx=idx,S_true=Sigma_true))
}




simulation_fun_gjam4<-function(data_set,Sp, Ntr, rval,nsamples=500, Ktrue,q=20, it=1000, burn=500){
  S<-Sp
  n<- nsamples
  r <- rval
  iterations<-it

  K=sum(S/(S+(1:S)-1)) #104, his prior number of clusters when alpha=S
  cat("Prior expected number of clusters : ",Ktrue,"\n")
  K_t= Ktrue
  xdata<-data_set$xdata
  Y<-data_set$Y
  idx<- data_set$idx
  Sigma_true<- data_set$S_true
  formula<-as.formula(~env1+env2)
   eps=0.1
  alp_sig<-as.data.frame(matrix(NA,nrow=20,ncol=3))
  colnames(alp_sig)<-c("alpha","sigma","is_less_150")
  alp_sig$sigma=seq(0.05,0.4,length.out = 20)
  #loop to run bisecetion on a grid for sigma
  for(i in 1:20){
    ####corrected added  -1
    func<-function(x) {(x/alp_sig[i,"sigma"])*(prod((x+alp_sig[i,"sigma"]+c(1:S)-1)/(x+c(1:S) -1))-1) - K_t}
    alp_sig[i,"alpha"]<-.bisec(func,0.0001,100)
    N_eps<-floor(.compute_tau_mean(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps) + 2*.compute_tau_var(alp_sig[i,"sigma"], alp_sig[i,"alpha"],eps))
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
  N_eps<-floor(.compute_tau_mean(sigma_py_max,alpha.max_val,eps) + 2*.compute_tau_var(sigma_py_max,alpha.max_val,eps))

  rl   <- list(r = r, N = N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
  ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
  fit<-.gjam_4(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
  alpha.chains<-fit$chains$alpha.PY_g
  sigma.chains<-fit$chains$discount.PY_g
  pk_chains<- fit$chains$pk_g
  Ntr<-N_eps+1
  alpha.DP<-alpha.PY
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  ind_trace<- seq(1,it,by=1)
  trace_short<- trace[ind_trace]
  df<-as.data.frame(trace)
  df$iter<-1:it
  #####Weights plot
     pk<- apply(fit$chains$pk_g[-c(1:burn),],2,mean)
    last_pk<- round(pk[Ntr-1],3)
    df_weights <- data.frame(matrix(NA, nrow = Ntr-1, ncol =1))
    df_weights$pw<-pk 
    df_weights$tr<-1:(Ntr-1)
    pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
      geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
      geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Weights for the case: S=",S," ,r=",r," true gr K=",K_t,"N=",Ntr, " pN=",last_pk), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
    # pl_weigths

  #####Alpha plot
     df_alpha <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
    df_alpha$alpha<- alpha.chains[-c(1:burn)]
    df_alpha$type<- "posterior"
    df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
    #df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
    alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
    df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)
    
    df_alpha_prior$type<- "prior"
    df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
    ###Compute mean
    mu <- ddply(df_alpha_all, "type", summarise, grp.mean=mean(alpha))
    mu$grp.mean[which(mu$type=="prior")]=alpha.DP
    p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
      geom_density(color="red")+labs(title=paste0("Posterior distribution for alpha"), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples," S=",S," ,r=",r," true gr K=",K_t, " ,N=",Ntr))+
      theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
      scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
    plot(p_alpha_2)
 
  #######Sigma plot

  df_sigma <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
  df_sigma$sigma<- sigma.chains[-c(1:burn)]
  df_sigma$type<- "posterior"
  ###Compute mean
  mu <- ddply(df_sigma, "type", summarise, grp.mean=mean(sigma))
  p_sigma<- ggplot(df_sigma, aes(x=sigma)) + geom_vline(data=mu, aes(xintercept=grp.mean),linetype="dashed")+
    geom_density()+labs(title=paste0("Distribution sigma: S=",S," ,r=",r," true gr K=",K_t,",N=",Ntr), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
  plot(p_sigma)
  
  dfs<-as.data.frame(sigma.chains)
  dfs$iter<-1:it
  #plot(apply(fit$chains$kgibbs,1,function(x) length(unique(x))))
  p_trace_sigma<-ggplot(dfs, aes(y=sigma.chains, x=iter)) + geom_point() + 
    labs(title=paste0("Trace plot for the sigma for S=",S," r=",r," true K=",K_t), caption=paste0("Number of iterations: ",it," burnin: ",burn,"number of samples: ",nsamples))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
  plot(p_trace_sigma)

  
  N_dim<-(it-burn)
  sigma<-array(dim=c(Sp,Sp,N_dim))
  for(j in 1:N_dim){
    K<-fit$chains$kgibbs[j,]
    Z  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)
    sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
  }
  sigma_mean<-apply(sigma,c(1,2),mean)
  err<-sum((sigma_mean-Sigma_true)^2)/(Sp*Sp)
  rmspe<-fit$fit$rmspeAll
  pN_chain<- pk_chains[-c(1:burn),(Ntr-1)]
  return(list(trace=trace_short,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,pk_val=pk, pkN=pN_chain, 
              coeff_t=Sigma_true,coeff_f=sigma_mean,
              err=err,fit=rmspe, sigma_chain=sigma.chains, sigma.PY=sigma_py))

}
####### Just one possible test case
#sim<-simulation_fun(Sp=50, Ntr=150, rval=3,nsamples=500, Ktrue=4,it=1000,burn=200)
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
# - K_t, the true number of clusters

###########Simulation for Continous data case :small S K=4###################################################



#####################################Simulation 2 K=10#######################################

####Small S, N==S, n=500


list5=list()
list0=list()
data_list=list()
lk<-list()
S_vec<-c(100,200)
r_vec<-5
#n_vec<-c(10)
n_samples<- 500
k<-1
it<-5000
burn<-2500
Ktr<-10
q<-20

path<- "~/Documents/GitHub/gjamed/sigma_post/"
for(i in 1:length(S_vec)){
  data_list=list()
  k=1
  for(l in (1:1)){
      data_list<- list.append(data_list,generate_data(Sp=S_vec[i],nsamples=n_samples,qval=q,Ktrue=Ktr))
      names(data_list)[[l]]<-paste0("S_",S_vec[i],"_q_",q,"n_",n_samples,"_K_",Ktr,"_l",l)
    }
    ########gjam 4  model list########################    
    l5<-list()
    l5<- lapply(data_list,simulation_fun_gjam4,Sp=S_vec[i], Ntr=150, q=q,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn)
    list5<-list.append(list5,assign(paste0("S_",S_vec[i],"_r_5_N_n_500_K",Ktr),l5))
    names(list5)[[k]]<-paste0("S_",S_vec[i],"_r_5_N_150_n_",n_samples,"_K",Ktr)
    save(list5, file = paste0(path,"Sigma_mod",S_vec[i],"K",Ktr,"_type4.Rda"))
    k=k+1
}




load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}




plot(list5$S_200_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$sigma_chain)


plot(density(list5$S_200_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$sigma_chain))



post1<- load_object(paste0(path,"Sigma_mod100K10_type4.Rda"))



df_sigma <- data.frame(matrix(NA, nrow =5000, ncol =1))
df_sigma$sigma<-post1$S_100_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$sigma_chain
df_sigma$type<- "posterior"
#df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
#df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
#alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
#df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)

#df_alpha_prior$type<- "prior"
#df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
###Compute mean
mu <- ddply(df_sigma, "type", summarise, grp.mean=mean(sigma))
mu1<- as.data.frame(post1$S_100_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$sigma.PY)
colnames(mu1)<- c("grp.mean")
mu1$type<- "prior"
mu<- rbind(mu, mu1)



pdf(paste0(path,"Posterior_density_sigmaT4.pdf"))
p_sigma<- ggplot(df_sigma, aes(x=sigma)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density(color="red",adjust = 2)+labs(title=paste0("Posterior distribution for sigma")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
p_sigma


it<-5000
burn<-2500
df_alpha <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
df_alpha$alpha<- post1$S_100_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$alpha.chains[-c(1:burn)]
df_alpha$type<- "posterior"
#df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
###Compute mean
mu <- ddply(df_alpha, "type", summarise, grp.mean=mean(alpha))
mu1<- as.data.frame(post1$S_100_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$alpha)
colnames(mu1)<- c("grp.mean")
mu1$type<- "prior"
mu<- rbind(mu, mu1)

p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density(color="red")+labs(title=paste0("Posterior distribution for alpha"))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
p_alpha_2

dev.off()


post2<- load_object(paste0(path,"Sigma_mod200K10_type4.Rda"))

df_sigma <- data.frame(matrix(NA, nrow =1000, ncol =1))
df_sigma$sigma<-post2$S_200_r_5_N_150_n_500_K10$S_100_q_20n_500_K_10_l1$sigma_chain
df_sigma$type<- "posterior"
#df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
#df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
#alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
#df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)

#df_alpha_prior$type<- "prior"
#df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
###Compute mean
mu <- ddply(df_sigma, "type", summarise, grp.mean=mean(sigma))
#mu1<- as.data.frame(LtT2$S_1000_r_5_N_150_n_500_K4$S_1000_q_20n_500_K_4_l2$alpha)
#colnames(mu1)<- c("grp.mean")
#mu1$type<- "prior"
#mu<- rbind(mu, mu1)



pdf(paste0(path,"Posterior_density_sigmaS200T4.pdf"))
p_sigma<- ggplot(df_sigma, aes(x=sigma)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density(color="red",adjust = 2)+labs(title=paste0("Posterior distribution for sigma")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
p_sigma
dev.off()



# 
# for(l in (1:2)){
#   data_list<- list.append(data_list,generate_data(Sp=S,nsamples=n_samples,qval=q,Ktrue=Ktr))
#   names(data_list)[[l]]<-paste0("S_",S,"_q_",q,"n_",n_samples,"_K_",Ktr,"_l",l)
# }
# ########gjam 4  model list########################    
# l5<-list()
# l5<- lapply(data_list,simulation_fun_oneDS,Sp=S_vec[i], Ntr=150, q=q,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn,type="4")
# list5<-list.append(list5,assign(paste0("S_",S_vec[i],"_r_5_N_150_n_500_K",Ktr),l5))
# names(list5)[[k]]<-paste0("S_",S_vec[i],"_r_5_N_150_n_",n_samples,"_K",Ktr)
# save(list5, file = paste0(path,"ODSim_smallS",S_vec[i],"K",Ktr,"_type4.Rda"))
# k=k+1
# 
# data_set<- data_list$S_10_q_20n_500_K_10_l1
# 
# f4<- simulation_fun_gjam4(data_set,Sp=10, Ntr=150, rval=5,nsamples=500, Ktrue=10,q=20, it=1000, burn=500)
# 
#########################################################################################



pdf(paste0(path,"Posterior_weights_all.pdf"))


path2<- "~/Documents/GitHub/gjamed/"


WtabS100K10<- load_object(paste0(path2,"weights_tablecompS100K10.Rds"))
WtabS100K10$Kchar<- "K=10"
WtabS300K20<- load_object(paste0(path2,"weights_tablecompS300K20.Rds"))
WtabS300K20$Schar<- "S=300"
WtabS300K20$Kchar<-"K=20"

WtabS500K50<- load_object(paste0(path2,"weights_tablecompS500K50.Rds"))
WtabS500K50$Kchar<-"K=50"
Weights_all<- rbind(WtabS100K10,WtabS300K20,WtabS500K50)

table_comp<-Weights_all
table_comp$Schar <- factor(table_comp$Schar, levels=c('S=100','S=300','S=500'))
#table_comp<-table_comp[,-1]
table_comp$Kchar <- factor(table_comp$K, levels=c('K=10','K=20','K=50'))


q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=Schar,y= as.numeric(res), fill = as.factor(mod)))+
  scale_y_continuous(name="Number of clusters",limits=c(0,60),breaks=seq(2,15,by=2))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  ylab("Posterior mean") +
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 10,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                             strip.background = element_blank(),
                                                             strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_clust







q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_weight



q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(sqrt(err)), fill = as.factor(mod)))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(sqrt(table_comp$err))))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_rmse_err



#dev.off()

Weights_to_table<- Weights_all[,c("Schar","mod", "res","lweight","err","Neps")]
Weights_to_table$mod<- as.factor(Weights_to_table$mod)

Weights_to_table_sum<-aggregate(Weights_to_table, by=list(Weights_to_table$mod,Weights_to_table$Schar),  FUN = "mean")

Weights_to_table_sum$lweight<- round(Weights_to_table_sum$lweight, 3)

Weights_to_table_sum_to_csv<- Weights_to_table_sum[which(Weights_to_table_sum$Group.1%in%c("gjam3","gjam4")),]




write.csv(Weights_to_table_sum_to_csv, file = "Weights.csv")



# caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples)


pdf(paste0(path2,"Posterior_weights_all_Kdiff.pdf"))



q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_x_discrete(name="Parameters", breaks=c("S=100", "S=300","S=500"),
                   labels=c("S=100, K=10", "S=300, K=20","S=500, K=50"))+
  scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  labs(title=expression(paste("Last ",p[N]," weight value for all models and different number of groups")))+
  theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_weight


dev.off()



 q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Kchar),y= as.numeric(lweight), fill = as.factor(mod))) +
   scale_x_discrete(name="Parameters", breaks=c("K=10", "K=20","K=50"),
                    labels=c("K=10", "K=20","K=50"),limits=c("K=10", "K=20","K=50"))+
   scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x") +
   labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
   theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
 q_weight

 
 dev.off()


