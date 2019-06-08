

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


library(ggplot2)
#setwd("/mnt/beegfs/mistis/dbystrov/Gjammod")
setwd("/Users/dariabystrova/Documents/GitHub/gjamed")
# LtGJ <- load_object("Sim_smallSK4_gjam.Rda")
# Ltgj0<- load_object( "Sim_smallSK4_gjam0.Rda")
# LtT1<-load_object("Sim_smallSK4_type1.Rda")
# LtT2<-load_object("Sim_smallSK4_type2.Rda")
# LtT3<-load_object("Sim_smallSK4_type3.Rda")
# LtT4<-load_object("Sim_smallSK4_type4.Rda")
#LtGJ <- load_object("sim_med/ODSim_smallS1000K10_gjam.Rda")
Ltgj0<- load_object( "sim_med/ODSim_smallS1000K10_gjam0.Rda")
LtT1<-load_object("sim_med/ODSim_smallS1000K10_type1.Rda")
LtT2<-load_object("sim_med/ODSim_smallS1000K10_type2.Rda")
LtT3<-load_object("sim_med/ODSim_smallS1000K10_type3.Rda")
LtT4<-load_object("sim_med/ODSim_smallS1000K10_type4.Rda")






S_vec<-c(1000)
r_vec<-c(5)
modn<-4
dsnum<-5


table_comp<-as.data.frame(matrix(NA, nrow=length(r_vec)*length(S_vec)*modn*dsnum, ncol=1))
table_comp$S<- rep(S_vec, each=modn*length(r_vec)*dsnum)
table_comp$mod<- rep(c("gjam1","gjam2","gjam3","gjam4"), length(r_vec)*length(S_vec), each=dsnum)
table_comp$num<- rep(1, each=dsnum)
table_comp$rv<- rep(c(5), each=dsnum, length(S_vec))
table_comp$avnum<- rep(1:5, 4)
burn<-500
it<- 1000
nsamples<-500
for(i in (1: nrow(table_comp))){
  j<-table_comp$num[i]
  jm<-table_comp$avnum[i]
  table_comp$Schar[i]<- paste0("S=",table_comp$S[i])
  # if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtGJ)[j]))>0)){ table_comp$res[i]<- mean(LtGJ[[j]]$trace[-c(1:burn)])}
  # if((table_comp$mod[i]=="gjam1")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT1)[j]))>0)){ table_comp$res[i]<-  mean(LtT1[[j]]$trace[-c(1:burn)])}
  # if((table_comp$mod[i]=="gjam2")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT2)[j]))>0)){ table_comp$res[i]<- mean(LtT2[[j]]$trace[-c(1:burn)])}
  # if((table_comp$mod[i]=="gjam3")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT3)[j]))>0)){ table_comp$res[i]<-  mean(LtT3[[j]]$trace[-c(1:burn)])}
  # if((table_comp$mod[i]=="gjam4")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT4)[j]))>0)){ table_comp$res[i]<- mean(LtT4[[j]]$trace[-c(1:burn)])}
  # 
  # if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(Ltgj0)[j]))>0)){
  #   table_comp$res[i]<- mean(LtGJ[[j]][[jm]]$trace)
  #   table_comp$lweight[i]<- mean(Ltgj0[[j]][[jm]]$pkN)
  #   table_comp$fit_err[i]<- Ltgj0[[j]][[jm]]$fit
  #   table_comp$err[i]<- Ltgj0[[j]][[jm]]$err
  # }
  if((table_comp$mod[i]=="gjam1")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT1)[j]))>0)){
    table_comp$res[i]<-  mean(LtT1[[j]][[jm]]$trace)
    table_comp$lweight[i]<- mean(LtT1[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT1[[j]][[jm]]$fit
    table_comp$err[i]<- LtT1[[j]][[jm]]$err
  }
  if((table_comp$mod[i]=="gjam2")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT2)[j]))>0)){
    table_comp$res[i]<- mean(LtT2[[j]][[jm]]$trace)
    table_comp$lweight[i]<- mean(LtT2[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT2[[j]][[jm]]$fit
    table_comp$err[i]<- LtT2[[j]][[jm]]$err
  }
  if((table_comp$mod[i]=="gjam3")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT3)[j]))>0)){
    table_comp$res[i]<-  mean(LtT3[[j]][[jm]]$trace)
    table_comp$lweight[i]<-  mean(LtT3[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT3[[j]][[jm]]$fit
    table_comp$err[i]<- LtT3[[j]][[jm]]$err
  }
  if((table_comp$mod[i]=="gjam4")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT4)[j]))>0)){
    table_comp$res[i]<- mean(LtT4[[j]][[jm]]$trace)
    table_comp$lweight[i]<-  mean(LtT4[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT4[[j]][[jm]]$fit
    table_comp$err[i]<- LtT4[[j]][[jm]]$err
  }
  
}



########combine tables
#S=1000
tab_1000<-load_object("tablecompS1000.rds")
load_object("tablecompS1000.rds")
#S=100
tab_100<-load_object("tablecompS100.rds")

#S=200
tab_200<-load_object("tablecompS200.rds")


#S=200
tab_500<-load_object("tablecompS500.rds")



table_all<-rbind(tab_1000,tab_100,tab_200,tab_500)
table_comp<-table_all
table_comp$Schar <- factor(table_comp$Schar, levels=c('S=100','S=200','S=500','S=1000'))
table_comp<-table_comp[,-1]

##########################################################################################
pdf("Clust_prop_gr4SmallS1000.pdf")
q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5"),
                   labels=c("r=5"),limits=c("5"))+
  scale_y_continuous(name="Number of clusters",limits=c(0,15),breaks=seq(2,15,by=2))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") +
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 10,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_clust

dev.off()
#######################################################################

pdf("Clust_prop_gr4SmallS100_1000.pdf")
q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=Schar,y= as.numeric(res), fill = as.factor(mod)))+
  scale_y_continuous(name="Number of clusters",limits=c(0,15),breaks=seq(2,15,by=2))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  ylab("Posterior mean") +
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 10,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                             strip.background = element_blank(),
                                                             strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_clust

dev.off()


######################################################################

####################
pdf("Weights_gr4SmallS1000.pdf")

q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_x_discrete(name="Parameters", breaks=c("5"),
                   labels=c("r=5"),limits=c("5"))+
  scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x") +
  labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_weight
dev.off()

##########################################################################################


pdf("Weights_gr4SmallS100_1000.pdf")

q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_weight
dev.off()


##########################################################################################
pdf("Fit_error_gr4SmallS100_1000.pdf")

q_fiterr<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(fit_err), fill = as.factor(mod))) +
  scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("r=5"),limits=c("5"))+
  scale_y_continuous(name="Error",limits=c(0,max(table_comp$fit_err)))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean")+
  labs(title="Fit error value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                              strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_fiterr

dev.off()

pdf("RMSE_gr4SmallS1000.pdf")


q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(err), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("5"),limits=c("5"))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(table_comp$err)))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_rmse_err

dev.off()



pdf("RMSE_gr4SmallS100_1000.pdf")


q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(err), fill = as.factor(mod)))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(table_comp$err)))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
   ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_rmse_err

dev.off()
########################################################################################################################





save(table_comp, file="tablecompS100_1000.Rds")














#####################Plot individual plots
list_mod<-c("GJAM","GJAM0","GJAM1","GJAM2","GJAM3","GJAM4")

for(k in (1:6)){
  pdf(paste0("Ind_",GJAM,"_average.pdf"))
  
}
pdf("Ind_GJAM_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(LtGJ[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
pdf("Ind_GJAM0_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(Ltgj0[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
pdf("Ind_GJAM1_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(LtT1[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
pdf("Ind_GJAM2_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(LtT2[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
pdf("Ind_GJAM3_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(LtT3[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
pdf("Ind_GJAM4_average.pdf")
for(i in 1:9){
  for(j in 1:10){
    lapply(LtT4[[i]][[j]]$pl_list, function(x) plot(x))
  }
}
dev.off()
