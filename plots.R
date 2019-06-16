

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
Ltgj0<- load_object( "ODSim_smallS1000K_10_gjam0.Rda")
LtT1<-load_object("ODSim_smallS1000K_10_type1.Rda")
LtT2<-load_object("ODSim_smallS1000K_10_type2.Rda")
LtT3<-load_object("OD2Sim_smallS1000K10_type3.Rda")
LtT4<-load_object("OD2Sim_smallS1000K10_type4.Rda")







S_vec<-c(1000)
r_vec<-c(5)
modn<-5
dsnum<-2
it<- 2000
burn<- 1000
nsamples<-500


table_comp<-as.data.frame(matrix(NA, nrow=length(r_vec)*length(S_vec)*modn*dsnum, ncol=1))
table_comp$S<- rep(S_vec, each=modn*length(r_vec)*dsnum)
table_comp$mod<- rep(c("gjam","gjam1","gjam2","gjam3","gjam4"), length(r_vec)*length(S_vec), each=dsnum)
table_comp$num<- rep(1, each=dsnum)
table_comp$rv<- rep(c(5), each=dsnum, length(S_vec))
table_comp$avnum<- rep(1:2, 5)
burn<-1000
it<- 2000
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
  if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(Ltgj0)[j]))>0)){
    table_comp$res[i]<- mean(Ltgj0[[j]][[jm]]$trace)
    table_comp$lweight[i]<- mean(Ltgj0[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- Ltgj0[[j]][[jm]]$fit
    table_comp$err[i]<- Ltgj0[[j]][[jm]]$err
  }
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

table_comp<-table_comp[,-1]

#save(table_comp, file = "tablecompS1000_all.rds")



########combine tables
#S=1000
tab_1000<-load_object("tablecompS1000_all.rds")

#S=100
tab_100<-load_object("tablecompS100.rds")

#S=200
tab_200<-load_object("tablecompS200.rds")


#S=500
tab_500<-load_object("tablecompS500.rds")





########combine tables
#S=1000
tab_1000<-load_object("tablecompS1000_cor.rds")

#S=100
tab_100<-load_object("tablecompS100_cor.rds")

#S=200
tab_300<-load_object("tablecompS300_cor.rds")


#S=500
tab_500<-load_object("tablecompS500_cor.rds")

tab_500<- tab_500[,-11]



nsamples=500
it<-2000
burn<-1000

table_all<-rbind(tab_100,tab_300,tab_500,tab_1000)
table_comp<-table_all
table_comp$Schar <- factor(table_comp$Schar, levels=c('S=100','S=300','S=500','S=1000'))
#table_comp<-table_comp[,-1]

##########################################################################################
# pdf("Clust_prop_gr4SmallS1000.pdf")
# q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
#   scale_x_discrete(name="Parameters", breaks=c("5"),
#                    labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name="Number of clusters",limits=c(0,15),breaks=seq(2,15,by=2))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") +
#   labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
#   geom_hline(yintercept = 10,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
#                                                             strip.background = element_blank(),
#                                                             strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_clust
# 
# dev.off()
#######################################################################

pdf("Clust_prop_gr4SmallS100_1000K10.pdf")
q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=Schar,y= as.numeric(res), fill = as.factor(mod)))+
  scale_y_continuous(name="Number of clusters",limits=c(0,15),breaks=seq(2,15,by=2))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  ylab("Posterior mean") +
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 4,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                             strip.background = element_blank(),
                                                             strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_clust

dev.off()




Weights_to_table10<- table_all[,c("Schar","mod", "res","lweight","err")]
Weights_to_table10$mod<- as.factor(Weights_to_table10$mod)

Weights_to_table_sum10<-aggregate(Weights_to_table10, by=list(Weights_to_table10$mod,Weights_to_table10$Schar),  FUN = "mean")

Weights_to_table_sum10$lweight<- round(Weights_to_table_sum10$lweight, 3)

write.csv(Weights_to_table_sum10, file = "WeightsK10.csv")

Wegihts_only<- Weights_to_table_sum10[,c(1,2,6)]
Wegihts_out<- dcast(data = WO,formula = Group.2~Group.1,value.var = "lweight")

write.csv(Wegihts_out, file = "WeightsK10.csv")

######################################################################

####################
# pdf("Weights_gr4SmallS1000.pdf")
# 
# q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(lweight), fill = as.factor(mod))) +
#   scale_x_discrete(name="Parameters", breaks=c("5"),
#                    labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x") +
#   labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
#   theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_weight
# dev.off()

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
# pdf("Fit_error_gr4SmallS100_1000.pdf")
# 
# q_fiterr<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(fit_err), fill = as.factor(mod))) +
#   scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name="Error",limits=c(0,max(table_comp$fit_err)))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean")+
#   labs(title="Fit error value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
#                                                                                                                                                               strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_fiterr
# 
# dev.off()
# 
# pdf("RMSE_gr4SmallS1000.pdf")
# 
# 
q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(sqrt(err)), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("5"),limits=c("5"))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(sqrt(table_comp$err))))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
 q_rmse_err
# 
# dev.off()



pdf("RMSE_gr4SmallS100_1000K10.pdf")


q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(sqrt(err)), fill = as.factor(mod)))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(sqrt(table_comp$err))))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_rmse_err

dev.off()









########combine tables Corrected
#S=1000
#tab_1000<-load_object("tablecompS300_cor.Rds")

#S=100
tab_100<-load_object("tablecompS100_cor.Rds")

#S=300
tab_300<-load_object("tablecompS300_cor.Rds")


#S=500
tab_500<-load_object("tablecompS500_cor.Rds")

tab_500<- tab_500[,c(1:10)]

#S=1000
tab_1000<-load_object("tablecompS1000_cor.Rds")


it<-2000
burn<-1000

table_all<-rbind(tab_100,tab_200,tab_500,tab_1000)
table_comp<-table_all
table_comp$Schar <- factor(table_comp$Schar, levels=c('S=100','S=200','S=500','S=1000'))
#table_comp<-table_comp[,-1]

##########################################################################################
# pdf("Clust_prop_gr4SmallS1000.pdf")
# q_clust<-  ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
#   scale_x_discrete(name="Parameters", breaks=c("5"),
#                    labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name="Number of clusters",limits=c(0,15),breaks=seq(2,15,by=2))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") +
#   labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
#   geom_hline(yintercept = 10,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
#                                                             strip.background = element_blank(),
#                                                             strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_clust
# 
# dev.off()
#######################################################################

pdf("Clust_prop_gr10SmallS100_1000.pdf")
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
# pdf("Weights_gr4SmallS1000.pdf")
# 
# q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(lweight), fill = as.factor(mod))) +
#   scale_x_discrete(name="Parameters", breaks=c("5"),
#                    labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x") +
#   labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
#   theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_weight
# dev.off()

##########################################################################################


pdf("Weights_gr10SmallS100_1000.pdf")

q_weight<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_y_continuous(name=expression(p[N]),limits=c(0,0.3))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  labs(title=expression(paste("Last ",p[N]," weight value for all models")), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  theme_bw()+theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_weight
dev.off()


##########################################################################################
# pdf("Fit_error_gr4SmallS100_1000.pdf")
# 
# q_fiterr<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(fit_err), fill = as.factor(mod))) +
#   scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("r=5"),limits=c("5"))+
#   scale_y_continuous(name="Error",limits=c(0,max(table_comp$fit_err)))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean")+
#   labs(title="Fit error value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
#                                                                                                                                                               strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_fiterr
# 
# dev.off()
# 
# pdf("RMSE_gr4SmallS1000.pdf")
# 
# 
# q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(sqrt(err)), fill = as.factor(mod)))+
#   scale_x_discrete(name="Parameters", breaks=c("5"),labels=c("5"),limits=c("5"))+
#   scale_y_continuous(name="RMSE error",limits=c(0,max(sqrt(table_comp$err))))+
#   scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
#   facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+
#   labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
#                                                                                                                                                                                                       strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
# q_rmse_err
# 
# dev.off()



pdf("RMSE_gr10SmallS100_1000.pdf")


q_rmse_err<- ggplot(data= table_comp) +geom_boxplot(aes(x=as.factor(Schar),y= as.numeric(sqrt(err)), fill = as.factor(mod)))+
  scale_y_continuous(name="RMSE error",limits=c(0,max(sqrt(table_comp$err))))+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
   ylab("Posterior mean") + xlab("S and r values")+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                                                      strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q_rmse_err

dev.off()
##############

df_alpha <- data.frame(matrix(NA, nrow =100, ncol =1))
df_alpha$alpha<-LtT2$S_1000_r_5_N_150_n_500_K4$S_1000_q_20n_500_K_4_l2$alpha.chains
df_alpha$type<- "posterior"
#df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
#df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
#alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
#df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)

#df_alpha_prior$type<- "prior"
#df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
###Compute mean
mu <- ddply(df_alpha, "type", summarise, grp.mean=mean(alpha))
mu1<- as.data.frame(LtT2$S_1000_r_5_N_150_n_500_K4$S_1000_q_20n_500_K_4_l2$alpha)
colnames(mu1)<- c("grp.mean")
mu1$type<- "prior"
mu<- rbind(mu, mu1)



pdf("Posterior_density_alphaT2.pdf")
p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density(color="red",adjust = 2)+labs(title=paste0("Posterior distribution for alpha")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
p_alpha_2
dev.off()





df_alpha <- data.frame(matrix(NA, nrow =100, ncol =1))
df_alpha$alpha<-LtT1$S_1000_r_5_N_150_n_500_K4$S_1000_q_20n_500_K_4_l1$alpha.chains
df_alpha$type<- "posterior"
#df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
#df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
#alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
#df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)

#df_alpha_prior$type<- "prior"
#df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
###Compute mean
mu <- ddply(df_alpha, "type", summarise, grp.mean=mean(alpha))
mu1<- as.data.frame(LtT1$S_1000_r_5_N_150_n_500_K4$S_1000_q_20n_500_K_4_l1$alpha)
colnames(mu1)<- c("grp.mean")
mu1$type<- "prior"
mu<- rbind(mu, mu1)



pdf("Posterior_density_alphaT1.pdf")
p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
  geom_density(color="red",adjust = 1.2)+labs(title=paste0("Posterior distribution for alpha")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
p_alpha_2
dev.off()



