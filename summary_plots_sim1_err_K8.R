


load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


setwd("~/Documents/GitHub/gjamed")
library(ggplot2)
# LtGJ <- load_object("Sim_smallSK4_gjam.Rda")
# Ltgj0<- load_object( "Sim_smallSK4_gjam0.Rda")
# LtT1<-load_object("Sim_smallSK4_type1.Rda")
# LtT2<-load_object("Sim_smallSK4_type2.Rda")
# LtT3<-load_object("Sim_smallSK4_type3.Rda")
# LtT4<-load_object("Sim_smallSK4_type4.Rda")
LtGJ <- load_object("ODSim_smallSK4_gjam.Rda")
Ltgj0<- load_object( "ODSim_smallSK4_gjam0.Rda")
LtT1<-load_object("ODSim_smallSK4_type1.Rda")
LtT2<-load_object("ODSim_smallSK4_type2.Rda")
LtT3<-load_object("ODSim_smallSK4_type3.Rda")
LtT4<-load_object("ODSim_smallSK4_type4.Rda")
 



S_vec<-c(20,50,80)
r_vec<-c(5,10,20)



table_comp<-as.data.frame(matrix(NA, nrow=length(r_vec)*length(S_vec)*5*10, ncol=1))
table_comp$S<- rep(S_vec, each=3*50)
table_comp$mod<- rep(c("gjam","gjam1","gjam2","gjam3","gjam4"), 3*3, each=10)
table_comp$num<- rep(1:9, each=50)
table_comp$rv<- rep(c(5,10,20), each=50, 3)
table_comp$avnum<- rep(1:10, 5*9)
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

  if((table_comp$mod[i]=="gjam")&(length(grep(paste0("r_",table_comp$rv[i]),names(Ltgj0)[j]))>0)){
    table_comp$res[i]<- mean(LtGJ[[j]][[jm]]$trace[-c(1:burn)])
    table_comp$lweight[i]<- mean(Ltgj0[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- Ltgj0[[j]][[jm]]$fit
    table_comp$err[i]<- Ltgj0[[j]][[jm]]$err
    }
  if((table_comp$mod[i]=="gjam1")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT1)[j]))>0)){
    table_comp$res[i]<-  mean(LtT1[[j]][[jm]]$trace[-c(1:burn)])
    table_comp$lweight[i]<- mean(LtT1[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT1[[j]][[jm]]$fit
    table_comp$err[i]<- LtT1[[j]][[jm]]$err
    }
  if((table_comp$mod[i]=="gjam2")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT2)[j]))>0)){
    table_comp$res[i]<- mean(LtT2[[j]][[jm]]$trace[-c(1:burn)])
    table_comp$lweight[i]<- mean(LtT2[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT2[[j]][[jm]]$fit
    table_comp$err[i]<- LtT2[[j]][[jm]]$err
    }
  if((table_comp$mod[i]=="gjam3")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT3)[j]))>0)){
    table_comp$res[i]<-  mean(LtT3[[j]][[jm]]$trace[-c(1:burn)])
    table_comp$lweight[i]<-  mean(LtT3[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT3[[j]][[jm]]$fit
    table_comp$err[i]<- LtT3[[j]][[jm]]$err
    }
  if((table_comp$mod[i]=="gjam4")&(length(grep(paste0("r_",table_comp$rv[i]),names(LtT4)[j]))>0)){
    table_comp$res[i]<- mean(LtT4[[j]][[jm]]$trace[-c(1:burn)])
    table_comp$lweight[i]<-  mean(LtT4[[j]][[jm]]$pkN)
    table_comp$fit_err[i]<- LtT4[[j]][[jm]]$fit
    table_comp$err[i]<- LtT4[[j]][[jm]]$err
    }

}





table_comp$Schar <- factor(table_comp$Schar, levels=c('S=20','S=50','S=80','S=100'))
table_comp<-table_comp[,-1]
p <- ggplot(data= table_comp)
q<- p +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5","10","20"),
                   labels=c("5","10","20"),limits=c("5","10","20"))+
  scale_y_continuous(name="Number of clusters")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,10))+
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 4,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q




q<- p +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5","10","20"),
                   labels=c("5","10","20"),limits=c("5","10","20"))+
  scale_y_continuous(name="Number of clusters")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,10))+
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 4,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q






q<- p +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(res), fill = as.factor(mod)))+
  scale_x_discrete(name="Parameters", breaks=c("5","10","20"),
                   labels=c("5","10","20"),limits=c("5","10","20"))+
  scale_y_continuous(name="Number of clusters")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,10))+
  labs(title="Ability to recover true number of groups for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+
  geom_hline(yintercept = 4,color = "red")+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q



p <- ggplot(data= table_comp)
q<- p +geom_boxplot(aes(x=as.factor(rv),y= as.numeric(lweight), fill = as.factor(mod))) +
  scale_x_discrete(name="Parameters", breaks=c("5","10","20"),
                   labels=c("5","10","20"),limits=c("5","10","20"))+
  scale_y_continuous(name="Number of clusters")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","GJAM2","GJAM3","GJAM4"))+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")+
  facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,0.3))+
  labs(title="Last p_k weight value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                            strip.background = element_blank(),
                                                            strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q






p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= fit_err, colour = mod) + geom_point( size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,max(table_comp$fit_err)))+
  labs(title="Fit error value for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                          strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q


p <- ggplot(data= table_comp)
q<- p +aes(x=rv,y= err, colour = mod) + geom_point( size = 2) +facet_wrap( ~ Schar, strip.position = "bottom", scales = "free_x")+ ylab("Posterior mean") + xlab("S and r values")+ylim(c(0,max(table_comp$err)))+
  labs(title="RMSE error between estimated and true covariance matrix for all models", caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples))+theme_bw()+theme(panel.spacing = unit(0, "lines"),
                                                                                                                                                                    strip.background = element_blank(),                                                                                                                                                                    strip.placement = "outside",legend.position = "top", plot.title = element_text(hjust = 0.5))
q
########################################################################################################################






