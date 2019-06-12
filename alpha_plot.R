
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



