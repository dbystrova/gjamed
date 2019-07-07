
library(cowplot)


p1 <- ggplot(data = data.frame(x = c(-10, 10)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 2)) + scale_fill_grey()+
  stat_function(fun = dnorm, n = 101, args = list(mean = 4, sd = 2)) +
  ylab("") +
  scale_y_continuous(breaks = NULL)
p1


x <- seq(-4, 4, length=100)

hx <- rnorm(10000, mean=1, sd=2)
df<- as.data.frame(hx)
df$type<- 1
df2<- as.data.frame(rnorm(10000, mean=4, sd=2))
names(df2)<-c("hx")
df2$type<- 2

df_all<- rbind(df, df2)
groupColors <- c(a="bisque1", b="blue")

p1 <- ggplot(data = df_all, aes(hx, fill=as.factor(type))) + geom_density(alpha=0.3, adjust=3)+scale_fill_manual(values=c("bisque1","darkblue"))+
  ylab("") + scale_y_continuous(breaks = NULL)+ 
p1
