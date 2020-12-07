
# Source functions
source("./R/plots.R")
source("./R/dgp.R")
source("./R/iden fcts.R")
library(PointFore)
library(xtable)



instruments_all <-list(
  c("lag(Y)","lag(lag(Y))"),
  c("lag(Y)","X")
)


rate <- .1

fun1 <- function(x) rejection_table(subset(s, instruments==paste0(x,collapse = "|")), alpha_level = rate,
                                    T=1000)


# Table for quantiles -----------------------------------------------------

file_name <- "./results/revision/specificationsY.Rda"
s <- readRDS(file_name)
s$instruments<-as.character(s$instruments)

for(T.cur in c(100,250,1000))
{
  tab_1 <- rejection_table(subset(s, instruments==paste0(instruments_all[[2]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  tab_2 <- rejection_table(subset(s, instruments==paste0(instruments_all[[1]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  print(xtable(cbind(tab_1,tab_2)))
}


# table for expectiles ----------------------------------------------------
file_name <-
s <- rbind(
  readRDS("./results/revision/specificationsYexpectile.Rda"),
  readRDS("./results/revision/specificationsY.Rda")
)

for(T.cur in c(100,250,1000))
{
  tab_1 <- rejection_table(subset(s,
                                  instruments==paste0(instruments_all[[2]],collapse = "|") &
                                  substr(s$test,1,1)=="e"),
                           alpha_level = rate,T=T.cur)
  tab_2 <- rejection_table(subset(s, instruments==paste0(instruments_all[[2]],collapse = "|")&
                                    substr(s$test,1,1)!="e"),
                           alpha_level = rate,T=T.cur)
  print(xtable(cbind(tab_1,tab_2)))
}


# table for state X -------------------------------------------------------

file_name <- "./results/revision/specificationsX.Rda"
s <- readRDS(file_name)
for(T.cur in c(100,250,1000))
{
  tab_1 <- rejection_table(subset(s, instruments==paste0(instruments_all[[2]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  tab_2 <- rejection_table(subset(s, instruments==paste0(instruments_all[[1]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  print(xtable(cbind(tab_1,tab_2)))
}



# # Scatter plot with true confidence ellipses ------------------------------
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)
library(mvtnorm)


num.draws.stationary <- 5*10^6
T.cur <- c(100,250,1000)

ellipse.level <- .9
ll <- 0.05
hl <- 0.95


parameter.vector <- c(.1,.25)

xlimits <- c(-.3,.6)
ylimits <- c(-.3,0.9)


breaksy <- seq(-.25,.75,by=0.25)
breaksx <- seq(-.25,.5,by=0.25)

inst.fore <- c("lag(Y)","X")

ss <- subset(s,T %in% T.cur &
               forecaster=="Xlinear" &
               test=="lin_Y" &
               instruments == paste0(inst.fore,collapse = "|"))


# assign true values
ss$theta1 <- ss$theta.1
ss$theta2 <- ss$theta.2



cov <- estimate.true.cov(parameter.vector = parameter.vector,
                         nn = num.draws.stationary,
                         instr = inst.fore,
                         forecast.type = "X",
                         drop.until = 100)
cov


normdat <- NULL
for(T in T.cur)
{
  normdat.T <- data.frame(mvtnorm::rmvnorm(10^4, mean = parameter.vector, sigma = cov/T))
  normdat.T$T <- T
  normdat <- rbind(normdat, normdat.T)
}


colnames(normdat)<-c("x","y","T")

joint.dat <- NULL
joint.dat$theta1 <- c(normdat$x, ss$thh.Theta.1.)
joint.dat$theta2 <- c(normdat$y, ss$thh.Theta.2.)
joint.dat$T <- c(normdat$T,ss$T)
joint.dat$theoretic <- c(rep(TRUE,length(normdat$x)), rep(FALSE,length(ss$thh.Theta.1.)))
joint.dat <- data.frame(joint.dat)
joint.dat$interaction <- interaction(joint.dat$T,joint.dat$theoretic)

library(dplyr)
summary.dat <- joint.dat %>%
  group_by(interaction) %>%
  summarize(low1 = quantile(theta1,ll,na.rm = TRUE),
            high1 = quantile(theta1,hl,na.rm = TRUE),
            low2 = quantile(theta2,ll,na.rm = TRUE),
            high2 = quantile(theta2,hl,na.rm = TRUE),
            T=mean(T),
            theoretic=unique(theoretic))

summary.dat$T <- as.factor(summary.dat$T)
head(summary.dat)





hist_top <- ggplot()+theme_classic(10)+
  geom_errorbar(data= subset(summary.dat, theoretic==TRUE),
                aes(
                  ymax=high1,
                  ymin=low1,
                  x=T,
                  color=T))+
  geom_segment(data= subset(summary.dat, theoretic==FALSE),
               aes(
                 y=high1,
                 yend=low1,
                 x=T,xend=T,
                 color=T),size=10, alpha=.2)+coord_flip(ylim=xlimits)+#ylab(TeX("$\\theta_1$"))+
  guides(color=FALSE)+theme(axis.title.x=element_blank())+
  scale_y_continuous(breaks = breaksx)




empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme_void()



scatter <- ggplot(ss[sample(1:nrow(ss),1000),],aes(x=thh.Theta.1.,y=thh.Theta.2., color=factor(T)))+
  geom_point(alpha=0.3)+
  scale_y_continuous(breaks = breaksy)+
  scale_x_continuous(breaks = breaksx)+
  stat_ellipse(data=normdat, aes(x=x,y=y,color=factor(T)),level = ellipse.level)+
  xlab(TeX("$\\theta_0$"))+ylab(TeX("$\\theta_1$"))+
  coord_cartesian(xlim=xlimits, ylim=ylimits)+
  geom_point(aes(x=theta1,y=theta2),shape=3, size=8, color="black",show.legend = FALSE)+
  #  geom_point(data=data.frame(theta1=c(0,qnorm(0.3)),theta2=c(0)),aes(x=theta1,y=theta2),shape=4, size=4, color="black",show.legend = FALSE)+
  #  geom_text(data=data.frame(theta1=c(0,qnorm(0.3)),theta2=c(0,0),name=c("median","0.3-quantile")),aes(x=theta1,y=theta2,label=name), size=4, color="black",show.legend = FALSE,vjust=-1)+
  guides(color=FALSE)+theme_classic(10)


#scatter <- scatter+coord_fixed(xlim=xlimits, ylim=ylimits)

hist_right <- ggplot()+theme_classic(10)+
  geom_errorbar(data= subset(summary.dat, theoretic==TRUE),
                aes(
                  ymax=high2,
                  ymin=low2,
                  x=T,
                  color=T))+
  geom_segment(data= subset(summary.dat, theoretic==FALSE),
               aes(
                 y=high2,
                 yend=low2,
                 x=T,xend=T,
                 color=T),size=10, alpha=.2)+
  coord_cartesian(ylim = ylimits)+
  theme(axis.title.y=element_blank())+
  scale_y_continuous(breaks = breaksy)

#empty <- get_legend(hist_right+theme(legend.position = c(.2, .5)))


hist_right2 <- hist_right+guides(color=FALSE)

g <- grid.arrange(hist_top, empty, scatter, hist_right2, ncol=2, nrow=2,
                  widths=c(8, 2.5),
                  heights=c(5,8))


ggsave("scatter_first.pdf",  arrangeGrob(g),width = 8,height = 4)


