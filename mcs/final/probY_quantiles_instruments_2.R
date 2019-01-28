### Replicate paper

set.seed(2091987)

# Source functions

source("./library.R")

reinbrowsern <- F
#reinbrowsern <- TRUE
### Macro Paramters

# data generating process
dgp <- "patton"

N <- 2000
file_name <- "./results/final_probY_quantiles_instruments_2.Rda"


########################## State dependent expectiless

# length of time series
T_int <- c(100,250,1000)


instruments_all <-list(
  c("lag(Y)","lag(lag(Y))"),
  c("lag(Y)","X")
)

parameter.vector <- c(.1,.5)
parameter.vector2 <- c(.1,.25)

#file.remove(file_name)
for(inst.fore in instruments_all)
{
  
simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "lin_Y", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)"),
           list(description = "break_Y", iden.fct = quantiles ,model = breakprobit, instruments= inst.fore, state = "lag(Y)"),
           list(description = "period_Y", iden.fct = quantiles ,model = periodprobit, instruments= inst.fore, state = "lag(Y)")
         ),forecasters = list(
           list(type = "Y", model = "linearprobit", parameter = parameter.vector2),
           list(type = "Y", model = "breakprobit", parameter = parameter.vector),
           list(type = "Y", model = "periodprobit", parameter = parameter.vector)
         ))
}


file_name <- "./results/final_probY_quantiles_instruments_2.Rda"

## show results
s <- readRDS(file_name)

rate <- .1
s$instruments<-as.character(s$instruments)

fun1 <- function(x) rejection_table(subset(s, instruments==paste0(x,collapse = "|")), alpha_level = rate,
                                    T=1000)
lapply(instruments_all, FUN = fun1)


for(T.cur in c(100,250,1000))
{
  tab_1 <- rejection_table(subset(s, instruments==paste0(instruments_all[[2]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  tab_2 <- rejection_table(subset(s, instruments==paste0(instruments_all[[1]],collapse = "|")),
                           alpha_level = rate,T=T.cur)
  
  print(xtable(cbind(tab_1,tab_2)))
}



1-rejection_table(results = subset(s, instruments==paste0(instruments_all[[1]],collapse = "|")), target = "hit.Theta.1.", T=100)

1-rejection_table(results = subset(s, instruments==paste0(instruments_all[[1]],collapse = "|")), target = "hit.Theta.1.", T=1000)



fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.1.", T=100)
lapply(instruments_all, FUN = fun2)
fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.2.", T=100)
lapply(instruments_all, FUN = fun2)


fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.1.", T=1000)
lapply(instruments_all, FUN = fun2)
fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.2.", T=1000)
lapply(instruments_all, FUN = fun2)


s <- readRDS(file_name)


ss <- s[
 s$test=="lin_Y" & s$forecaster=="Ylinearprobit" |
   s$test=="period_Y" & s$forecaster=="Yperiodprobit" |
   s$test=="break_Y" & s$forecaster=="Ybreakprobit" 
,]

# assign true values
ss$theta1 <- .1
ss$theta2 <- .5
ss$theta2[ss$test=="lin_Y" & ss$forecaster=="Ylinearprobit"]<-.25

#ss <- subset(ss,T==1000)

levels(ss$instruments) <- c(expression("w"["t"]*"=(1,Y"["t-1"]*",Y"["t-2"]*")"),
                            expression("w"["t"]*"=(1,Y"["t-1"]*",X"["t"]*")"))
levels(ss$forecaster) <- c("linear.", "break.", "periodic.")

ss$T <- as.factor(ss$T)

library(latex2exp)
ggplot(ss,aes(x=thh.Theta.1.,y=thh.Theta.2.,color=T,shape=T))+geom_point(alpha=.1)+
  facet_wrap(instruments~forecaster, labeller = label_parsed)+
  xlab(TeX("$\\theta_1$"))+ylab(TeX("$\\theta_2$"))+
  xlim(c(-1,1))+ylim(c(-1,1))+
  geom_point(aes(x=theta1,y=theta2),shape=3, size=12, color="black",show.legend = FALSE)+
  scale_color_manual(values=c("red","green","blue"))+
  scale_shape_manual(values=c(1,2,5))+
  guides(colour = guide_legend(override.aes = list(alpha = 1),title="      T"),
         shape = guide_legend(title="      T"))+theme_classic()
  

#ggsave("scatter_facets.pdf", height=4, width=7)


#### Scatter plot with true confidence ellipses
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)

num.draws.stationary <- 1000000
T.cur <- c(100,1000)

ellipse.level <- .9
ll <- 0.05
hl <- 0.95

file_name <- "./results/final_probY_quantiles_instruments.Rda"

s <- readRDS(file_name)
ss <- subset(s,T %in% T.cur &
               forecaster=="Ylinearprobit" & 
               test=="lin_Y" &
               instruments == paste0(c("lag(Y)","X"),collapse = "|"))


# assign true values
ss$theta1 <- .1
ss$theta2 <- .25
parameter.vector2 <- c(.1,.25)


cov <- estimate.true.cov(parameter.vector = parameter.vector2,nn = num.draws.stationary)
cov

normdat <- NULL
for(T in T.cur)
{
  normdat.T <- data.frame(rmvnorm(num.draws.stationary, mean = parameter.vector2, sigma = cov/T))
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

summary.dat <- joint.dat %>%
  group_by(interaction) %>%
  summarize(low1 = quantile(theta1,ll),
            high1 = quantile(theta1,hl),
            low2 = quantile(theta2,ll),
            high2 = quantile(theta2,hl),
            T=mean(T),
            theoretic=unique(theoretic))

summary.dat$T <- as.factor(summary.dat$T)
head(summary.dat)



xlimits <- c(-0.3,0.5)
ylimits <- c(-0.1,0.6)


hist_top <- ggplot()+
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
  guides(color=FALSE)+theme(axis.title.x=element_blank())




empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme_void()

scatter <- ggplot(ss,aes(x=thh.Theta.1.,y=thh.Theta.2., color=factor(T)))+geom_point(alpha=.1)+
  stat_ellipse(data=normdat, aes(x=x,y=y,color=factor(T)),level = ellipse.level)+
  xlab(TeX("$\\theta_0$"))+ylab(TeX("$\\theta_1$"))+
  coord_cartesian(xlim=xlimits, ylim=ylimits)+
  geom_point(aes(x=theta1,y=theta2),shape=3, size=12, color="black",show.legend = FALSE)+
  guides(color=FALSE)



hist_right <- ggplot()+
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
  coord_cartesian(ylim = ylimits)+theme(axis.title.y=element_blank())

empty <- get_legend(hist_right+theme(legend.position = c(.2, .5)))
plot(empty)
hist_right <- hist_right+guides(color=FALSE)

g <- grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, 
             widths=c(4, 1.5), 
             heights=c(2, 4))

#ggsave("scatter_margins1.pdf",  arrangeGrob(g))

