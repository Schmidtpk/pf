### Replicate paper

set.seed(2091987)

# Source functions

source("./library.R")

reinbrowsern <- F
#reinbrowsern <- TRUE
### Macro Paramters

# data generating process
dgp <- "patton"



#fast version
N <- 20
T_int <- c(100,250)
file_name <- "./results/temp_spline_fullinfoX_2.Rda"
#file.remove(file_name)


# full version
N <- 2000
T_int <- c(50,100,250,500,1000,2000,4000)


file_name <- "./results/final_spline_fullinfoX_2.Rda"



############################ Define the spline identification function

spline.flex.median <- function(x,y,stateVariable,theta,model)
{
  
  e <- y-x
  
  model_variable<-y
  parameter <- theta
  
  node <-  c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  a <- cbind( gamm(1,parameter,model_variable),
              gamm(2,parameter,model_variable),
              gamm(3,parameter,model_variable),
              gamm(4,parameter,model_variable))
  
  b1 <- 2* a[,2] * node[1]-2*a[,1]*node[1] 
  c1 <-  a[,2]*node[1]^2-a[,1]*node[1]^2-b1 * node[1]
  
  b4 <- 2* a[,3] * node[3]-2*a[,4]*node[3]
  c4 <-  a[,3]*node[3]^2-a[,4]*node[3]^2-b4 * node[3]
  
  ff <- (
    (2* a[,1] * e + b1)       *(e <= node[1]) 
    + 2* a[,2] * e       *(e > node[1])*(e <= node[2]) 
    + 2* a[,3] * e       *(e > node[2])*(e <= node[3]) 
    + (2* a[,4] * e + b4)       *(e > node[3])
  )   
  
  return(ff)
}

gamm <- function(i,teta,y)
{
  ifelse(i==1 , f <- Tau(teta[1] + teta[2]*y - log (3)),
         ifelse(i==2 , f <- (1- gamm(1,teta,y)) * Tau(teta[3]+teta[4]*y-log(2)) ,
                ifelse(i==3 , f <- (1- gamm(1,teta,y)-gamm(2,teta,y)) * Tau(teta[5]+teta[6]*y-log(1)) ,
                       f <- 1-gamm(1,teta,y)-gamm(2,teta,y)-gamm(3,teta,y) )))
  return(f)
}

########################## State dependent expectiless




inst.fore <- c("X", "lag(X-Y)", "lag(X-Y)^2", "lag(lag(X-Y))", "lag(X)", "lag(lag(X-Y)^2)")

parameter.vector <- c(-0.25,0)


simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "spline", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "quantileslin_X", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "expectileslin_X", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "X")
           #list(description = "break_X", iden.fct = quantiles ,model = breakprobit, instruments= inst.fore, state = "X"),
           #list(description = "period_X", iden.fct = quantiles ,model = periodprobit, instruments= inst.fore, state = "X")
         ),forecasters = list(
           list(type = "Y", model = "linearprobit", parameter = parameter.vector)
         ))



source("./library.R")
file_name <- "./results/final_spline_fullinfoX_2.Rda"
## show results
plot.on.T(file_name,nice.names = TRUE)+theme_classic(20)
#ggsave("sizemcs.pdf", height = 5, width=8)


### generate table of rejection rate
 s <- readRDS(file_name)
 rate <- .1
 s$instruments<-as.character(s$instruments)
 rejection_table(s, alpha_level = rate)
 
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.1.", T=100)
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.1.", T=250)
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.1.", T=1000)
 
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.2.", T=100)
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.2.", T=250)
 1-rejection_table(s, alpha_level = rate,target = "hit.Theta.2.", T=1000)


# Scatter plot with true confidence ellipses ------------------------------
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)

num.draws.stationary <- 5*10^6
T.cur <- c(100,250,1000)

ellipse.level <- .9
ll <- 0.05
hl <- 0.95

inst.fore <- c("X", "lag(X-Y)", "lag(X-Y)^2", "lag(lag(X-Y))", "lag(X)", "lag(lag(X-Y)^2)")

parameter.vector <- c(-0.25,0)

xlimits <- c(-1.3,.8)
ylimits <- c(-.6,0.6)


breaksy <- seq(-.5,.5,by=0.25)
breaksx <- seq(-1,.5,by=0.25)


s <- readRDS(file_name)
ss <- subset(s,T %in% T.cur &
               forecaster=="Ylinearprobit" & 
               test=="quantileslin_X" &
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
  normdat.T <- data.frame(rmvnorm(10^4, mean = parameter.vector, sigma = cov/T))
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
  summarize(low1 = quantile(theta1,ll,na.rm = TRUE),
            high1 = quantile(theta1,hl,na.rm = TRUE),
            low2 = quantile(theta2,ll,na.rm = TRUE),
            high2 = quantile(theta2,hl,na.rm = TRUE),
            T=mean(T),
            theoretic=unique(theoretic))

summary.dat$T <- as.factor(summary.dat$T)
head(summary.dat)





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
  guides(color=FALSE)


#scatter <- scatter+coord_fixed(xlim=xlimits, ylim=ylimits)

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
  coord_cartesian(ylim = ylimits)+
  theme(axis.title.y=element_blank())+
  scale_y_continuous(breaks = breaksy)

#empty <- get_legend(hist_right+theme(legend.position = c(.2, .5)))


hist_right2 <- hist_right+guides(color=FALSE)

g <- grid.arrange(hist_top, empty, scatter, hist_right2, ncol=2, nrow=2, 
                  widths=c(8, 2.5), 
                  heights=c(5,8))


#ggsave("scatter_marginsX_ratiofixed.pdf",  arrangeGrob(g),width = 8,height = 4)



# Interpretation ----------------------------------------------------------

#simulate time series
dat_sim <- generate_data_y(DGP = "patton", T = 100000)
dat_sim <- generate_forecast(dat_sim, list(type = "X", model = "linear", parameter = parameter.vector))



#compute quantiles of state variable X
extreme.states <- quantile(dat_sim$X,probs=c(0.05,0.95))

parameter.cur <- c(-0.25,0.3)
linearprobit(extreme.states[1],parameter.cur)
#0.27
linearprobit(extreme.states[2],parameter.cur)
#0.48

# -> state-dependence such that forecasted quantile levels are between 0.27 and 0.48 with 90%
