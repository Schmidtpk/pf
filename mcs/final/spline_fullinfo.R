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
file_name <- "./results/temp_spline_fullinfo.Rda"
file.remove(file_name)


# full version
N <- 2000
T_int <- c(50,100,250,500,1000,2000,4000)
file_name <- "./results/final_spline_fullinfo.Rda"



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

parameter.vector <- c(qnorm(0.25),0)


simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "spline", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore, state = "lag(Y)"),
           list(description = "quantileslin_Y", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)"),
           list(description = "expectileslin_Y", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)")
           #list(description = "break_Y", iden.fct = quantiles ,model = breakprobit, instruments= inst.fore, state = "lag(Y)"),
           #list(description = "period_Y", iden.fct = quantiles ,model = periodprobit, instruments= inst.fore, state = "lag(Y)")
         ),forecasters = list(
           list(type = "Y", model = "linearprobit", parameter = parameter.vector)
         ))



source("./library.R")
file_name <- "./results/final_spline_fullinfo.Rda"
## show results
plot.on.T(file_name)
plot.on.T(file_name,only.complete.cases = TRUE)


### generate table of rejection rate
 s <- readRDS(file_name)
 rate <- .4
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

num.draws.stationary <- 1000000
T.cur <- c(100,250,1000)

ellipse.level <- .9
ll <- 0.05
hl <- 0.95

#true value for covariance estimation
parameter.vector2 <- c(-0.6744898,0)

inst.fore <- c("X", "lag(X-Y)", "lag(X-Y)^2", "lag(lag(X-Y))", "lag(X)", "lag(lag(X-Y)^2)")


xlimits <- c(-1.2,-0.4)
ylimits <- c(-.3,.3)



s <- readRDS(file_name)
ss <- subset(s,T %in% T.cur &
               forecaster=="Ylinearprobit" & 
               test=="quantileslin_Y" &
               instruments == paste0(inst.fore,collapse = "|"))


# assign true values
ss$theta1 <- ss$theta.1
ss$theta2 <- ss$theta.2


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

#ggsave("scatter_margins_patton.pdf",  arrangeGrob(g))


