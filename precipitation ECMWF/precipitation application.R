
rm(list = ls())
source("library.R")
library(scales)

# load data
load("./precipitation ECMWF/precipitation data.RData")

#prec <- prec_Brussels
#Temp <- temp_Brussels

prec <- prec_London
Temp <- temp_London



N <- dim(prec)[1]


reinbrowsern<- F

#cur.forecast <- 'ENS'
#cur.forecast <- 'CNT'
cur.forecast <- 'HRES'

interval.time <- (dim(prec)[1]-N):dim(prec)[1]
Y <- prec[interval.time,"24h", c('Y')]
X <- prec[interval.time,"24h", c(paste0('X_',cur.forecast))]
Temp <- Temp[interval.time, "24h", c(paste0('Z_',cur.forecast))]


########## Helpers
extra <- 2
TT <- length(Y)-extra

#Uncomment line for robustness check. Uses only the data points applied in Patton and Timmermann 2007.
#TT <- 123-extra

T <- TT+extra


# Interval
int  <- (extra+1):T
int1 <- (extra):(T-1)
int2 <- (extra-1):(T-2)


####################################### Paper

###### MZ

# standard
res <- lm(Y~X);summary(res)$coef
test <- linearHypothesis(res, vcov = vcovHAC(res), c("(Intercept) = 0", "X = 1"))
test
signif(test$`Pr(>F)`[2],digits=2)

# coef 0.29 0.79
# test H_0 p_value < 10^-7


# plot standard gmm

plot.dat <- data.frame(X=X,Y=Y)
ggplot(plot.dat, aes(x=X,y=Y),ylim = c(0,25))+
  geom_point(alpha=.3)+geom_smooth(method=lm, se=F, size = 1)+
  geom_abline(slope=1,intercept=0,col='red',linetype=2, size = 1)+
  theme_classic(20)+xlim(c(0,25))+ylim(c(0,25))
#ggsave('pre_scatter.pdf', height = 5, width=5)



# standard with gmm

res <- estimate.functional(iden.fct = l1_MZ, model_type = 'MZ',
                           instruments = 'forecastonly',
                           state_variable = NULL, Y = Y, X=X)

summary(res$gmm)
linearHypothesis(res$gmm, c("Theta[1] = 0", "Theta[2] = 1"))


# with median regression

# standard
n.steps <- 30
for(cur.tau in seq(0.5,0.7,length.out = n.steps))
{
  res <- Rq(Y~X,tau = cur.tau)
  res
  test <- linearHypothesis(res, vcov. = vcov(res), c("Intercept = 0", "X = 1"))
  test
  cat("\n tau = ", cur.tau, 
      ", p-value = ", signif(test$`Pr(>Chisq)`[2],digits=2))
  if(signif(test$`Pr(>Chisq)`[2],digits=2)>0.01){cat("(* not rejected)")}
}

# best regression
res <- Rq(Y~X,tau = 0.624)
res

res <- rq(Y~X,tau = 0.624)
summary.rq(res, se='nid')
linearHypothesis(res, vcov. = summary(res)$vcov, c("Intercept = 0", "X = 1"))

##### Constant expectile

res <- estimate.functional(iden.fct = expectile_model, model_type = 'logit_const',
                           Y = Y, X=X)
summary(res$gmm)

#p-value < 10^-8
#coef

inv.logit(0.174)
#0.54-expectile

##### linear on forecast

# Expectile

res <- estimate.functional(iden.fct =   expectile_model ,model_type = 'linear',
                           state_variable = X, Y = Y, X=X)
summary(res$gmm)

#p-value 0.49

threshold <- median(X)
median(X)
#0.5433756

p <- plot.levels(res, show.p.value = F,plot.hist = FALSE, limits = c(0,15))
p <- p+ geom_density(data = data.frame(state.variable = X[X>=threshold]), 
                     aes(x=state.variable,y=..scaled..),  adjust=2/5,fill="green", alpha=.2)
p <- p+scale_x_continuous("predicted precipitation",limits = c(0,15))+
  scale_y_continuous("expectile level", limits = c(0,1))
p 

#ggsave('prec.pdf', height = 5, width=5)


sum(res$state.variable>10)
sum(res$state.variable>5)
max(res$state.variable)


# Quantile

res <- estimate.functional(iden.fct =   quantile_model ,model_type = 'linear',
                           state_variable = X, Y = Y, X=X)
summary(res$gmm)

#p-value 0.71


p <- plot.levels(res, show.p.value = F, limits = c(0,25))
p <- p+scale_x_continuous("precipitation")+scale_y_continuous("quantile level", limits = c(0,1))
p

#ggsave('prec_quantile.pdf', height = 5, width=5)


################## Descriptive Plots

### Time series

## Only 2013 

#one row

upper.limit <- 20
plot.list <- list()
i <- 1
plot.dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

plot.dat <- subset(plot.dat, year(Date)==2013)


#cut upper limit for plot
if(sum(plot.dat$Y>upper.limit)>0) {warning(paste('In total', sum(plot.dat$Y>upper.limit), 'points above ylim'))}
plot.dat$Y[plot.dat$Y>upper.limit]<- upper.limit
if(sum(plot.dat$X>upper.limit)>0) {warning(paste('In total', sum(plot.dat$X>upper.limit), 'points above ylim'))}
plot.dat$X[plot.dat$X>upper.limit]<- upper.limit



    sset.plot <- subset(plot.dat, month(Date)< 7 & year(Date)==2013)
      
      p <- ggplot(sset.plot)+
        geom_line(aes(x=Date,y=Y))+
        geom_point(aes(x=Date,y=X), color = 'red', size = 2, shape=4)+
        ylab('precipitation')+
        xlab('')+
        ylim(0,upper.limit)+
        theme_classic(15)+ coord_fixed(ratio=2)+scale_x_date(labels = date_format("%b %Y"))
      
p
ggsave(p, filename = 'prec2013b.pdf', width = 10, height = 3)



#two rows

upper.limit <- 20
plot.list <- list()
i <- 1
plot.dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

plot.dat <- subset(plot.dat, year(Date)==2013)


#cut upper limit for plot
if(sum(plot.dat$Y>upper.limit)>0) {warning(paste('In total', sum(plot.dat$Y>upper.limit), 'points above ylim'))}
plot.dat$Y[plot.dat$Y>upper.limit]<- upper.limit
if(sum(plot.dat$X>upper.limit)>0) {warning(paste('In total', sum(plot.dat$X>upper.limit), 'points above ylim'))}
plot.dat$X[plot.dat$X>upper.limit]<- upper.limit



for(cur.year in unique(year(plot.dat$Date)))
{
  for(cur.half in c(1,2))
  {
    if(cur.half == 1)
    {
      sset.plot <- subset(plot.dat, month(Date)< 7 & year(Date)==cur.year)
    } else
    {
      sset.plot <- subset(plot.dat, month(Date)>= 7 & year(Date)==cur.year)
    }
    
    if(dim(sset.plot)[1]>0)
    {  
      p <- ggplot(sset.plot)+
        geom_line(aes(x=Date,y=Y))+
        geom_point(aes(x=Date,y=X), color = 'red', size = 2, shape=4)+
        ylab('precipitation')+
        ylim(0,upper.limit)+
        xlab('')+
        theme_classic(15)+ coord_fixed(ratio=2)
      
      plot.list[[i]]<-p
      i<-i+1
    }
  }
}

p <- plot_grid(plotlist=plot.list, ncol=1)
p
ggsave(p, filename = 'prec2013.pdf', width = 10, height = 6)


# one plot only

upper.limit <- 20
plot.dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

plot.dat <- subset(plot.dat, year(Date)==2013)


#cut upper limit for plot
if(sum(plot.dat$Y>upper.limit)>0) {warning(paste('In total', sum(plot.dat$Y>upper.limit), 'points above ylim'))}
plot.dat$Y[plot.dat$Y>upper.limit]<- upper.limit
if(sum(plot.dat$X>upper.limit)>0) {warning(paste('In total', sum(plot.dat$X>upper.limit), 'points above ylim'))}
plot.dat$X[plot.dat$X>upper.limit]<- upper.limit



      p <- ggplot(plot.dat)+
        geom_line(aes(x=Date,y=Y))+
        geom_point(aes(x=Date,y=X), color = 'red', size = 2, shape=4)+
        ylab('')+
        ylim(0,upper.limit)+
        xlab('')+
        theme_classic(30)

ggsave(p, filename = 'prec2013b.pdf', width = 20, height = 8)


#Din A4

upper.limit <- 20
plot.list <- list()
i <- 1
plot.dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

#cut upper limit for plot
if(sum(plot.dat$Y>upper.limit)>0) {warning(paste('In total', sum(plot.dat$Y>upper.limit), 'points above ylim'))}
plot.dat$Y[plot.dat$Y>upper.limit]<- upper.limit
if(sum(plot.dat$X>upper.limit)>0) {warning(paste('In total', sum(plot.dat$X>upper.limit), 'points above ylim'))}
plot.dat$X[plot.dat$X>upper.limit]<- upper.limit



for(cur.year in unique(year(plot.dat$Date)))
{
  for(cur.half in c(1,2))
  {
    if(cur.half == 1)
    {
      sset.plot <- subset(plot.dat, month(Date)< 7 & year(Date)==cur.year)
    } else
    {
      sset.plot <- subset(plot.dat, month(Date)>= 7 & year(Date)==cur.year)
    }
    
    if(dim(sset.plot)[1]>0)
    {  
      p <- ggplot(sset.plot)+
      geom_line(aes(x=Date,y=Y))+
      geom_point(aes(x=Date,y=X), color = 'red', size = .5)+
        ylab('')+
        ylim(0,upper.limit)+
        xlab('')
      plot.list[[i]]<-p
      i<-i+1
    }
  }
}

p <- plot_grid(plotlist=plot.list, ncol=1)
ggsave(p, filename = 'prec_ts.pdf', width = 210, height = 297, units = "mm")


# 2 pages of Din A4

upper.limit <- 40
plot.list <- list()
i <- 1
plot.dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

#cut upper limit for plot
if(sum(plot.dat$Y>upper.limit)>0) {warning(paste('In total', sum(plot.dat$Y>upper.limit), 'points above ylim'))}
plot.dat$Y[plot.dat$Y>upper.limit]<- upper.limit
if(sum(plot.dat$X>upper.limit)>0) {warning(paste('In total', sum(plot.dat$X>upper.limit), 'points above ylim'))}
plot.dat$X[plot.dat$X>upper.limit]<- upper.limit



for(cur.year in unique(year(plot.dat$Date)))
{
  for(cur.half in c(1,2))
  {
    if(cur.half == 1)
    {
      sset.plot <- subset(plot.dat, month(Date)< 7 & year(Date)==cur.year)
    } else
    {
      sset.plot <- subset(plot.dat, month(Date)>= 7 & year(Date)==cur.year)
    }
    
    if(dim(sset.plot)[1]>0)
    {  
      p <- ggplot(sset.plot)+
        geom_line(aes(x=Date,y=Y))+
        geom_point(aes(x=Date,y=X), color = 'red', size = .5)+
        ylab('')+
        ylim(0,upper.limit)+
        xlab('')
      plot.list[[i]]<-p
      i<-i+1
    }
  }
}

p <- plot_grid(plotlist=plot.list, ncol=1)
ggsave(p, filename = 'prec_ts2.pdf', width = 210, height = 2*297, units = "mm")


#### Bias

# Plot bias vs. month. 
dat <- data.frame(X=X,Y=Y, Date = as.Date(names(Y), "%d-%m-%Y"))

bias.month <- numeric(12)
for(cur.month in 1:12)
{
  
  sset <- subset(dat, month(Date)==cur.month)
  bias.month[cur.month] <-  mean(sset$X)-mean(sset$Y)
    
}


plot.dat <- data.frame(month = as.factor(1:12), bias=bias.month)

ggplot(plot.dat, aes(x=month, y=bias))+geom_point(size=2)+geom_abline(slope = 0)
ggsave( filename = 'prec_biasmonth.pdf', height = 5, width=7)


# Plot bias vs. fcst bin (0-1, ..., 9-10, 10-15, 15-20, 20-25). 

dat <- data.frame(X=X,Y=Y)

dat$bin <-  1*data.table::between(dat$X, 0,1)+
            2*data.table::between(dat$X, 1,5)+
            3*data.table::between(dat$X, 5,10)+
            4*data.table::between(dat$X, 10,20)+
            5*data.table::between(dat$X, 20,40)
  

bias.bin <- numeric(max(dat$bin))
for(cur.bin in 1:length(bias.bin))
{
  
  sset <- subset(dat, bin==cur.bin)
  bias.bin[cur.bin] <-  mean(sset$X)-mean(sset$Y)
  
}


plot.dat <- data.frame(bin = as.character(1:length(bias.bin)), bias=bias.bin)

ggplot(plot.dat, aes(x=bin, y=bias))+geom_point(size=2)+geom_abline(slope = 0)+
  scale_x_discrete(labels=c("1" = "0-1", "2" = "1-5",
                            "3" = "5-10","4"="10-20","5"="20-40"),name="precipitation forecast")

ggsave( filename = 'prec_biasforecast.pdf', height = 5, width=7)



############################################ Test different specifications
# 
# 
# 
# ################# Constant models 
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model,t_model_type = 'logit_const',
#                             Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model,t_model_type = 'logit_const',
#                            instruments = 'forecast.sq',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model,t_model_type = 'logit_const',
#                            instruments = 'er.fo.ou',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model,t_model_type = 'logit_const',
#                            instruments = 'er.fo.ou.sq',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model,t_model_type = 'logit_const',
#                            instruments = 'forecast',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model,t_model_type = 'logit_const',
#                            instruments = 'forecast.sq',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model,t_model_type = 'logit_const',
#                            instruments = 'er.fo.ou.sq',
#                            state_variable = NULL, Y = Y, X=X)
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0.5"))
# 
# 
# 
# ################# Mincer Zarnowitz
# 
# 
# res <- estimate.functional(iden.fct = l1_MZ, t_model_type = 'MZ',
#                            instruments = 'er.fo.ou',
#                            state_variable = NULL, Y = Y, X=X)
# 
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0", "Theta[2] = 1"))
# 
# 
# res <- estimate.functional(iden.fct = l1_MZ, t_model_type = 'MZ',
#                            instruments = 'forecastonly',
#                            state_variable = NULL, Y = Y, X=X)
# 
# summary(res$gmm)
# linearHypothesis(res$gmm, c("Theta[1] = 0", "Theta[2] = 1"))
# 
# res <- lm(Y~X);summary(res)$coef
# linearHypothesis(res, vcov = vcovHAC(res), c("(Intercept) = 0", "X = 1"))
# 
# ######################### Break models ##########################
# 
# 
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break',
#                            instruments = "forecast",
#                            state_variable = seq(-1,1,length.out = length(Y[,1])), Y = Y, X=X)
# summary(res$gmm)
# 
# 
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break',
#                            instruments = "ff66",
#                            state_variable = seq(-1,1,length.out = length(Y[,1])), Y = Y, X=X)
# summary(res$gmm)
# 
# 
# res <- estimate.functional(iden.fct =  l1_quad_quad_model ,t_model_type = 'break',
#                            instruments = "ff66",
#                            state_variable = seq(-1,1,length.out = length(Y[,1])), Y = Y, X=X)
# summary(res$gmm)
# 
# 
# res <- estimate.functional(iden.fct =  l1_quad_quad_model ,t_model_type = 'break',
#                            instruments = "er.fo.ou",
#                            state_variable = seq(-1,1,length.out = length(Y[,1])), Y = Y, X=X)
# summary(res$gmm)
# 
# 
# 
# 
# winter.half <- month(as.Date(dimnames(Y)[[1]],format= "%d-%m-%Y")) %in% c(10,11,12,1,2,3)
# spring.only <- month(as.Date(dimnames(Y)[[1]],format= "%d-%m-%Y")) %in% c(3,4,5)
# summer.only <- month(as.Date(dimnames(Y)[[1]],format= "%d-%m-%Y")) %in% c(6,7,8)
# automn.only <- month(as.Date(dimnames(Y)[[1]],format= "%d-%m-%Y")) %in% c(9,10,11)
# winter.only <- month(as.Date(dimnames(Y)[[1]],format= "%d-%m-%Y")) %in% c(12,1,2)
# 
# # Break at half
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = int<(length(Y)/2), Y = Y, X=X)
# summary(res$gmm)
# 
# 
# # Break at winter
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = winter.half, Y = Y, X=X)
# summary(res$gmm)
# 
# 
# # only one term
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = winter.only, Y = Y, X=X)
# summary(res$gmm)
# 
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = spring.only, Y = Y, X=X)
# summary(res$gmm)
# 
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = summer.only, Y = Y, X=X)
# summary(res$gmm)
# 
# res <- estimate.functional(iden.fct =  l1_lin_lin_model ,t_model_type = 'break_logit_const',
#                            state_variable = automn.only, Y = Y, X=X)
# summary(res$gmm)
# 
# 
# 
# ####################### Lagged outcome as state variable
# res <- estimate.functional(iden.fct =   quantile_model ,t_model_type = 'current_level',
#                            state_variable = Y[-length(Y)], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res)
# 
# res <- estimate.functional(iden.fct =   l1_lin_lin_model ,t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Y[-length(Y[,1]),1], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res)
# 
# 
# res <- estimate.functional(iden.fct =   l1_quad_quad_model ,t_model_type = 'current_level',
#                            state_variable = Y[-length(Y[,1]),1], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1])
# 
# 
# res <- estimate.functional(iden.fct =   l1_quad_quad_model ,t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Y[-length(Y[,1]),1], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1])
# 
# 
# 
# ######################## Forecast as state variable
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2])
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'forecast.sq',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2])
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'er.fo.ou.sq',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2], t_model_type = "current_level")
# 
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2])
# linearHypothesis(res$gmm, c("Theta[2] = 0"))
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'forecast.sq',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2])
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'er.fo.ou.sq',
#                            state_variable = Y[,2], Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[,2])
# 
# ######################## Temperature as state variable
# 
# # lin lin
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'forecast',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'er.fo.ou.sq',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'lags2',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# # quad quad 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'forecast',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'lags2',
#                            state_variable = Temp, Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Temp)
# 
# 
# 
# ######################### Forecast error as state variable
# 
# # lin lin
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'forecast',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# res <- estimate.functional(iden.fct = l1_lin_lin_model, t_model_type = 'current_level',
#                            instruments = 'lags2',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# # quad quad
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'forecast',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'ff66',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# res <- estimate.functional(iden.fct = l1_quad_quad_model, t_model_type = 'current_level',
#                            instruments = 'lags2',
#                            state_variable = Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2] , Y = Y, X=X)
# summary(res$gmm)
# plot.levels(res,Y[-length(Y[,1]),1]-Y[-length(Y[,1]),2])
# 
# 
# 
# 
# 

