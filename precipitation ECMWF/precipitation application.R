source("library.R")
library(scales)

# load data
load("./precipitation ECMWF/precipitation data.RData")

prec <- prec_London
Temp <- temp_London



N <- dim(prec)[1]

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





