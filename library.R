rm(list=ls())
######################## Libraries

# new version of gmm had trouble with iterative and cue weighting matrix
# require(devtools)
# install_version("gmm", version = "1.5-2", repos = "http://cran.us.r-project.org" )

# install.packages("devtools")
#devtools::install_github("Schmidtpk/PointFore")
library(gmm)
library(PointFore)
library(aod)
library(lubridate)
library(ggplot2)
library(boot)
library(xtable)
library(car)
library(cowplot)
library(quantreg)
#library(rms)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(grid)
library(mvtnorm)



source("analyze.R")
source("simulate_single.R")
source("MCS.R")


breaklogit <- function (stateVariable, theta) 
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for break model")
  }
  return(boot::inv.logit((stateVariable>0)*theta[2] + theta[1]))
}

periodlogit <- function (stateVariable, theta) 
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for periodic model")
  }
  return(boot::inv.logit((sin(stateVariable*2*pi/2))*theta[2] + theta[1]))
}

breakprobit <- function (stateVariable, theta) 
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for break model")
  }
  return(pnorm((stateVariable>0)*theta[2] + theta[1]))
}

periodprobit <- function (stateVariable, theta) 
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for periodic model")
  }
  return(pnorm((sin(stateVariable*2*pi/4))*theta[2] + theta[1]))
}


linearprobit <- function (stateVariable, theta) 
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for linear probit model")
  }
  return(pnorm(stateVariable*theta[2] + theta[1]))
}


splines <- function(x,y,stateVariable, theta, model) #(parameter, x,y,model_variable, ...)
{
 
  if(length(theta)!=6){stop("Wrong number of parameters")}
  
  e <- y-x
  
  node <-  c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  a1 <- - inv.logit(theta[1]+theta[2]*stateVariable-log(3))
  
  
  a2 <- -(1-a1)* inv.logit(theta[3]+theta[4]*stateVariable-log(3))
  
  a3 <-  (1-a1-a2)* inv.logit(theta[5]+theta[6]*stateVariable-log(3))
  
  a4 <-   1-a1-a2-a3
  
  return(
    ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
           ifelse(e<node[2], abs(e)*a2,
                  ifelse(e<node[3], abs(e)*a3 ,
                         (node[3]-node[1])*a3 + abs(e)*a4
                  )))
  )
  
}






find.model.var <- function(test.model, data,int,int1)
{
  if(test.model=="linear" | test.model=="current_level" | test.model=="current_level2"){model_variable_t<-data$Y[int1]}
  else if(test.model %in% c("logistic0","logisticx","linearx")){model_variable_t<-data$X[int]}
  else if(test.model=="reallinear"){model_variable_t<-cbind(rep(1,length(data$Y[int1])),data$Y[int1])}
  else if (test.model == "break_logit_const"){model_variable_t<-(data$Y[int1]<0)}#(int<(T_max/2))}
  else if (test.model == "break_logit_constx"){model_variable_t<-(data$X[int]<0)}#(int<(T_max/2))}
  else if (test.model =="logit_const"){model_variable_t<-rep(0,length(int))}
  else if (test.model =="periodic" | test.model =="periodic_simple"){model_variable_t<-sin(data$Y[int1]*2*pi/2)}#sin(int*2*pi/period_t)}
  else if (test.model =="periodicx" | test.model =="periodic_simple"){model_variable_t<-sin(data$X[int]*2*pi/2)}#sin(int*2*pi/period_t)}
  else if (test.model =="periodic3"){model_variable_t<-int}
  else if (test.model =="spline"){model_variable_t<-data$Y[int]}
  else if (test.model =="splinex"){model_variable_t<-data$X[int]}
  else if (test.model =="spliney"){model_variable_t<-data$Y[int]}
  else if (test.model =="const.spline"){model_variable_t<-rep(0,length(int))}
  else {stop("test.model unknown")}

  return(model_variable_t)
}




estimate.functional.old <- function(iden.fct = quantile_model, model_type = 'current_level', 
                                state_variable=NULL, instruments = 'forecast', 
                                Y, X, 
                                demean.state = F,
                                wmatrix='optimal',
                                bandwidth="bwAndrews",
                                node = NULL,
                                gmm.type = "twoStep")
{
  ##### helper variables
  extra <- 3
  TT <- length(Y)-extra
  T <- TT+extra
  # Interval
  int  <- (extra+1):T
  int1 <- (extra):(T-1)
  int2 <- (extra-1):(T-2)
  int3 <- (extra-2):(T-3)
  
  
  #### demean state variable
  if(demean.state == T)
  {
    state_variable <- state_variable-mean(state_variable)
  }
  
  
  ############ Instruments
  if(class(instruments)=='character')
  {
     
    if(instruments=="constant")
    {
      w <- cbind(rep(1,TT))
                 
    } else if(instruments=='forecast')
    {
      
      w <- cbind(rep(1,TT), # constant
                 Y[int1], # lagged observation
                 X[int])  # forecast
      
    } else if(instruments=='forecastabs')
    {
      
      w <- cbind(rep(1,TT), # constant
                 Y[int1], # lagged observation
                 X[int],
                 abs(X[int]))  # forecast
      
    } else if (instruments=='outcome2')##############################
    {
      
      w <- cbind( rep(1,length(int)),                    #constant
                  Y[int1],
                  Y[int2])
      
      
    } else if (instruments=='ff66')##############################
    {
      
      w <- cbind( rep(1,length(int)),                    #constant
                  X[int],                   #current value of forecast
                  Y[int1]-X[int1],     #lagged value of forecast error       
                  (Y[int1]-X[int1])^2 ,      #lagged squared forecast error       
                  Y[int2]-X[int2],    # two lags forecast error
                  X[int1],     #lagged value of forecast 
                  (Y[int2]-X[int2])^2            #two lags of squared forecast error
      )
      
   
    }  else if (instruments=='patton.correct2')##############################
    {
      
      w <- cbind( rep(1,length(int)),                    #constant
                  X[int],                   #current value of forecast
                  Y[int1]-X[int1],     #lagged value of forecast error       
                  Y[int2]-X[int2],      #lagged squared forecast error       
                  Y[int1],    
                  Y[int2],
                  rep(3, length(int)),
                  rep(4, length(int))# substituted in gmm with lagged V(x,y,theta)
      )
      
      
    } else if (instruments=='patton.short')##############################
    {
      
      w <- cbind( rep(1,length(int)),                    #constant
                  X[int],                   #current value of forecast
                  Y[int1]-X[int1],     #lagged value of forecast error       
                  Y[int1],    
                  rep(3, length(int))
                )
      
      
    } else if (instruments=='lagff66')##############################
    {
      
      w <- cbind( rep(1,length(int)),                    #constant
                  X[int1],                   #current value of forecast
                  Y[int2]-X[int2],     #lagged value of forecast error       
                  (Y[int2]-X[int2])^2 ,      #lagged squared forecast error       
                  Y[int3]-X[int3],    # two lags forecast error
                  X[int2],     #lagged value of forecast 
                  (Y[int3]-X[int3])^2            #two lags of squared forecast error
                  
      )
      
      
      
      
    } else {stop('Instrument description unknown')}
 
    # If instruments are specified via character variables, use interval definition
    Y <- Y[int]
    X <- X[int]
    

    #if state vector use appropriate subvector
    if(is.null(dim(state_variable)))
      state_variable<- state_variable[(1+length(state_variable)-length(Y)) :length(state_variable)]
    else #use appropriate submatrix
      state_variable<- state_variable[(1+dim(state_variable)[1]-length(Y)) :dim(state_variable)[1],]
      
    
   } else
  {
    w <- instruments
  }
  
  
  ########### Identification function
  
  V <- function(theta,x,y,model_variable)
  {
    return(iden.fct(x=x,y=y,model_variable=model_variable, type=model_type , parameter=theta, node=node))
  }
  
  ######## Starting value
  
  # defines apropriate starting values for each model
  theta_0 <- find.theta0(model_type,state_variable)
  
  
  ######### Checks
  
  if(dim(w)[1]!=length(Y) | dim(w)[1]!=length(X)){stop('Wrong dimensions')}
  if(dim(w)[1]<length(theta_0)){stop('Not enough moment conditions')}
  
  
  
  # Determines the algorithm used in the GMM estimation (optim for multidimensional, nlminb for one-dimensional paramter space)
  if (length(theta_0)>1){optfct <- 'optim'} else { optfct <- 'nlminb'}
    

  

  res <- estimate.gmm(forecast = X,
                         observation = Y, 
                         instruments = w,
                         V = V,
                         theta_0 = theta_0,
                         model_variable = state_variable,
                         optfct=optfct,
                         wmatrix=wmatrix,
                         bandwidth = bandwidth,
                        gmm_type = gmm.type)


  
  return(list(gmm = res, 
              iden.fct = iden.fct, 
              t_model_type = model_type, 
              instruments = instruments,
              state.variable = state_variable,
              V=V))
}



find.theta0 <- function(t_model_type,state_variable=NULL)
{
  if (t_model_type == "const"){ theta_0 <-.5}
  else if (t_model_type == "logit_const"){theta_0 <- 0}
  else if (t_model_type %in% c( "current_level" ,"linear" ,"logisticx","logistic0","linearx","lineary")){theta_0 <- c(0,0)}
  else if (t_model_type == "reallinear"){theta_0 <- c(0.5,rep(0,dim(state_variable)[2]-1))}
  else if (t_model_type == "current_level2"){theta_0 <- c(0,0,0)}
  else if (t_model_type == "break"){theta_0 <- c(0,0,0)}
  else if (substr(t_model_type,1,17)== "break_logit_const"){theta_0 <- c(0,0)}
  else if (substr(t_model_type,1,8)== "periodic"){theta_0 <- c(0,0)}
  #else if (t_model_type == "periodic3"){theta_0 <- c(0,0,0)}
  else if (t_model_type == "MZ"){theta_0 <- c(0,0)}
  else if (substr(t_model_type,1,6)== "spline"){theta_0 <- rep(0,6)}
  else if (t_model_type =="const.spline"){theta_0 <- rep(0,3)}
  else {stop("unkown t_model_type")}

  return(theta_0)
}



# Plots the specification model given the GMM result res for the support of the state variable state.variable
plot.levels <- function(res, state.variable=NULL, conf.levels = c(0.6,0.9), show.p.value = T, limits = NULL,plot.hist=T, adjust.factor=1/5)
{
  
  if(is.null(state.variable)){state.variable<-res$state.variable}
  t_model_type <- res$t_model_type
  
  # safe coefficients
  tetastern <- c(coef(res$gmm))
  var_theta <- res$gmm$vcov
  
  # define function for level
  alpha <- function(y){return(level_model(model_variable = y, type=t_model_type, parameter = tetastern))}
  
  if(t_model_type=="current_level" | t_model_type == "linear" | t_model_type =="logisticx" | t_model_type =="logistic0")
  {
    mue_tau <- function(y){return(tetastern[1]+tetastern[2]*y)}
    var_tau <- function(y){return(var_theta[1,1]+var_theta[2,2]*y^2+2*var_theta[1,2]*y)}
  } else  if(t_model_type=="current_level2")
  {
    mue_tau <- function(y){return(level_model(model_variable=y, type="current_level2", parameter=tetastern))}
    var_tau <- function(y){return(var_theta[1,1]+var_theta[2,2]*y^2+2*var_theta[1,2]*y)}
  } else {stop("Can't handle this t_model_type")}
  
  con_interval <- function(y,quan){
    return(c(
      qnorm((1-quan)/2, mue_tau(y), sqrt(var_tau(y))),
      qnorm(1-(1-quan)/2,mue_tau(y), sqrt(var_tau(y)))
    ))
  }
  
  
  if(is.null(limits)){
    interval_y <- seq(quantile(state.variable, probs = 0.01),quantile(state.variable, probs = 0.99), length.out=20)
  } else {
    if(length(limits)!=2) {stop('Limits not well-defined')}
    interval_y <- seq(limits[1],limits[2], length.out=200)
}
  
    
  alpha_int <- numeric(length(interval_y))
  alpha_low <- numeric(length(interval_y))
  alpha_high <- numeric(length(interval_y))
  alpha_low2 <- numeric(length(interval_y))
  alpha_high2 <- numeric(length(interval_y))
  
  
  for ( i in 1:length(interval_y))
  {
    alpha_int[i] <- alpha(interval_y[i])
    alpha_low[i] <- inv.logit(con_interval(interval_y[i],conf.levels[2])[1])
    alpha_high[i] <- inv.logit(con_interval(interval_y[i],conf.levels[2])[2])
    alpha_low2[i] <- inv.logit(con_interval(interval_y[i],conf.levels[1])[1])
    alpha_high2[i] <- inv.logit(con_interval(interval_y[i],conf.levels[1])[2])
  }
  
  
  ######## Create plot of quantile levels
  
  plot_data <- data.frame(cbind(interval_y,alpha_int, alpha_low,alpha_high, alpha_low2,alpha_high2))
 
   p.quantile <- ggplot()+
    geom_line(data=plot_data, aes(x=interval_y, y=alpha_int), size=1.2)+
    geom_ribbon(data=plot_data, aes(x=interval_y,ymin=alpha_low,ymax=alpha_high), alpha=0.2)+
    geom_ribbon(data=plot_data, aes(x=interval_y,ymin=alpha_low2,ymax=alpha_high2), alpha=0.4)
  
  if(show.p.value==T){p.quantile <- p.quantile+ggtitle(paste0("p-value J-test = ",
                                                              signif(specTest(res$gmm)$test[2],digits = 2)))}
  

  if(!is.null(plot.hist))
    if(plot.hist==TRUE)
    p.quantile <- p.quantile + geom_density(data = data.frame(state.variable), aes(x=state.variable,y=..scaled..), adjust=adjust.factor, fill="green", alpha=.2)+xlim(limits)
  
  p.quantile <- p.quantile +scale_y_continuous("quantile level", limits=c(0,1))+
    theme_classic(20)
  #+  scale_x_continuous("state variable")
  
  p.quantile
}




######################## Flexible Fucntionals

level_model <- function(model_variable=NULL, type="const", parameter=NULL)
{
  if (type=="const")
  {
    if(length(parameter)!=1){stop("Wrong parameter")}
    if(parameter>=1 | parameter <=0){stop("Wrong parameter value for constant forecasting model")}
    
    return(parameter)
  } else if (type=="logit_const")
  {
    if(length(parameter)!=1){stop("Wrong parameter")}
    
    return(inv.logit(parameter))
  }else if (type=="reallinear")
  {

    if(is.null(model_variable)){stop("Model variable missing")}
    
    if(is.null(dim(model_variable)))#if model_variable vector of 2
    {
      if(length(model_variable)!=length(parameter)){stop("Wrong number of paramters or model_variables")}
      level.loc <- sum(model_variable*parameter)
    } else #model_variable is matrix
    {
      level.loc <- model_variable%*%parameter
      if(dim(model_variable)[1]==1){stop("wrong interpretation")}
      if(dim(model_variable)[2]!=length(parameter)){stop("Wrong number of paramters or model_variables")}
    }
          
    level.loc <- ifelse(level.loc<0,
           0.001,
           ifelse(level.loc>1,
           0.999,
           level.loc))
           
    return(level.loc)
  }else if (type=="current_level"| type == "linear" |type == "linearx" | type=="logisticx")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(model_variable*parameter[2]+parameter[1]))
  }else if (type=="logistic0")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(model_variable*parameter[2]+parameter[1])*(model_variable>0))
  }else if (type=="current_level2")
  {
    if(length(parameter)!=3){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(model_variable^2*parameter[3]+model_variable*parameter[2]+parameter[1]))
  } else if (type=="break")
  {
    if(length(parameter)!=3){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit( parameter[1]+ parameter[2]*(model_variable>parameter[3])))
    
  }else if (type=="break_logit_const" | type=="break_logit_constx")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit( parameter[1]+ parameter[2]*(model_variable)))
    
  }else if (type=="break_logit_const0")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit( parameter[1]+ parameter[2]*(model_variable>0)))
    
  }else if (type=="periodic" | type=="periodicx")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(parameter[1]+ parameter[2]*(model_variable)))
    
  }else if (type=="periodic2")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(parameter[1]+ parameter[2]*sin(model_variable*2*pi/2)))
    
  }else if (type=="periodic_simple")
  {
    if(length(parameter)!=2){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(parameter[1]+ parameter[2]*(model_variable))
    
  }else if (type=="periodic3")
  {
    if(length(parameter)!=3){stop("Wrong parameter")}
    if(is.null(model_variable)){stop("Model variable missing")}
    
    return(inv.logit(parameter[1]+parameter[2]*(sin(model_variable*2*pi/(1+abs(parameter[3]))))))
    
  } else {stop("Unknown model forecasting type")}
}

### expectile_flexible
l1_quad_quad <- function(theta,x,y,z,...)
{
  e <- y-x
  ff <- (alpha(y,theta)-(e<0)) *abs(e)
  return(ff)
}

l1_lin_lin_model<- function(x,y,model_variable , type,parameter,...)
{
  (y<=x)-level_model(model_variable=model_variable, type=type,parameter=parameter)
}

#for easier notation
quantile_model <- l1_lin_lin_model

l1_quad_quad_model<- function(x,y,model_variable , type,parameter,...)
{
abs((y<=x)-level_model(model_variable=model_variable, type=type,parameter=parameter))*(x-y)
}

#for easier notation
expectile_model <- l1_quad_quad_model



### SPLINE 1
sp.const <- function(parameter, x,y,...)
{
  if(length(parameter)!=3){stop("Wrong number of paramters")}
  
  node <- c(-2,0,2)
  
  a1 <- - abs(parameter[1])
  a2 <- - abs(parameter[2])
  a3 <-   1
  a4 <-   abs(parameter[3])  
  
  
  e <- y-x
  
  return(
  ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
  ifelse(e<node[2], abs(e)*a2,
  ifelse(e<node[3], abs(e)*a3 ,
                    (node[3]-node[1])*a3 + abs(e)*a4
         )))
  )
}

sp.const.med <- function(parameter, x,y,...)
{
  if(length(parameter)!=3){stop("Wrong number of paramters")}
  
  e <- y-x
  node <- c(median(e[which(e<0)]),0,median(e[which(e>0)])) 
  
  # For plotting without data
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  a1 <- - abs(parameter[1])
  a2 <- - abs(parameter[2])
  a3 <-   1
  a4 <-   abs(parameter[3])  

  
  return(
    ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
           ifelse(e<node[2], abs(e)*a2,
                  ifelse(e<node[3], abs(e)*a3 ,
                         (node[3]-node[1])*a3 + abs(e)*a4
                  )))
  )
}


sp.flex <- function(parameter, x,y,model_variable, node=NULL, ...)
{
  
  if(length(parameter)!=6){stop("Wrong number of paramters")}
  
  e <- y-x
  
  if(is.null(node))
    node <- c(-2,0,2)
  
  
  a1 <- - inv.logit(parameter[1]+parameter[2]*model_variable-log(3))
  
  
  a2 <- -(1-a1)* inv.logit(parameter[3]+parameter[4]*model_variable-log(3))
  
  a3 <-  (1-a1-a2)* inv.logit(parameter[5]+parameter[6]*model_variable-log(3))
  
  a4 <-   1-a1-a2-a3
  
  
  return(
    ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
           ifelse(e<node[2], abs(e)*a2,
                  ifelse(e<node[3], abs(e)*a3 ,
                         (node[3]-node[1])*a3 + abs(e)*a4
                  )))
  )
  
}



sp.flex.new <- function(parameter, x,y,model_variable, node=NULL, ...)
{
  
  if(length(parameter)!=6){stop("Wrong number of paramters")}
  
  e <- y-x
  
  if(is.null(node))
    node <- c(-2,0,2)
  
  
  a1 <- - inv.logit(parameter[1]+parameter[2]*model_variable)
  
  
  a2 <- - inv.logit(parameter[3]+parameter[4]*model_variable)
  
  a3 <-   inv.logit(parameter[5]+parameter[6]*model_variable)
  
  a4 <-   1
  
  
  return(
    ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
           ifelse(e<node[2], abs(e)*a2,
                  ifelse(e<node[3], abs(e)*a3 ,
                         (node[3]-node[1])*a3 + abs(e)*a4
                  )))
  )
  
}

sp.flex.new.med <- function(parameter, x,y,model_variable, ...)
{
  
  if(length(parameter)!=6){stop("Wrong number of paramters")}
  
  e <- y-x
  
  node <-  c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  
  a1 <- - inv.logit(parameter[1]+parameter[2]*model_variable)
  
  
  a2 <- - inv.logit(parameter[3]+parameter[4]*model_variable)
  
  a3 <-   inv.logit(parameter[5]+parameter[6]*model_variable)
  
  a4 <-   1
  
  
  return(
    ifelse(e<node[1], (node[2]-node[1])*a2 + abs(e)*a1,
           ifelse(e<node[2], abs(e)*a2,
                  ifelse(e<node[3], abs(e)*a3 ,
                         (node[3]-node[1])*a3 + abs(e)*a4
                  )))
  )
  
}




##### SPLINE 2
Tau <- function(y)
{
  return((1+exp(-y))^(-1))
}


gamm <- function(i,teta,y)
{
  ifelse(i==1 , f <- Tau(teta[1] + teta[2]*y - log (3)),
         ifelse(i==2 , f <- (1- gamm(1,teta,y)) * Tau(teta[3]+teta[4]*y-log(2)) , 
                ifelse(i==3 , f <- (1- gamm(1,teta,y)-gamm(2,teta,y)) * Tau(teta[5]+teta[6]*y-log(1)) ,
                       f <- 1-gamm(1,teta,y)-gamm(2,teta,y)-gamm(3,teta,y) )))
  return(f)
}

spline <- function(parameter,x,y, ...)
{
  e <- y-x

  #node <- c(-1,0,1)
  node <- c(mean(e[which(e<0)]),0,mean(e[which(e>0)])) 
  #node <- c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  # For plotting without data
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  a <- cbind( gamm(1,parameter,y),
              gamm(2,parameter,y),
              gamm(3,parameter,y),
              gamm(4,parameter,y))
  
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


splinex <- function(parameter,x,y, ...)
{
  e <- y-x
  
  
  #node <- c(-1,0,1)
  node <- c(mean(e[which(e<0)]),0,mean(e[which(e>0)])) 
  #node <- c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  # For plotting without data
  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }
  
  a <- cbind( gamm(1,parameter,x),
              gamm(2,parameter,x),
              gamm(3,parameter,x),
              gamm(4,parameter,x))
  
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

spline.flex <- function(parameter,x,y, model_variable,...)
{
  e <- y-x
  
  
  #node <- c(-1,0,1)
  node <- c(mean(e[which(e<0)]),0,mean(e[which(e>0)])) 
  #node <- c(median(e[which(e<0)]),0,median(e[which(e>0)]))  
  
  # For plotting without data
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

const.spline <- function(parameter,x,y,model_variable,...)
{
  if(sum(abs(model_variable))!=0){stop("const.spline does not have constant model variable")}
  
  return(spline.flex.median(c(parameter[1],0,parameter[2],0,parameter[3],0),x,y,model_variable,...))
}

spline.flex.median <- function(parameter,x,y, model_variable,node=NULL,...)
{
  e <- y-x
  
    #node <- c(-1,0,1)
    #node <- c(mean(e[which(e<0)]),0,mean(e[which(e>0)])) 
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


spline2 <- function(teta,x,y, node,...)
{
  
  e <- y-x

  
  a <- cbind( gamm(1,teta,y),
              gamm(2,teta,y),
              gamm(3,teta,y),
              gamm(4,teta,y))
  
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

iden <- function(theta, z)
{
  return(z)
}


transformation_paper <- function(theta,z)
{  
  res <- z
  res[,5]<- derivative(theta_0, res[,5])
  if (length(res[1,])>7)
  {res[8,]<- derivative(theta_0, res[8,])}
  
  return(res)
}

###################### Estimation with GMM
estimate.gmm<-function(forecast,observation,instruments,V, theta_0, model_variable=NULL, gmm_type= "twoStep"
                       ,bandwidth="bwAndrews",optfct='optim', wmatrix='optimal')
{
  #n:                        length of dataset
  #forecast:                (n*1)-vector
  #obs:                     (n*1)-vector
  #instruments:             (n*k)-matrix , where k is number of instruments
  #V:                        identification function with parameters(theta, forecast, observation, instruments). Theta is to be estimated.
  #t0                        starting value for optimization
  
  #gmm_type= c("cue","iterative","twoStep")
  #         "cue" continuously updated gmm, slowest
  #          "iterative" iterative gmm
  #          "twoStep" two-Step gmm, fastest, sometimes more robust
  
  #bandwidth= c("bwAndrews","bwNeweyWest")
  #          Different methods for the estimation of the covariance matrix
  
  T_sample <- length(forecast)
  
  #safe length of model variable
  if(is.null(dim(model_variable)))
    model.dim <- 1
  else
    model.dim<-dim(model_variable)[2]
  
  # Generates function g
  g <- function(theta, m_data)
  {
    x <- m_data[,1]
    y <- m_data[,2]
    
    z <- m_data[,3:(dim(m_data)[2]-model.dim)]
    
    model_variable <- m_data[,(dim(m_data)[2]-model.dim+1):(dim(m_data)[2])]
 
    n <- length(x)
    
    #substitute constant 2,3 rows with V and lagged V
    if(!is.null(dim(z)))
    {
      for (i in 1:dim(z)[2])
      {
        if(identical(z[,i],rep(2,n)))
          z[,i]<-V(theta=theta,x=x,y=y,model_variable=model_variable)
        if(identical(z[,i],rep(3,n)))
          z[,i]<-c(V(theta=theta,x=x,y=y,model_variable=model_variable)[-1],0)
        if(identical(z[,i],rep(4,n)))
          z[,i]<-c(V(theta=theta,x=x,y=y,model_variable=model_variable)[-c(1,2)],0,0)
        if(identical(z[,i],rep(5,n)))
          z[,i]<-c(V(theta=theta,x=x,y=y,model_variable=model_variable)[-c(1,2,3)],0,0,0)
      }
    }  
    
    f <- diag(as.vector(V(theta=theta,x=x,y=y,model_variable=model_variable)))%*%cbind(z)
    return(f)
  }
  
  
  if (gmm_type!="cue" && gmm_type!="iterative" && gmm_type!="twoStep")
  {
    print("unkonwn gmm_type")
    gmm_type="twoStep"
  }
  
  if(is.null(model_variable)){model_variable <- rep(0,length(forecast))}
  
  matrix_data <-cbind(forecast, observation, instruments, model_variable)

  if(exists("reinbrowsern")){if(reinbrowsern==T){browser()}}
  
  # executes the gmm method
  res <- gmm(g, x=matrix_data,t0=theta_0, bw=get(bandwidth) , type=gmm_type, optfct=optfct,wmatrix = wmatrix, prewhite = F)

  return(res)
}


###################### Data generation
generate_data_y <- function(DGP="patton", T=1000, extra = 100, forecasting_model_type = NULL, forecasting_parameter = NULL,...)
{
  Y <- numeric(T+extra)
  
  var <- numeric(T+extra)  
  mean <- numeric(T+extra)
  
  
  if(DGP=="gauss")
  {

    Y <- gen.gauss(T+extra,...)
    var <- NA
    mean <- NA
  } else 
  {
    
    
    a <- .5
    b <- 0.1
    c <- 0.8
    d <- 0.1
    
    Y[1]<-0
    
    var[1]<-1
    mean[1]<-0
    epsilon <- rnorm(T+extra,0,1)
  
    for (t in 2:(T+extra))
    {
      if (DGP=="patton")
      {
        if(t==2){
          Y[t]<-0
            } else {
            
            var[t] <- b+c*var[t-1]+d*var[t-1]*epsilon[t-1]^2
            mean[t] <- a * Y[t-1]
            Y[t] <- a * Y[t-1]+sqrt(var[t])*epsilon[t]
          }
      } else if (DGP=="AR1")
      {
        if(t==2){
          Y[t]<-0
        } else {
          
          var[t] <- 1
          mean[t] <- a * Y[t-1]
          
          Y[t] <- a * Y[t-1]+sqrt(var[t])*epsilon[t]
        }
      } else if (DGP=="AR2")
      {
        if(t==2){
          Y[t]<-0
        } else {
          
          var[t] <- 1
          mean[t] <-  a * Y[t-1]+ b * Y[t-2]
          Y[t] <- a * Y[t-1]+ b * Y[t-2] + sqrt(var[t])*epsilon[t]
        }
      } else {stop("DGP unknown")}
    }
  }

  return(data.frame(Y=Y[(extra+1):(T+extra)], var=var[(extra+1):(T+extra)],mean=mean[(extra+1):(T+extra)]))
}






generate_forecast <- function(data, forecaster)
#(Y, var=NULL,  forecasting_type="model", forecasting_model_type = "const", forecasting_parameter = .5,...)
{
  #generate forecaster
  T <- nrow(data)
  X <- numeric(T)
  Y <- data$Y
  var <- data$var
  mean <- data$mean 
  
  
  
  if (class(forecaster$model)=="character" & grepl("Y",forecaster$type))
    model.fct <- get(forecaster$model)
  
  # if(forecasting_model_type=="reallinear")
  #   model_variable_loc<-matrix(0,nrow = T,ncol=length(forecasting_parameter))
  
  
  # if (is.null(var)) {var <- numeric(T)}
  
  a <- .5
  b <- 0.1
  c <- 0.8
  d <- 0.1
  
  X[1]<-0
  X[2]<-0
  
  for (t in 3:(T))
  {
    
    if (forecaster$type =="X")
    {
      

      mean.now <- mean[t]
      sd.now <- sqrt(var[t])
      
      if(forecaster$model =="linear")
      {
        X[t] <- (mean.now + sd.now*forecaster$parameter[1]) /
              (1-sd.now*forecaster$parameter[2])    
      } else if (forecaster$model=="break")
      {
        stop("not ready")
        X[t] <- (mean.now + sd.now*forecaster$parameter[1])
        if(X[t]<0)
          X[t] <- X[t]-forecaster$parameter[2]*sd.now
      } else if (forecaster$model=="period")
      {
        stop("not ready")
        X[t] <- (mean.now + sd.now*forecaster$parameter[1])
        if(X[t]<0)
          X[t] <- X[t]-forecaster$parameter[2]*sd.now
      }
        
    } else if (forecaster$type =="Y")
    {
      
      
      model.var.now <- Y[t-1]
      
      level_now <- model.fct(model.var.now, forecaster$parameter)
      
      X[t] <- mean[t] +  qnorm(level_now) * sqrt(var[t])     
      
    } else if (forecaster$type =="X2")
    {
      
      mean.now <- a^2*Y[t-2]
      sd.now <- sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1]) 
      
      X[t] <- (mean.now + sd.now*forecaster$parameter[1]) /
        (1-sd.now*forecaster$parameter[2])    
      
    } else if (forecaster$type =="Y2")
    {
      
      model.var.now <- Y[t-2]
      level_now <- model.fct(model.var.now, forecaster$parameter)
      
      mean.now <- a^2*Y[t-2]
      sd.now <- sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1]) 
      
      X[t] <- mean.now +  qnorm(level_now) * sd.now     
      
    } else if (forecaster$type =="patton")
    {
      
      
      X[t] <- .5*Y[t-1]+sqrt(var[t])*.25
      
      
    } else if (forecaster$type =="patton2")
    {

      
      X[t] <- .25*Y[t-2]+sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1])*.25

      
    } else {stop("forecasting type unknown")}
  }
  
  
  return(data.frame(Y=Y,X=X))
}



#old version
generate_data_x <- function(Y, var=NULL,  forecasting_type="model", forecasting_model_type = "const", forecasting_parameter = .5,...)
{
  
  T <- length(Y)
  X <- numeric(T)
  
  model_variable_loc=numeric(T)
  if(forecasting_model_type=="reallinear")
    model_variable_loc<-matrix(0,nrow = T,ncol=length(forecasting_parameter))
  
  
  if (is.null(var)) {var <- numeric(T)}
  
  a <- .5
  b <- 0.1
  c <- 0.8
  d <- 0.1
  
  X[1]<-0
  
  for (t in 2:(T))
  {
    if (forecasting_type =="patton")
    {
      
      
      X[t] <- .5*Y[t-1]+sqrt(var[t])*.25
      
      
    } else if (forecasting_type =="patton2")
    {
      
      if(t>2)
      {
        X[t] <- .25*Y[t-2]+sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1])*.25
      } else {X[t] <- .5*Y[t-1]+sqrt(var[t])*.25}
      
    } else if (forecasting_type=="model")
    {
      if(forecasting_model_type=="current_level"| forecasting_model_type == "linear" |forecasting_model_type=="current_level2"){model_variable_loc[t]<-Y[t-1]}
      else if (forecasting_model_type == "reallinear"){model_variable_loc[t,]<-c(1,Y[t-1])}#(t<(T/2))}
      else if (forecasting_model_type == "break_logit_const"){model_variable_loc[t]<-(Y[t-1]<0)}#(t<(T/2))}
      else if (forecasting_model_type =="logit_const"){model_variable_loc[t]<-0}
      else if (forecasting_model_type =="periodic" | forecasting_model_type == "periodic_simple"){model_variable_loc[t]<-sin(Y[t-1]*2*pi/2)}#sin(t*2*pi/period_f)}
      else if (forecasting_model_type =="periodic3"){model_variable_loc[t]<-t}
      else {stop("forecasting_model_type unknown")}
      
      model.var.now <- model_variable_loc[t]
      if(forecasting_model_type=="reallinear")
        model.var.now <- model_variable_loc[t,]
      
      level_now <- level_model(model_variable=model.var.now, type=forecasting_model_type, parameter=forecasting_parameter)
      
      X[t] <- a*Y[t-1] +  qnorm(level_now) * sqrt(var[t])
      
      
      
    }else if (forecasting_type=="modelx")
    {
      if(forecasting_model_type=="current_level"| forecasting_model_type == "linear" |forecasting_model_type=="current_level2"){model_variable<-function(x) x}
      #else if (forecasting_model_type == "reallinear"){model_variable_loc[t,]<-c(1,X[t])}#(t<(T/2))}
      else if (forecasting_model_type == "break_logit_const"){model_variable<-function(x) x<0}#(t<(T/2))}
      else if (forecasting_model_type =="logit_const"){model_variable<-function(x) 0}
      else if (forecasting_model_type =="periodic" | forecasting_model_type == "periodic_simple"){model_variable<-function(x) sin(x*2*pi/2)}#sin(t*2*pi/period_f)}
      # else if (forecasting_model_type =="periodic3"){model_variable_loc[t]<-t}
      else {stop("forecasting_model_type unknown")}
      
      target_fun <- function(x)
      {
        return(abs(x - a*Y[t-1] +  qnorm(
          level_model(model_variable=model_variable(x), type=forecasting_model_type, parameter=forecasting_parameter)
        ) * sqrt(var[t])))
      }
      
      X[t] <- optimize(f = target_fun,interval = c(-100,100))$minimum
      
      model_variable_loc[t] <- X[t]
      
    } else if (forecasting_type=="model2")
    {
      if(t==2){X[t]<-a*Y[t-1]
      model_variable_loc[t]<-0}
      else{
        
        if(forecasting_model_type=="linear" | forecasting_model_type=="current_level"|forecasting_model_type=="current_level2"){model_variable_loc[t]<-Y[t-2]}
        else if (forecasting_model_type == "break_logit_const"){model_variable_loc[t]<-((t-1)<(T/2))}
        else if (forecasting_model_type =="logit_const"){model_variable_loc[t]<-0}
        else if (forecasting_model_type =="periodic" | forecasting_model_type == "periodic_simple"){model_variable_loc[t]<-sin((t-1)*2*pi/period_f)}
        else if (forecasting_model_type =="periodic3"){model_variable_loc[t]<-(t-1)}
        else {stop("forecasting_model_type unknown")}
        
        level_now <- level_model(model_variable=model_variable_loc[t], type=forecasting_model_type, parameter=forecasting_parameter)
        
        sd.cur <- sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1]) 
        
        mean.cur <- a^2*Y[t-2]
        
        X[t] <- mean.cur +  qnorm(level_now) * sd.cur
      }
      
    } else {stop("forecasting type unknown")}
  }
  
  
  return(data.frame(Y=Y,X=X,model_variable=model_variable_loc))
}










