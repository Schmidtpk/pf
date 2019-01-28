### Replicate paper



# Source functions

source("./library.R")

reinbrowsern <- F
#reinbrowsern <- TRUE
### Macro Paramters

# data generating process
dgp <- "patton"



#fast version
N <- 10
T_int <- c(100,250,500)
file_name <- "./results/temp_spline_infofrictionX.Rda"
file.remove(file_name)


# full version
N <- 1000
T_int <- c(50,100,250,500,1000,2000,4000)
file_name <- "./results/final_spline_infofrictionX.Rda"



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


# Simulate

inst.fore <- c("X", "lag(X-Y)", "lag(X-Y)^2", "lag(lag(X-Y))", "lag(X)", "lag(lag(X-Y)^2)")
inst.fore2 <- paste0("lag(",inst.fore,")")

parameter.vector <- c(-0.25,0)

set.seed(2091987)
simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "spline", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "quantileslin_X", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "expectileslin_X", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "spline_lag", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore2, state = "lag(X)"),
           list(description = "quantileslin_X_lag", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore2, state = "lag(X)"),
           list(description = "expectileslin_X_lag", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore2, state = "lag(X)")
         ),forecasters = list(
           list(type = "Y2", model = "linearprobit", parameter = parameter.vector)
         ))

set.seed(3091987)
simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "spline", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "quantileslin_X", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "expectileslin_X", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "X"),
           list(description = "spline_lag", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.fore2, state = "lag(X)"),
           list(description = "quantileslin_X_lag", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore2, state = "lag(X)"),
           list(description = "expectileslin_X_lag", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore2, state = "lag(X)")
         ),forecasters = list(
           list(type = "Y2", model = "linearprobit", parameter = parameter.vector)
         ))

source("./library.R")
file_name <- "./results/final_spline_infofrictionX.Rda"

## show results
plot.on.T(file_name,nice.names = TRUE)+theme_classic(20)
#ggsave("mcs.pdf", height = 5, width=8)

### generate table of rejection rate
# s <- readRDS(file_name)
# rate <- .1
# s$instruments<-as.character(s$instruments)
# rejection_table(s, alpha_level = rate)



