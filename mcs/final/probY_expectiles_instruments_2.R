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
N <- 2000

file_name <- "./results/final_probY_expectiles_instruments_2.Rda"


########################## State dependent expectiless

# length of time series
T_int <- c(100,250,1000)


instruments_all <-list(
  c("lag(Y)","X")
)

parameter.vector <- c(.1,.5)
parameter.vector2 <- c(.1,.25)

#file.remove(file_name)
for(inst.fore in instruments_all)
{
  
simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,
         tests = list(
           list(description = "elin_Y", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)"),
           list(description = "ebreak_Y", iden.fct = expectiles ,model = breakprobit, instruments= inst.fore, state = "lag(Y)"),
           list(description = "eperiod_Y", iden.fct = expectiles ,model = periodprobit, instruments= inst.fore, state = "lag(Y)")
         ),forecasters = list(
           list(type = "Y", model = "linearprobit", parameter = parameter.vector2),
           list(type = "Y", model = "breakprobit", parameter = parameter.vector),
           list(type = "Y", model = "periodprobit", parameter = parameter.vector)
         ))
}



## show results
s <- readRDS(file_name)

s <- rbind(s,readRDS("./results/final_probY_quantiles_instruments_2.Rda"))

rate <- .1
s$instruments<-as.character(s$instruments)

s$test<-as.character(s$test)


for(T.cur in c(100,250,1000))
{
  tab_1 <- rejection_table(subset(s, instruments==paste0(instruments_all[[1]],collapse = "|") & substr(test,1,1)=="e"),
                           alpha_level = rate,T=T.cur)
  tab_2 <- rejection_table(subset(s, instruments==paste0(instruments_all[[1]],collapse = "|") & substr(test,1,1)!="e"),
                           alpha_level = rate,T=T.cur)
  
  print(xtable(cbind(tab_1,tab_2)))
}




# fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.1.", T=250)
# lapply(instruments_all, FUN = fun2)
# #
# fun2 <- function(x) 1-rejection_table(results = subset(s, instruments==paste0(x,collapse = "|")), target = "hit.Theta.2.", T=250)
# lapply(instruments_all, FUN = fun2)

# 
# s$test<-as.character(s$test)
# s$forecaster<-as.character(s$forecaster)
# T.cur <- 1000
# test.cur <- "lin_Y"
# fore.cur <- "Ylinearprobit"
# 
# test.cur <- "break_Y"
# fore.cur <- "Ybreakprobit"
# 
# test.cur <- "period_Y"
# fore.cur <- "Yperiodprobit"
# fun3 <- function(x,...) plot(subset(s, instruments==paste0(x,collapse = "|") &
#                                   T==T.cur &
#                                   test==test.cur &
#                                   forecaster == fore.cur)$thh.Theta.1.,
#                          subset(s, instruments==paste0(x,collapse = "|") &
#                                   T==T.cur &
#                                   test==test.cur &
#                                   forecaster == fore.cur)$thh.Theta.2.,...)
# lapply(instruments_all, FUN = fun3, xlim=c(-2,2),ylim=c(-2,2))

