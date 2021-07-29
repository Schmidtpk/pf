### Replicate paper

set.seed(2091987)

# Source functions
library(MCSPointFore)
library(PointFore)

reinbrowsern <- F
#reinbrowsern <- TRUE
### Macro Paramters

# data generating process
dgp <- "patton"

N <- 2000
file_name <- "/home/quamet/patrick/PF/revision/specificationsY.Rda"

########################## State dependent expectiless
# length of time series
T_int <- c(100,250,1000)


instruments_all <-list(
  c("lag(Y)","lag(lag(Y))"),
  c("lag(Y)","X")
)

parameter.vector <- c(.1,.5)
parameter.vector2 <- c(.1,.25)

for(inst.fore in instruments_all)
{
  simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,cov.est = "HAC",
           tests = list(
             list(description = "lin_Y", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE),
             list(description = "break_Y", iden.fct = quantiles ,model = breakprobit, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE),
             list(description = "period_Y", iden.fct = quantiles ,model = periodprobit, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE)
           ),forecasters = list(
             list(type = "Y", model = "probit_linear", parameter = parameter.vector2),
             list(type = "Y", model = "breakprobit", parameter = parameter.vector),
             list(type = "Y", model = "periodprobit", parameter = parameter.vector)
           ))
}


file_name <- "/home/quamet/patrick/PF/revision/specificationsX.Rda"

for(inst.fore in instruments_all)
{
  simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,cov.est = "HAC",
           tests = list(
             list(description = "lin_Y", iden.fct = quantiles ,model = probit_linear, instruments= inst.fore, state = "X",centeredVcov = FALSE),
             list(description = "break_Y", iden.fct = quantiles ,model = breakprobit, instruments= inst.fore, state = "X",centeredVcov = FALSE),
             list(description = "period_Y", iden.fct = quantiles ,model = periodprobit, instruments= inst.fore, state = "X",centeredVcov = FALSE)
           ),forecasters = list(
             list(type = "X", model = "linear", parameter = parameter.vector2)
           ))
}
