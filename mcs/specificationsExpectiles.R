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
file_name <- "/home/quamet/patrick/PF/revision/specificationsYexpectile.Rda"

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
             list(description = "elin_Y", iden.fct = expectiles ,model = probit_linear, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE),
             list(description = "ebreak_Y", iden.fct = expectiles ,model = breakprobit, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE),
             list(description = "eperiod_Y", iden.fct = expectiles ,model = periodprobit, instruments= inst.fore, state = "lag(Y)",centeredVcov = FALSE)
           ),forecasters = list(
             list(type = "Y", model = "probit_linear", parameter = parameter.vector2),
             list(type = "Y", model = "breakprobit", parameter = parameter.vector),
             list(type = "Y", model = "periodprobit", parameter = parameter.vector)
           ))
}

