### Replicate paper
library(PointFore)
source("./R/MCS.R")
source("./R/plots.R")
source("./R/dgp.R")
source("./R/iden fcts.R")

#set seed
set.seed(2091987)

### Macro Paramters

# data generating process
dgp <- "patton"

# fast version
N <- 50
T_int <- c(50,100,250)
file_name <- "./results/standard.Rda"

inst.patton <- c("X", "lag(Y)","X^2", "lag(Y)^2","X^3", "lag(Y)^3")
inst.small <- c("X","lag(Y)")
parameter.vector <- c(-0.25,0)

simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,cov.est = iid,
         tests = list(
           list(description = "spline", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.patton, state = "X",
                cov.est ="HAC",centeredVcov = FALSE),
           list(description = "quantiles", iden.fct = quantiles ,model = probit_linear, instruments= inst.patton, state = "X",
                cov.est ="HAC",centeredVcov = FALSE),
           list(description = "expectiles", iden.fct = expectiles ,model = probit_linear, instruments= inst.patton, state = "X",
                cov.est ="HAC",centeredVcov = FALSE),
           list(description = "quantiles_inst", iden.fct = quantiles ,model = probit_linear, instruments= inst.small, state = "X",
                cov.est ="HAC",centeredVcov = FALSE),
           list(description = "expectiles_inst", iden.fct = expectiles ,model = probit_linear, instruments= inst.small, state = "X",
                cov.est ="HAC",centeredVcov = FALSE),
           list(description = "quantile", iden.fct = quantiles ,model = probit_constant, instruments= inst.small, state = NULL,
                cov.est ="iid",centeredVcov = FALSE),
           list(description = "expectile", iden.fct = expectiles ,model = probit_constant, instruments= inst.small, state = NULL,
                cov.est ="iid",centeredVcov = FALSE),
           list(description = "fixed_exp", iden.fct = fixed_exp, model=NULL, instruments= inst.small, state = NULL,
                cov.est =iid),
           list(description = "fixed", iden.fct = fixed, model=NULL, instruments= inst.small, state = NULL,
                cov.est =iid)
         ),forecasters = list(
           list(type = "Y", model = "linearprobit", parameter = parameter.vector)
         ))
