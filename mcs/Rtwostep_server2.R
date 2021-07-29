### Replicate paper
library(PointFore)


source("./R/MCS.R")
source("./R/plots.R")
source("./R/dgp.R")
source("./R/iden fcts.R")

#set seed
#set.seed(2091987)
 set.seed(3091987)
# set.seed(4091987)
# set.seed(5091987)

### Macro Paramters
# data generating process
dgp <- "patton"
center.cur <- FALSE

# fast version
N <- 500
T_int <- c(50,100,250,500,1000,2000,4000)
#file_name <- "./results/revision/friction_ser1.Rda"
file_name <- "./results/revision/friction_ser2.Rda"
# file_name <- "./results/revision/friction_ser3.Rda"
# file_name <- "./results/revision/friction_ser4.Rda"

inst.patton <-  c("X",      "lag(Y)",  "X^2",      "lag(Y)^2",  "X^3",      "lag(Y)^3")
inst.patton2 <- c("lag(X)", "lag(Y,2)","lag(X)^2", "lag(Y,2)^2","lag(X)^3", "lag(Y,2)^3")
inst.small <-  c("X",     "lag(Y)")
inst.small2 <- c("lag(X)","lag(Y,2)")
parameter.vector <- c(0,0)

simulate(save_to = file_name, T_int=T_int, N=N, dgp = dgp,cov.est = HAC,
         tests = list(
           list(description = "splineP", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.patton,cov.est ="HAC",centeredVcov = center.cur),
           list(description = "splineP_lag", iden.fct = spline.flex.median ,model = probit_linear, instruments= inst.patton2, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectilesP", iden.fct = expectiles ,model = probit_linear, instruments= inst.patton, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectilesP_lag", iden.fct = expectiles ,model = probit_linear, instruments= inst.patton2, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantilesP", iden.fct = quantiles ,model = probit_linear, instruments= inst.patton, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantilesP_lag", iden.fct = quantiles ,model = probit_linear, instruments= inst.patton2, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectiles", iden.fct = expectiles ,model = probit_linear, instruments= inst.small, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectiles_lag", iden.fct = expectiles ,model = probit_linear, instruments= inst.small2, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantiles", iden.fct = quantiles ,model = probit_linear, instruments= inst.small, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantiles_lag", iden.fct = quantiles ,model = probit_linear, instruments= inst.small2, state = "X",cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantile", iden.fct = quantiles ,model = probit_constant, instruments= inst.small, state = NULL,cov.est ="HAC",centeredVcov = center.cur),
           list(description = "quantile_lag", iden.fct = quantiles ,model = probit_constant, instruments= inst.small2, state = NULL,cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectile", iden.fct = expectiles ,model = probit_constant, instruments= inst.small, state = NULL,cov.est ="HAC",centeredVcov = center.cur),
           list(description = "expectile_lag", iden.fct = expectiles ,model = probit_constant, instruments= inst.small2, state = NULL,cov.est ="HAC",centeredVcov = center.cur),
           list(description = "mean", iden.fct = fixed_mean, model=NULL, instruments= inst.small, state = NULL),
           list(description = "mean_lag", iden.fct = fixed_mean, model=NULL, instruments= inst.small2, state = NULL),
           list(description = "median", iden.fct = fixed_median, model=NULL, instruments= inst.small, state = NULL),
           list(description = "median_lag", iden.fct = fixed_median, model=NULL, instruments= inst.small2, state = NULL)
         ),
         forecasters = list(
           list(type = "Y2_mean", model = "linearprobit", parameter = parameter.vector)
         ))

