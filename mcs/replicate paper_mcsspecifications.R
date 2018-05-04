### Replicate paper

set.seed(2091987)

# Source functions

source("./library.R")

reinbrowsern <- F
### Macro Paramters

# data generating process
dgp <- "patton"


# full version
N <- 2000
file_name <- "./results/save_quantiles.Rda"


#fast version
# N <- 100
# file_name <- "./results/temp_quantiles.Rda"


########################## State dependent quantiles

# length of time series
T_int <- c(100, 250,1000)






    
    
    
simulate(save_to = file_name,
         T_int=T_int,
         N=N,
         tests = list(
           list(iden.fct = "quantile_model",model = "linear", instruments="forecast"),
           list(iden.fct = "quantile_model",model = "periodic", instruments="forecast"),
           list(iden.fct = "quantile_model",model = "break_logit_const", instruments ="forecast")
         ),forecasters = list(
           list(forecaster = "model", model_type = "linear", model_parameter = c(1,-1))
         ))
simulate(save_to = file_name,
         T_int=T_int,
         N=N,
         tests = list(
           list(iden.fct = "quantile_model",model = "linear", instruments="forecast"),
           list(iden.fct = "quantile_model",model = "periodic", instruments="forecast"),
           list(iden.fct = "quantile_model",model = "break_logit_const", instruments ="forecast")
         ),forecasters = list(
           list(forecaster = "model", model_type = "periodic", model_parameter = c(1,-1) )
         ))

    simulate(save_to = file_name,
             T_int=T_int,
             N=N,
             tests = list(
               list(iden.fct = "quantile_model",model = "linear", instruments="forecast"),
               list(iden.fct = "quantile_model",model = "periodic", instruments="forecast"),
               list(iden.fct = "quantile_model",model = "break_logit_const", instruments ="forecast")
             ),forecasters = list(
                 list(forecaster = "model", model_type = "break_logit_const", model_parameter=c(1,-1))
               ))


## show results

# read results
file_name <- "./results/save_quantiles.Rda"
#file_name <- "./results/temp_quantiles.Rda"
s <- readRDS(file_name)
head(s)

tab <- rejection_table(s, alpha_level = .05)

t(rejection_table(s, alpha_level = .05,T=100))
t(rejection_table(s, alpha_level = .05,T=250))
t(rejection_table(s, alpha_level = .05,T=1000))
