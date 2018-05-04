### Replicate paper

set.seed(2091987)

# Source functions

source("./library.R")

### Macro Paramters

# data generating process
dgp <- "patton"

# Full version
N <- 2000
file_name <- "./results/save_spline_friction.Rda"
T_int <- c(50,100,250,500,1000,2000,4000)


# Fast version
# N <- 50
# file_name <- "./results/temp_spline_friction.Rda"
# T_int <- c(100,250,1000)



#################################### Replicate MCS ##############################

########################## Spline comparison with constant expectile - full information



    simulate_single(forecasters = list(
                              list(forecaster="model2",model_type="linear",model_parameter=c(logit(pnorm(0.25)),0))
                                      ),
                    tests = list(
                              #  list(model = "linear", instruments="forecast",iden.fct="expectile_model"),
                              #  list(model = "linear", instruments="forecast",iden.fct="quantile_model"),
                              list(model = "linearx", instruments="ff66",iden.fct="expectile_model"),
                              list(model = "linearx", instruments="ff66",iden.fct="quantile_model"),
                              list(model = "spline", instruments="ff66", iden.fct ="spline.flex.median"),
                              list(model = "linearx", instruments="lagff66",iden.fct="expectile_model"),
                              list(model = "linearx", instruments="lagff66",iden.fct="quantile_model"),
                              list(model = "spline", instruments="lagff66", iden.fct ="spline.flex.median")
                                ),
                    
             save_to = file_name,  
             T.int=T_int,
             N=N,
             dgp=dgp)




### show results

#load results
#file_name <- "./results/save_spline_friction.Rda"
file_name <- "./results/temp_spline_friction.Rda"

plot.on.T(file_name,0.05,nice.names = TRUE)+theme(legend.position="right")+theme_classic(20)

#ggsave("mcs.pdf", height = 5, width=8)

