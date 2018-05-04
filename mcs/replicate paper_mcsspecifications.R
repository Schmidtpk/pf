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




forecasters <- list(
list(forecaster = "model", model_type = "break_logit_const0", model_parameter=c(1,-1)),
list(forecaster = "model", model_type = "linear", model_parameter = c(1,-.3)),
list(forecaster = "model", model_type = "periodic2", model_parameter = c(1,-1) )
     )

p <- ggplot(data = data.frame(state = 0), mapping = aes(x = x))+
  geom_density(data=data.frame(x=generate_data_y(T = 100000)$Y),fill="green", alpha=.1)

interval <- seq(-3,3, length.out = 1000)


count <- 0
for (forecaster.cur in forecasters)
{
  #Plot specification models
  plot.fun <-  function(x) level_model(model_variable = x, 
                                      type = forecaster.cur$model_type,
                                       parameter = forecaster.cur$model_parameter)

  
  
  plot.y <- sapply(interval,plot.fun)
  
  
  
  count <- count+1
  p <- p+geom_line(data = data.frame(state=interval,quantile.level = plot.y, specification = factor(count)),size=1, aes(x=state,y=quantile.level,color=specification))
}
plot(p+xlab("state variable")+ylab("quantile level")+
       theme_classic(20)+
       xlim(-3,3)+ylim(0,1)+
        scale_colour_manual("Specification",
                              labels = c("break","logistic", "periodic"), 
                              values = c("green","red", "blue"))
     )
     
     ggsave("quantilelevel.pdf", height=4, width=8)
