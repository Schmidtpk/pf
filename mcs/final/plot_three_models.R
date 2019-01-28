library(ggplot2)
library(PointFore)

##### Plot the three forecasts
forecasters <- list(
  list(forecaster = "model", model_type = "breakprobit", model_parameter=c(.1,.5)),
  list(forecaster = "model", model_type = "probit_linear", model_parameter = c(.1,.25)),
  list(forecaster = "model", model_type = "periodprobit", model_parameter = c(.1,.5))
)

# state is lag(Y)
p <- ggplot(data = data.frame(state = 0), mapping = aes(x = x))+
  geom_density(data=data.frame(x=generate_data_y(T = 100000)$Y),fill="green", alpha=.1)



interval <- seq(-3,3, length.out = 1000)


count <- 0
for (forecaster.cur in forecasters)
{
  #Plot specification models
  plot.fun <-  forecaster.cur$model_type
  
  
  
  plot.y <- sapply(interval,plot.fun,theta=forecaster.cur$model_parameter)
  
  
  
  count <- count+1
  p <- p+geom_line(data = data.frame(state=interval,quantile.level = plot.y, specification = factor(count)),size=1, aes(x=state,y=quantile.level,color=specification))
}

plot(p+xlab("state variable")+ylab("quantile level")+
       theme_classic(20)+
       xlim(-3,3)+ylim(0,1)+
       scale_colour_manual("Specification",
                           labels = c("break","linear", "periodic"), 
                           values = c("green","red", "blue"))
)

#ggsave("quantilelevel.pdf", height=3, width=8)


#### For forecast

##### Plot the three forecasts
forecasters <- list(
  list(forecaster = "model", model_type = "probit_linear", model_parameter = c(.1,.25))
)

# state is X
p <- ggplot(data = data.frame(state = 0), mapping = aes(x = x))+
  geom_density(data = data.frame(x=
                                   generate_forecast(
                                     data = generate_data_y(T = 100000),
                                     forecaster=list(type = "X", model = "linear", parameter = forecasters[[1]]$model_parameter))$X),
               fill="green", alpha=.1)



interval <- seq(-3,3, length.out = 1000)


count <- 0
for (forecaster.cur in forecasters)
{
  #Plot specification models
  plot.fun <-  forecaster.cur$model_type
  
  
  
  plot.y <- sapply(interval,plot.fun,theta=forecaster.cur$model_parameter)
  
  
  
  count <- count+1
  p <- p+geom_line(data = data.frame(state=interval,quantile.level = plot.y, specification = factor(count)),size=1, aes(x=state,y=quantile.level,color=specification))
}
plot(p+xlab("state variable")+ylab("quantile level")+
       theme_classic(20)+
       xlim(-3,3)+ylim(0,1)+
       theme(legend.position="top")+
       scale_colour_manual("Specification",
                           labels = c("linear probit"), 
                           values = c("red"))
)

#ggsave("quantilelevel_X.pdf", height=4, width=6)
