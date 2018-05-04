source("library.R")


# load observations
Y<- read.csv("./gdp/observations.csv", sep=";", dec=",", header=TRUE)
# format variable 'second' 
Y$second <- as.numeric(gsub(",",".",as.character(Y$second)))

# read forecasts
X<- read.csv("./gdp/forecasts.csv", sep=";", dec=",", header=TRUE)
# format date
X$date <- as.Date(as.character(X$GBdate),"%Y%m%d")


####### Choose forecast closest to middle of the respective quarter
year(X$date)<-2000
X$comparison <-as.Date("02152000","%m%d%Y")
X$comparison[(month(X$date)>3)] <- as.Date("05152000","%m%d%Y")
X$comparison[(month(X$date)>6)] <- as.Date("08152000","%m%d%Y")
X$comparison[(month(X$date)>9)] <- as.Date("11152000","%m%d%Y")

X$diff<-abs(X$comparison-X$date)

X<-transform(X, 
             date.rank = ave(diff, DATE, 
                             FUN = function(x) rank(x, ties.method = "first")))
#drop other forecasts
X.new <- X[X$date.rank==1,]


####### merge observations and forecasts
Y$X <-gsub(":Q",".",Y$X)
X.new$t <- 0
X.new$t[1:(dim(X.new)[1]-1)] <- X.new$DATE[-1]
Y<-merge(X.new,Y,by.x = 't',by.y = 'X')

# drop unnecessary variables
Y<-Y[,c('t','DATE','gRGDPF1','first','second','recent')]


###### Choose second or most recent vintage as robustness check
#Y<-Y[,c('second','gRGDPF1')]
#Y<-Y[,c('recent','gRGDPF1')]
Y<-Y[,c('first','gRGDPF1')]

#estimate functional without state-dependence
res <- estimate.functional(iden.fct = quantile_model,model_type = 'logit_const',
                           state_variable = NULL, Y = Y[,1], X=Y[,2])
summary(res$gmm)


#estimate functional with state-dependence on y_{t-1}
res <- estimate.functional(iden.fct = quantile_model, model_type = "linear",
                           state_variable = Y[1:(length(Y[,1])-1),1], Y=Y[,1],X=Y[,2])

plot.levels(res, show.p.value = T,limits = c(-5,10))
summary(res$gmm)

plot.levels(res, show.p.value = F,limits = c(-5,10))+  scale_x_continuous("predicted growth rate", limits = c(-5,10))+
  theme_classic(20)



#estimate functional with state-dependence on x_{t}
res <- estimate.functional(iden.fct = quantile_model, model_type = "linear",
                           state_variable = Y[,2], Y=Y[,1],X=Y[,2])
summary(res$gmm)
plot.levels(res, show.p.value = F,limits = c(-5,10),Y = Y[,2])+  scale_x_continuous("predicted growth rate", limits = c(-5,10))+
  theme_classic(20)

#ggsave('greenbook_confidence.pdf', height = 5, width=5)


res$gmm$coefficients
res$gmm$vcov

confint(res$gmm, level = 0.9, lambda = FALSE)     

summary(res$gmm)


# Show wald test
linearHypothesis(res$gmm, c('Theta[2] = 0'), test = 'Chisq')



