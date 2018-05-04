########### Setup

source("library.R")

library(reshape)
library('VGAM')

test_names <- list(
  'spline.flex.medianspliney'="Spline y",
  'spline.flex.mediansplinex'="Spline x",
  'quantile_modellogisticx'="Quantile",
  'expectile_modellogisticx'="Expectile"
)

test_names_const <- list(
  'const.splinesplinex'="Spline",
  'quantile'="Quantile",
  'expectile'="Expectile"
)

  my_labeller <- function(variable,value){
    return(test_names[value])
  }

  my_labeller_const <- function(variable,value){
    return(test_names_const[value])
  }
  


  ############################################################################
  #######################      Constant spline     ###########################
  ############################################################################

  
  ######## ff66 ###########################################
  
  
file <- "./results/supp_constspline.Rda"

cur.instruments <- "ff66"

# simulate_single(
#   
#   T.int = c(500,1000,2000,3000,4000,5000),
#   tests = list(
# 
#            list(iden.fct = "const.spline", model="const.spline", instruments=cur.instruments),
#            list(iden.fct = "quantile_model", model="logit_const", instruments=cur.instruments),
#            list(iden.fct = "expectile_model", model="logit_const", instruments=cur.instruments)
# 
#          ), N=10, save_to = file)


res <- readRDS(file)

res <- res[res$T>=1000,]

limits <- c(-1,1)

#copy results to same variable as spline
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$Theta.1.[grepl("expectile",res$test)]


#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#move from theta1 to theta2
res$thh.Theta.2.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("quantile",res$test)]<-NA

#rename
res$test <- as.character(res$test)
res$test[grepl("expectile",res$test)]<-"expectile"
res$test[grepl("quantile",res$test)]<-"quantile"
res$test[grepl("spline",res$test)]<-"spline"
res$test <- as.factor(res$test)

limits <- c(-1,1)

#plot
plot.single(save=res,theta0 = c(0,0,0),y.var = c("thh.Theta.1.","thh.Theta.2.","thh.Theta.3."),
            limits = limits,ncol=3,show.legend = "bottom",
            delete.y.lab=c(2,3))
ggsave("supp_splines_const_T.pdf",width = 5,height = 3)




## Histogram for T = 1000

file <- "./results/supp_constspline1000.Rda"

cur.instruments <- "patton.short"
# simulate_single(delete.old.results = TRUE,
# 
#   T.int = c(1000),
#   tests = list(
# 
#                list(iden.fct = "const.spline", model="const.spline", instruments=cur.instruments),
#                list(iden.fct = "quantile_model", model="logit_const", instruments=cur.instruments),
#                list(iden.fct = "expectile_model", model="logit_const", instruments=cur.instruments)
# 
#   ), N=1000, save_to = file)


#### boxplots

res <-readRDS(file)


#copy results to same variable as spline
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$Theta.1.[grepl("expectile",res$test)]

#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#safe in long format
ress <- melt(res,measure.vars = grep("thh.Theta",names(res)))

#drop NAs
ress <- ress[!is.na(ress$value),]

limits <- c(-2,2)
ggplot(ress,aes(x=variable,y=value,fill=test))+geom_boxplot(outlier.shape = 1)+
  theme(
    axis.title.x = element_blank(),legend.position = "none",
    axis.text.x = element_blank())+ 
  coord_cartesian(ylim = limits)+
  facet_grid(.~test,labeller=my_labeller)+
  ylab("estimation error")
ggsave(filename = "supp_const_boxplot1000.pdf")



######## patton instruments ###########################################


file <- "./results/supp_constspline_T_p.Rda"

cur.instruments <- "patton.short"

# simulate_single(
# 
#   T.int = c(500,1000,2000,3000,4000,5000),
#   tests = list(
# 
#            list(iden.fct = "const.spline", model="const.spline", instruments=cur.instruments),
#            list(iden.fct = "quantile_model", model="logit_const", instruments=cur.instruments),
#            list(iden.fct = "expectile_model", model="logit_const", instruments=cur.instruments)
# 
#          ), N=10, save_to = file)


res <- readRDS(file)

res <- res[res$T>=1000,]

limits <- c(-1,1)

#copy results to same variable as spline
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$Theta.1.[grepl("expectile",res$test)]


#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#move from theta1 to theta2
res$thh.Theta.2.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("quantile",res$test)]<-NA

#rename
res$test <- as.character(res$test)
res$test[grepl("expectile",res$test)]<-"expectile"
res$test[grepl("quantile",res$test)]<-"quantile"
res$test[grepl("spline",res$test)]<-"spline"
res$test <- as.factor(res$test)

limits <- c(-1,1)

#plot
plot.single(save=res,theta0 = c(0,0,0),y.var = c("thh.Theta.1.","thh.Theta.2.","thh.Theta.3."),
            limits = limits,ncol=3,show.legend = "bottom",
            delete.y.lab=c(2,3))
#ggsave("supp_splines_const_T_p.pdf",width = 5,height = 3)




## Histogram for T = 1000

file <- "./results/supp_constspline1000_p.Rda"

cur.instruments <- "patton.short"
# simulate_single(delete.old.results = TRUE,
# 
#   T.int = c(1000),
#   tests = list(
# 
#                list(iden.fct = "const.spline", model="const.spline", instruments=cur.instruments),
#                list(iden.fct = "quantile_model", model="logit_const", instruments=cur.instruments),
#                list(iden.fct = "expectile_model", model="logit_const", instruments=cur.instruments)
# 
#   ), N=1000, save_to = file)


#### boxplots

res <-readRDS(file)


#copy results to same variable as spline
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$Theta.1.[grepl("expectile",res$test)]

#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#safe in long format
ress <- melt(res,measure.vars = grep("thh.Theta",names(res)))

#drop NAs
ress <- ress[!is.na(ress$value),]

limits <- c(-5,5)
ggplot(ress,aes(x=variable,y=value,fill=test))+geom_boxplot(outlier.shape = 1)+
  theme(
    axis.title.x = element_blank(),legend.position = "none",
    axis.text.x = element_blank())+ 
  coord_cartesian(ylim = limits)+
  facet_grid(.~test,labeller=my_labeller_const)+
  ylab("estimation error")
#ggsave(filename = "supp_const_boxplot1000_p.pdf")



############################################################################
#######################      Flexible spline     ###########################
############################################################################

file <- "./results/supp_flexspline_p.Rda"

cur.instruments <- "patton.correct2"

# simulate_single(
# 
#   T.int = c(1000,2000,3000,4000,5000),
#   tests = list(
# 
#            list(iden.fct = "spline.flex.median", model="spliney", instruments=cur.instruments),
#            list(iden.fct = "spline.flex.median", model="splinex", instruments=cur.instruments),
#            list(iden.fct = "quantile_model", model="logisticx", instruments=cur.instruments),
#            list(iden.fct = "expectile_model", model="logisticx", instruments=cur.instruments)
# 
#          ), N=10, save_to = file)


res <- readRDS(file)

limits <- c(-2,2)

plot.single(save=subset(res,test=="spline.flex.medianspliney"),theta0 = c(0,0,0,0,0,0),
            limits = limits,nrow=2,ncol=3,show.legend = FALSE)
#ggsave("supp_spliney_T.pdf")
# 
# plot.single(save=subset(res,test=="spline.flex.mediansplinex"),theta0 = c(0,0,0,0,0,0),
#             limits = limits,nrow=2,ncol=3,show.legend = FALSE)
# 
# plot.single(save=subset(res,test=="quantile_modellogisticx"),theta0 = c(-logit(pnorm(0.25)),0,NA,NA,NA,NA),
#             limits=limits,show.legend = FALSE)
# 
# plot.single(save=subset(res,test=="expectile_modellogisticx"),theta0 = c(-logit(penorm(0.25)),0,NA,NA,NA,NA),
#             limits = limits,show.legend = FALSE)
# 




#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#move from theta12 to theta45
res$thh.Theta.4.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]
res$thh.Theta.5.[grepl("quantile",res$test)]<-res$thh.Theta.2.[grepl("quantile",res$test)]
res$thh.Theta.1.[grepl("quantile",res$test)]<-NA
res$thh.Theta.2.[grepl("quantile",res$test)]<-NA

#rename
res$test <- as.character(res$test)
res$test[grepl("expectile",res$test)]<-"expectile"
res$test[grepl("quantile",res$test)]<-"quantile"
res$test[grepl("splinex",res$test)]<-"spline"
res$test <- as.factor(res$test)

limits <- c(-1,1)
#plot all but inadmissible spline
plot.single(save=subset(res,test!="spline.flex.medianspliney"),theta0 = c(0,0,0,0,0,0),
            limits = limits,nrow=2,ncol=3,show.legend = "bottom")
#ggsave("supp_splines_T_p.pdf")


## Histogram for T = 1000

file <- "./results/supp_flexspline1000_p.Rda"

cur.instruments <- "patton.correct2"

# simulate_single(
# 
#   T.int = c(1000),
#   tests = list(
# 
#     list(iden.fct = "spline.flex.median", model="spliney", instruments=cur.instruments),
#     list(iden.fct = "spline.flex.median", model="splinex", instruments=cur.instruments),
#     list(iden.fct = "quantile_model", model="logisticx", instruments=cur.instruments),
#     list(iden.fct = "expectile_model", model="logisticx", instruments=cur.instruments)
# 
#   ), N=1000, save_to = file)


#### boxplots

res <-readRDS(file)



#adjust to estimation error
res$thh.Theta.1.[grepl("quantile",res$test)]<-res$thh.Theta.1.[grepl("quantile",res$test)]+logit(pnorm(0.25))
res$thh.Theta.1.[grepl("expectile",res$test)]<-res$thh.Theta.1.[grepl("expectile",res$test)]+logit(penorm(0.25))

#safe in long format
ress <- melt(res,measure.vars = grep("thh.Theta",names(res)))

#drop NAs
ress <- ress[!is.na(ress$value),]

#without facets
# ggplot(ress,aes(x=interaction(variable,test),y=value,fill=test))+geom_boxplot(outlier.shape = NA)+
#   theme(
#   axis.title.x = element_blank(),
#   axis.text.x = element_blank())+
#   coord_cartesian(ylim = limits)#+
#   #scale_fill_grey()

limits <- c(-5,5)
ggplot(ress,aes(x=variable,y=value,fill=test))+geom_boxplot(outlier.shape = 1)+
  theme(
    axis.title.x = element_blank(),legend.position = "none",
    axis.text.x = element_blank())+ 
  coord_cartesian(ylim = limits)+
  facet_grid(.~test,labeller=my_labeller)+
  ylab("estimation error")
#ggsave(filename = "supp_boxplot1000_p.pdf")

