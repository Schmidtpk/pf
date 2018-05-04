#analyze dependent forecasters




#################### Plot graphs
plot.single <- function(file_name,  save=NULL, y.var.estimate = NULL, 
                        y.var = NULL, 
                        print=FALSE,
                        limits=NULL,
                        show.sqrt.T=FALSE,
                        theta0=NULL,
                        show.t.test=FALSE,
                        show.legend= "bottom",
                        delete.y.lab=NULL,...) {
  
  if(is.null(save))
    save <- readRDS(file_name)
  

  if (is.null(y.var))
    estimates.names <- names(save)[grepl(pattern = "thh.",x= names(save)) | 
                                     grepl(pattern = "Theta.",x= names(save)) &
                                     !grepl(pattern = "pval",x= names(save))]
  else
    estimates.names <- y.var

  #if true, show t-test p-values instead of estimates
  if(show.t.test)
  {
    estimates.names <- names(save)[grepl(pattern = "pval.t",x= names(save))]
    #doesn't make sense here
    show.sqrt.T<- FALSE
    y.lab.cur <- "t-test p-value"
  } else
  {
    y.lab.cur <- 'estimator'
  }
  
  
  #create plot data
  array_T <- sort(unique(save$T))
  pos_T <- 1:length(array_T)
  
  #save$T <- findInterval(save$T,array_T)
  
  ##### Some formatting
  i <- sapply(save, is.factor)
  save[i] <- lapply(save[i], as.character)
  
  colnames(save)[4]<-"Test"
  
  
  if(length(unique(save$instruments))>1)
    save$Test <- paste(save$Test,save$instruments,sep = "|")
  
  
  if(print==T)
    print(save)

  plot.list <- list()
  i <- 1
  
  for(y.var.estimate in estimates.names)
  {
    #save version that can be manipulated
    ssave <- save
    
    
    #skip if no data
    if(sum(!is.na(ssave[,y.var.estimate]))==0)
      next
    

  #if true, show the output times square route the sample size
  if(show.sqrt.T)
  {
    y.lab.cur <- "estimator x sqrt(T)"
    if(is.null(theta0))
      ssave[,y.var.estimate] <- ssave[,y.var.estimate] * sqrt(save$T)
    else
      ssave[,y.var.estimate] <- (ssave[,y.var.estimate]-theta0[which(estimates.names==y.var.estimate)]) * sqrt(ssave$T)
  } else
  {
    if(!is.null(theta0))
    {
      message("substract theta0")
      y.lab.cur <- "estimation error"
      ssave[,y.var.estimate] <- (ssave[,y.var.estimate]-theta0[which(estimates.names==y.var.estimate)])
          
    }
  }
          
  # generate plot
  qpp <- ggplot(data=ssave, 
                aes_string(x="factor(T)", y=y.var.estimate, colour="Test", shape="Test", group="interaction(Test,data.key)"))
  
  
  qpp <- qpp +geom_point(size=1) +
    geom_line(size=1.2, alpha=.2, aes(colour=Test))
  qpp <- qpp + 
   ylab(y.lab.cur)+xlab('sample size')+scale_linetype_discrete(guide=FALSE)+theme_classic(8)
  
  if(!is.null(limits))
    qpp<-qpp+coord_cartesian(ylim=limits)
  

  if(i %in% delete.y.lab)
    qpp<-qpp+theme(axis.title.y = element_blank())
  
  plot.list[[i]]<-qpp
  i <- i+1
  }
  
  ggarrange(plotlist = plot.list,common.legend =TRUE, legend = show.legend,...)
}




plot.hist <- function(file, var="pval.t")
{
  #choose highest sample size
  dat.cur <- readRDS(file)
  max.T <- max(dat.cur$T)
  dat.cur <- subset(dat.cur, T==max.T)
  

  
  #plot histrogram
  ggplot(data = dat.cur, aes_string(x=var,fill="test"))+
    geom_histogram(position = "dodge",alpha=.5,bins = 20)+
    xlim(c(-.1,1.1))+
    ylab(paste("T =",max.T))
}  
  
  



#################### Plot graphs
plot.on.T <- function(file_name, alpha.level = 0.1,only.complete.cases=F, average=F,limits=c(0,1),show.instruments =F, nice.names=F) {
  
  save <- readRDS(file_name)
  
  if(length(unique(save$bw))!=1)
    warning("Different bandwidth selection in results")
  
  ###Create plot data
  pos_T <-   c( 1,  2,  3,  4,  5,6,7)
  array_T <- c(50,100,250,500,1000,2000,4000)
  
  array_T <- sort(unique(save$T))
  pos_T <- 1:length(array_T)
  
  ##### Some formatting
  i <- sapply(save, is.factor)
  save[i] <- lapply(save[i], as.character)
  
  colnames(save)[4]<-"Test"
  
  

  if(show.instruments)
    {
      save$Test <- paste(save$Test,save$instruments,sep = "|")
    }
  
  #if some lags, show power with dashed line
  if(sum(substr(save$instruments,1,3)=="lag"))
  {
    save$power <- "dashed"
    save$power[substr(save$instruments,1,3)=="lag"] <- "solid"
  } else
  {
    save$power <- "solid"
  }

  
  res <- data.frame(T=numeric(), Test=character(), rejection.rate=numeric(),completed=numeric(),total=numeric(),power=character())
  
  for (i in 1:length(array_T))
  {
    for (Test_cur in unique(save$Test))
    {
      for(cur.power in c("solid","dashed"))
      {
        T_cur <- array_T[i]
        T_curpos <- pos_T[i]
        
        sset_s <- subset(save, T==T_cur & Test==Test_cur & power==cur.power)
        
        
        
        compl_s <- sum(!is.na(sset_s$p_value))/length(sset_s[,1])
        
        #ignore non-complete cases
        if(only.complete.cases)
          sset_s <- sset_s[!is.na(sset_s$p_value),]
        else
          sset_s$p_value[is.na(sset_s$p_value)]<-0
        
        
        res<- rbind(res, 
                    data.frame(T=T_curpos, 
                               Test=Test_cur, 
                               rejection.rate = ifelse(average==F, mean(sset_s$p_value<alpha.level), mean(sset_s$p_value)), 
                               completed=compl_s, total=length(sset_s[,1]),
                               power=cur.power)
        )
      }
    }
    
    
  }
  
  print(res[complete.cases(res),])
  

  if(nice.names)
   {
    res$Test<-as.character(res$Test)

    
    res$Test[grepl(pattern ="expectile" ,x =  res$Test)]<-"Expectile"
    res$Test[grepl(pattern ="quantile" ,x =  res$Test)]<-"Quantile"
    res$Test[grepl(pattern ="spline" ,x =  res$Test)]<-"Spline"
    
    res$Test<-as.factor(res$Test)
  }

  # generate plot
  qpp <- ggplot(data=subset(res), 
                aes(x=T, y=rejection.rate, colour=Test, shape=Test))
  
  
  qpp <- qpp +geom_point(size=3, aes(shape=Test, colour=Test)) +
    geom_line(size=1.2, alpha=.5, aes(colour=Test,linetype=power))

  qpp <- qpp + scale_x_continuous(breaks = pos_T, labels=array_T)+geom_abline(intercept=alpha.level, slope=0)+
    ylim(limits)+ylab('rejection rate')+xlab('sample size')+scale_linetype_discrete(guide=FALSE)+theme_classic(10)+
    theme(legend.position="bottom")
  
  return(qpp)
}


rejection_table<-function(results, alpha_level = .1, T=NULL)
{
  if(is.null(T))
  {
    T_int <- unique(results$T)
  } else {T_int <- T}
  
  
  
  for(T_cur in T_int)
  {
    min.observation.number <-99999

    resT <- subset(results, T==T_cur)
    resT <- resT[complete.cases(resT$p_value),]
    tests <- unique(resT[c("test","instruments")])
    n_tests <- dim(tests)[1]
    forecasts <-unique(resT[c("forecaster")])
    n_forecasts <- dim(forecasts)[1]
    
    output <- matrix(ncol=n_tests,nrow=n_forecasts)
    colnames(output)<-paste(tests$test[1:n_tests],tests$instruments)
    rownames(output)<-forecasts$forecaster[1:n_forecasts]
    
    
    for(test_number in 1:n_tests)
    {
      test_cur<-tests[test_number,]
      restest <- subset(resT, test==test_cur$test & instruments==test_cur$instruments)
      
      
      for(forecast_number in 1:n_forecasts)
      {
        forecast_cur <- forecasts[forecast_number,]
        resforecast <- subset(restest, forecaster==forecast_cur)
        
        min.observation.number <- min(length(resforecast$p_value),min.observation.number)
        
        output[forecast_number,test_number]<- round(mean(resforecast$p_value < alpha_level), digits=3)
        
      }
    }
  
    cat("\n \n sample size = ", T_cur , "with minimal number of observations ", min.observation.number, "\n")
    print(output)
    
  }
  
  return(output)
}

