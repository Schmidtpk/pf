grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}
#estiamtes the true covariance for linear model with state y_t-1 and patton gdp
estimate.true.cov <- function(nn=10000, drop.until = 1, parameter.vector, instr=c("lag(Y)","X"),forecast.type="Y") {
  
  #simulate time series
  dat_sim <- generate_data_y(DGP = "patton", T = nn)
  if(forecast.type=="Y")
    dat_sim <- generate_forecast(dat_sim, list(type = "Y", model = "linearprobit", parameter = parameter.vector))
  else if(forecast.type=="X")
    dat_sim <- generate_forecast(dat_sim, list(type = "X", model = "linear", parameter = parameter.vector))
  else
    stop("unknown forecaster")
  
  ### compute covariance matrix
  #define m
  m <- function(s) probit_linear(s,parameter.vector)
  
  
  #define derivative of m
  md_all <- function(s,theta) rbind(dnorm(s*theta[2]+theta[1]),dnorm(s*theta[2]+theta[1])*s)
  md <- function(s) md_all(s,parameter.vector)
  
  #define instruments
  w <- rep(1,nn)
  for(inst_cur in instr)
  {
    #add instrument
    w <- cbind(w,eval(parse(text=inst_cur),dat_sim))
  }
  
  q <- ncol(w)
  
  if(forecast.type=="Y")
    mm <- md(lag(dat_sim$Y))
  else
    mm <- md(dat_sim$X)
  
  g<- array(dim = c(nn/drop.until,2,q))
  for(i in seq(1,nn/drop.until))
  {
    g[i, ,]<-matrix(mm[,i*drop.until],ncol=1) %*% matrix(w[i*drop.until,] ,nrow=1)
  }
  G <- apply(g, c(2,3), mean, na.rm=TRUE)
  
  
  #compute first part (vectorized)
  if(forecast.type=="Y")
    diff.cur <- ((dat_sim$Y<dat_sim$X)-m(lag(dat_sim$Y)))^2
  else
    diff.cur <- ((dat_sim$Y<dat_sim$X)-m(dat_sim$X))^2  
  
  sigma<- array(dim = c(nn/drop.until,q,q))
  for(i in seq(1,nn/drop.until))
  {
    sigma[i, ,]<-diff.cur[i*drop.until]*w[i*drop.until,]%*%t(w[i*drop.until,])
  }
  Sigma <- apply(sigma, c(2,3), mean, na.rm=TRUE)
  
  cov <- solve((G) %*% solve(Sigma) %*% t(G))
  cov
}




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




plot.hist <- function(file, var="p_value", xlim=c(-2,2))
{
  if(is.character(file))
    dat.cur <- readRDS(file)
  else
    dat.cur <- file

  #choose highest sample size
  max.T <- max(dat.cur$T)
  dat.cur <- subset(dat.cur, T==max.T)
  

  
  #plot histrogram
  ggplot(data = dat.cur, aes_string(x=var,fill="test"))+
    geom_histogram(position = "dodge",alpha=.5,bins = 20)+
    xlim(xlim)+
    ylab(paste("T =",max.T))
}  
  
  



#################### Plot graphs
plot.on.T <- function(file_name, alpha.level = 0.1,only.complete.cases=F, average=F,limits=c(0,1),show.instruments =F, nice.names=F) {
  
  save <- readRDS(file_name)
  
  if(!is.null(save$bw))
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
  
  colnames(save)[5]<-"Test"
  
  

  if(show.instruments)
    {
      save$Test <- paste(save$Test,save$instruments,sep = "|")
    }
  
  #if some lags, show power with dashed line
  if(sum(grepl("lag",save$Test))>0)
  {
    save$power <- "dashed"
    save$power[grepl("lag",save$Test)] <- "solid"
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
  
  res <- res[complete.cases(res),]
  print(res)
  

  if(nice.names)
   {
    res$Test<-as.character(res$Test)

    res$Test[grepl(pattern ="expectile" ,x =  res$Test)]<-"expectile"
    res$Test[grepl(pattern ="quantile" ,x =  res$Test)]<-"quantile"
    res$Test[grepl(pattern ="spline" ,x =  res$Test)]<-"spline"
    
    res$Test<-as.factor(res$Test)
  }
  

  # generate plot
  qpp <- ggplot(data=res, 
                aes(x=T, y=rejection.rate, colour=Test, shape=Test))
  
  
  qpp <- qpp +geom_point(size=3, aes(shape=Test, colour=Test)) +
    geom_line(size=1.2, alpha=.5, aes(colour=Test,linetype=power))

  qpp <- qpp + scale_x_continuous(breaks = pos_T, labels=array_T)+geom_abline(intercept=alpha.level, slope=0)+
    ylim(limits)+ylab('rejection rate')+xlab('sample size')+scale_linetype_discrete(guide=FALSE)+theme_classic(10)+
    theme(legend.position="bottom")
  
  return(qpp)
}


rejection_table<-function(results, alpha_level = .1, T=NULL, target = "p_value")
{
  if(is.null(T))
  {
    T_int <- unique(results$T)
  } else {T_int <- T}
  
  
  
  for(T_cur in T_int)
  {
    min.observation.number <-99999

    
    resT <- subset(results, T==T_cur)
    message("Drop ", sum(!complete.cases(resT[,target])))
    resT <- resT[complete.cases(resT[,target]),]
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
        
        min.observation.number <- min(length(resforecast[,target]),min.observation.number)
        
        output[forecast_number,test_number]<- round(mean(resforecast[,target] < alpha_level), digits=3)
        
      }
    }
  
    cat("\n \n sample size = ", T_cur , "with minimal number of observations ", min.observation.number, "\n")

  }
  
  
  # don't show instruments if unique
  if(length(unique(resT$instruments))==1)
    colnames(output)<-tests$test
  
  
  # don't show test type if unique
  if(length(unique(resT$test))==1)
    colnames(output)<-tests$instruments
  
  return(output)
}

