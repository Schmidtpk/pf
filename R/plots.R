
#' @export
plot_on_T <- function(file_name,
                      alpha.level = 0.1,
                      only.complete.cases=F,
                      average=F,
                      limits=c(0,1),
                      show.instruments =F,
                      nice.names=F,
                      keep=NULL,
                      drop=NULL,
                      drop2=NULL,
                      size.adjusted=F,
                      revision.levels=F) {

  if(is.character(file_name))
    save <- readRDS(file_name)
  else
    save <- file_name

  if(!is.null(save$bw))
    if(length(unique(save$bw))!=1)
      warning("Different bandwidth selection in results")

  ###Create plot data
  array_T <- sort(unique(save$T))
  pos_T <- 1:length(array_T)

  ##### Some formatting
  i <- sapply(save, is.factor)
  save[i] <- lapply(save[i], as.character)

  colnames(save)[colnames(save)=="test"]<-"Test"


  if(show.instruments)
  {
    save$Test <- paste(save$Test,save$instruments,sep = "|")
  }

  if(!is.null(keep))
  {
    save <- subset(save, grepl(keep,save$Test))
  }

  if(!is.null(drop))
  {
    save <- subset(save, !grepl(drop,save$Test))
  }
  if(!is.null(drop2))
  {
    save <- subset(save, !grepl(drop2,save$Test))
  }

  #if some lags, show power with dashed line
  if(sum(grepl("lag",save$Test))>0)
  {
    save$power <- "dashed"
    save$power[grepl("lag",save$Test)] <- "solid"
    save$Test <- gsub("_lag","",save$Test)
  } else
  {
    save$power <- "solid"
  }


  res <- data.frame(T=numeric(), Test=character(), rejection.rate=numeric(),completed=numeric(),total=numeric(),power=character())

  for (i in 1:length(array_T))
  {
    for (Test_cur in unique(save$Test))
    {
      cur.level <- NULL
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

        #skip if no observations
        if(nrow(sset_s)==0)
          next

        if(size.adjusted & cur.power=="solid")
        {
          res.cur <- mean(sset_s$p_value<alpha.level)
          cur.level <- quantile(sset_s$p_value,prob=alpha.level)

          if(cur.level==0){
            ratio.zero <- mean(sset_s$p_value==cur.level)
          } else {
            ratio.zero <- alpha.level
          }

        } else if(size.adjusted & cur.power=="dashed")
        {
          res.cur <- mean(sset_s$p_value<=cur.level)*alpha.level/ratio.zero
        } else
        {
          res.cur <- ifelse(average==F, mean(sset_s$p_value<alpha.level), mean(sset_s$p_value))
        }


        res<- rbind(res,
                    data.frame(T=T_curpos,
                               Test=Test_cur,
                               rejection.rate = res.cur,
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

    res$Test<-gsub("_tiny", "",res$Test)
    res$Test<-gsub("_c", "",res$Test)

    res$Test<-gsub("_inst", "",res$Test)

    #res$Test[grepl(pattern ="expectiles" ,x =  res$Test)]<-"expectiles"
    #res$Test[grepl(pattern ="quantiles" ,x =  res$Test)]<-"quantiles"
    res$Test[grepl(pattern ="spline" ,x =  res$Test)]<-"spline"





  }

  if(revision.levels)
  {
    if(revision.levels==2){
      res$Test[res$Test=="quantiles"] <- "linear (quantiles)"
      res$Test[res$Test=="expectiles"] <- "linear (expectiles)"

      res$Test[res$Test=="quantilesP"] <- "linear (quantiles)"
      res$Test[res$Test=="expectilesP"] <- "linear (expectiles)"

      res$Test <- factor(res$Test,levels=c("spline",
                                           "linear (quantiles)",
                                           "linear (expectiles)"))
    } else {
      res$Test[res$Test=="quantile"] <- "constant"
      res$Test[res$Test=="expectile"] <- "constant"

      res$Test[res$Test=="fixed_exp"] <- "fixed"
      res$Test[res$Test=="fixed"] <- "fixed"
      res$Test[res$Test=="mean"] <- "fixed"
      res$Test[res$Test=="median"] <- "fixed"

      res$Test[res$Test=="quantiles"] <- "linear"
      res$Test[res$Test=="expectiles"] <- "linear"

      res$Test <- factor(res$Test,levels=c("spline",
                                           "linear",
                                           "constant",
                                           "fixed"))
    }

  }

  # generate plot
  qpp <- ggplot(data=res,
                aes(x=T, y=rejection.rate, colour=Test, shape=Test))+geom_abline(intercept=alpha.level, slope=0)


  qpp <- qpp +geom_point(size=3, aes(shape=Test, colour=Test)) +
    geom_line(size=1.2, alpha=.5, aes(colour=Test,linetype=power))

  qpp <- qpp + scale_x_continuous(breaks = pos_T, labels=array_T)+
    ylim(limits)+ylab('rejection rate')+xlab('sample size')+scale_linetype_discrete(guide=FALSE)+theme_classic(10)+
    theme(legend.position="bottom")

  return(qpp)
}

#' Title
#'
#' @param results
#' @param alpha_level
#' @param target
#'
#' @return
#' @export
#'
#' @examples
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

# estiamtes the true covariance for linear model with state y_t-1 and patton gdp
#' @export
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
  #m <- function(s) probit_linear(s,parameter.vector)
  m <- function(s) linearprobit(s,parameter.vector)


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

