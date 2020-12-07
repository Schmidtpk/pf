#' @export
generate_data_y <- function(DGP="patton", T=1000, extra = 100, forecasting_model_type = NULL, forecasting_parameter = NULL,...)
{
  Y <- numeric(T+extra)

  var <- numeric(T+extra)
  mean <- numeric(T+extra)


  if(DGP=="gauss")
  {

    Y <- gen.gauss(T+extra,...)
    var <- NA
    mean <- NA
  } else if(DGP=="EKT")
  {
    W1 <- rnorm(T+extra,1,1)
    W2 <- rnorm(T+extra,-1,1)
    #W1 <- arima.sim(list(ar=.5),n=T+extra)
    #W2 <- arima.sim(list(ar=.5),n=T+extra)

    Y <- 1+0.5*W1+0.5*W2+rnorm(T+extra,0,.5)

    return(data.frame(Y=Y[(extra+1):(T+extra)],
                      W1=W1[(extra+1):(T+extra)],
                      W2=W2[(extra+1):(T+extra)],
                      var=rep(.5,T),
                      mean=1+0.5*W1[(extra+1):(T+extra)]+0.5*W2[(extra+1):(T+extra)]))
  } else
  {


    a <- .5
    b <- 0.1
    c <- 0.8
    d <- 0.1

    Y[1]<-0

    var[1]<-1
    mean[1]<-0
    epsilon <- rnorm(T+extra,0,1)

    for (t in 2:(T+extra))
    {
      if (DGP=="patton")
      {
        var[t] <- b+c*var[t-1]+d*var[t-1]*epsilon[t-1]^2
        mean[t] <- a * Y[t-1]
        Y[t] <- a * Y[t-1]+sqrt(var[t])*epsilon[t]
      } else if (DGP=="AR1")
      {
        if(t==2){
          Y[t]<-0
        } else {

          var[t] <- 1
          mean[t] <- a * Y[t-1]

          Y[t] <- a * Y[t-1]+sqrt(var[t])*epsilon[t]
        }
      } else if (DGP=="AR2")
      {
        if(t==2){
          Y[t]<-0
        } else {
          var[t] <- 1
          mean[t] <-  a * Y[t-1]+ b * Y[t-2]
          Y[t] <- a * Y[t-1]+ b * Y[t-2] + sqrt(var[t])*epsilon[t]
        }
      } else {stop("DGP unknown")}
    }
  }

  #return(data.frame(Y=Y[(extra+1):(T+extra)], var=var[(extra+1):(T+extra)],mean=mean[(extra+1):(T+extra)]))
  return(data.frame(Y=Y, var=var,mean=mean))
}





#' @export
generate_forecast <- function(data, forecaster, extra = 100)
  #(Y, var=NULL,  forecasting_type="model", forecasting_model_type = "const", forecasting_parameter = .5,...)
{
  #generate forecaster
  T <- nrow(data)
  X <- numeric(T)
  Y <- data$Y
  var <- data$var
  mean <- data$mean

  if(!is.null(forecaster$model)){
    if (class(forecaster$model)=="character" & grepl("Y",forecaster$type))
      model.fct <- get(forecaster$model)}

  # if(forecasting_model_type=="reallinear")
  #   model_variable_loc<-matrix(0,nrow = T,ncol=length(forecasting_parameter))


  # if (is.null(var)) {var <- numeric(T)}

  a <- .5
  b <- 0.1
  c <- 0.8
  d <- 0.1

  X[1]<-0
  X[2]<-0

  for (t in 3:(T))
  {

    if (forecaster$type=="model")
    {
      if(t<forecaster$T_model)
      {
        X[t]<-NA
        next
      }

      #interval.cur <- (t-forecaster$T_model):(t-1)
      interval.cur <- (t-1):1
      #qm <- quantreg::rq(y~x+W1+W2, data=data.frame(y=Y[interval.cur],x=lag(Y[interval.cur]),W1 = data$W1[interval.cur],W2 = data$W2[interval.cur]),tau=pnorm(-.25))
      qm <- quantreg::rq(y~W1, data=data.frame(y=Y[interval.cur],lagy=lag(Y[interval.cur]),W1 = data$W1[interval.cur],W2 = data$W2[interval.cur]),tau=pnorm(-.25))
      X[t]<-predict(qm,data.frame(
        #lagy=Y[t-1],
        W1=data$W1[t]
        ))

    } else if (forecaster$type =="X")
    {


      mean.now <- mean[t]
      sd.now <- sqrt(var[t])

      if(forecaster$model =="linear")
      {
        X[t] <- (mean.now + sd.now*forecaster$parameter[1]) /
          (1-sd.now*forecaster$parameter[2])
      } else if (forecaster$model=="break")
      {
        stop("not ready")
        X[t] <- (mean.now + sd.now*forecaster$parameter[1])
        if(X[t]<0)
          X[t] <- X[t]-forecaster$parameter[2]*sd.now
      } else if (forecaster$model=="period")
      {
        stop("not ready")
        X[t] <- (mean.now + sd.now*forecaster$parameter[1])
        if(X[t]<0)
          X[t] <- X[t]-forecaster$parameter[2]*sd.now
      }

    } else if (forecaster$type =="Y")
    {


      model.var.now <- Y[t-1]

      level_now <- model.fct(model.var.now, forecaster$parameter)

      X[t] <- mean[t] +  qnorm(level_now) * sqrt(var[t])

    } else if (forecaster$type =="X2")
    {

      mean.now <- a^2*Y[t-2]
      sd.now <- sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1])

      X[t] <- (mean.now + sd.now*forecaster$parameter[1]) /
        (1-sd.now*forecaster$parameter[2])

    } else if (forecaster$type =="Y2")
    {

      model.var.now <- Y[t-2]
      level_now <- model.fct(model.var.now, forecaster$parameter)

      mean.now <- a^2*Y[t-2]
      sd.now <- sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1])

      X[t] <- mean.now +  qnorm(level_now) * sd.now

    } else if (forecaster$type =="Y2_mean")
    {

      mean.now <- a^2*Y[t-2]
      X[t] <- mean.now

    } else if (forecaster$type =="Y_mean")
    {

      mean.now <- a*Y[t-1]
      X[t] <- mean.now

    } else if (forecaster$type =="patton")
    {


      X[t] <- .5*Y[t-1]+sqrt(var[t])*.25


    } else if (forecaster$type =="Y_mean")
    {

      mean.now <- a*Y[t-1]
      X[t] <- mean.now

    } else if (forecaster$type =="patton2")
    {


      X[t] <- .25*Y[t-2]+sqrt(a^2 * var[t-1] + .1 + .9 * var[t-1])*.25


    } else {stop("forecasting type unknown")}
  }

  indices <- (extra+1):T

  if(forecaster$type!="model")
    return(data.frame(Y=Y[indices],X=X[indices]))
  else
    return(data.frame(Y=Y[indices],X=X[indices],W1=data$W1[indices],W2=data$W2[indices]))
}
