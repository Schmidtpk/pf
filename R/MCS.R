# The simulate function is a wrapper for the simulation
# iden.fct defines the identification function
# instr the instruments used in the tests
# f_model_type/parameter defines forecaster type and paramter
# t_model type defines test type
# N_cur the number of simulation runs
# save_to defines path where results are safed
# T_int is a vector of sample sizes which the simulation applies

#' @export
simulate <- function(
  forecasters = list(
    list(type="X2",model="logistic",parameter=c(logit(pnorm(0.25)),0)),
    list(type="Y2",model="probit",parameter=c(logit(pnorm(0.25)),0))),
  tests = list(
    list(iden.fct = expectiles, model=logistic, state="lag(X)",
         instruments=c("X","lag(X-Y)","lag(X-Y)^2",
                       "lag(X)","lag(lag(X-Y))","lag(lag(X-Y))^2"))),
  T_int = c(200),
  dgp="patton",
  N=2,
  save_to ="./results/test.Rda",
  delete.old.results=F,
  cov.est=NW_HAC,...){

  #delete old file
  if(delete.old.results)
    file.remove(save_to)

  old.save <- NULL


  for(i in 1:N)
  {

    for (T_cur in T_int)
    {

      datay <- generate_data_y(T=T_cur, DGP=dgp)
      data.key <- sample(1:1000000,1)

      for(forecaster in forecasters)
      {



        data <- generate_forecast(datay,forecaster)
        data <- data[complete.cases(data),]

        for(test in tests)
        {

          #define statevar
          if(!is.null(test$state))
          {
            if (is.na(test$state))
              statevar <- NULL
            else
              statevar <- eval(parse(text=test$state),data)
          } else
          {
            statevar <- data$X
          }


          if(identical(test$model,breakprobit))
            theta0 <- c(0,0)
          else if(identical(test$model,periodprobit))
            theta0 <- c(0,0)
          else if(identical(test$iden.fct,spline.flex.median) | grepl(pattern = "spline",test$description))
            theta0 <- rep(0,6)
          else if(identical(test$model, fixed))
             theta0 <- vector()
          else if(identical(test$model, probit_constant))
            theta0 <- 0
          else if(identical(test$model, probit_constant2))
            theta0 <- .5
          else
            theta0 <- NULL


          if(is.null(test$description))
            test.char <- "NA"
          else
            test.char <- test$description


          if(is.null(test$state))
            test$state <- NA

          if(is.null(test$cov.est))
            test$cov.est <- cov.est

          if(exists("reinbrowsern"))
            if(reinbrowsern)
              browser()


          # test fifxed moment -------------------------------------------------------
          if(!is.null(test$moment)|
             identical(test$iden.fct,fixed)|
             identical(test$iden.fct,fixed_exp)|
             identical(test$iden.fct,fixed_mean)|
             identical(test$iden.fct,fixed_median)|
             identical(test$iden.fct,MZ_test)|
             identical(test$iden.fct,MZ)|
             identical(test$iden.fct,MZ_PointFore))
          {
            if(
              !is.null(test$moment)|
              identical(test$iden.fct,fixed)|
              identical(test$iden.fct,fixed_exp)|
              identical(test$iden.fct,fixed_mean)|
              identical(test$iden.fct,fixed_median))
              {

                temp_save <- try(suppressMessages(
                  test_moment(X=data$X,Y=data$Y,
                                  instruments = test$instruments,
                                  iden.fct=test$iden.fct,
                                  other_data = data,
                                  cov.est = test$cov.est
                                  ))
                  ,silent=T)
              }else if(identical(test$iden.fct,MZ_test))
              {
                temp_save <- try(suppressMessages(
                  MZ_test(X=data$X,
                          Y=data$Y,
                          instruments = test$instruments,
                          cov.est = test$cov.est
                  )),silent=T)
              }else if(identical(test$iden.fct,MZ_PointFore))
              {
                temp_save <- try(suppressMessages(
                  test_MZ_PointFore(X=data$X,
                          Y=data$Y,
                          test=test,
                          stateVariable = statevar,
                          ,...
                  )),silent=T)
              }


            if (class(temp_save)=="try-error"|any(is.na(temp_save))){

              #safe as NA
              new.save <- data.frame(T=T_cur, p_value=NA,
                                     state_test = test$state,
                                     forecaster=paste0(forecaster$type,forecaster$model),
                                     test=test.char,
                                     instruments=paste(test$instruments,collapse = '|'),
                                     theta=t(forecaster$parameter))
            } else{

            if(length(temp_save)>1)
              stop("wrong output in simulation for pvalue")

            new.save <- data.frame(T=T_cur,
                                   p_value=temp_save,
                                   forecaster=paste0(forecaster$type,forecaster$model),
                                   test=test.char,
                                   instruments=paste(test$instruments,collapse = '|'),
                                   theta=t(forecaster$parameter))
            }

          } else # test flexible moment ----------------------------------------------------
          {

            # use HAC in gmm if not specified otherwise
            if(is.null(test$cov.est))
              test$cov.est <- "HAC"

            if(is.null(test$prewhite))
              test$prewhite <- FALSE

            if(is.null(test$centeredVcov))
              test$centeredVcov <- TRUE

            temp_save <- try(suppressMessages(estimate.functional(...,
                                                                  vcov = test$cov.est,
                                                                  iden.fct = test$iden.fct,
                                                                  model = test$model,
                                                                  theta0 = theta0,
                                                                  prewhite = test$prewhite,
                                                                  centeredVcov = test$centeredVcov,
                                                                  stateVariable = statevar,
                                                                  instruments = test$instruments,
                                                                  Y = data$Y, X = data$X,
                                                                  other_data = data))
                             ,silent=T)


            if (class(temp_save)=="try-error"){

              #safe as NA
              new.save <- data.frame(T=T_cur, p_value=NA,
                                     state_test = test$state,
                                     forecaster=paste0(forecaster$type,forecaster$model),
                                     test=test.char,
                                     instruments=paste0(test$instruments,collapse = '|'),
                                     theta=t(forecaster$parameter),
                                     data.key=data.key)

              #message(temp_save)
              #warning(paste0("NA in result for test ",test$description,"\n"))

            } else {

              conf <- confint(temp_save$gmm,level=.8)$test
              hit <- forecaster$parameter > conf[,1] & forecaster$parameter < conf[,2]

              p_temp <- try(gmm::specTest(temp_save$gmm)$test[2],silent=T)

              if(class(p_temp)=="try-error")
                p_temp <- NA

              new.save <- data.frame(T=T_cur,
                                     p_value=p_temp,
                                     state_test = test$state,
                                     forecaster=paste0(forecaster$type,forecaster$model),
                                     test=test.char,
                                     instruments=paste(test$instruments,collapse = '|'),
                                     theta=t(forecaster$parameter),
                                     thh=t(temp_save$gmm$coef),
                                     sdh=matrix(diag(temp_save$gmm$vcov),nrow=1),
                                     hit=t(hit)
              )
            }
          }



          old.save <- plyr::rbind.fill(old.save,new.save)

        }
      }
    }




    ######### Safe results

    if((i*100/N)%%10==0)
    {
      cat('\r',(i*100/N), "% finished  ")
      flush.console()


      #append if already existing results
      if(file.exists(save_to))
      {
        old.save <- plyr::rbind.fill(readRDS(file=save_to), old.save)
      }
      #safe
      saveRDS(old.save, file=save_to)

      old.save <- NULL
    }
  }

  readRDS(file=save_to)
}




# test_MZ_PointFore <- function(X,Y,stateVariable,test, other_data=NULL)
# {
#   theta0 <- c(1)
#   res <- estimate.functional(iden.fct = test$iden.fct,
#                       model = test$model,
#                       theta0 = theta0,
#                       stateVariable = stateVariable,
#                       instruments = test$instruments,
#                       Y = Y, X = X, wmatrix="ident",
#                       other_data = other_data)
#
#   cur.lm <- res$gmm
#   browser()
#   res <- car::linearHypothesis(cur.lm,c("Theta[1]=1"))#,"Theta[2]=1"))
#
#   pval <- res$`Pr(>F)`[2]
# }


# test_MZ <- function(X,Y,w,instruments,iden.fct=fixed, other_data=NULL)
# {
#   if(is.null(w))
#   {
#     cur.lm <- lm(Y~X,data= data.frame(Y=Y,X=X))
#     res <- car::linearHypothesis(cur.lm,c("(Intercept) = 0", "X = 1"),vcov=sandwich::vcovHAC,bw=PointFore:::bwNeweyWest1987)
#     pval <- res$`Pr(>F)`[2]
#   } else if(is.vector(w)){
#     cur.lm <- lm(Y~X+w,data= data.frame(Y=Y,X=X))
#     if(sum(is.na(cur.lm$coefficients))>0)
#       stop("NA in lm estimation indicating colinear relationship between X and instruments")
#     res <- car::linearHypothesis(cur.lm,c("(Intercept) = 0", "X = 1","w = 0"),vcov=sandwich::vcovHAC,bw=PointFore:::bwNeweyWest1987)
#     pval <- res$`Pr(>F)`[2]
#   } else {
#     stop("MZ can't handle multiple instruments")
#   }
# }


# test_MZ <- function(X,Y,w,instruments,iden.fct=fixed, other_data=NULL)
# {
#   if(is.null(w))
#   {
#     cur.lm <- rms::Rq(Y~X,data= data.frame(Y=Y,X=X),tau = 0.4012937)
#     #res <- car::linearHypothesis(cur.lm,c("Intercept = 0", "X = 1"))#,vcov=sandwich::vcovHAC,bw=PointFore:::bwNeweyWest1987)
#     res <- car::linearHypothesis(cur.lm,c("Intercept = 0", "X = 1"),vcov.=quantreg::summary.rq(cur.lm,covariance = TRUE,se="boot",method="wxy")$cov)
#     pval <- res$`Pr(>Chisq)`[2]
#   } else if(is.vector(w)){
#     cur.lm <- rms::Rq(Y~X+w,data= data.frame(Y=Y,X=X),tau = 0.4012937)
#     if(sum(is.na(cur.lm$coefficients))>0)
#       stop("NA in lm estimation indicating colinear relationship between X and instruments")
#     #res <- car::linearHypothesis(cur.lm,c("Intercept = 0", "X = 1","w = 0"))#,vcov=sandwich::vcovHAC,bw=PointFore:::bwNeweyWest1987)
#     res <- car::linearHypothesis(cur.lm,c("Intercept = 0", "X = 1","w = 0"),vcov.=quantreg::summary.rq(cur.lm,covariance = TRUE,se="boot",method="wxy")$cov)
#     pval <- res$`Pr(>Chisq)`[2]
#   } else {
#     stop("MZ can't handle multiple instruments")
#   }
# }


test_moment <- function(X,Y,instruments,iden.fct=fixed,
                        cov.est = NW_HAC,
                        other_data=NULL)
{
  # instrument creation from estimate.functional ------------------------------------------------
  w <- rep(1, length(Y))
  if (!(length(instruments) == 1 & grepl("const",
                                         instruments[1]))) {
    if (is.null(other_data))
      other_data <- data.frame(X = X, Y = Y)
    else other_data <- cbind(data.frame(X = X, Y = Y),
                             other_data)
    for (inst_cur in instruments) {
      if (grepl("Y", inst_cur) & !grepl("lag",
                                        inst_cur))
        warning("Y without lags is not a valid instrument as it is not in the information set of the forecaster.")
      w <- cbind(w, eval(parse(text = inst_cur), other_data))
    }
    compl <- complete.cases(w)
    message(paste("Drop ", length(Y) - sum(compl),
                  "case(s) because of chosen instruments"))
    w <- w[compl, ]
    Y <- Y[compl]
    X <- X[compl]
  }

  instruments <- w

  # test from meanmedianmode ------------------------------------------------

  n <- length(Y)

  res <- iden.fct(X,Y)*instruments


  if(is.matrix(res))
  {
    #compute average of identification fct
    mean_id <- colMeans(res)

    #compute dof
    J_DoF <- dim(instruments)[2]
  } else
  {
    mean_id <- mean(res)
    J_DoF <- 1
  }

  # safe covariance info
  var_info <- tryCatch(cov.est(res), error=function(x) NA)

  # safe covariance estimation result
  var_id <- if(!any(is.na(var_info))) var_info$cov else NA

  #compute pvalue
  J <- tryCatch(n * mean_id %*% solve(var_id) %*% mean_id, error=function(e) NA)
  pval <- tryCatch(1-pchisq(J,J_DoF), error=function(e) NA)

  return(pval)
}


bwNeweyWest1987 <- function(x,...) {
  sandwich::bwNeweyWest(x,lag=nrow(gmm::estfun.gmmFct(x))^(0.2),...)
}

#' @export
NW_HAC <- function(ids)
{
  if(is.matrix(ids))
  {
    n <- nrow(ids)
  } else
  {
    n <- length(ids)
  }
  list(cov= n*sandwich::kernHAC(lm(ids~1),bw=bwNeweyWest1987,kernel="Bartlett",prewhite = FALSE), name = "NW_HAC")
}


#' @export
HAC <- function(ids)
{
  if(is.matrix(ids))
  {
    n <- nrow(ids)
  } else
  {
    n <- length(ids)
  }
  list(cov= n*sandwich::kernHAC(lm(ids~1),prewhite = FALSE), name = "HAC")
}

#' HAC estimator with prewhite for 1 lag
#'
#' @param ids values
#' @param ...
#'
#' @return list of estimation result cov and name (for saving purpose)
#' @export
HAC_prewhite1 <- function(ids,...)
{
  if(is.matrix(ids))
  {
    n <- nrow(ids)
  } else
  {
    n <- length(ids)
  }
  res <- n*sandwich::vcovHAC(lm(ids~1),prewhite=1,...)
  list(cov = res, name = "HAC_prewhite1")
}

#' iid variance estimator
#'
#' @param ids values
#' @param ...
#'
#' @return list of estimation result cov and name (for saving purpose)
#' @export
iid <- function(ids,...)
{
  list(cov= var(ids), name = "iid")
}


#' HC estimator
#'
#' @param ids values
#' @param ...
#'
#' @return list of estimation result cov and name (for saving purpose)
#' @export
HC <- function(ids)
{
  if(is.matrix(ids))
  {
    n <- nrow(ids)
  } else
  {
    n <- length(ids)
  }
  res <- n*sandwich::vcovHC(lm(ids~1))
  list(cov = res, name = "HC")
}

#' @export
HC3 <- function(ids)
{
  if(is.matrix(ids))
  {
    n <- nrow(ids)
  } else
  {
    n <- length(ids)
  }
  res <- n*sandwich::vcovHC(lm(ids~1),type="HC3")
  list(cov = res, name = "HC3")
}
