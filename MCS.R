# The simulate function is a wrapper for the simulation
# iden.fct defines the identification function
# instr the instruments used in the tests
# f_model_type/parameter defines forecaster type and paramter
# t_model type defines test type
# N_cur the number of simulation runs 
# save_to defines path where results are safed
# T_int is a vector of sample sizes which the simulation applies
simulate <- function(forecasters = list(
                        list(type="X2",model="logistic",parameter=c(logit(pnorm(0.25)),0)),
                        list(type="Y2",model="probit",parameter=c(logit(pnorm(0.25)),0))
                        ),
                    tests = list(
                       list(description = "quantile_logistic", iden.fct = quantiles, model=logistic,state="lag(lag(Y))", 
                            instruments=c("X","lag(X-Y)","lag(X-Y)^2",
                                          "lag(X)","lag(lag(X-Y))","lag(lag(X-Y))^2")),
                       list(iden.fct = expectiles, model=logistic, state="lag(X)",
                            instruments=c("X","lag(X-Y)","lag(X-Y)^2", 
                                          "lag(X)","lag(lag(X-Y))","lag(lag(X-Y))^2"))
                                ),
                    T_int = c(200),
                    dgp="patton",
                    N=2,
                    save_to ="./results/test.Rda",
                    delete.old.results=F,...
                    )
  
{
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
          else if(identical(test$iden.fct,spline.flex.median))
            theta0 <- rep(0,6)
          else
            theta0 <- NULL
          
          
          if(exists("reinbrowsern"))
            if(reinbrowsern)
              browser()
          
          temp_save <- try(suppressMessages(estimate.functional(...,
                                               iden.fct = test$iden.fct,
                                               model = test$model,
                                               theta0 = theta0,
                                               stateVariable = statevar, 
                                               instruments = test$instruments,
                                               Y = data$Y, X = data$X))
                           ,silent=T)
          
         
          if(is.null(test$description))
            test.char <- "NA"
          else
            test.char <- test$description
          
          
          if(is.null(test$state))
            test$state <- NA
          
          
          if (class(temp_save)=="try-error"){
            
            #safe as NA
            new.save <- data.frame(T=T_cur, p_value=NA,
                                          state_test = test$state,
                                          forecaster=paste0(forecaster$type,forecaster$model), 
                                          test=test.char, 
                                          instruments=paste0(test$instruments,collapse = '|'), 
                                          theta=t(forecaster$parameter),
                                          data.key=data.key) 
            
            message(temp_save)
            warning(paste0("NA in result for test ",test$description,"\n"))
            
          } else {
            
            conf <- confint(temp_save$gmm,level=.8)$test
            hit <- forecaster$parameter > conf[,1] & forecaster$parameter < conf[,2]
            
            p_temp <- try(specTest(temp_save$gmm)$test[2],silent=T)

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
          
     
          
          old.save <- rbind.fill(old.save,new.save)
       
          
            

          
          
          
        }
      }   
    }

    
    
    
    ######### Safe results
    
    
    #append if already existing results
    if(file.exists(save_to))
    {
      old.save <- rbind.fill(readRDS(file=save_to), old.save)
    }
    #safe
    saveRDS(old.save, file=save_to) 
    
    old.save <- NULL
    
    if((i*100/N)%%10==0)
    {
      cat('\r',(i*100/N), "% finished \n")
      flush.console()
    }
  }
  
  readRDS(file=save_to)
}  




  


