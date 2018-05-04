# The simulate function is a wrapper for the simulation
# iden.fct defines the identification function
# instr the instruments used in the tests
# f_model_type/parameter defines forecaster type and paramter
# t_model type defines test type
# N_cur the number of simulation runs 
# save_to defines path where results are safed
# T_int is a vector of sample sizes which the simulation applies
simulate <- function(forecasters= list(
                        list(forecaster="model",model_type="linear",model_parameter=c(logit(pnorm(0.25)),0))
                        ),
                    tests = list(
                       list(iden.fct = "spline.flex", model="spline0", instruments="ff66"),
                       list(iden.fct = "spline.flex", model="spliney", instruments="ff66"),
                       list(iden.fct = "splinex", model="spline", instruments="ff66"),
                       list(iden.fct = "spline.flex", model="splinex", instruments="ff66"),
                       list(iden.fct = "quantile_model", model="linear", instruments="forecast"),
                       list(iden.fct = "expectile_model", model="linear", instruments="forecast")
                  
                     ),
                    T_int = c(100),
                    dgp="patton",wmatrix="optimal",
                    N=1,
                    save_to ="./results/test.Rda",
                    delete.old.results=F,
                    bandwidth="bwAndrews",
                    gmm.type = "twoStep")
  
{
  #delete old file
  if(delete.old.results)
    file.remove(save_to)
  
  old.save <- NULL

  
  for(i in 1:N)
  {

    for (T_cur in T_int)
    {
      int <- 5:(T_cur)
      int1 <- 4:(T_cur-1)
      int2 <- 3:(T_cur-2)
      int3 <- 2:(T_cur-3)
      int4 <- 1:(T_cur-4)
      
      data <- generate_data_y(T=T_cur, DGP=dgp)
      data.key <- sample(1:1000000,1)
      
      for(forecaster in forecasters)
      {

     
        
        
        
        data <- generate_data_x(Y=data$Y, var=data$var, forecasting_type=forecaster$forecaster,
                                forecasting_model_type = forecaster$model_type,
                                forecasting_parameter = forecaster$model_parameter)
        
       
      
        
        for(test in tests)
        {
     
          # connect real identification function to string
          iden.fct.real <- get(test$iden.fct)
          
          


          #define model variable
          model_variable_t <- find.model.var(test$model, data,int,int1)
          
          
          if(test$iden.fct=="sp.flex")
            node.cur <- c(-1,0,1)
          else
            node.cur <- NULL
          
          ####################### Estimate functional and save results
   

          temp_save <- try(estimate.functional(iden.fct = iden.fct.real,
                                               state_variable = model_variable_t, 
                                               model_type = test$model,
                                               instruments = test$instruments,
                                               Y = data$Y[int], X = data$X[int],
                                               bandwidth = bandwidth,
                                               wmatrix = wmatrix,
                                               node = node.cur,
                                               gmm.type=gmm.type), 
                           silent=T)
          
          
          
          if (class(temp_save)=="try-error"){
            
            #safe as NA
            new.save <- data.frame(T=T_cur, p_value=NA,
                                          forecaster=paste0(forecaster$forecaster,forecaster$model_type), 
                                          test=paste0(test$iden.fct,test$model), 
                                          instruments=test$instruments, 
                                          theta=t(forecaster$model_parameter),
                                          data.key=data.key,
                                          bw=bandwidth,
                                          wmatrix=wmatrix) 
            
            message(temp_save)
            warning(paste0("NA in result for ",test$iden.fct,"\n"))
            
          } else {
            
            p_temp <- try(specTest(temp_save$gmm)$test[2],silent=T)
            
            if(class(p_temp)=="try-error")
              p_temp <- NA
            
            
            
            
            new.save <- data.frame(T=T_cur, 
                                          p_value=p_temp, 
                                          forecaster=paste0(forecaster$forecaster,forecaster$model_type), 
                                          test=paste0(test$iden.fct,test$model), 
                                          instruments=test$instruments, 
                                          theta=t(forecaster$model_parameter),
                                          thh=t(temp_save$gmm$coef),                                          data.key=data.key,
                                          bw=bandwidth,
                                          wmatrix=wmatrix)
            
      
            
            names(new.save)[names(new.save)=="Theta.1."]<- "thh.Theta.1."

            
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
      cat((i*100/N), "% finished \n")
    }
  }
  
  readRDS(file=save_to)
}  




  


