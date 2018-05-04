
# The simulate function is a wrapper for the simulation
# iden.fct defines the identification function
# instr the instruments used in the tests
# f_model_type/parameter defines forecaster type and paramter
# t_model type defines test type
# N_cur the number of simulation runs 
# save_to defines path where results are safed
# T_int is a vector of sample sizes which the simulation applies
simulate_single <- function(forecasters= list(
                        list(forecaster="model",model_type="linear",model_parameter=c(logit(pnorm(0.25)),0))
                        ),
                    tests = list(
                     
                       list(iden.fct = "spline.flex", model="splinex", instruments="ff66"),
                       list(iden.fct = "quantile_model", model="linear", instruments="forecast"),
                       list(iden.fct = "expectile_model", model="linear", instruments="forecast")
                  
                     ),
                    T_max = 1000,
                    T.int = NULL,
                    dgp="patton",wmatrix="optimal",
                    N=1,
                    save_to ="./results/test.Rda",
                    delete.old.results=F,
                    bandwidth="bwAndrews",
                    gmm.type = "twoStep",...)
  
{
  
  if(delete.old.results & file.exists(save_to))
    file.remove(save_to)
  
  #define interval
  if(is.null(T.int))
    T.int <- seq.int(0,T_max,length.out = 11)[-1]
  else
    T_max <- max(T.int)
  
  old.save <- NULL
  
  for(i in 1:N)
  {
    data.key <- sample(1:1000000,1)
    data.all <- generate_data_y(T=T_max, DGP=dgp,...)
    
      for(forecaster in forecasters)
      {

        
        
       
        
        data.all <- generate_data_x(Y=data.all$Y, var=data.all$var, forecasting_type=forecaster$forecaster,
                                forecasting_model_type = forecaster$model_type,
                                forecasting_parameter = forecaster$model_parameter,...)
        
        
        for(T.cur in T.int)
        {
          
          #select data
          data <- data.all[1:T.cur,]
        
          int <- 5:(T.cur)
          int1 <- 4:(T.cur-1)
          int2 <- 3:(T.cur-2)
          int3 <- 2:(T.cur-3)
          int4 <- 1:(T.cur-4)
          
          for(test in tests)
          {
       
            # connect real identification function to string
            iden.fct.real <- get(test$iden.fct)
            
            model_variable_t <- find.model.var(test$model, data,int,int1)
                      
            
            ####################### Estimate functional and save results
         
     
            temp_save <- try(estimate.functional(iden.fct = iden.fct.real,
                                                 state_variable = model_variable_t, 
                                                 model_type = test$model,
                                                 instruments = test$instruments,
                                                 Y = data$Y[int], X = data$X[int],
                                                 bandwidth=bandwidth,
                                                 wmatrix = wmatrix,
                                                 gmm.type = gmm.type), 
                             silent=T)
            
            
            
            if (class(temp_save)=="try-error"){
              
              #safe as NA
              new.save <- data.frame(T=T.cur, p_value=NA,
                                            forecaster=paste0(forecaster$forecaster,forecaster$model_type), 
                                            test=paste0(test$iden.fct,test$model), 
                                            instruments=test$instruments, 
                                            theta=t(forecaster$model_parameter), 
                                            thh.Theta.1.=NA,
                                            data.key=data.key,
                                            bw=bandwidth,
                                            wmatrix=wmatrix) 
              message(temp_save)
              warning(paste0("NA in result for ",test$iden.fct,"\n"))
              
            } else {
              
              p_temp <- try(specTest(temp_save$gmm)$test[2],silent=T)
              
              if(class(p_temp)=="try-error")
                p_temp <- NA
              
          
            
              
              
              
              
              new.save <- data.frame(T=T.cur, 
                                            p_value=p_temp, 
                                            forecaster=paste0(forecaster$forecaster,forecaster$model_type), 
                                            test=paste0(test$iden.fct,test$model), 
                                            instruments=test$instruments, 
                                            theta=t(forecaster$model_parameter), 
                                            thh=t(temp_save$gmm$coefficients),
                                            pval.t = t(summary(temp_save$gmm)$coef[,4]),
                                            data.key=data.key,
                                            bw=bandwidth,
                                            wmatrix=wmatrix)
              
              
            }
          
     
          
            old.save <- rbind.fill(old.save,new.save)
       
          
            

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
          
          
          flush.console()
            cat('\r',(T.cur), " finished and ", signif(i/N*100,digits = 2) ,"%")
            
          
          
        }
 
    }

    
    
    
    
  }
  
  readRDS(file=save_to)
}  




  


