## load and install multiple packages
  install.packages(c('dplyr', 'broom', 'ggplot2','ggpubr', 'glmnet', 
                     'caret', 'hablar', 'doMC', 'msaenet', 'tidyr'))
  packages <- c("dplyr", "broom", "ggplot2", "ggpubr", "glmnet", 
                "caret", "hablar", "doMC", "msaenet", "tidyr")
  lapply(packages, require, character.only = TRUE)
 
## set working directory
  setwd("/p01-hdd/dsb/asever/For_Aleksandar/SCRIPTS")

## load the data frame
  df <- read.csv("Aggregated_HCPMMP_manual_edits_Time_01.csv")
  ICV <- read.csv("ICV_subjects_manual_edits_Time_01")
  df <- cbind(df[,1:5], ICV, df[,6:length(df)])

## create sub data frames only containing the variables of interest
  df_thick <- select(df, contains("thickness")) %>%
    scale()
  df_thick <- na.replace(df_thick, colMeans(df_thick, na.rm = TRUE)) %>%
    as.data.frame()
  
  df_area <- select(df, contains("area"))
  df_area$ICV <- df$ICV
  

##### ICV normalization for area parcel values #####
  coef_area <- list()
  resid_area <- list()
  resi_norm_area <- list()
  
 ##### Loop through each variable (parcel) and create new data frame containing only the normalized residuals
  
  for (c in 1:(ncol(df_area)-1)){
    model_area <- lm(as.formula(paste0(colnames(df_area)[c] , "~ ICV")), data = df_area)
    coef_area[[c]] <- model_area$coefficients[2]    #the 2 stands for B1 --> the slope. we need this one.  
    resid_area[[c]] <- residuals(model_area)
    resi_norm_area[[c]] <- residuals(model_area) + (model_area$coefficients[2]*mean(df_area$ICV))   #residual method as stated in voevodskaya paper.  
    #fig <- paste0('Resinorm_', colnames(df_area)[c], 'jpg')    #for fun - plot them and save in the working directory.  
    #jpeg(file=fig)
    #plot(resi_norm_area[[c]])
    #dev.off() 
  }
  coef_area <- as.data.frame(do.call(rbind,coef_area))
  resi_norm_area <- as.data.frame(t(do.call(rbind,resi_norm_area)))
  resi_area <- as.data.frame(t(do.call(rbind,resid_area)))
  
  
  colnames(resi_norm_area) <- colnames(df_area)[1:(ncol(df_area)-1)]   #reassign colnames 
  
  df_area <- resi_norm_area #df_area gets reassigned the values of the normalized residuals
  
  df_area <- df_area %>%
    scale()
  df_area <- na.replace(df_area, colMeans(df_area, na.rm = TRUE)) %>%
    as.data.frame()
  
##################################### 
#multistep adaptive elastic net 

  
#why do a msaenet? --> the p values of a normal elastic net regression are biased. so if we want to be precise we can use the multistep adaptive elastic net
#because of the oracly property this regression includes, the p values won't be biased anymore. 
#only problem is, how do we extract them? --> still looking for a solution there.  
  
#also, not sure if it would be much simpler to just calculate elastic net and then proceed with a simple multiple regression using all non zero variables.  

#first: clean the data and remove subjects that have missing data on one of the two behavior measures
#5 Subjects have missing data on Aksara_Recognition and Word_Reading. There is no subj. missing data on
#one, but not the other, so cleaning for one is sufficient.
  y_aksara <- df$Aksara_Recogntion     
  y_word <- df$Word_Reading
  
  #remove subjects from df_thick/area/curv that have missing data in literacy measures and convert them to matrices
  x_thick <- df_thick %>%
    slice(-which(is.na(y_aksara))) %>%
    droplevels() %>%
    as.matrix()
  
  x_area <- df_area %>%
    slice(-which(is.na(y_aksara))) %>%
    droplevels() %>%
    as.matrix()
  
  #update y_aksara and y_word to contain only subjects without missing data
  y_aksara <- y_aksara[!is.na(y_aksara)]
  y_word <- y_word[!is.na(y_word)]
  
####### models for variable thickness #######################################################

##### aksara
  msaenet.fit_thick_aksara <- msaenet(x_thick, y_aksara,
                                                   alphas = seq(0.1, 0.9,0.1),  #alphas from zero to 1 in 0.1 steps
                                                   nsteps = 10L, tune.nsteps = "ebic",
                                                   seed = 1005) 
  
  msaenet.fit_thick_aksara
  plot(msaenet.fit_thick_aksara, label = TRUE)
  
  plot(msaenet.fit_thick_aksara, type = "dotplot", label = TRUE, label.cex = 1)
 
  summary(msaenet.fit_thick_aksara)                  
  
  print(msaenet.fit_thick_aksara)
  #Df      %Dev       Lambda
  #5 0.3327544 4.195124e+16 Whi is lambda so high?
  #model prediction 
  
  msaenet.pred_thick_aksara <- predict(msaenet.fit_thick_aksara, x_thick)
  msaenet.pred_thick_aksara
  msaenet.rmse(y_aksara, msaenet.pred_thick_aksara) #14.14126  
  msaenet.rmsle(y_aksara, msaenet.pred_thick_aksara)  #logarithmic:   1.436366 
  plot(msaenet.fit_thick_aksara)
  
  
  
  # Multiple R-squared
  rsq_msaenet <- cor(y_aksara, msaenet.pred_thick_aksara)^2
  rsq_msaenet #0.1887373
  
  
  #extract model coefficients:
  coef(msaenet.fit_thick_aksara)
  var_thick_aksara <- which(coef(msaenet.fit_thick_aksara) != 0)
  colnames(df_thick[var_thick_aksara])
  #"average_cortical_thickness_L_MST"  "average_cortical_thickness_L_8Ad"  "average_cortical_thickness_L_p10p"
  #"average_cortical_thickness_R_STV"  "average_cortical_thickness_R_MI"
  
  msaenet.fn(msaenet.fit_thick_aksara, 1:5)#number of positive values   
  msaenet.fp(msaenet.fit_thick_aksara, 1:5)
  
##### word
  msaenet.fit_thick_word <- msaenet(x_thick, y_word,
                                    alphas = seq(0.1, 0.9,0.1),  #alphas from zero to 1 in 0.1 steps
                                    nsteps = 10L, tune.nsteps = "ebic",
                                    seed = 1005) 
  
  msaenet.fit_thick_word
  # Df      %Dev       Lambda
  # 17 0.7568534 6.348886e+16
  plot(msaenet.fit_thick_word, label = TRUE)
  
  plot(msaenet.fit_thick_word, type = "dotplot", label = TRUE, label.cex = 1)
 
  summary(msaenet.fit_thick_word)                  
  
  print(msaenet.fit_thick_word)
  
  #model prediction 
  
  msaenet.pred_thick_word <- predict(msaenet.fit_thick_word, x_thick)
  msaenet.pred_thick_word
  msaenet.rmse(y_word, msaenet.pred_thick_word) #16.52394 
  msaenet.rmsle(y_word, msaenet.pred_thick_word)  #[1] NaN Warning message: In log(ypred + 1) : NaNs produced ??????????
  plot(msaenet.fit_thick_word)
  
  
  
  # Multiple R-squared
  rsq_msaenet <- cor(y_word, msaenet.pred_thick_word)^2
  rsq_msaenet  
  # 0.7627011
  
  #extract model coefficients:
  coef(msaenet.fit_thick_word)
  var_thick_word <- which(coef(msaenet.fit_thick_word) != 0)
  colnames(df_thick[var_thick_word])
  #"average_cortical_thickness_L_V1"   "average_cortical_thickness_L_MST"  "average_cortical_thickness_L_55b" 
  #"average_cortical_thickness_L_LIPv" "average_cortical_thickness_L_8BL"  "average_cortical_thickness_L_9p"  
  #"average_cortical_thickness_L_LIPd" "average_cortical_thickness_L_FOP1" "average_cortical_thickness_L_p10p"
  #"average_cortical_thickness_L_A4"   "average_cortical_thickness_R_V4"   "average_cortical_thickness_R_23d" 
  #"average_cortical_thickness_R_10v"  "average_cortical_thickness_R_LIPd" "average_cortical_thickness_R_52"  
  #"average_cortical_thickness_R_FOP3" "average_cortical_thickness_R_A4"  
  
  msaenet.fn(msaenet.fit_thick_word, 1:5)#number of positive values   
  msaenet.fp(msaenet.fit_thick_word, 1:5)
  
#######models for variable area#######################################################

##### aksara
  msaenet.fit_area_aksara <- msaenet(x_area, y_aksara,
                                     alphas = seq(0.1, 0.9,0.1),  #alphas from zero to 1 in 0.1 steps
                                     nsteps = 10L, tune.nsteps = "ebic",
                                     seed = 1005) 
  
  msaenet.fit_area_aksara
  plot(msaenet.fit_area_aksara, label = TRUE)
  
  plot(msaenet.fit_area_aksara, type = "dotplot", label = TRUE, label.cex = 1)
  #"total_surface_area_0x28mm0x5E20x29_L_8Av"  "total_surface_area_0x28mm0x5E20x29_L_STGa"
  summary(msaenet.fit_area_aksara)                  
  
  print(msaenet.fit_area_aksara)
  
  #model prediction 
  
  msaenet.pred_area_aksara <- predict(msaenet.fit_area_aksara, x_area)
  msaenet.pred_area_aksara
  msaenet.rmse(y_aksara, msaenet.pred_area_aksara) #15.35709  
  msaenet.rmsle(y_aksara, msaenet.pred_area_aksara)  #logarithmic:   1.427725 
  plot(msaenet.fit_area_aksara)
  
  
  
  # Multiple R-squared
  rsq_msaenet <- cor(y_aksara, msaenet.pred_area_aksara)^2
  rsq_msaenet  
  #0.2131158
  
  #extract model coefficients:
  coef(msaenet.fit_area_aksara)
  var_area_aksara <- which(coef(msaenet.fit_area_aksara) != 0)
  colnames(df_area[var_area_aksara])
  #"total_surface_area_0x28mm0x5E20x29_L_8Av"  "total_surface_area_0x28mm0x5E20x29_L_STGa"
  msaenet.fn(msaenet.fit_area_aksara, 1:5)#number of positive values   
  msaenet.fp(msaenet.fit_area_aksara, 1:5)
  
##### word 
  msaenet.fit_area_word <- msaenet(x_area, y_word,
                                   alphas = seq(0.1, 0.9,0.1),  #alphas from zero to 1 in 0.1 steps
                                   nsteps = 10L, tune.nsteps = "ebic",
                                   seed = 1005) 
  
  msaenet.fit_area_word
  #Df       %Dev       Lambda
  # 1 0.08445994 2.645446e+16
  plot(msaenet.fit_area_word, label = TRUE)
  
  plot(msaenet.fit_area_word, type = "dotplot", label = TRUE, label.cex = 1)
 
  summary(msaenet.fit_area_word)                  
  
  print(msaenet.fit_area_word)
  
  #model prediction 
  
  msaenet.pred_area_word <- predict(msaenet.fit_area_word, x_area)
  msaenet.pred_area_word
  msaenet.rmse(y_word, msaenet.pred_area_word) #32.06404  
  msaenet.rmsle(y_word, msaenet.pred_area_word)  #logarithmic:   2.215139 
  plot(msaenet.fit_area_word)
  
  
  
  # Multiple R-squared
  rsq_msaenet <- cor(y_word, msaenet.pred_area_word)^2
  rsq_msaenet  
  #0.08456501
  
  #extract model coefficients:
  coef(msaenet.fit_area_word)
  var_area_word <- which(coef(msaenet.fit_area_word) != 0)
  colnames(df_area[var_area_word])
  #total_surface_area_0x28mm0x5E20x29_R_RSC
  msaenet.fn(msaenet.fit_area_word, 1:5)#number of positive values   
  msaenet.fp(msaenet.fit_area_word, 1:5)
  
###############################
# elastic net regression 
##############################

  #include behavior data in the measures data frame separately and update them to only contain the subjects with behavior data
  
  df_thick_aksara <- df_thick %>% 
    cbind(df$Aksara_Recogntion) %>% 
    rename("aksara" = "df$Aksara_Recogntion") %>%
    slice(-which(is.na(aksara))) %>%
    droplevels()
  df_thick_word <- df_thick %>% 
    cbind(df$Word_Reading) %>% 
    rename("word" = "df$Word_Reading") %>%
    slice(-which(is.na(word))) %>%
    droplevels()
  
  df_area_aksara <- df_area %>% 
    cbind(df$Aksara_Recogntion) %>% 
    rename("aksara" = "df$Aksara_Recogntion") %>%
    slice(-which(is.na(aksara))) %>%
    droplevels()
  df_area_word <- df_area %>% 
    cbind(df$Word_Reading) %>% 
    rename("word" = "df$Word_Reading") %>%
    slice(-which(is.na(word))) %>%
    droplevels()

  
##### thickness ################################################################
##### aksara 
  control <- trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 5,
                          search = "random",
                          verboseIter = TRUE)
  
  
  elastic_modelelastic_model <- train(aksara ~ .,
                         data = df_thick_aksara,
                         method = "glmnet",
                         preProcess = c("center", "scale"),
                         tuneGrid =expand.grid(alpha=seq(0,1,length=100),
                                               lambda = seq(0.0001,10,length=100)), #if I substitute the tunegrid parameter with just using tunelength = 25 i get sensible values for alpha and lambda
                         trControl = control)
  
  #If not using the tunegrid specification, i get the following. Have to find out the rationale behind it.
  #Aggregating results
  #Selecting tuning parameters
  #Fitting alpha = 0.538, lambda = 5.03 on full training set
  mean(elastic_model$resample$RMSE)
  
  
  
  #plotting the model
  
  plot(elastic_model, main = "Elastic Net Regression")
  
  #plotting important variables
  
  plot(varImp(elastic_model,scale=TRUE))
  
  
   varImp(elastic_model)
   
  # average_cortical_thickness_L_8Ad   100.000
  # average_cortical_thickness_R_MI     93.169
  # average_cortical_thickness_L_MST    79.412
  # average_cortical_thickness_L_p10p   76.551
  # average_cortical_thickness_R_STV    63.134
  # average_cortical_thickness_L_LIPv   40.089
  # average_cortical_thickness_R_LIPd   37.970
  # average_cortical_thickness_L_FOP5   28.772
  # average_cortical_thickness_R_23c    28.240
  # average_cortical_thickness_L_IFJp   26.163
  # average_cortical_thickness_R_LO2     7.376
  # average_cortical_thickness_R_TA2     2.223

##### word
  control <- trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 5,
                          search = "random",
                          verboseIter = TRUE)
  
  elastic_model <- train(word ~ .,
                         data = df_thick_word,
                         method = "glmnet",
                         preProcess = c("center", "scale"),
                         tuneGrid =expand.grid(alpha=seq(0,1,length=100),
                                               lambda = seq(0.0001,10,length=100)),
                         trControl = control)
  #Aggregating results, using tunelength and ditching the tunegrid parameter. Else i get lambda =10 and alpha =0
  #Selecting tuning parameters
  # Fitting alpha = 0.052, lambda = 0.0317 on full training set
  # elastic_model
  
  
  
  mean(elastic_model$resample$RMSE)
  #32.48116
  
  
  #plotting the model
  
  plot(elastic_model, main = "Elastic Net Regression")
  
  #plotting important variables
  
  plot(varImp(elastic_model,scale=TRUE))
  
  
  varImp(elastic_model)
  # average_cortical_thickness_R_52    100.00
  # average_cortical_thickness_L_LIPd   89.30
  # average_cortical_thickness_R_VIP    77.66
  # average_cortical_thickness_R_LIPd   76.74
  # average_cortical_thickness_L_5mv    74.62
  # average_cortical_thickness_L_MST    67.33
  # average_cortical_thickness_R_1      62.48
  # average_cortical_thickness_R_23d    60.44
  # average_cortical_thickness_R_TE2p   57.04
  # average_cortical_thickness_R_10v    56.30
  # average_cortical_thickness_L_LIPv   53.34
  # average_cortical_thickness_L_V1     52.31
  # average_cortical_thickness_L_p10p   51.61
  # average_cortical_thickness_L_V4t    51.50
  # average_cortical_thickness_L_IFSp   50.33
  # average_cortical_thickness_L_A4     50.31
  # average_cortical_thickness_L_31a    47.73
  # average_cortical_thickness_R_TGd    47.66
  # average_cortical_thickness_R_IP1    46.79
  # average_cortical_thickness_R_STV    45.49
##### area #####################################################################
##### aksara #####
  control <- trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 5,
                          search = "random",
                          verboseIter = TRUE)
  
  elastic_model <- train(aksara ~ .,
                         data = df_area_aksara,
                         method = "glmnet",
                         preProcess = c("center", "scale"),
                         tuneGrid =expand.grid(alpha=seq(0,1,length=100),
                                               lambda = seq(0.0001,10,length=100)),
                         trControl = control)
  # Aggregating results
  # Selecting tuning parameters
  # Fitting alpha = 1, lambda = 4.65 on full training set
  # Warning message:
  # In nominalTrainWorkflow(x = x, y = y, wts = weights, info = trainInfo,  :
  #                           There were missing values in resampled performance measures.
  elastic_model
  
  
  
  mean(elastic_model$resample$RMSE)
  #17.05239
  
  
  #plotting the model
  
  plot(elastic_model, main = "Elastic Net Regression")
  
  #plotting important variables
  
  plot(varImp(elastic_model,scale=TRUE))
  
  
  varImp(elastic_model)
  # total_surface_area_0x28mm0x5E20x29_L_8Av       100.00
  # total_surface_area_0x28mm0x5E20x29_L_STGa       55.11
##### word #####
  control <- trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 5,
                          search = "random",
                          verboseIter = TRUE)
  
  elastic_model <- train(word ~ .,
                         data = df_area_word,
                         method = "glmnet",
                         preProcess = c("center", "scale"),
                         tuneGrid =expand.grid(alpha=seq(0,1,length=100),
                                               lambda = seq(0.0001,10,length=100)),
                         trControl = control)
  #Fitting alpha = 1, lambda = 10 on full training set (lambda = 10???)
  #Fitting alpha = 0.994, lambda = 7.43 on full training set if not using tuneGrid**
  elastic_model
  
  
  
  mean(elastic_model$resample$RMSE)
  #34.25432
  #35.09228**
  
  #plotting the model
  
  plot(elastic_model, main = "Elastic Net Regression")
  
  #plotting important variables
  
  plot(varImp(elastic_model,scale=TRUE))
  
  
  varImp(elastic_model)
  #All NaN
  # If **
  # total_surface_area_0x28mm0x5E20x29_R_TE1p   100.00
  # total_surface_area_0x28mm0x5E20x29_R_RSC     98.66
  # total_surface_area_0x28mm0x5E20x29_L_PFm     96.59
  # total_surface_area_0x28mm0x5E20x29_L_STGa    63.19
  # total_surface_area_0x28mm0x5E20x29_L_8Av     59.22
  # total_surface_area_0x28mm0x5E20x29_R_TGd     51.18
  # total_surface_area_0x28mm0x5E20x29_L_TA2     48.13
  # total_surface_area_0x28mm0x5E20x29_L_FOP3    46.63
  # total_surface_area_0x28mm0x5E20x29_L_10d      0.00
  # total_surface_area_0x28mm0x5E20x29_R_SCEF     0.00
  

  
  
  
  