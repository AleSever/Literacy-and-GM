## load and install multiple packages
  install.packages(c('dplyr', 'broom', 'ggplot2','ggpubr', 'glmnet', 
                     'caret', 'hablar', 'doMC', 'msaenet', 'tidyr', 'apaTables'))
  packages <- c("dplyr", "broom", "ggplot2", "ggpubr", "glmnet", 
                "caret", "hablar", "doMC", "msaenet", "tidyr", "apaTables")
  lapply(packages, require, character.only = TRUE)
## set working directory
  setwd("/p01-hdd/dsb/asever/For_Aleksandar/SCRIPTS")

## load the data frame
  df <- read.csv("Aggregated_HCPMMP_manual_edits_Time_01.csv")
  ICV <- read.csv("ICV_subjects_manual_edits_Time_01")
  df <- cbind(df[,1:5], ICV, df[,6:length(df)])

## create sub data frames only containing the variables of interest, extract the behavior variables aksara and word.
## remove subjects that have missing data on the behavior measures from the data frames
  y_aksara <- df$Aksara_Recogntion 
  y_word <- df$Word_Reading
  
  #thickness
  df_thick <- select(df, contains("thickness")) %>%
    slice(-which(is.na(y_aksara))) %>%
    droplevels()
  
  df_thick <- na.replace(df_thick, colMeans(df_thick, na.rm = TRUE)) %>%
    scale() %>%
    as.data.frame()
  
  #correlation among the variables in df_thick
  cor_matrix <- cor(df_thick)
  cor_vector <- cor_matrix[lower.tri(cor_matrix)]
  mean(cor_vector)#0.3811514
  sd(cor_vector) #0.1650645
  
  # update the behavior measures, removing subjects with missing data
  y_aksara <- y_aksara[!is.na(y_aksara)]
  y_word <- y_word[!is.na(y_word)]
  
  # plot them to reason binary coding
  hist(y_aksara)
  hist(y_word, ylab = "N. of Subjects", xlab = "N. of words correct", main = "")
  
  # produce binary coding. If you get less than half right you get assigned
  # Illiterate(0) and at least half right Literate(1)
  y_aksara <- ifelse(y_aksara>=23, 1, 0)
  y_word <- ifelse(y_word>=43, 1, 0)
  
  # get the sex variable to show females and males in each group
  sex <- df$Sex
  sex <- sex[!is.na(sex)]
  
  df_sex_word <- data.frame(y_word, sex)
  table(df_sex_word$y_word)
  # 0  1  number of participants in the illiterate 0 and literate 1 group after the binary coding
  # 52 24 
  table(df_sex_word$sex)
  # 0  1  number of females 0 and males 1
  # 46 30 
  sum(df_sex_word$sex[df_sex_word$y_word == 1] == 0, na.rm = TRUE)
  # 8 females in the literate group
  
  
##################################### 
#multistep adaptive elastic net 

#convert df's to matrices
  x_thick <- as.matrix(df_thick)
  x_thick_ICV <- as.matrix(df_thick_ICV)


  
####### models for variable thickness #######################################################
##### aksara
  msaenet.fit_thick_aksara <- msaenet(x_thick, y_aksara, family = "binomial",
                                      alphas = seq(0.2, 0.8, 0.2),  
                                      nsteps = 10L, 
                                      tune.nsteps = "ebic", 
                                      nfolds = 10,
                                      seed = 1005) 
                                                     
  
  msaenet.fit_thick_aksara
  # Df     %Dev       Lambda
  # 1  4 0.264749 7.804543e+13
  plot(msaenet.fit_thick_aksara, label = TRUE)
  plot(msaenet.fit_thick_aksara, type = "dotplot", label = TRUE, label.cex = 1)
  
  #extract model coefficients and the corresponding variable names:
  coef_thick_aksara <- coef(msaenet.fit_thick_aksara)%>%
    as.data.frame()
  write.csv(coef_thick_aksara, file = "Results_T1/coef_thick_aksara_full.csv")
  var_thick_aksara <- msaenet.nzv(msaenet.fit_thick_aksara)
  colnames(df_thick[var_thick_aksara])
  #"average_cortical_thickness_L_LIPv" "average_cortical_thickness_L_IFJp" "average_cortical_thickness_R_LO2"  "average_cortical_thickness_R_STV" 
  
  ## train model on training data set and test on test set 100 times
  ## output: vector predictions (cor of predicted and actual values), pseudo-R² for each run,
  ## and list of coefficients selected in each run
  
  predictions_thick_aksara <- numeric(100)
  r_sqr_thick_aksara <- numeric(100)
  coef_list_thick_aksara <- list()
  
  # train .8 and test .2 of the full data 100 times
  for (i in 1:100) {
    skip_iteration <- FALSE
    
    tryCatch({
      set.seed(i*277)
      # create train and test indices
      train_index <- sample(nrow(x_thick), floor(0.8 * nrow(x_thick)))
      test_index <- setdiff(1:nrow(x_thick), train_index)
      
      # split data into training and test sets
      x_thick_tr <- x_thick[train_index, ]
      x_thick_te <- x_thick[test_index, ]
      y_aksara_tr <- y_aksara[train_index]
      y_aksara_te <- y_aksara[test_index]
      
      # train model
      msaenet.fit_thick_aksara_tr <- msaenet(x_thick_tr, y_aksara_tr,
                                             alphas = seq(0.2, 0.8, 0.2),  
                                             nsteps = 10L, 
                                             tune.nsteps = "ebic", 
                                             nfolds = 10,
                                             seed = 1005)
      
      # predict y-values for test data
      msaenet.pred_thick_aksara <- predict(msaenet.fit_thick_aksara_tr, x_thick_te)
      msaenet_sqr <- cor(predict(msaenet.fit_thick_aksara_tr, x_thick_tr), y_aksara_tr)^2
      coef <- msaenet.nzv(msaenet.fit_thick_aksara_tr)
      # calculate correlation between predicted and actual y-values of the test set
      pred_thick_aksara <- cor(y_aksara_te, msaenet.pred_thick_aksara)
      
      # store the result in the results vectors
      predictions_thick_aksara[i] <- pred_thick_aksara
      r_sqr_thick_aksara[i] <- msaenet_sqr
      coef_list_thick_aksara <- append(coef_list_thick_aksara, list(coef))
    }, error = function(e) {
      # Handle the specific error
      if (grepl("Null model produced by the full fit", e$message)) {
        # Set the flag to skip the iteration
        skip_iteration <- TRUE
      }
      
      # Handle other errors if needed
      # ...
    })
    
    if (skip_iteration) {
      next
    }
  }
  
  # Here we separate the variables that were selected at each iteration. 
  # The total count of selections is then exported for visualization.
  var_frequency_thick_aksara <- table(unlist(coef_list_thick_aksara)) %>%
    as.data.frame()%>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric)
  c <- numeric(358)
  for (i in 1:length(var_frequency_thick_aksara$Var1)){
    a = var_frequency_thick_aksara[i,1]
    c[a] = var_frequency_thick_aksara[i,2]
  }
 
  write.csv(c, file = "Results_T1/coef_thick_aksara.csv")    

##### word
  msaenet.fit_thick_word <- msaenet(x_thick, y_word, family = "binomial",
                                    alphas = seq(0.2, 0.8, 0.2),  
                                    nsteps = 10L, 
                                    tune.nsteps = "ebic", 
                                    nfolds = 10,
                                    seed = 1005) 
  
  
  msaenet.fit_thick_word
  # Df      %Dev       Lambda
  # 1 11 0.7204229 5.039106e+13
  plot(msaenet.fit_thick_word, label = TRUE)
  plot(msaenet.fit_thick_word, type = "dotplot", label = TRUE, label.cex = 1)
  #extract model coefficients and the corresponding variable names:
  coef_thick_word <- coef(msaenet.fit_thick_word)%>%
    as.data.frame()
  write.csv(coef_thick_word, file = "Results_T1/coef_thick_word_full.csv")
  var_thick_word <- msaenet.nzv(msaenet.fit_thick_word)
  colnames(df_thick[var_thick_word])
  # [1] "average_cortical_thickness_L_V1"    "average_cortical_thickness_L_MST"   "average_cortical_thickness_L_8BM"   "average_cortical_thickness_L_9p"   
  # [5] "average_cortical_thickness_R_LO2"   "average_cortical_thickness_R_STV"   "average_cortical_thickness_R_LIPd"  "average_cortical_thickness_R_52"   
  # [9] "average_cortical_thickness_R_FOP3"  "average_cortical_thickness_R_PBelt" "average_cortical_thickness_R_PHA2" 
  
  ## train model on training data set and test on test set with 100 iterations.
  predictions_thick_word <- numeric(100)
  r_sqr_thick_word <- numeric(100)
  coef_list_thick_word <- list()
  
  # train .8 and test .2 of the full data 100 times
  for (i in 1:100) {
    set.seed(i*277)
    # create train and test indices
    train_index <- sample(nrow(x_thick), floor(0.8 * nrow(x_thick)))
    test_index <- setdiff(1:nrow(x_thick), train_index)
    
    # split data into training and test sets
    x_thick_tr <- x_thick[train_index, ]
    x_thick_te <- x_thick[test_index, ]
    y_word_tr <- y_word[train_index]
    y_word_te <- y_word[test_index]
    
    # train model
    msaenet.fit_thick_word_tr <- msaenet(x_thick_tr, y_word_tr, family = "binomial",
                                         alphas = seq(0.2, 0.8, 0.2),  
                                         nsteps = 10L, 
                                         tune.nsteps = "ebic", 
                                         nfolds = 10,
                                         seed = 1005)
    
    
    # predict y-values for test data
    msaenet.pred_thick_word <- predict(msaenet.fit_thick_word_tr, x_thick_te)
    msaenet_sqr <- cor(predict(msaenet.fit_thick_word_tr, x_thick_tr), y_word_tr)^2
    coef <- msaenet.nzv(msaenet.fit_thick_word_tr)
    # calculate correlation between predicted and actual y-values
    pred_thick_word <- cor(y_word_te, msaenet.pred_thick_word)
    
    # store the result in the results vector
    predictions_thick_word[i] <- pred_thick_word
    r_sqr_thick_word[i] <- msaenet_sqr
    coef_list_thick_word <- append(coef_list_thick_word, list(coef))
  }
  
  # get the mean and standard dev. of the R² and the predictions and their correlation
  mean(predictions_thick_word) # 0.251567
  sd(predictions_thick_word) # 0.227788
  mean(r_sqr_thick_word)# 0.5577387
  sd(r_sqr_thick_word) #0.1291793
  cor.test(predictions_thick_word, r_sqr_thick_word)
  # -0.4079132
  
  # Separate the variables that were selected at each iteration. 
  # The total count of selections is then exported for visualization.
  var_frequency_thick_word <- table(unlist(coef_list_thick_word)) %>%
    as.data.frame()%>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric)
  c <- numeric(358)
  for (i in 1:length(var_frequency_thick_word$Var1)){
    a = var_frequency_thick_word[i,1]
    c[a] = var_frequency_thick_word[i,2]
  }
  
  write.csv(c, file = "Results_T1/coef_thick_word.csv")
 
  # the mean and standard deviation of the 20 iterations with the highest prediction correlation
  mean(head(sort(predictions_thick_word, decreasing = TRUE), 20))
  # 0.5644173
  sd(head(sort(predictions_thick_word, decreasing = TRUE), 20))
  # 0.10661
  
  # extract the variables that go selected in the 20 iterations that had he highest predictive accuracy for visualization
  top20_predictive_iterations <- coef_list_thick_word[order(predictions_thick_word, decreasing = TRUE)[1:20]]
  top20_predictive_iterations <- table(unlist(top20_predictive_iterations)) %>%
    as.data.frame()%>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric)
  c2 <- numeric(358)
  for (i in 1:length(top20_predictive_iterations$Var1)){
    a = top20_predictive_iterations[i,1]
    c2[a] = top20_predictive_iterations[i,2]
  }
  write.csv(c2, file = "Results_T1/coef_thick_word_top20.csv")
  
  top20_predictive_iterations$Var1 <- colnames(df_thick)[top20_predictive_iterations$Var1]
  top20_predictive_iterations <- top20_predictive_iterations %>%
    filter(Freq > 5) %>%
    arrange(desc(Freq))
  colnames(top20_predictive_iterations)[1] <- "Region"
  colnames(top20_predictive_iterations)[2] <- "Frequency"
  top20_predictive_iterations
  # Region                                    Frequency
  # 1  average_cortical_thickness_L_MST        15
  # 2  average_cortical_thickness_R_LO2        13
  # 3 average_cortical_thickness_R_FOP3        11
  # 4   average_cortical_thickness_R_52         9
  # 5  average_cortical_thickness_L_8BM         7
  # 6 average_cortical_thickness_R_LIPd         7
  # 7 average_cortical_thickness_R_PHA2         6
  
  # get names of the variables that were selected 20 or more times and their corresponding selection count
  var_frequency_thick_word$Var1 <- colnames(df_thick)[var_frequency_thick_word$Var1]
  var_frequency_thick_word <- var_frequency_thick_word %>%
    filter(Freq > 20) %>%
    arrange(desc(Freq))
  colnames(var_frequency_thick_word)[1] <- "Region"
  colnames(var_frequency_thick_word)[2] <- "Frequency"
  var_frequency_thick_word
  # Region Frequency
  # 1  average_cortical_thickness_L_MST                      69
  # 2 average_cortical_thickness_R_FOP3                      55
  # 3   average_cortical_thickness_R_52                      52
  # 4  average_cortical_thickness_R_LO2                      49
  # 5 average_cortical_thickness_R_LIPd                      37
  # 6 average_cortical_thickness_R_PHA2                      24

  
  
  # create boxplots with the regions that were selected more than 5 times in the 20 iterations with the highest predictive accuracy
  y_word <- df$Word_Reading
  df_thick <- select(df, contains("thickness")) %>%
    slice(-which(is.na(y_word))) %>%
    droplevels()
  y_word <- y_word[!is.na(y_word)]
  y_word <- ifelse(y_word>=43, 1, 0)
  df_thick$word <- y_word
  df_thick$word <- ifelse(df_thick$word == 1, "Literate", "Illiterate")

  t.test(df_thick$average_cortical_thickness_L_MST ~ word, data = df_thick)
  # t = -5.1071, df = 46.665, p-value = 5.942e-06
  boxplot(df_thick$average_cortical_thickness_L_MST ~ word, data = df_thick, main = "L_MST", xlab = "", ylab = "Cortical Thickness")
  
  t.test(df_thick$average_cortical_thickness_R_FOP3 ~ word, data = df_thick)
  # t = -3.6279, df = 45.323, p-value = 0.0007223
  boxplot(df_thick$average_cortical_thickness_R_FOP3 ~ word, data = df_thick, main = "R_FOP3", xlab = "", ylab = "Cortical Thickness")

  t.test(df_thick$average_cortical_thickness_R_52 ~ word, data = df_thick)
  # t = 0.855, df = 41.864, p-value = 0.3974
  boxplot(df_thick$average_cortical_thickness_R_52 ~ word, data = df_thick, main = "R_52", xlab = "", ylab = "Cortical Thickness")
  
  t.test(df_thick$average_cortical_thickness_R_LO2 ~ word, data = df_thick)
  # t = -4.2134, df = 37.192, p-value = 0.0001537
  boxplot(df_thick$average_cortical_thickness_R_LO2 ~ word, data = df_thick, main = "R_LO2", xlab = "", ylab = "Cortical Thickness")
  
  t.test(df_thick$average_cortical_thickness_R_LIPd ~ word, data = df_thick)
  # t = 0.9631, df = 34.004, p-value = 0.3423
  boxplot(df_thick$average_cortical_thickness_R_LIPd ~ word, data = df_thick, main = "R_LIPd", xlab = "", ylab = "Cortical Thickness")
  
  t.test(df_thick$average_cortical_thickness_R_PHA2 ~ word, data = df_thick)
  # t = 0.36291, df = 34.528, p-value = 0.7189
  boxplot(df_thick$average_cortical_thickness_R_PHA2 ~ word, data = df_thick, main = "R_PHA2", xlab = "", ylab = "Cortical Thickness")
  
  
  boxplot(df_thick$average_cortical_thickness_L_8BL ~ word, data = df_thick, main = "L_8BM", xlab = "", ylab = "Cortical Thickness")
  # excluding the outliers from 8BM
  literate_data <- subset(df_thick, word == "Literate")
  illiterate_data <- subset(df_thick, word == "Illiterate")
  max_literate <- max(literate_data$average_cortical_thickness_L_8BM)
  max_illiterate <- max(illiterate_data$average_cortical_thickness_L_8BM)
  literate_data <- subset(literate_data, average_cortical_thickness_L_8BM < max_literate)
  illiterate_data <- subset(illiterate_data, average_cortical_thickness_L_8BM < max_illiterate)
  t.test(average_cortical_thickness_L_8BM ~ word, data = rbind(literate_data, illiterate_data))
  #t = -4.9482, df = 37.114, p-value = 1.644e-05
  
  
  
  
  