library(ggplot2)
library(pROC)
library(DescTools)
library(Hmisc)
library(plyr)
library(boot)
library(MASS)
library(rms)
library(readr)
library(table1)

set.seed(2313)


#View(imputed_df)
hist <- #data before the shift (historical)
target <- #data after the shift (target population)

hist_Time <- hist; hist_Time$Time <- 0 #add dummy variable for the target population membership as feature
target_Time <- target; target_Time$Time <- 1  #add dummy variable for the target population membership as feature


######################################################
#membership model weighting functions; original weighting, weights limited to one, 
#and weights scaled with scaling factor

propensity_weighting_original <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist$set <- 0
  target$set <- 1
  all<- rbind(hist, target)
  
  membership_mod <- glm(set~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=all, family = 'binomial')
  preds <- predict(membership_mod,type="response", newdata = hist)
  propensity_weight <- preds / (1 - preds)
  
  model_weights <- c(propensity_weight, target$set)

  return(model_weights)
}

propensity_weighting_limit_one <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist$set <- 0
  target$set <- 1
  all<- rbind(hist, target)
  
  membership_mod <- glm(set~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=all, family = 'binomial')
  preds <- predict(membership_mod,type="response", newdata = hist)
  
  propensity_weight <- preds / (1 - preds)
  propensity_weight[propensity_weight > 1]= 1
  
  model_weights <- c(propensity_weight, target$set)

  
  return(model_weights)
}



propensity_weighting_scaling <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist$set <- 0
  target$set <- 1
  all<- rbind(hist, target)
  
  membership_mod <- glm(set~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=all, family = 'binomial')
  preds <- predict(membership_mod,type="response", newdata = hist)
  
  propensity_weight <- preds / (1 - preds)
  
  propensity_weight_scaled <- 0
  if(sum(propensity_weight) > nrow(hist)){
    print('scaled')
    scale_factor <- nrow(hist)/sum(propensity_weight)
    propensity_weight_scaled <- propensity_weight * scale_factor
  }else{
    propensity_weight_scaled <- propensity_weight
  }
  
  model_weights <- c(propensity_weight_scaled, target$set) 

  
  return(model_weights)
  
}
#########################
#mahalanobis distance model weighting functions; original weighting, weights limited to one, 
#and weights scaled with scaling factor

mahalanobis_weighting_original <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist_num <- data.frame(lapply(hist,as.numeric))
  target_num <- data.frame(lapply(target,as.numeric))
  
  
  #Mahalanobis distance weights
  distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num)) #get mahalanobis distances
  hist_num$mah_distances <- distances
  
  fx <-  replicate(nrow(hist), 1) #fx uniform distribution
  pvalues <- pchisq(distances, df=7, lower.tail=FALSE) #get pvalues of distances from chi-square distribution
  
  gx<- density(x=pvalues) #estimate gx
  gx_prob <- approxfun(gx$x, gx$y) #get gx probabilities
  
  hist_weights <- fx/gx_prob(pvalues)
  target_weights <- replicate(nrow(target), 1)
  
  model_weights <- c(hist_weights, target_weights)

  return(model_weights)
}



mahalanobis_weighting_limit_one <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist_num <- data.frame(lapply(hist,as.numeric))
  target_num <- data.frame(lapply(target,as.numeric))
  
  
  #Mahalanobis distance weights
  distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num)) #get mahalanobis distanes
  hist_num$mah_distances <- distances
  
  fx <-  replicate(nrow(hist), 1) #fx uniform disrtibution
  pvalues <- pchisq(distances, df=7, lower.tail=FALSE) #get pvalues of disnces from chi-square distribution
  
  gx<- density(x=pvalues) #estimate gx
  gx_prob <- approxfun(gx$x, gx$y) #get gx probabilities
  
  hist_weights <- fx/gx_prob(pvalues)
  hist_weights[hist_weights > 1]= 1
  
  target_weights <- replicate(nrow(target), 1)
  
  model_weights <- c(hist_weights, target_weights)
  
  return(model_weights)
}


mahalanobis_weighting_scaling <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist_num <- data.frame(lapply(hist,as.numeric))
  target_num <- data.frame(lapply(target,as.numeric))
  
  
  #Mahalanobis distance weights
  distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num)) #get mahalanobis distanes
  hist_num$mah_distances <- distances
  
  fx <-  replicate(nrow(hist), 1) #fx uniform disrtibution
  pvalues <- pchisq(distances, df=7, lower.tail=FALSE) #get pvalues of disnces from chi-square distribution
  
  gx<- density(x=pvalues) #estimate gx
  gx_prob <- approxfun(gx$x, gx$y) #get gx probabilities
  
  hist_weights <- fx/gx_prob(pvalues)
  target_weights <- replicate(nrow(target), 1)
  
  hist_weights_scaled <- 0
  
  if(sum(hist_weights) > nrow(hist)){
    print('scaled')
    MD_scale_factor <- nrow(hist)/sum(hist_weights)
    hist_weights_scaled <- hist_weights * MD_scale_factor
  }else{
    hist_weights_scaled <- hist_weights 
  }
  
  model_weights <- c(hist_weights_scaled, target_weights)
  
  
  return(model_weights)
}


#####################

#appent validation results
apparent_results_fun <- function(model, original_samp, model_name){
  val_results <- matrix(nrow = 1,ncol = 5)
  
  pr_test <- predict(model,type="response", newdata = original_samp) # predict probabilities from the bootstrap model in the original sample
  lp_test <- predict(model, newdata = original_samp) # predict lp from the bootstrap model in the original sample
  
  
  # calculate the apparent performance of the model in the original sample
  test_cstat_model <- roc(outcome~pr_test,data=original_samp)
  val_results[1,1] <- test_cstat_model$auc
  test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=original_samp)
  val_results[1,2] <- summary(test_citl_model)$coefficients[1,1]
  test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=original_samp)
  val_results[1,3] <- summary(test_cslope_model)$coefficients[2,1]
  val_results[1,4] <- BrierScore(model)
  val_results[1,5] <- model_name
  
  return(val_results)
}


#bootstrap validation results
validation_results <- function(model, boot_samp, original_samp){
  val_results <- matrix(nrow = 1,ncol = 7)
  
  pr_bs <- predict(model,type="response", newdata = boot_samp) # predict probabilities from the bootstrap model in the bs sample
  lp_bs <- predict(model,newdata = boot_samp ) # predict lp from the bootstrap model in the bs sample
  
  pr_test <- predict(model,type="response", newdata = original_samp) # predict probabilities from the bootstrap model in the original sample
  lp_test <- predict(model, newdata = original_samp) # predict lp from the bootstrap model in the original sample
  
  
  # calculate the apparent performance of the bootstrap model in the bs sample
  app_cstat_model <- roc(outcome~pr_bs,data=boot_samp)
  val_results[1,1] <- app_cstat_model$auc
  app_citl_model <- glm(outcome ~ offset(lp_bs),family=binomial, data=boot_samp)
  val_results[1,2] <- summary(app_citl_model)$coefficients[1,1]
  app_cslope_model <- glm(outcome ~ lp_bs,family=binomial(link='logit'), data=boot_samp)
  val_results[1,3] <- summary(app_cslope_model)$coefficients[2,1]
  
  # calculate the test performance of the bootstrap model in the original sample
  test_cstat_model <- roc(outcome~pr_test,data=original_samp)
  val_results[1,4] <- test_cstat_model$auc
  test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=original_samp)
  val_results[1,5] <- summary(test_citl_model)$coefficients[1,1]
  test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=original_samp)
  val_results[1,6] <- summary(test_cslope_model)$coefficients[2,1]
  val_results[1,7] <- BrierScore(model)
  return(val_results)
}



#propensity models
propensity_weight_original <- propensity_weighting_original(hist, target) #propensity weights
propensity_weight_limit_one <- propensity_weighting_limit_one(hist, target) #propensity weights
propensity_weight_scaling <- propensity_weighting_scaling(hist, target) #propensity weights


propensity_model_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= propensity_weight_original, family='binomial')
propensity_model_Time_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= propensity_weight_original, family='binomial')

propensity_model_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= propensity_weight_limit_one, family='binomial')
propensity_model_Time_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= propensity_weight_limit_one, family='binomial')

propensity_model_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= propensity_weight_scaling, family='binomial')
propensity_model_Time_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= propensity_weight_scaling, family='binomial')


#mahalanobis models

mahalanobis_weight_original <- mahalanobis_weighting_original(hist, target) # mahalanobis weights
mahalanobis_weight_limit_one <- mahalanobis_weighting_limit_one(hist, target) # mahalanobis weights
mahalanobis_weight_scaling <- mahalanobis_weighting_scaling(hist, target) # mahalanobis weights


mahalanobis_model_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= mahalanobis_weight_original, family='binomial')
mahalanobis_model_Time_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= mahalanobis_weight_original, family='binomial')

mahalanobis_model_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= mahalanobis_weight_limit_one, family='binomial')
mahalanobis_model_Time_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= mahalanobis_weight_limit_one, family='binomial')

mahalanobis_model_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), weights= mahalanobis_weight_scaling, family='binomial')
mahalanobis_model_Time_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), weights= mahalanobis_weight_scaling, family='binomial')


#unweighted models on all data
unweighted_model_on_AllData <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), family='binomial')
unweighted_model_on_AllData_Time <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), family='binomial')

#unweighted model on recent data
recent_data_model <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=target, family='binomial')



#get the apparent validation results on target dataset only

apprent_results <- data.frame(matrix(ncol = 5, nrow = 0))
names <- c("app_c_stat","app_citl","app_c_slope","BrierScore","model_name"); colnames(apprent_results) <- names

apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_original, target, 'weighted_propensity_model_original')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_original, target, 'weighted_mahalanobis_model_original')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Time_original, target_Time, 'weighted_propensity_model_Time_original')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Time_original, target_Time, 'weighted_mahalanobis_model_Time_original')

apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_limit_one, target, 'weighted_propensity_model_limit_one')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_limit_one, target, 'weighted_mahalanobis_model_limit_one')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Time_limit_one, target_Time, 'weighted_propensity_model_Time_limit_one')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Time_limit_one, target_Time, 'weighted_mahalanobis_model_Time_limit_one')

apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_scaling, target, 'weighted_propensity_model_scaling')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_scaling, target, 'weighted_mahalanobis_model_scaling')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Time_scaling, target_Time, 'weighted_propensity_model_Time_scaling')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Time_scaling, target_Time, 'weighted_mahalanobis_model_Time_scaling')

apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(unweighted_model_on_AllData, target, 'unweighted_model_on_AllData')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(unweighted_model_on_AllData_Time, target_Time, 'unweighted_model_on_AllData_Time')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(recent_data_model, target, 'unweighted_model_on_recent_data')





########################
#bootstrap Validation
#



manual_boot <- function(hist_set, target_set, target_set_Time,samples){
  results <- data.frame(matrix(ncol = 9, nrow = 0))
  names <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope","BrierScore","booting_round","model_name")
  colnames(results) <- names
  

  set.seed(231398)
  for (i in 1:samples) {
    
    #create bootstrap sample
    samp_index_hist <- sample(1:nrow(hist_set), nrow(hist_set), rep=TRUE) # create a sampling index vector
    samp_index_target <- sample(1:nrow(target_set), nrow(target_set), rep=TRUE) # create a sampling index vector
    
    bs_samp_hist <- hist_set[samp_index_hist,] # index the original dataset using the sampling vector to give the bs sample
    bs_samp_target <- target_set[samp_index_target,]
    
    bs_samp_hist_Time <- bs_samp_hist; bs_samp_hist_Time$Time <- 0
    bs_samp_target_Time <- bs_samp_target; bs_samp_target_Time$Time <- 1
    
    combined_bs_samp <- rbind(bs_samp_hist, bs_samp_target)
    combined_bs_samp_Time <- rbind(bs_samp_hist_Time, bs_samp_target_Time)
    
    #propensity models

    prop_mod_weight_original <<- propensity_weighting_original(bs_samp_hist, bs_samp_target) #propensity weights
    prop_mod_weight_limit_one <<- propensity_weighting_limit_one(bs_samp_hist, bs_samp_target) #propensity weights
    prop_mod_weight_scaling <<- propensity_weighting_scaling(bs_samp_hist, bs_samp_target) #propensity weights
    
    #original propensity weights
    weighted_propensity_model_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= prop_mod_weight_original, family='binomial')
    prop_results_original <- validation_results(weighted_propensity_model_original,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_original, i, 'weighted_propensity_model_original')

    weighted_propensity_model_Time_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= prop_mod_weight_original, family='binomial')
    prop_results_Time_original <- validation_results(weighted_propensity_model_Time_original,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_original, i, 'weighted_propensity_model_Time_original')
    
    #limited to 1 propensity weights
    weighted_propensity_model_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= prop_mod_weight_limit_one, family='binomial')
    prop_results_limit_one  <- validation_results(weighted_propensity_model_limit_one ,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_limit_one, i, 'weighted_propensity_model_limit_one')
    
    weighted_propensity_model_Time_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= prop_mod_weight_limit_one, family='binomial')
    prop_results_Time_limit_one <- validation_results(weighted_propensity_model_Time_limit_one,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_limit_one, i, 'weighted_propensity_model_Time_limit_one')
    
    #scaled propensity weights
    weighted_propensity_model_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= prop_mod_weight_scaling, family='binomial')
    prop_results_scaling  <- validation_results(weighted_propensity_model_scaling ,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_scaling, i, 'weighted_propensity_model_scaling')
    
    weighted_propensity_model_Time_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= prop_mod_weight_scaling, family='binomial')
    prop_results_Time_scaling <- validation_results(weighted_propensity_model_Time_scaling,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_scaling, i, 'weighted_propensity_model_Time_scaling')
    
    

    
    #mahalanobis models
    mah_mod_weight_original <- mahalanobis_weighting_original(bs_samp_hist, bs_samp_target) # mahalanobis weights
    mah_mod_weight_limit_one <- mahalanobis_weighting_limit_one(bs_samp_hist, bs_samp_target) # mahalanobis weights
    mah_mod_weight_scaling <- mahalanobis_weighting_scaling(bs_samp_hist, bs_samp_target) # mahalanobis weights
    
    #original mahalanobis weights
    weighted_mah_model_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= mah_mod_weight_original, family='binomial')
    mah_results_original <- validation_results(weighted_mah_model_original,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_original, i, 'weighted_mahalanobis_model_original')
  
    weighted_mah_model_Time_original <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= mah_mod_weight_original, family='binomial')
    mah_results_Time_original <- validation_results(weighted_mah_model_Time_original,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_original, i, 'weighted_mahalanobis_model_Time_original')
    
    #limited to 1 mahalanobis weights
    weighted_mah_model_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= mah_mod_weight_limit_one, family='binomial')
    mah_results_limit_one <- validation_results(weighted_mah_model_limit_one,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_limit_one, i, 'weighted_mahalanobis_model_limit_one')
    
    weighted_mah_model_Time_limit_one <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= mah_mod_weight_limit_one, family='binomial')
    mah_results_Time_limit_one <- validation_results(weighted_mah_model_Time_limit_one,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_limit_one, i, 'weighted_mahalanobis_model_Time_limit_one')
    
    #scaled mahalanobis weights
    weighted_mah_model_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, weights= mah_mod_weight_scaling, family='binomial')
    mah_results_scaling <- validation_results(weighted_mah_model_scaling,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_scaling, i, 'weighted_mahalanobis_model_scaling')
    
    weighted_mah_model_Time_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, weights= mah_mod_weight_scaling, family='binomial')
    mah_results_Time_scaling <- validation_results(weighted_mah_model_Time_scaling,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_scaling, i, 'weighted_mahalanobis_model_Time_scaling')
    
    
    
    #unweighted models on all data
    unweighted_model <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, family='binomial')
    unweighted_results <- validation_results(unweighted_model,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(unweighted_results, i, 'unweighted_model_on_AllData')
  
    
    unweighted_model_Time <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, family='binomial')
    unweighted_results_Time <- validation_results(unweighted_model_Time,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(unweighted_results_Time, i, 'unweighted_model_on_AllData_Time')
    
    #unweighted models on recent data
    recent_data_model <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=bs_samp_target, family='binomial')
    recent_results <- validation_results(recent_data_model,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(recent_results, i, 'unweighted_model_on_recent_data')
    
  }
  return(results)
}



boot_results<- manual_boot(hist, target, target_Time, 200)






