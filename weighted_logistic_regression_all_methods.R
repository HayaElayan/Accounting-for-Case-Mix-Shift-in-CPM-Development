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
library(PRROC)

#View(imputed_df)
hist <- #data before the shift (historical)

target <- #data after the shift (target population)
    
hist_Time <- hist; hist_Time$Time <- 0 #add dummy variable for the target population membership as feature
target_Time <- target; target_Time$Time <- 1  #add dummy variable for the target population membership as feature

  


propensity_weighting_original <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist$set <- 0
  target$set <- 1
  all<- rbind(hist, target)
  
  membership_mod <- glm(set~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=all, family = 'binomial')
  
  preds <- predict(membership_mod,type="response", newdata = hist) #predict p(set=1|X) target
  
  propensity_weight <- preds / (1 - preds) #p(p(set=1|X)/ p(set=0|X)) target/source
  propensity_weight <- propensity_weight * (nrow(hist)/ nrow(target)) #target/source * s/n
  
  model_weights <- c(propensity_weight, rep(1, nrow(target)))
  
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
  propensity_weight <- propensity_weight * (nrow(hist)/ nrow(target))
  
  propensity_weight[propensity_weight > 1]= 1
  
  
  model_weights <- c(propensity_weight, rep(1, nrow(target)))
  
  
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
  propensity_weight <- propensity_weight * (nrow(hist)/ nrow(target))
  
  propensity_weight_scaled <- 0
  if(sum(propensity_weight) > nrow(hist)){
    print('scaled')
    scale_factor <- nrow(hist)/sum(propensity_weight)
    propensity_weight_scaled <- propensity_weight * scale_factor
  }else{
    propensity_weight_scaled <- propensity_weight
  }
  
  model_weights <- c(propensity_weight, rep(1, nrow(target)))
  
  
  return(model_weights)
  
}
#########################
#mahalanobis distance model weighting functions; original weighting, weights limited to one, 
#and weights scaled with scaling factor

mahalanobis_weighting_original <- function(hist, target){
  tryCatch({
    drops <- c('outcome')
    hist <- hist[ , !(names(hist) %in% drops)]
    target <- target[ , !(names(target) %in% drops)]
    
    hist_num <- data.frame(lapply(hist,as.numeric))
    target_num <- data.frame(lapply(target,as.numeric))
    
    
    #Mahalanobis distance weights
    distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num), tol=1e-20) #get mahalanobis distances
    hist_num$mah_distances <- distances
    
    fx <-  replicate(nrow(hist), 1) #fx uniform distribution
    pvalues <- pchisq(distances, df=7, lower.tail=FALSE) #get pvalues of distances from chi-square distribution
    
    gx<- density(x=pvalues) #estimate gx
    gx_prob <- approxfun(gx$x, gx$y) #get gx probabilities
    
    hist_weights <- fx/gx_prob(pvalues)
    target_weights <- replicate(nrow(target), 1)
    
    model_weights <- c(hist_weights, target_weights)
    
    return(model_weights)
  }, error = function(e) {
    print("mahalanobis_weighting_original")
  }
  )
}



mahalanobis_weighting_limit_one <- function(hist, target){
  drops <- c('outcome')
  hist <- hist[ , !(names(hist) %in% drops)]
  target <- target[ , !(names(target) %in% drops)]
  
  hist_num <- data.frame(lapply(hist,as.numeric))
  target_num <- data.frame(lapply(target,as.numeric))
  
  
  #Mahalanobis distance weights
  distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num), tol=1e-20) #get mahalanobis distanes
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
  distances <- mahalanobis(hist_num, colMeans(target_num), cov(target_num), tol=1e-20) #get mahalanobis distanes
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



apparent_results_fun <- function(model, original_samp, model_name){
  val_results <- matrix(nrow = 1,ncol = 6)
  
  pr_test <- predict(model,type="response", newdata = original_samp) # predict probabilities from the bootstrap model in the original sample
  lp_test <- predict(model, newdata = original_samp) # predict lp from the bootstrap model in the original sample
  
  
  # calculate the apparent performance of the model in the original sample
  test_cstat_model <- roc(outcome~pr_test,data=original_samp)
  val_results[1,1] <- test_cstat_model$auc
  
  labels <- original_samp$outcome
  test_PRAUC_model <- pr.curve(scores.class0 = pr_test[labels == 1], scores.class1 = pr_test[labels == 0], curve = TRUE)
  val_results[1,2] <- test_PRAUC_model$auc.integral
  
  test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=original_samp)
  val_results[1,3] <- summary(test_citl_model)$coefficients[1,1]
  test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=original_samp)
  val_results[1,4] <- summary(test_cslope_model)$coefficients[2,1]
  val_results[1,5] <- mean((original_samp$outcome - pr_test)^2)
  val_results[1,6] <- model_name
  
  return(val_results)
}



validation_results <- function(model, boot_samp, original_samp){
  val_results <- matrix(nrow = 1,ncol = 10)
  
  pr_bs <- predict(model,type="response", newdata = boot_samp) # predict probabilities from the bootstrap model in the bs sample
  lp_bs <- predict(model,newdata = boot_samp ) # predict lp from the bootstrap model in the bs sample
  
  pr_test <- predict(model,type="response", newdata = original_samp) # predict probabilities from the bootstrap model in the original sample
  lp_test <- predict(model, newdata = original_samp) # predict lp from the bootstrap model in the original sample
  
  
  # calculate the apparent performance of the bootstrap model in the bs sample
  app_cstat_model <- roc(outcome~pr_bs,data=boot_samp)
  val_results[1,1] <- app_cstat_model$auc
  
  
  labels_boost <- boot_samp$outcome
  app_PRAUC_model <- pr.curve(scores.class0 = pr_bs[labels_boost == 1], scores.class1 = pr_bs[labels_boost == 0], curve = TRUE)
  val_results[1,2] <- app_PRAUC_model$auc.integral
  
  app_citl_model <- glm(outcome ~ offset(lp_bs),family=binomial, data=boot_samp)
  val_results[1,3] <- summary(app_citl_model)$coefficients[1,1]
  app_cslope_model <- glm(outcome ~ lp_bs,family=binomial(link='logit'), data=boot_samp)
  val_results[1,4] <- summary(app_cslope_model)$coefficients[2,1]
  
  val_results[1,5] <- mean((boot_samp$outcome - pr_bs)^2)
  
  # calculate the test performance of the bootstrap model in the original sample
  test_cstat_model <- roc(outcome~pr_test,data=original_samp)
  val_results[1,6] <- test_cstat_model$auc
  
  labels_org <- original_samp$outcome
  test_PRAUC_model <- pr.curve(scores.class0 = pr_test[labels_org == 1], scores.class1 = pr_test[labels_org == 0], curve = TRUE)
  val_results[1,7] <- test_PRAUC_model$auc.integral
  
  test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=original_samp)
  val_results[1,8] <- summary(test_citl_model)$coefficients[1,1]
  test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=original_samp)
  val_results[1,9] <- summary(test_cslope_model)$coefficients[2,1]
  val_results[1,10] <- mean((original_samp$outcome - pr_test)^2)
  
  return(val_results)
  
}


manual_boot <- function(hist_set, target_set, target_set_Time,samples){
  results <- data.frame(matrix(ncol = 12, nrow = 0))
  names <- c("app_c_stat", 'app_PRAUC',"app_citl","app_c_slope",'app_BrierScore',
             "test_c_stat",'test_PRAUC',"test_citl","test_c_slope",
             "test_BrierScore","booting_round","model_name")
  colnames(results) <- names
  
  stability_df <- data.frame()
  
  
  for (i in 1:samples) {
    set.seed(23139*i)
    samp_index_hist <- sample(1:nrow(hist_set), nrow(hist_set), rep=TRUE) # create a sampling index vector
    samp_index_target <- sample(1:nrow(target_set), nrow(target_set), rep=TRUE) # create a sampling index vector
    
    bs_samp_hist <- hist_set[samp_index_hist,] # index the orignal dataset using the sampling vector to give the bs sample
    bs_samp_target <- target_set[samp_index_target,]
    
    bs_samp_hist_Time <- bs_samp_hist; bs_samp_hist_Time$Time <- 0
    bs_samp_target_Time <- bs_samp_target; bs_samp_target_Time$Time <- 1
    
    combined_bs_samp <- rbind(bs_samp_hist, bs_samp_target)
    combined_bs_samp_Time <- rbind(bs_samp_hist_Time, bs_samp_target_Time)
    
    #propensity models
    boot_propensity_weight_original <<- propensity_weighting_original(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    boot_propensity_weight_limit_one <<- propensity_weighting_limit_one(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    boot_propensity_weight_scaling <<- propensity_weighting_scaling(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    
    
    propensity_model_NoRecalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, 
                                                           data=combined_bs_samp,
                                                           weights= boot_propensity_weight_original, family='binomial')
    
    prop_results_or <- validation_results(propensity_model_NoRecalibrate_original_weights,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_or, i, 'propensity_model_NoRecalibrate_original_weights')
    
    
    propensity_model_Recalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                                         data=combined_bs_samp_Time,
                                                         weights= boot_propensity_weight_original, family='binomial')
    
    prop_results_Time_or <- validation_results(propensity_model_Recalibrate_original_weights,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_or, i, 'propensity_model_Recalibrate_original_weights')
    
    #####
    
    propensity_model_NoRecalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, 
                                                 data=combined_bs_samp,
                                                 weights= boot_propensity_weight_limit_one, family='binomial')
    
    prop_results_lim <- validation_results(propensity_model_NoRecalibrate_limit1,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_lim, i, 'propensity_model_NoRecalibrate_limit1')
    
    
    propensity_model_Recalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                               data=combined_bs_samp_Time,
                                               weights= boot_propensity_weight_limit_one, family='binomial')
    
    prop_results_Time_lim <- validation_results(propensity_model_Recalibrate_limit1,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_lim, i, 'propensity_model_Recalibrate_limit1')
    
    
    #####
    
    propensity_model_NoRecalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, 
                                                  data=combined_bs_samp,
                                                  weights= boot_propensity_weight_scaling, family='binomial')
    
    prop_results_scale <- validation_results(propensity_model_NoRecalibrate_scaling,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(prop_results_scale, i, 'propensity_model_NoRecalibrate_scaling')
    
    
    propensity_model_Recalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                                data=combined_bs_samp_Time,
                                                weights= boot_propensity_weight_scaling, family='binomial')
    
    prop_results_Time_scale <- validation_results(propensity_model_Recalibrate_scaling,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(prop_results_Time_scale, i, 'propensity_model_Recalibrate_scaling')
    
    
    #mahalanobis models
    
    boost_mah_weight_original <<- mahalanobis_weighting_original(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    boost_mah_weight_limit_one <<- mahalanobis_weighting_limit_one(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    boost_mah_weight_scaling <<- mahalanobis_weighting_scaling(bs_samp_hist, bs_samp_target) #propensity_weighting or mahalanobis_weighting
    
    
    mahalanobis_model_NoRecalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR,
                                                            data=combined_bs_samp,
                                                            weights= boost_mah_weight_original, family='binomial')
    
    mah_results_or <- validation_results(mahalanobis_model_NoRecalibrate_original_weights,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_or, i, 'mahalanobis_model_NoRecalibrate_original_weights')
    
    mahalanobis_model_Recalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time,
                                                          data=combined_bs_samp_Time,
                                                          weights= boost_mah_weight_original, family='binomial')
    
    mah_results_Time_or <- validation_results(mahalanobis_model_Recalibrate_original_weights,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_or, i, 'mahalanobis_model_Recalibrate_original_weights')
    
    #################
    mahalanobis_model_NoRecalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR,
                                                  data=combined_bs_samp,
                                                  weights= boost_mah_weight_limit_one, family='binomial')
    
    mah_results_lim <- validation_results(mahalanobis_model_NoRecalibrate_limit1,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_lim, i, 'mahalanobis_model_NoRecalibrate_limit1')
    
    mahalanobis_model_Recalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time,
                                                data=combined_bs_samp_Time,
                                                weights= boost_mah_weight_limit_one, family='binomial')
    
    mah_results_Time_lim <- validation_results(mahalanobis_model_Recalibrate_limit1,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_lim, i, 'mahalanobis_model_Recalibrate_limit1')
    
    ##############
    
    mahalanobis_model_NoRecalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR,
                                                   data=combined_bs_samp,
                                                   weights= boost_mah_weight_scaling, family='binomial')
    
    mah_results_scale <- validation_results(mahalanobis_model_NoRecalibrate_scaling,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(mah_results_scale, i, 'mahalanobis_model_NoRecalibrate_scaling')
    
    mahalanobis_model_Recalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time,
                                                 data=combined_bs_samp_Time,
                                                 weights= boost_mah_weight_scaling, family='binomial')
    
    mah_results_Time_scale <- validation_results(mahalanobis_model_Recalibrate_scaling,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(mah_results_Time_scale, i, 'mahalanobis_model_Recalibrate_scaling')
    
    #unweighted models
    unweighted_model_allData_NoRecalibrate <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=combined_bs_samp, family='binomial')
    unweighted_results <- validation_results(unweighted_model_allData_NoRecalibrate,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(unweighted_results, i, 'unweighted_model_allData_NoRecalibrate')
    
    
    unweighted_model_allData_Recalibrate <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=combined_bs_samp_Time, family='binomial')
    unweighted_results_Time <- validation_results(unweighted_model_allData_Recalibrate,bs_samp_target_Time, target_set_Time)
    results[nrow(results) + 1,] <- cbind(unweighted_results_Time, i, 'unweighted_model_allData_Recalibrate')
    
    #recent model
    unweighted_model_target_only <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=bs_samp_target, family='binomial')
    recent_results <- validation_results(unweighted_model_target_only,bs_samp_target, target_set)
    results[nrow(results) + 1,] <- cbind(recent_results, i, 'unweighted_model_target_only')
  }
  return(results)
}




#WEIGHTS
#propensity models

propensity_weight_original <- propensity_weighting_original(hist, target) #propensity_weighting or mahalanobis_weighting
propensity_weight_limit_one <- propensity_weighting_limit_one(hist, target) #propensity_weighting or mahalanobis_weighting
propensity_weight_scaling <- propensity_weighting_scaling(hist, target) #propensity_weighting or mahalanobis_weighting


#mahalanobis models

mahalanobis_weight_original <- mahalanobis_weighting_original(hist, target) #propensity_weighting or mahalanobis_weighting
mahalanobis_weight_limit_one <- mahalanobis_weighting_limit_one(hist, target) #propensity_weighting or mahalanobis_weighting
mahalanobis_weight_scaling <- mahalanobis_weighting_scaling(hist, target) #propensity_weighting or mahalanobis_weighting


  
#MODELS
#PROPENSITY_original weights
propensity_model_NoRecalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                                         weights= propensity_weight_original, family='binomial')

propensity_model_Recalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                                       data=rbind(hist_Time, target_Time), weights= propensity_weight_original, family='binomial')
  

  #PROPENSITY limit one weights
propensity_model_NoRecalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                               weights= propensity_weight_limit_one, family='binomial')

  
propensity_model_Recalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                             data=rbind(hist_Time, target_Time), weights= propensity_weight_limit_one, family='binomial')
  

  #
  #PROPENSITY scaling factor weights
propensity_model_NoRecalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                                weights= propensity_weight_scaling, family='binomial')

propensity_model_Recalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, 
                                              data=rbind(hist_Time, target_Time), weights= propensity_weight_scaling, family='binomial')

  
  
  
#DISTANCE
#original weights
mahalanobis_model_NoRecalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                                          weights= mahalanobis_weight_original, family='binomial')

  
mahalanobis_model_Recalibrate_original_weights <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time),
                                                        weights= mahalanobis_weight_original, family='binomial')

  
  #
#limit one weights
mahalanobis_model_NoRecalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                                weights= mahalanobis_weight_limit_one, family='binomial')

mahalanobis_model_Recalibrate_limit1 <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time),
                                              weights= mahalanobis_weight_limit_one, family='binomial')
  
  
#scaling weights
mahalanobis_model_NoRecalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target),
                                                 weights= mahalanobis_weight_limit_one, family='binomial')

  
  
mahalanobis_model_Recalibrate_scaling <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time),
                                               weights= mahalanobis_weight_limit_one, family='binomial')

  
  
#unweighted models
unweighted_model_allData_NoRecalibrate <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=rbind(hist, target), family='binomial')

  
  
unweighted_model_allData_Recalibrate <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR+Time, data=rbind(hist_Time, target_Time), family='binomial')

  
#Recent model
unweighted_model_target_only <- glm(outcome~Age+Sex+AF_atrial_flutter+Diabetes+BMI+LVEF35+eGFR, data=target, family='binomial')

#get the apparent validation results
apprent_results <- data.frame(matrix(ncol = 6, nrow = 0))
names <- c("app_c_stat","app_PRAUC","app_citl","app_c_slope","BrierScore","model_name")
colnames(apprent_results) <- names

#prop
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_NoRecalibrate_original_weights, target, 'propensity_model_NoRecalibrate_original_weights')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Recalibrate_original_weights, target_Time, 'propensity_model_Recalibrate_original_weights')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_NoRecalibrate_limit1, target,'propensity_model_NoRecalibrate_limit1')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Recalibrate_limit1, target_Time, 'propensity_model_Recalibrate_limit1')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_NoRecalibrate_scaling,  target, 'propensity_model_NoRecalibrate_scaling')    
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(propensity_model_Recalibrate_scaling, target_Time, 'propensity_model_Recalibrate_scaling')

# #mah  
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_NoRecalibrate_original_weights, target, 'mahalanobis_model_NoRecalibrate_original_weights')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Recalibrate_original_weights,  target_Time, 'mahalanobis_model_Recalibrate_original_weights')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_NoRecalibrate_limit1,target, 'mahalanobis_model_NoRecalibrate_limit1')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Recalibrate_limit1, target_Time, 'mahalanobis_model_Recalibrate_limit1')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_NoRecalibrate_scaling, target, 'mahalanobis_model_NoRecalibrate_scaling')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(mahalanobis_model_Recalibrate_scaling, target_Time,'mahalanobis_model_Recalibrate_scaling')

#All data  
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(unweighted_model_allData_NoRecalibrate, target, 'unweighted_model_allData_NoRecalibrate')
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(unweighted_model_allData_Recalibrate, target_Time, 'unweighted_model_allData_Recalibrate')       

#target only  
apprent_results[nrow(apprent_results) + 1,] <- apparent_results_fun(unweighted_model_target_only,target, 'unweighted_model_target_only') 


  ########################
#bootstrap Validation

boot_results<- manual_boot(hist, target, target_Time, 200)



