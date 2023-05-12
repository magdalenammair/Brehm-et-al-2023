# Aim: optimize predictions from xgboost by nested CV
# plus: use optimized model for predictive performance and inferences
library(dplyr) # for pipes
library(xgboost) #for xgboost
library(Matrix) #for sparse matrix
library(utiml) #for mldr

dat <- readRDS("Data/imputed_dat.rds")
str(dat)

# prepare response:--------
response = dat$yi

# prepare predictors-------
# remove response
preds = dat %>% select(-yi)
str(preds)

# one-hot-encoding------
# prepare sparse matrix: -1 = without intercept!
preds_matrix = sparse.model.matrix(~ .-1, data = preds)
colnames(preds_matrix)
# Have a look at the final dataset:
matrix(preds_matrix, ncol =65, byrow = F)[10:50,1:3]

# Tuning-------
## define CV parameters------
outer_n = 5 # number of outer CV steps
inner_n = 5 # number of inner CV steps
tune_xgb = 50 # number of tuning parameter sets to be tested in hyper parameter tuning

## evaluation functions --------
# for calculating RMSE; o: observed values, p: predicted values
rmse = function(o, p) {sqrt(mean((o-p)**2))}
#define function that returns RMSE values and spearman's correlation coefficient:
get_stats = function(obs, pred) { 
  return(c(rmse(obs, pred), cor(obs, pred, method = "spearman"))) 
}

## CV split setup --------

n = 1:nrow(dat) # number of samples
outer_mldr = mldr_from_dataframe(dat, labelIndices = c(length(dat),length(dat))) # labelIndices = yi values
outer_folds = utiml::create_kfold_partition(outer_mldr, k =  outer_n ,method = "stratified") 
outer = lapply(outer_folds$fold, function(f) list(train = n[-f], test = n[f])) # samples in train and test datasets for each of the outer_n CV splits
inner = # now for each of the outer splits, create an inner CV split
  lapply(1:outer_n, function(i) {
    inner_mldr = mldr_from_dataframe(dat[outer[[i]]$train,], labelIndices = c(length(dat),length(dat)))
    inner = utiml::create_kfold_partition(inner_mldr, k =  inner_n ,method = "stratified") 
    return(lapply(inner$fold, function(f) list(train = outer[[i]]$train[-f], test = outer[[i]]$train[f])))
  })


## training function--------

### xgboost tuning

# define function that trains BRT on training dataset and returns predictions on test dataset:
fit_xgb = function(Y_train, P_train, P_test, ...) { #Y_train = response, P_train = features, P_test = features of test dataset from CV splits
  r = xgboost::xgboost(data = P_train, label = Y_train, ...)
  return(predict(r, newdata = P_test))
}

## tuning parameters--------
# Define tuning parameter combinations from which random samples will be taken for hyper parameter tuning
max_depth = sample.int(20, 14)+1 
eta = seq(0.1,1,0.1) 
n.rounds = sample.int(20, 14)+1 
gamma = c(0, 0.1, 0.2) 
tuning_pars = expand.grid(max_depth, eta, n.rounds, gamma)
# sample from parameter grid
tune = tuning_pars[sample.int(nrow(tuning_pars), tune_xgb),]
# change column names
colnames(tune) = c( "max_depth", "eta", "n.rounds", "gamma")

## output data frame------
results = data.frame(matrix(NA, outer_n, 6))
inner_preds = vector("list", outer_n)
colnames(results) = c(colnames(tune), "rmse", "spear")

## run nested CV------
for(i in 1:outer_n) {
  
  results_inner = 
    lapply(1:inner_n, function(j) {
      train = inner[[i]][[j]]$train
      test = inner[[i]][[j]]$test
      train_feat = preds_matrix[train,]
      test_feat = preds_matrix[test,] 
      cl = parallel::makeCluster(4L)
      parallel::clusterExport(cl, varlist = list("response", "preds_matrix", "tune", "get_stats", "train", "test", "train_feat", "test_feat", "rmse", "fit_xgb"), 
                              envir = environment())
      results_inner = 
        parSapply(cl, 1:nrow(tune), function(k) {
          tt = tune[k,]
          pred = fit_xgb(response[train], train_feat, test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
          res = get_stats(response[test], pred) #rmse and cor
          gc()
          return(res) 
        })
      parallel::stopCluster(cl)
      return(t(results_inner))
    })
  
  mean_results_inner = apply(abind::abind(results_inner, along = 0L), 2:3, mean) #
  #inner_preds[[i]] = cbind(tune, data.frame(ll = abind::abind(results_inner, along = 0L)[1,,2]), pred_set = i)
  best_par = which.min(mean_results_inner[,1])
  
  ### train with best on outer
  results_outer = 
    sapply(1:outer_n, function(k) {
      train_outer = outer[[k]]$train
      test_outer = outer[[k]]$test
      outer_train_feat = preds_matrix[train_outer,]
      outer_test_feat = preds_matrix[test_outer,]
      tt = tune[best_par, ]
      pred = fit_xgb(response[train_outer], outer_train_feat, outer_test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
      res = get_stats(response[test_outer], pred)
      return(res)
    })
  results_outer = t(results_outer)
  res = cbind(tune[best_par, ], matrix(apply(results_outer, 2, mean), nrow = 1L))
  results[i,] = res
}

results

## save results -------
saveRDS(results, "Results/xgb_tuning.RDS")


# Best model training -----
# Train model with best parameter combination from outer CV on full dataset:
# results = readRDS("Results/xgb_tuning.RDS")
final_par = results[which.min(results[,5]),]
# train best model
best <- xgboost(data = preds_matrix, label = response, max_depth = final_par$max_depth,
               eta = final_par$eta, gamma = final_par$gamma, nrounds = final_par$n.rounds)

# get R^2:
pred = predict(best, newdata = preds_matrix)
plot(response, pred)
abline(0,1)
sstot <- sum((pred - mean(response))^2)
ssresid <- sum((pred-response)^2)
boxplot((pred - mean(response))^2, (pred-response)^2)
100 * (1 - ssresid / sstot)

# get feature importances
importance <- xgb.importance(feature_names = NULL, model = best)

# create column with grouping variable "predictor"
predlist = colnames(dat)
for(i in predlist) {
  importance$predictor[grepl(i,importance$Feature)] = i
}

# combine different uncontrolled surface modifications into one predictor group:
# DOC_added, biotic_env and other_added := Uncontrolled SM yes/no
importance$predictor[importance$predictor %in% c("DOC_added", "biotic_env", "other_added")] = "USM"

# sum up gain, cover, frequency and importance per predictor group
pred_importance = importance %>%  group_by(predictor) %>% 
  summarise(gain = sum(Gain), cover = sum(Cover), frequency = sum(Frequency))
# order: decreasing
pred_importance = pred_importance[order(pred_importance$gain, decreasing = T),]

# prepare variable groups to assign colors for plotting later (according to mindmap in Figure 1 in Brehm et al. 2023):
grouping = data.frame(
  predictor = pred_importance$predictor,
  linegroup = c(2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2), #1 = yellow, 2 = blue, 3 = green
  boxgroup = c(2, 3, 3, 2, 2, 2, 2, 3, 1, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2)
)
# attach color codes to importances dataframe:
pred_importance = merge(pred_importance, grouping, by = "predictor", sort = FALSE)
# define labels
pred_importance$label = c("NMP modification type", "Concentration mg/ml", "Concentration particles/ml",  
                          "NMP size", "NMP modification yes/no", "Density (buoyancy)","Polymer type",
                          "Temperature", "Age", "Food", "NMP shape", "Exposure duration", "Uncontrolled SM (biotic)",
                          "Type of DOC during exposure", "MP color", "Fluorescence","Uncontrolled SM (UV weathered)",
                          "NMP color yes/no", "Surfactant", "Species", "Surface charge")

# save importances
saveRDS(list(importance = importance, grouped_importances = pred_importance), file = "Results/importances.RDS")


# Baseline model---------

# shuffle response values:
base_response = sample(response, length(response), replace = FALSE)

# Tune baseline model
## output data frame------
base_results = data.frame(matrix(NA, outer_n, 6))
inner_preds = vector("list", outer_n)
colnames(base_results) = c(colnames(tune), "rmse", "spear")

## run nested CV------
for(i in 1:outer_n) {
  results_inner = 
    lapply(1:inner_n, function(j) {
      train = inner[[i]][[j]]$train
      test = inner[[i]][[j]]$test
      train_feat = preds_matrix[train,]
      test_feat = preds_matrix[test,] 
      cl = parallel::makeCluster(4L)
      parallel::clusterExport(cl, varlist = list("base_response", "preds_matrix", "tune", "get_stats", "train", "test", "train_feat", "test_feat", "rmse", "fit_xgb"), 
                              envir = environment())
      results_inner = 
        parSapply(cl, 1:nrow(tune), function(k) {
          tt = tune[k,]
          pred = fit_xgb(base_response[train], train_feat, test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
          res = get_stats(base_response[test], pred) 
          gc()
          return(res) 
        })
      parallel::stopCluster(cl)
      return(t(results_inner))
    })
  
  mean_results_inner = apply(abind::abind(results_inner, along = 0L), 2:3, mean) 
  best_par = which.min(mean_results_inner[,1])
  
  ### train with best on outer
  results_outer = 
    sapply(1:outer_n, function(k) {
      train_outer = outer[[k]]$train
      test_outer = outer[[k]]$test
      outer_train_feat = preds_matrix[train_outer,]
      outer_test_feat = preds_matrix[test_outer,]
      tt = tune[best_par, ]
      pred = fit_xgb(base_response[train_outer], outer_train_feat, outer_test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
      res = get_stats(base_response[test_outer], pred)
      return(res)
    })
  results_outer = t(results_outer)
  res = cbind(tune[best_par, ], matrix(apply(results_outer, 2, mean), nrow = 1L))
  base_results[i,] = res
}

base_results


## save results -------
saveRDS(base_results, "Results/base_xgb_tuning.RDS")

## Best base model training -----
# Train model with best parameter combination from outer CV on full dataset:

base_final_par = base_results[which.min(base_results[,5]),]
base_best <- xgboost(data = preds_matrix, label = base_response, max_depth = base_final_par$max_depth,
                eta = base_final_par$eta, gamma = base_final_par$gamma, nrounds = base_final_par$n.rounds)

# get R^2:
base_pred = predict(base_best, newdata = preds_matrix)
plot(base_response, base_pred)
abline(0,1)
base_sstot <- sum((base_pred - mean(base_response))^2)
base_ssresid <- sum((base_pred-base_response)^2)
boxplot((base_pred - mean(base_response))^2, (base_pred-base_response)^2)
100 * (1 - base_ssresid / base_sstot)

# Comparison of true data model with baseline model---------
## define 10fold CV data splits----
outer_n2 = 10
n = 1:nrow(dat) #number of samples
outer_mldr = mldr_from_dataframe(dat, labelIndices = c(length(dat),length(dat))) #labelIndices = yi values here
outer_folds = utiml::create_kfold_partition(outer_mldr, k =  outer_n2 ,method = "stratified") 
outer = lapply(outer_folds$fold, function(f) list(train = n[-f], test = n[f])) 

# run CV to derive predictive performance (rmse, spearman rank correlation) on the 10 different test datasets using defined data splits:

## True data model-----
## output data frame
performance_results = data.frame(matrix(NA, outer_n, 3))
colnames(performance_results) = c("rmse", "spear", "which.model")

## run 
for(i in 1:outer_n2) {
  train_outer = outer[[i]]$train
  test_outer = outer[[i]]$test
  outer_train_feat = preds_matrix[train_outer,]
  outer_test_feat = preds_matrix[test_outer,]
  tt = final_par
  pred = fit_xgb(response[train_outer], outer_train_feat, outer_test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
  res = get_stats(response[test_outer], pred)
  performance_results[i,] = c(res, "data")
}
performance_results

## baseline model------
## base output data frame
base_performance_results = data.frame(matrix(NA, outer_n, 3))
colnames(base_performance_results) = c("rmse", "spear", "which.model")

## base run nested CV
for(i in 1:outer_n2) {
  train_outer = outer[[i]]$train
  test_outer = outer[[i]]$test
  outer_train_feat = preds_matrix[train_outer,]
  outer_test_feat = preds_matrix[test_outer,]
  tt = base_final_par
  pred = fit_xgb(base_response[train_outer], outer_train_feat, outer_test_feat, max_depth = tt$max_depth, eta = tt$eta, nrounds = tt$n.rounds, gamma = tt$gamma)
  res = get_stats(base_response[test_outer], pred)
  #results_outer = t(results_outer)
  base_performance_results[i,] = c(res, "baseline")
}
base_performance_results

(all_performance = rbind(performance_results, base_performance_results))

all_performance$rmse = as.numeric(all_performance$rmse)
all_performance$which.model = as.factor(all_performance$which.model)

## save performance data-----
saveRDS(all_performance, file = "Results/performance_xgboost.RDS")

## END ---------

