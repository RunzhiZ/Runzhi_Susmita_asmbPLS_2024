###############################################################################
# load R package and function
###############################################################################
library(glmnet)
library(blockForest)
library(ipflasso)
library(compositions)
library(SGL)
library(readxl)
library(dplyr)
library(survminer)
library(survival)
library(asmbPLS)
library(prioritylasso)

###############################################################################
# read real data
###############################################################################
setwd("/blue/datta/runzhi.zhang/asmbPLS/TCGA/LUSC")
load("TCGA_LUSC_gene_mirna.RData")

###############################################################################
# quantile combination setting
###############################################################################
quantile_1 <- c(0.9997, 0.9998, 0.9999, 0.99999)
quantile_2 <- c(0.99, 0.995, 0.999)
quantilelist <- list(quantile_1, quantile_2)
quantile.comb <- quantileComb(quantilelist)
quantile.comb.mbPLS <- quantileComb(list(1, 1))

PLS_term_selected <- 3
k = 5

###############################################################################
# prepare X and Y matrix for model fitting
###############################################################################
survival_time <- LUSC_clinical_TP[, c("days_to_last_follow_up", "vital_status")]
Non_NA_index <- complete.cases(survival_time$days_to_last_follow_up) & complete.cases(survival_time$vital_status) & (survival_time$days_to_last_follow_up > 0)
X_dim <- c(ncol(data_gene_common_TP), ncol(data_mirna_common_TP))
X <- cbind(data_gene_common_TP[Non_NA_index, ], data_mirna_common_TP[Non_NA_index, ])
## imputation method
survival_time_non_NA <- survival_time[Non_NA_index,]
survival_time_non_NA$vital_status <- as.numeric(survival_time_non_NA$vital_status == "Dead")
survival_time_non_NA <- survival_time_non_NA %>% mutate_at(c("days_to_last_follow_up"), as.numeric)
survival_time_non_NA$days_to_last_follow_up_log <- log(survival_time_non_NA$days_to_last_follow_up)
survival_time_imp <- meanimp(survival_time_non_NA[,c("days_to_last_follow_up_log", "vital_status")])$imputed_table
Y <- as.matrix(survival_time_imp[, 4])
Y_indicator <- survival_time_imp[, 2]

## blocks setting ##
blocks_vector <- rep(1:length(X_dim), times = X_dim)
blocks <- lapply(1:length(X_dim), function(x) which(blocks_vector == x))
blocks_1 <- lapply(1:2, function(x) which(blocks_vector == x))
blocks_2 <- lapply(2:1, function(x) which(blocks_vector == x))

#### compare function: compare predicted and true Y ####
result_compare <- function (prediction, true) {
  result_abs = sum(abs(prediction-true))/length(true)
  result_sq = sum((prediction-true)^2)/length(true)
  return(list(result_sq = result_sq,
              result_abs = result_abs))
}

###############################################################################
# comparison between different methods:
# 1) LOO, MSE_prediction

# methods include:
# 1) asmbPLS
# 2) mbPLS
# 3) BF
# 4) IPF-lasso
# 5) priority_lasso
# 6) SGL

###############################################################################
observed_index <- which(Y_indicator == 1)

results_LOO <- matrix(NA, nrow = length(observed_index), ncol = 5)
colnames(results_LOO) <- c("index", "true", "asmbPLS", 
                           "mbPLS", "IPF")
results_LOO_MSE <- matrix(NA, nrow = length(observed_index), ncol = 4)
colnames(results_LOO_MSE) <- c("index", "asmbPLS", 
                               "mbPLS", "IPF")

for(i in 1:length(observed_index)){
  print(paste0("****************************** index ", i, " ******************************"))
  #### leave-one-out ####
  i_index <- observed_index[i]
  X_test <- matrix(X[i_index, ], ncol = ncol(X))
  Y_test <- matrix(Y[i_index,], ncol = ncol(Y))
  colnames(X_test) <- colnames(X)
  X_train <- X[-i_index, ]
  Y_train <- matrix(Y[-i_index,], ncol = ncol(Y))
  Y_indicator_train <- Y_indicator[-i_index]
  
  results_LOO[i, "index"] <- i_index
  results_LOO[i, "true"] <- Y_test
  
  results_LOO_MSE[i, "index"] <- i_index
  
  ##========================= asmbPLS =========================##
  ## cv ##
  asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train,
                                   Y.matrix = Y_train,
                                   PLS.comp = PLS_term_selected,
                                   X.dim = X_dim,
                                   quantile.comb.table = quantile.comb[,1:length(X_dim)],
                                   Y.indicator = Y_indicator_train ,
                                   k = k)
  n_optimal <- asmbPLS_results_cv$optimal_nPLS
  PLS_table <- asmbPLS_results_cv$quantile_table_CV[, 1:length(X_dim)]
  
  ## fit ##
  asmbPLS_fit <- asmbPLS.fit(X.matrix = X_train, 
                             Y.matrix = Y_train,
                             PLS.comp = n_optimal, 
                             X.dim = X_dim,
                             quantile.comb = PLS_table)
  
  ## predict ##
  Y_predict_asmbPLS <- asmbPLS.predict(asmbPLS_fit, X_test, n_optimal)$Y_pred
  
  results_LOO[i, "asmbPLS"] <- Y_predict_asmbPLS
  results_LOO_MSE[i, "asmbPLS"] <- result_compare(Y_predict_asmbPLS, Y_test)[[1]]
  
  print("asmbPLS done")
  
  ##========================= mbPLS =========================##
  ## cv ##
  mbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                 Y.matrix = Y_train, 
                                 PLS.comp = PLS_term_selected, 
                                 X.dim = X_dim,
                                 quantile.comb.table = quantile.comb.mbPLS,
                                 Y.indicator = Y_indicator_train,
                                 k = k)
  n_optimal <- mbPLS_results_cv$optimal_nPLS
  
  ## fit ##
  mbPLS_fit <- mbPLS.fit(X.matrix = X_train, 
                         Y.matrix = Y_train,
                         PLS.comp = n_optimal, 
                         X.dim = X_dim)
  
  ## predict ##
  Y_predict_mbPLS <- asmbPLS.predict(mbPLS_fit, X_test, n_optimal)$Y_pred
  
  results_LOO[i, "mbPLS"] <- Y_predict_mbPLS
  results_LOO_MSE[i, "mbPLS"] <- result_compare(Y_predict_mbPLS, Y_test)[[1]]
  
  print("mbPLS done")
  
  # ##========================= Block Forest =========================##
  # Y_train_vector <- as.vector(Y_train)
  # blockforobj <- blockfor(X_train, 
  #                         Y_train_vector, 
  #                         replace = TRUE, 
  #                         blocks = blocks, 
  #                         block.method = "BlockForest", 
  #                         importance = "impurity", 
  #                         seed = i)
  # 
  # ## predict ##
  # Y_predict_BF <- predict(blockforobj$forest, data = X_test)$predictions
  # results_LOO[i, "BF"] <- Y_predict_BF
  # results_LOO_MSE[i, "BF"] <- result_compare(Y_predict_BF, Y_test)[[1]]
  
  ##========================= IPF =========================##
  set.seed(i)
  Y_train_vector <- as.vector(Y_train)
  Y_vector <- as.vector(Y)
  pflist <- list(c(1, 1), c(2, 1), c(1, 2), c(3, 1), c(1, 3), c(4, 1), c(1, 4), c(5, 1), c(1, 5))
  ipf_fit <- cvr2.ipflasso(X = X_train, 
                           Y = Y_train_vector, 
                           family = "gaussian", 
                           type.measure = "mse",
                           blocks = blocks,
                           pflist = pflist, 
                           nfolds = 5, 
                           ncv = 5)
  
  ## predict ##
  Y_predict_IPF <- ipflasso.predict(ipf_fit, X_test)$linpredtest
  results_LOO[i, "IPF"] <- Y_predict_IPF
  results_LOO_MSE[i, "IPF"] <- result_compare(Y_predict_IPF, Y_test)[[1]]
  print("IPF done")
  
  # ##========================= Prioirity Lasso =========================##
  # priority_results_cv <- cvm_prioritylasso(X = X_train, 
  #                                          Y = Y_train, 
  #                                          family = "gaussian", 
  #                                          type.measure = "mse", 
  #                                          standardize = T, 
  #                                          blocks.list = list(blocks_1, blocks_2), 
  #                                          block1.penalization = TRUE, 
  #                                          lambda.type = "lambda.1se")
  # priority_fit <- priority_results_cv$best.model
  # Y_predict_priority_lasso <- predict(priority_fit, X_test, type = "response")
  # 
  # ## predict ##
  # Y_predict_priority_lasso <- predict(priority_fit, X_test, type = "response")
  # results_LOO[i, "priority_lasso"] <- Y_predict_priority_lasso
  # results_LOO_MSE[i, "priority_lasso"] <- result_compare(Y_predict_priority_lasso, Y_test)[[1]]
  # print("priority_lasso done")
  
  # ##========================= SGL =========================##
  # data_SGL = list(x = X_train, y = Y_train)
  # index_SGL <- c(rep(1, X_dim[1]), rep(2, X_dim[2]))
  # cvFit = cvSGL(data_SGL, index_SGL, type = "linear", nfold = 5)
  # Fit = SGL(data_SGL, index_SGL, type = "linear", lambdas = cvFit$lambdas)
  # Y_predict_SGL <- predictSGL(Fit, X_test, which.min(cvFit$lldiff))
}

###############################################################################
# comparison between different methods:
# 1) Fit
# 2) Feature Selection

# methods include:
# 1) asmbPLS
# 2) IPF-lasso

###############################################################################
FS_table_1 <- matrix(NA, nrow = X_dim[1], ncol = 2)
FS_table_2 <- matrix(NA, nrow = X_dim[2], ncol = 2)
row.names(FS_table_1) <- colnames(X)[1:X_dim[1]]
row.names(FS_table_2) <- colnames(X)[(X_dim[1]+1):(sum(X_dim))]
colnames(FS_table_1) <- colnames(FS_table_2) <- c("asmbPLS", "IPF")

asmbPLS_results_cv_FS <- asmbPLS.cv(X.matrix = X,
                                    Y.matrix = Y,
                                    PLS.comp = PLS_term_selected,
                                    X.dim = X_dim,
                                    quantile.comb.table = quantile.comb[,1:length(X_dim)],
                                    Y.indicator = Y_indicator,
                                    k = k)

n_optimal_FS <- asmbPLS_results_cv_FS$optimal_nPLS
PLS_table_FS <- asmbPLS_results_cv_FS$quantile_table_CV[, 1:length(X_dim)]

## fit ##
asmbPLS_fit_FS <- asmbPLS.fit(X.matrix = X, 
                              Y.matrix = Y,
                              PLS.comp = n_optimal_FS, 
                              X.dim = X_dim,
                              quantile.comb = PLS_table_FS)

FS_table_1[, "asmbPLS"] <- asmbPLS_fit_FS$X_weight[[1]]
FS_table_2[, "asmbPLS"] <- asmbPLS_fit_FS$X_weight[[2]]

#### IPF ####
set.seed(123)
Y_vector <- as.vector(Y)
pflist <- list(c(1, 1), c(2, 1), c(1, 2), c(3, 1), c(1, 3), c(4, 1), c(1, 4), c(5, 1), c(1, 5))
ipf_fit_FS <- cvr2.ipflasso(X = X, Y = Y_vector, family = "gaussian", type.measure = "mse",
                            blocks = blocks,
                            pflist = pflist, nfolds = 5, ncv = 5)
FS_table_1[, "IPF"] <- ipf_fit_FS$coeff[,ipf_fit_FS$ind.bestlambda][-1][blocks[[1]]]
FS_table_2[, "IPF"] <- ipf_fit_FS$coeff[,ipf_fit_FS$ind.bestlambda][-1][blocks[[2]]]

save(list = c('results_LOO', 
              'results_LOO_MSE',
              "FS_table_1",
              "FS_table_2",
              "asmbPLS_fit_FS",
              "ipf_fit_FS"), file = "Results_TCGA_LOO_FS.RData")