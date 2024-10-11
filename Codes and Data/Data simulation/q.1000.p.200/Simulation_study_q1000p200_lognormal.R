rm(list = ls()) # clears the global environment
arr_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(arr_id)

library(asmbPLS)
library(blockForest)
library(ipflasso)
library(SGL)
library(prioritylasso)
library(compositions)

setwd("/blue/datta/runzhi.zhang/asmbPLS/Sim202312/")
load("Parameter_list.RData")
Parameter_lognormal <- Parameter_comb[Parameter_comb$Dist == "Lognormal" & Parameter_comb$Censoring_rate != 0,]

##==========================##
##    Results Comparison
##==========================##
result_compare <- function (prediction, true) {
  result_abs = round(sum(abs(prediction-true))/length(true), 3)
  result_sq = round(sum((prediction-true)^2)/length(true), 3)
  return(list(result_sq = result_sq,
              result_abs = result_abs))
}


###############################################################################
# parameters setting
###############################################################################
Censoring_rate = Parameter_lognormal$Censoring_rate[arr_id]
Noise_level = Parameter_lognormal$Noise_level[arr_id]
Effect_scale = Parameter_lognormal$Effect_scale[arr_id]

c <- c("00", "01", "03", "05", "07")[which(c(0, 0.1, 0.3, 0.5, 0.7) == Censoring_rate)]
r <- c("00", "01", "02", "05", "10")[which(c(0, 0.1, 0.2, 0.5, 1) == Noise_level)]
e <- c("05", "20")[which(c(0.5, 2) == Effect_scale)]

niter = 100
#### q = 1000 and p = 200 ####
q = 1000
p = 200
k = 5
PLS_term_selected = 3

print(paste0("**** q = ", q, "; p = ", p, "; c = ", Censoring_rate, "; r = ", Noise_level, "; e = ", Effect_scale, " ****"))

X_dim = c(q, p)
##==== block setting ====##
blocks_vector <- rep(1:length(X_dim), times = X_dim)
blocks <- lapply(1:2, function(x) which(blocks_vector == x))
blocks_1 <- lapply(1:2, function(x) which(blocks_vector == x))
blocks_2 <- lapply(2:1, function(x) which(blocks_vector == x))

##==== lambda combination ====##
quantile_1 <- c(0.9, 0.95, 0.975, 0.99, 0.995)
quantile_2 <- c(0.7, 0.8, 0.9, 0.95, 0.975)
quantilelist <- list(quantile_1, quantile_2)
quantile.comb <- quantileComb(quantilelist)

quantile.comb.mbPLS <- quantileComb(list(1, 1))


setwd(paste0("/blue/datta/runzhi.zhang/asmbPLS/Sim202312/q.", q, ".p.", p, "/Data_sim/"))
file_name <- paste0("SimData_c", c, "_r", r, "_e", e, ".RData")
load(file_name)

train_index <- 1:100
test_index <- 101:200

results_FS_list <- list()
predict_table <- fit_table <- time_table <- matrix(NA, nrow = 600, ncol = 7)
colnames(predict_table) <- colnames(fit_table) <- colnames(time_table) <- c("beta", "asmbPLS", "mbPLS", "BF", "IPF", "priority_lasso", "SGL")
predict_table[, "beta"] <- fit_table[, "beta"] <- time_table[, "beta"] <- rep(paste0("beta_", 1:6), each = 100)

results_simulation_sigma_sqr <- matrix(NA, nrow = 100, ncol = 6)
colnames(results_simulation_sigma_sqr) <- paste0("beta_", 1:6)

for(b in 1:6) {
  print(paste0("****************************** Beta ", b, " ******************************"))
  
  results_FS_list[[b]] <- list()
  
  for(l in 1:niter) {
    print(paste0("************************** Iteration ", l, " **************************"))
    X_data <- X_data_list[[l]]
    
    ### clr transformation for X ####
    X_clr <- X_data
    X_clr[, 1:X_dim[1]] <- clr(X_data[, 1:X_dim[1]])
    X_train <- X_clr[train_index, ]
    X_test <- X_clr[test_index, ]
    Y_data <- Y_data_lognormal_list[[b]][[l]]
    Y_test <- matrix(Y_data[test_index, "log_right_censored_time"], nrow = length(test_index))
    
    results_simulation_sigma_sqr[l, b] <- Sigma_sqr_list[[b]][[l]]
    
    ## Mean imputation
    Y_train <- matrix(meanimp(Y_data[train_index, 3:4])$imputed_table[, "Imputed_time"])
    Y_indicator <- Y_data[train_index, 4]
    
    ##======================== features significance ========================##
    ## Lognormal
    results_FS_list[[b]][[l]] <- list()
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]] <- matrix(NA, nrow = X_dim[i], ncol = 8)
      row.names(results_FS_list[[b]][[l]][[i]]) <- colnames(X_train)[blocks[[i]]]
      colnames(results_FS_list[[b]][[l]][[i]]) <- c("p_value", "p_value_adjusted", "asmbPLS", "mbPLS", "BF", "IPF", "priority_lasso", "SGL")
    }
    for(i in 1:length(X_dim)) {
      for(j in 1:X_dim[i]) {
        if(sum(X_train[,blocks[[i]][j]]!=0)) {results_FS_list[[b]][[l]][[i]][j, "p_value"] <- summary(lm(Y_train ~ X_train[,blocks[[i]][j]]))$coefficients[2, 4]}
      }
      results_FS_list[[b]][[l]][[i]][, "p_value_adjusted"] <- p.adjust(results_FS_list[[b]][[l]][[i]][, "p_value"], method = "BH")
    }
    
    ##========================================== asmbPLS ==========================================##
    time_start <- Sys.time()
    ## cv ##
    asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                     Y.matrix = Y_train, 
                                     PLS.comp = PLS_term_selected, 
                                     X.dim = X_dim,
                                     quantile.comb.table = quantile.comb,
                                     Y.indicator = Y_indicator)
    n_optimal <- asmbPLS_results_cv$optimal_nPLS
    PLS_table <- matrix(asmbPLS_results_cv$quantile_table_CV[1:n_optimal, 1:2], nrow = n_optimal)
    
    ## fit ##
    asmbPLS_fit <- asmbPLS.fit(X.matrix = X_train, 
                               Y.matrix = Y_train,
                               PLS.comp = n_optimal, 
                               X.dim = X_dim,
                               quantile.comb = PLS_table)
    time_end <- Sys.time()
    time_asmbPLS <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for asmbPLS model fit: ", round(time_asmbPLS, 4)))
    
    ## predict ##
    Y_predict_asmbPLS <- asmbPLS.predict(asmbPLS_fit, X_test, n_optimal)$Y_pred
    Y_fit_asmbPLS <- asmbPLS.predict(asmbPLS_fit, X_train, n_optimal)$Y_pred
    
    predict_table[(b - 1)*100 + l, "asmbPLS"] <- result_compare(Y_predict_asmbPLS, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "asmbPLS"] <- result_compare(Y_fit_asmbPLS, Y_train)[[1]]
    
    ## feature selection ##
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "asmbPLS"] <- asmbPLS_fit$X_weight[[i]][, 1]
    }
    
    ##========================================== mbPLS ==========================================##
    time_start <- Sys.time()
    ## cv ##
    mbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                   Y.matrix = Y_train, 
                                   PLS.comp = PLS_term_selected, 
                                   X.dim = X_dim,
                                   quantile.comb.table = quantile.comb.mbPLS,
                                   Y.indicator = Y_indicator)
    n_optimal <- mbPLS_results_cv$optimal_nPLS
    
    ## fit ##
    mbPLS_fit <- mbPLS.fit(X.matrix = X_train, 
                           Y.matrix = Y_train,
                           PLS.comp = n_optimal, 
                           X.dim = X_dim)
    time_end <- Sys.time()
    time_mbPLS <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for mbPLS model fit: ", round(time_mbPLS, 4)))
    
    ## predict ##
    Y_predict_mbPLS <- asmbPLS.predict(mbPLS_fit, X_test, n_optimal)$Y_pred
    Y_fit_mbPLS <- asmbPLS.predict(mbPLS_fit, X_train, n_optimal)$Y_pred
    
    predict_table[(b - 1)*100 + l, "mbPLS"] <- result_compare(Y_predict_mbPLS, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "mbPLS"] <- result_compare(Y_fit_mbPLS, Y_train)[[1]]
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "mbPLS"] <- mbPLS_fit$X_weight[[i]][,1]
    }
    
    ##========================================== Block Forest ==========================================##
    Y_train_vector <- as.vector(Y_train)
    time_start <- Sys.time()
    blockforobj <- blockfor(X_train, 
                            Y_train_vector, 
                            replace = TRUE, 
                            blocks = blocks, 
                            block.method = "BlockForest", 
                            importance = "impurity", 
                            seed = l)
    time_end <- Sys.time()
    time_BF <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for BF model fit: ", round(time_BF, 4)))
    Y_predict_BF <- predict(blockforobj$forest, data = X_test)$predictions
    Y_fit_BF <- predict(blockforobj$forest, data = X_train)$predictions
    
    predict_table[(b - 1)*100 + l, "BF"] <- result_compare(Y_predict_BF, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "BF"] <- result_compare(Y_fit_BF, Y_train)[[1]]
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "BF"] <- blockforobj$forest$variable.importance[blocks[[i]]]
    }
    
    ##========================================== IPF ==========================================##
    set.seed(l)
    ## cv & fit##
    pflist <- list(c(1, 1), c(2, 1), c(1, 2), c(3, 1), c(1, 3), c(4, 1), c(1, 4), c(5, 1), c(1, 5))
    time_start <- Sys.time()
    ipf_fit <- cvr2.ipflasso(X = X_train, 
                             Y = Y_train_vector, 
                             family = "gaussian", 
                             type.measure = "mse",
                             blocks = blocks,
                             pflist = pflist, 
                             nfolds = 5, 
                             ncv = 5)
    time_end <- Sys.time()
    time_IPF <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for IPF model fit: ", round(time_IPF, 4)))
    
    ## predict ##
    Y_predict_IPF <- ipflasso.predict(ipf_fit, X_test)$linpredtest
    Y_fit_IPF <- ipflasso.predict(ipf_fit, X_train)$linpredtest
    
    predict_table[(b - 1)*100 + l, "IPF"] <- result_compare(Y_predict_IPF, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "IPF"] <- result_compare(Y_fit_IPF, Y_train)[[1]]
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "IPF"] <- ipf_fit$coeff[,ipf_fit$ind.bestlambda][-1][blocks[[i]]]
    }
    
    ##========================================== Prioirity Lasso ==========================================##
    time_start <- Sys.time()
    priority_results_cv <- cvm_prioritylasso(X = X_train, 
                                             Y = Y_train, 
                                             family = "gaussian", 
                                             type.measure = "mse", 
                                             standardize = T, 
                                             blocks.list = list(blocks_1, blocks_2), 
                                             block1.penalization = TRUE, 
                                             lambda.type = "lambda.1se")
    priority_fit <- priority_results_cv$best.model
    
    time_end <- Sys.time()
    time_priority_lasso <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for priority lasso model fit: ", round(time_priority_lasso, 4)))
    
    Y_predict_priority_lasso <- predict(priority_fit, X_test, type = "response")
    Y_fit_priority_lasso <- predict(priority_fit, X_train, type = "response")
    
    predict_table[(b - 1)*100 + l, "priority_lasso"] <- result_compare(Y_predict_priority_lasso, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "priority_lasso"] <- result_compare(Y_fit_priority_lasso, Y_train)[[1]]
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "priority_lasso"] <- priority_fit$coefficients[blocks[[i]]]
    }
    
    ##========================================== SGL ==========================================##
    data_SGL = list(x = X_train, y = Y_train)
    index_SGL <- c(rep(1, X_dim[1]), rep(2, X_dim[2]))
    time_start <- Sys.time()
    cvFit = cvSGL(data_SGL, index_SGL, type = "linear", nfold = 5)
    Fit = SGL(data_SGL, index_SGL, type = "linear", lambdas = cvFit$lambdas)
    time_end <- Sys.time()
    time_SGL <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for SGL model fit: ", round(time_SGL,4)))
    Y_predict_SGL <- predictSGL(Fit, X_test, which.min(cvFit$lldiff))
    Y_fit_SGL <- predictSGL(Fit, X_train, which.min(cvFit$lldiff))
    
    predict_table[(b - 1)*100 + l, "SGL"] <- result_compare(Y_predict_SGL, Y_test)[[1]]
    fit_table[(b - 1)*100 + l, "SGL"] <- result_compare(Y_fit_SGL, Y_train)[[1]]
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_FS_list[[b]][[l]][[i]][, "SGL"] <- Fit$beta[, which.min(cvFit$lldiff)][blocks[[i]]]
    }
    
    ## Time ##
    time_table[(b - 1)*100 + l, "asmbPLS"] <- time_asmbPLS
    time_table[(b - 1)*100 + l, "mbPLS"] <- time_mbPLS
    time_table[(b - 1)*100 + l, "BF"] <- time_BF
    time_table[(b - 1)*100 + l, "IPF"] <- time_IPF
    time_table[(b - 1)*100 + l, "priority_lasso"] <- time_priority_lasso
    time_table[(b - 1)*100 + l, "SGL"] <- time_SGL
  }
}

###############################################################################
# Save results
###############################################################################
path <- paste0("/blue/datta/runzhi.zhang/asmbPLS/Sim202312/q.", q, ".p.", p,"/Sim_results_lognormal/")
setwd(path)
file_name <- paste0("Results_q", q, "p", p, "_c", c, "_r", r, "_e", e, ".RData")
save(list = c('predict_table',
              'fit_table',
              "time_table",
              "results_FS_list",
              "results_simulation_sigma_sqr"),
     file = file_name)


