library(asmbPLS)
library(dplyr)
##============================================================##
##
## asmbPLS: q = 1000, p = 200, data simulation
##
##============================================================##

setwd("/blue/datta/runzhi.zhang/asmbPLS/Sim202312/")
load("Data Simulation function.RData")
load("Parameter_list.RData")

## fixed part ##
iter = 100
n = 200
n_censored = 100
q = 1000
p = 200
m = 100000
n_effect = 500

## dimension for blocks ##
X_dim <- c(q, p)

#### Different strategies for beta ####
## First 5 - ,1:5,  remaining 0
f <- 5
beta_1_microbiome <- c(1:5, rep(0, q - f))
beta_1_metabolome <- c(1:5, rep(0, p - f))
beta_1 <- c(beta_1_microbiome, beta_1_metabolome)
beta_1 = beta_1/c(sqrt(t(beta_1)%*%beta_1))
## beta = exp(-j)
beta_2_microbiome <- exp(- (1:q))
beta_2_metabolome <- exp(- (1:p))
beta_2 <- c(beta_2_microbiome, beta_2_metabolome)
beta_2 = beta_2/c(sqrt(t(beta_2)%*%beta_2))
## beta = 1/j
beta_3_microbiome <- 1 / (1:q)
beta_3_metabolome <- 1 / (1:p)
beta_3 <- c(beta_3_microbiome, beta_3_metabolome)
beta_3 = beta_3/c(sqrt(t(beta_3)%*%beta_3))
## beta = 1 for all
beta_4_microbiome <- rep(1, q)
beta_4_metabolome <- rep(1, p)
beta_4 <- c(beta_4_microbiome, beta_4_metabolome)
beta_4 = beta_4/c(sqrt(t(beta_4)%*%beta_4))
## different proportion for different blocks
beta_5_microbiome <- c(1:5, 1:5, rep(0, q - 10))
beta_5_metabolome <- c(1:5, rep(0, p - 5))
beta_5 <- c(beta_5_microbiome, beta_5_metabolome)
beta_5 = beta_5/c(sqrt(t(beta_5)%*%beta_5))
## different proportion for different blocks
beta_6_microbiome <- c(1:5, rep(0, q - 5))
beta_6_metabolome <- c(1:5, 1:5, rep(0, p - 10))
beta_6 <- c(beta_6_microbiome, beta_6_metabolome)
beta_6 = beta_6/c(sqrt(t(beta_6)%*%beta_6))

setwd("/blue/datta/runzhi.zhang/asmbPLS/Sim202312/q.1000.p.200/Data_sim/")
for(i in 1:nrow(Parameter_comb)) {
  Dist <- Parameter_comb$Dist[i]
  Censoring_rate <- Parameter_comb$Censoring_rate[i]
  Noise_level <- Parameter_comb$Noise_level[i]
  Effect_scale <- Parameter_comb$Effect_scale[i]
  
  X_data_list <- list()
  Y_data_lognormal_list <- list()
  Y_data_Weibull_list <- list()
  Sigma_sqr_list <- list()
  
  for(b in 1:6) {
    beta <- eval(parse(text = paste0("beta_", b)))
    Y_data_lognormal_list[[b]] <- list()
    Y_data_Weibull_list[[b]] <- list()
    Sigma_sqr_list[[b]] <- list()
    for(l in 1:iter) {
      Survival_data <- one_step_simulation(n = n, 
                                           n_censored = n_censored, 
                                           q = q, 
                                           p = p, 
                                           m = m,
                                           n_effect = n_effect, 
                                           effect_scale = Effect_scale,
                                           r = Noise_level, 
                                           c = Censoring_rate, 
                                           beta = beta,
                                           seed = l)
      if(b == 1) {X_data_list[[l]] <- Survival_data$predictor_data}
      Y_data_lognormal_list[[b]][[l]] <- Survival_data$survival_time_lognormal
      Y_data_Weibull_list[[b]][[l]] <- Survival_data$survival_time_weibull
      Sigma_sqr_list[[b]][[l]] <- Survival_data$sigma_sqr
    }
  }
  c <- c("00", "01", "03", "05", "07")[which(c(0, 0.1, 0.3, 0.5, 0.7) == Censoring_rate)]
  r <- c("00", "01", "02", "05", "10")[which(c(0, 0.1, 0.2, 0.5, 1) == Noise_level)]
  e <- c("05", "20")[which(c(0.5, 2) == Effect_scale)]
  file_name <- paste0("SimData_c", c, "_r", r, "_e", e, ".RData")
  save(list = c("X_data_list",
                "Y_data_lognormal_list",
                "Y_data_Weibull_list",
                "Sigma_sqr_list"),
       file = file_name)
  print(paste0("************** Complete parameter combination ", i, " **************"))
}

