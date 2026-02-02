########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
seed<- as.integer(args[1])
if(is.na(seed)){seed <- 1}


########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(grf))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(SuperLearner))
suppressPackageStartupMessages(library(survSuperLearner))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(randomForestSRC))

suppressPackageStartupMessages(library(flexsurv))

library(data.table)
library(binom)

# Check the number of threads in use
getDTthreads()

# Set the number of threads to the maximum available
setDTthreads(20)

## suppressPackageStartupMessages(library(caret))

########################################
### source code
########################################
source("./source_code.R")
source("./data_generation.R")
source("./simu.R")


########################################
### run simulations
########################################
## configurations
setting_list = c("ld_setting1", "ld_setting2", "ld_setting3", "ld_setting4", "hd_homosc", "hd_heterosc") 
## parameters
alpha <- .1    # target level 1-alpha
conf_level <- 0.95
n <- 1000
n_test <- n
n_train <- n
n_calib <- n
xmin <- 0 
xmax <- 4
beta <- 20 / sqrt(n)
exp_rate <- .1
mode <- 'tc'   # tc - training-set conditional, marg - marginal

for(setting in setting_list){
  write(paste(setting, n, alpha, exp_rate, sep = ", "), file = "output.txt", append = TRUE)
  res <- list()
  reslen <- list()

  res_marg <- list()
  reslen_marg <- list()
  for (i in 1:100) {
    if(setting %in% c("hd_homosc","hd_heterosc")){
      p <- 10
    }else{
      p <- 1
    }
    if (mode == "marg") {
      simures <- simu_marg(seed + i, setting, n, p, 
                    n_train, n_calib, n_test, 
                    beta, xmin, xmax, 
                    exp_rate, alpha, "survSuperLearner") 
    } else if (mode == "tc") {
      simures <- simu(seed + i, setting, n, p, 
                    n_train, n_calib, n_test, 
                    beta, xmin, xmax, 
                    exp_rate, alpha, "survSuperLearner") 
    }
    
    res<- c(res, simures[1])
    reslen <- c(reslen, simures[2]) 

    res_marg<- c(res_marg, simures[3])
    reslen_marg <- c(reslen_marg, simures[4]) 
    save_dir <- sprintf("../results/%s_seed_n%d_%d_alpha_%.2f_exp_rate%.2f.csv", setting, n_train, seed+i, alpha, exp_rate)
    write.csv(simures, save_dir)
  }

  res <- res[!sapply(res, function(x) is.na(x) || is.nan(x))]
  reslen <- reslen[!sapply(reslen, function(x) is.na(x) || is.nan(x))]

  res_marg <- res_marg[!sapply(res_marg, function(x) is.na(x) || is.nan(x))]
  reslen_marg <- reslen_marg[!sapply(reslen_marg, function(x) is.na(x) || is.nan(x))]
  # Total number of observations
  n <- length(res)

  # Count the number of successes where X > 0.9
  successes <- sum(res >= 1 - alpha)

  # Calculate the sample proportion (point estimate)
  phat <- successes / n

  cat("Average LPB:", mean(unlist(reslen)), "\n", "Average coverage:", mean(unlist(res)), "\n")
  cat("Marginal average LPB:", mean(unlist(reslen_marg)), "\n", "Marginal average coverage:", mean(unlist(res_marg)), "\n")

  # Display the point estimate
  cat("Estimated P(X > 0.9):", phat, "\n")

  # Calculate the 95% Wilson score confidence interval
  conf_interval <- binom.confint(successes, n, conf.level = conf_level, methods = "wilson")

  # Display the confidence interval
  cat("95% Wilson Score Confidence Interval:\n")
  print(conf_interval)

  filename <- paste0("coverage_n", n_train, "_setting_", setting,"_alpha_", alpha, exp_rate, ".RData")

  # Save both res1 and res2 into the same file
  save(res, reslen, res_marg, reslen_marg, file = filename)

  write(paste(as.character(unlist(conf_interval)), collapse = ", "), file = "output.txt", append = TRUE, sep = "\n")
}

