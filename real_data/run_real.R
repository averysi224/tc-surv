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
# suppressPackageStartupMessages(library(randomForestSRC))

# library(SuperLearner)
# library(flexsurv) 

# ls("SurvSuperLearner")
# browser()
# library(SuperLearner)
# library(randomForestSRC)
# library(glmnet)
# library(xgboost)
# library(survival)

## suppressPackageStartupMessages(library(caret))

########################################
### source code
########################################
# source("./source_code.R")
# source("./data_generation.R")
source("./simu.R")
library(data.table)

# Check the number of threads in use
getDTthreads()

# Set the number of threads to the maximum available
setDTthreads(32)


########################################
### run simulations
########################################
## configurations
setting_list = c("real9_12") #,"hd_heterosc")
                #  "ld_setting1","ld_setting2","ld_setting3","ld_setting4","hd_homosc","hd_heterosc")

## parameters
# Use superlearner, we make pred mse down from 44 to 12
# but even superlearner is not good enough
# the rmse is going down until 2000, when n = 5000, the rmse goes up
alpha <- .1    # target level 1-alpha
n <- 1000
n_test <- n
n_train <- n
n_calib <- n
xmin <- 0 
xmax <- 4
beta <- 20 / sqrt(n)
exp_rate <- .1

for(setting in setting_list){
  res <- list()
  reslen <- list()

  res_marg <- list()
  reslen_marg <- list()
  for (i in 2:100) {
    if(setting %in% c("hd_homosc","hd_heterosc")){
      p <- 10
    }else{
      p <- 1
    }
    p <- 3
    simures <- simu(seed + i, setting, n, p, 
                    n_train, n_calib, n_test, 
                    beta, xmin, xmax, 
                    exp_rate, alpha, "survSuperLearner") 
    res<- c(res, simures[1])
    reslen <- c(reslen, simures[2]) # reslen + simures[2]

    res_marg<- c(res_marg, simures[3])
    reslen_marg <- c(reslen_marg, simures[4]) # reslen + simures[2]
    save_dir <- sprintf("../results/%s_seed_%d_marg.csv", setting, seed+i)
    write.csv(simures, save_dir)
  }

  # browser()
  # res <- na.omit(res)
  # reslen <- na.omit(reslen)
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

  # Load the binom package for confidence interval calculations
  # Install the package if it's not already installed
  if (!require(binom)) {
    install.packages("binom", dependencies = TRUE)
    library(binom)
  } else {
    library(binom)
  }

  # Calculate the 95% Wilson score confidence interval
  conf_interval <- binom.confint(successes, n, conf.level = 0.95, methods = "wilson")

  # Display the confidence interval
  cat("95% Wilson Score Confidence Interval:\n")
  print(conf_interval)

  filename <- paste0("coverage_n", n_train, "_setting_", setting, "marg_GP.RData")

  # Save both res1 and res2 into the same file
  save(res, reslen, res_marg, reslen_marg, file = filename)
}
# To load the objects back into R:
# load(filename)

