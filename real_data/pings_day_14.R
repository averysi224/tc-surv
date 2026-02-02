########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(grf))
## library(ggExtra)
source("pings_run.R")

m0 = 14
L = 21
K = 16 # number of dataset splitting


data <- read.csv(paste0("data_", m0,".csv"))
data = data[,-1]


N0 = nrow(data)
m = 0
seed = 1234 + 7
set.seed(seed)
samp = NULL
n_test = 500
n = 1000+n_test
n_train = (n-n_test)/2
tr = 1:n_train
dat_cp = setdiff(1:N0,tr)
for(j in 1:K){
  samp = rbind(samp,c(tr,sample(dat_cp,(n+n_test)/2)))
}
T1 = data$censored_T
trueT = T1
c_list = NULL



alpha = 0.1
nam = c("adaptive LPB (T)","adaptive LPB (CT)",
        "Fix-LPB","Naive-CQR","Cox","RF")
num_method = length(nam)
mod_list = c("cox")
sim_res = NULL 
for(i in 1:K){
  sim_res <- real.simu(data[samp[i,],],trueT[samp[i,]],alpha,n_test,c_list)
  save_dir = sprintf("../results/real_data_seed_%d.csv", i)
  write.csv(sim_res, save_dir)
}
