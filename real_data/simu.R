simu <- function(seed, setting, n, p, 
                 n_train, n_calib, n_test, 
                 beta, xmin, xmax, 
                 exp_rate, alpha, mod_list){
  
  # ## Initialization
  set.seed(seed)
  # xnames <- paste0("X",1:p) 

  # ## Generate data according to the setting
  m0 = 9
  L = 21
  K = 16 # number of dataset splitting

  data <- read.csv(paste0("data_", m0,"_12.csv"))
  data = data[,-1]
  
  N0 = nrow(data)
  tr <- 1:n_train
  dat_cp = setdiff(1:N0,tr)
  samp = sample(dat_cp,3000,replace = TRUE)
  # trueT = data$censored_T
  
  df0 = data[samp,]
  # trueT = trueT[samp,]
  source("../simulation/source_code.R")

  N = nrow(df0)
  p = 3 # ncol(df0) - 3
  xnames <- paste0('X', 1:p)
  data = df0
  calib <- (setdiff(1:N,tr))[1:n_calib]
  test <- setdiff(1:N,c(tr,calib))

  data_fit <- data[tr,]
  data_calib <- data[calib, ]
  data_test <- data[test,]
  # browser()
  T_test <- data_test$T
  

  data_fit <- data_fit %>% select(-T)
  data_calib <- data_calib %>% select(-T)
  data_test <- data_test %>% select(-T)

  simures = NULL
  simulen = NULL
  for(mod in mod_list){
    # Step 2: Prepare covariates and follow-up time data frames
    covariates_fit <- data_fit[, c(xnames), drop=FALSE]
    follow.up.time_fit <- data_fit[, c("censored_T", "event", "C")]

    covariates_calib <- data_calib[, c(xnames), drop=FALSE]
    follow.up.time_calib <- data_calib[, c("censored_T", "event", "C")]

    # Step 3: Set variable names for CovDRsurv function
    time.var <- "censored_T"
    event.var <- "event"
    # Q.formula <- ~ .  # Use all covariates for estimating P(T > t | W = w)

    taus <- seq(0.01, 0.2, by = 0.002)  # Adjust as needed

    # Step 5: Call the CovDRsurv function with appropriate arguments
    # tau <- CovDRsurv(
    #   covariates = covariates_fit,
    #   follow.up.time = follow.up.time_fit,
    #   cal.covariates = covariates_calib,
    #   cal.follow.up.time = follow.up.time_calib,
    #   candidate.taus = taus,
    #   time.var = time.var,
    #   event.var = event.var,
    #   event.formula = ~ .,       # Use all covariates for event model
    #   censor.formula = ~ .,      # Use all covariates for censoring model
    #   # Q.formula = ~ .,           # Use all covariates for estimating P(T > t | W = w)
    #   event.method = mod,    # You can choose any supported method
    #   censor.method = mod,   # Should match event.method if using 'survSuperLearner'
    #   conf.level = 0.95,
    # )

    # tau_marg <- CovDRSurvMarg(
    #   covariates = covariates_fit,
    #   follow.up.time = follow.up.time_fit,
    #   cal.covariates = covariates_calib,
    #   cal.follow.up.time = follow.up.time_calib,
    #   candidate.taus = taus,
    #   time.var = time.var,
    #   event.var = event.var,
    #   event.formula = ~ .,       # Use all covariates for event model
    #   censor.formula = ~ .,      # Use all covariates for censoring model
    #   # Q.formula = ~ .,           # Use all covariates for estimating P(T > t | W = w)
    #   event.method = mod,    # You can choose any supported method
    #   censor.method = mod,   # Should match event.method if using 'survSuperLearner'
    #   conf.level = 0.95,
    # )

    tau_marg <- CovDRSurvMargC(
      covariates = covariates_fit,
      follow.up.time = follow.up.time_fit,
      cal.covariates = covariates_calib,
      cal.follow.up.time = follow.up.time_calib,
      candidate.taus = taus,
      time.var = time.var,
      event.var = event.var,
      event.formula = ~ .,       # Use all covariates for event model
      censor.formula = ~ .,      # Use all covariates for censoring model
      # Q.formula = ~ .,           # Use all covariates for estimating P(T > t | W = w)
      event.method = mod,    # You can choose any supported method
      censor.method = "coxph",   # Should match event.method if using 'survSuperLearner'
      conf.level = 0.95,
      p=p
    )

    # Step 6: Prepare test covariates for prediction
    covariates_test <- data_test[, c(xnames), drop=FALSE]
    event.surv.data <- cbind(covariates_fit, follow.up.time_fit)
    event.surv.test.data <- cbind(covariates_test, follow.up.time_fit)

    fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                           paste(as.character("~ ."),collapse=""),
                           collapse=""))
    event.control <- fit_surv_option(option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc")))
    # event.control=fit_surv_option(
    #     option=list(event.SL.library=c("survSL.weibreg"),
    #                 cens.SL.library=c("survSL.weibreg")))
    
    fit_surv_arg<-c(
        list(method=mod,formula=fmla_event,data=event.surv.data,newdata=event.surv.test.data,time.var=time.var,event.var=event.var),
        event.control
    )
    pred_event_obj<-do.call(fit_surv, fit_surv_arg)
    
    find_quantile_time_final <- function(surv_probs, times, quantile) {
        # browser()
        index <- which.min(surv_probs >= quantile)
        return(times[index])
    }
    # L_tau <- apply(pred_event_obj$surv, 1, find_quantile_time_final, times = pred_event_obj$time, quantile = 1-tau)
    L_tau_marg <- apply(pred_event_obj$surv, 1, find_quantile_time_final, times = pred_event_obj$time, quantile = 1-tau_marg)

    # L_tau <- matrix(L_tau, nrow = length(L_tau), ncol = 1)  # Convert vector to matrix
    L_tau_marg <- matrix(L_tau_marg, nrow = length(L_tau_marg), ncol = 1) 

    # Now, 'predictions_df' contains the estimated P(T > t | W = w) for each individual in test data at each time point in 'tvals'
    cov_rt = function(lb_res, T_test){sum(T_test >= lb_res)/length(lb_res)}

    simulen <- 0
    simures <- 0

    simures_marg <- 0
    simulen_marg <- 0
    
    # simures = apply(L_tau, 2, cov_rt, T_test=T_test)
    # simulen = apply(L_tau, 2, mean)

    simures_marg = apply(L_tau_marg, 2, cov_rt, T_test=T_test)
    simulen_marg = apply(L_tau_marg, 2, mean)

    # Print the coverage
    cat(seed, "Coverage probability at alpha =", alpha, "is:", simures, "lpb: ", simulen, "\n")
    cat(seed, "Marginal coverage probability at alpha =", alpha, "is:", simures_marg, "lpb: ", simulen_marg, "\n")
  }
  simu_out = c(simures, simulen, simures_marg, simulen_marg)
  return(simu_out=simu_out)
}


