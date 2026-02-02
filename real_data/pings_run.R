cov_rt = function(lb_res,T_test){sum(T_test >= lb_res)/length(T_test)}
cov_rt_comp = function(lb_res,T0,T_test){1 - sum(T_test < lb_res)/length(T0)}


real.simu = function(df0,trueT,alpha,n_test,c_list){
  
  ########################################
  ### source code
  ########################################
  source("../simulations/source_code.R")
  N = nrow(df0)
  p = ncol(df0) - 3
  xnames <- paste0('X', 1:p)
  data = df0
  
  n_train <- (N-n_test)/2
  n_calib <- n_train
  n <- n_calib
  tr <- 1:n_train
  calib <- (setdiff(1:N,tr))[1:n_calib]
  test <- setdiff(1:N,c(tr,calib))
  simures <- NULL
  
  xmin <- 0; xmax <- 4
  exp_rate <- 0.4
  alpha_list <- alpha
  mod <- "cox"
  
  
  data_fit <- data[tr,]
  data_calib <- data[calib, ]
  data_test <- data[test,]
  T0 <- data_test$censored_T
  trueT0 <- trueT[test]
  data <- rbind(data_fit, data_calib)
  

  
  ## Arguments to be passed to the subsequent functions
  x=data_test[,names(data_test) %in% xnames]
  Xtrain = data[,names(data) %in% xnames]
  C = data$C
  event = data$event
  time=data$censored_T
  fit = data_fit
  # fit$C <- -fit$C
  # mdl0 <- GauPro(X = as.matrix(fit[,names(fit) %in% xnames]),
  #                Z = fit$C, D = p,
  #                type = "Gauss")
  
  
  ## Obtaining lower bounds
  simures = NULL
  simulen = NULL
  for(mod in mod_list){
    cat("Compute the result of qt...\n")
    cat("Compute the result of qct...\n")
    lb_res = cfsurv_c(x=x,
                      Xtrain = Xtrain,
                      C = C,
                      event = event,
                      time=time,
                      alpha=alpha,
                      mdl0=mdl0)
    
    cat("Compute the result of qc0...\n")
    
    res0 = cfsurv(x,c_list=c_list,
                  pr_list=NULL,
                  pr_new_list=NULL,
                  Xtrain,C,event,time,
                  alpha=alpha,
                  type="quantile",
                  model = mod,
                  dist= "weibull",
                  I_fit = 1:n_train,
                  ftol=.1,tol=.1,
                  n.tree=100)
    
    res1 = lb_res$lower_bnd_qtl
    res2 = lb_res$lower_bnd_qctl
    output <- data.frame(qtl = res1)
    output$qctl <- res2
    output$qc0 <- res0$res
    
    # parametric methods
    ########################################
    ## vanilla CQR
    ########################################
    cat(" - Computing the result of vanilla-CQR...\n")
    res <- lapply(alpha_list,
                  cqr,
                  x=data_test[,names(data_test) %in% xnames],
                  Xtrain = data[,names(data) %in% xnames],
                  Ytrain = data$censored_T,
                  I_fit = 1:n_train,
                  seed = seed +7)
    res <- do.call(rbind,lapply(res,as.data.frame))
    output$cqr.bnd <- res[,1]
    cat("done.\n")
    
    
    ########################################
    ### utility functinos
    ########################################
    ## A function to get result
    extract_res_univariate <- function(x){
      res <- mdl_coef[1] + mdl_coef[-1]*x
      return(res)
    }
    
    ## A utility function to extract quantiles from a coxph object
    extract_quant <- function(mdl, x, alpha){
      res <- summary(survfit(mdl, newdata = x))
      time_point <- res$time
      survcdf <- 1 - res$surv
      if(sum(survcdf >= alpha)==0){
        quant = max(time_point)
      }else{
        quant <- time_point[min(which(survcdf >= alpha))]
      }
      # quant <- time_point[min(which(survcdf >= alpha))]
      return(quant)
    }
    
    ########################################
    ## Quantile regression
    ########################################
    ## Cox
    cat("Compute the result of Cox...\n")
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                             paste(xnames, collapse= "+")))
    mdl <- coxph(fmla, data = data)
    res <- c()
    for(alpha in alpha_list){
      res <- c(res, apply(data_test, 1, extract_quant, mdl = mdl,  alpha = alpha))
    }
    output$cox.bnd <- res
    cat("Done.\n")
 
    
    ## random forest
    cat("Compute the result of random forest...")
    ntree <- 1000
    nodesize <- 80
    res <- c()
    fmla <- as.formula(paste("censored_T ~ ",
                             paste(xnames, collapse= "+")))
    for(alpha in alpha_list){
      if(p==1){
        mdl <- crf.km(as.formula("censorerd_T~X1"),
                      ntree = ntree, 
                      nodesize = nodesize,
                      data_train = data[,colnames(data) %in% c("X1","censored_T","event")], 
                      data_test = data.frame(X1 = data_test$X1), 
                      yname = 'censored_T', 
                      iname = 'event',
                      tau = alpha,
                      method = "grf")
      }else{
        mdl <- crf.km(fmla,
                      ntree = ntree, 
                      nodesize = nodesize,
                      data_train = data[,colnames(data) %in% c(xnames,"censored_T","event")], 
                      data_test = data_test[,colnames(data_test) %in% xnames], 
                      yname = 'censored_T', 
                      iname = 'event',
                      tau = alpha,
                      method = "grf")
      }
      res <- c(res,mdl$predicted)
    }
    output$rf.bnd  <- res
    cat("Done.\n")
    
    id2 = (data_test$event)
    simures2 = apply(output, 2, cov_rt, T_test=trueT0)
    simures1 = apply(output, 2, cov_rt, T_test=T0)
    simures = apply(output[data_test$event,], 2, cov_rt_comp, T0, T_test=T0[data_test$event])
    simulen = apply(output, 2, mean)
  }
  simu_out = c(simures1,simures,simures2,simulen)
  return(simu_out)
}
