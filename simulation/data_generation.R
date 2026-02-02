model_generating_fun <- function(n_train, n_calib, n_test, 
                                 setting, beta, xnames, xmin, xmax, exp_rate){
  
  if(setting == "ld_setting0"){
    p <- 1
    ########################################
    ## Generate training data
    ########################################
    X <- runif(n_train, xmin, xmax)
    # T <- rexp(n_train, rate = 0.1)  # Survival times from an exponential distribution
    T <- rexp(n_train, rate = X ) 
    # Simulate censoring times (C)
    C <- rexp(n_train, rate = 0.5)
    event <- (T < C)
    censored_T <- pmin(T, C)
    # remove all C for data generation
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X <- runif(n_calib + n_test, xmin, xmax)     
    # Generate survival times (T), affected by age
    T <- rexp(n_calib + n_test, rate = X )  # Hazard rate increases with age

    # Generate censoring times (C)
    C <- rexp(n_calib + n_test, rate = 0.5)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, event = event, C = C, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting1"){
    p <- 1
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(beta * x +  2 * rnorm(length(x)))
    gen_c <- function(x) rexp(rate = exp_rate, n = length(x))
    
    ########################################
    ## Generate training data
    ########################################
    X <- runif(n_train, xmin, xmax)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    # remove all C for data generation
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X <- runif(n_calib + n_test, xmin, xmax)     
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, event = event, C = C, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting2"){
    p <- 1

    ########################################
    ## Data generating models
    ########################################
    sigma_x <- function(x) (5 + x)/10
    gen_t <- function(x) exp(3*(x>2) + 1*x*(x<=2) + 0.5 * rnorm(length(x)))
    gen_c <- function(x) rexp(rate = exp_rate, n = length(x))
    
    ########################################
    ## Generate training data
    ########################################
    X <- runif(n_train, xmin, xmax)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X <- runif(n_calib + n_test, xmin, xmax)     
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, event = event, C = C, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting3"){
    p <- 1
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(2 * (x>2) + 1 * x *(x<=2) + 0.5 * rnorm(length(x)))
    gen_c <- function(x) rexp(rate = exp_rate * (2.5 + (6+x)/10), n = length(x))
    
    ########################################
    ## Generate training data
    ########################################
    # set.seed(24601)
    X <- runif(n_train, xmin, xmax)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    # set.seed(seed)
    X <- runif(n_calib + n_test, xmin, xmax)     
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, event = event, C = C, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  if(setting == "ld_setting4"){
    
    p <- 1
    
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(3 * (x>2) + 1.5 * x *(x<=2) + 0.5 * rnorm(length(x)))
    gen_c <- function(x) exp(1 + (2-x) / 50 + 0.5 * rnorm(length(x)))
    
    ########################################
    ## Generate training data
    ########################################
    # set.seed(24601)
    X <- runif(n_train, xmin, xmax)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    # set.seed(seed)
    X <- runif(n_calib + n_test, xmin, xmax)     
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, event = event, C = C, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  

  if(setting == "hd_homosc"){
    
    p <- length(xnames)

    ########################################
    ## Data generating models
    ########################################
    mu_x <- function(x) (beta * x[,1] + beta * sqrt(x[,3] * x[,5])) / 5 + 1
    gen_t <- function(x) exp(mu_x(x) + rnorm(dim(x)[1]))
    gen_c <- function(x) rexp(rate = exp_rate * (x[,10] + 0.5), n = dim(x)[1])
    
    ## Generate training data
    X <- matrix(runif(n_train * p, min = xmin, max = xmax), n_train)
    T <- gen_t(X)
    C <- gen_c(X) 
    event <- (T<C)
    censored_T <- pmin(T,C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    colnames(data_fit) <- c(xnames, "C", "censored_T", "event")
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X <- matrix(runif((n_calib + n_test) * p, min = xmin, max = xmax), n_calib + n_test)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T<C)
    censored_T <- pmin(T,C)
    data <- data.frame(X, C = C, censored_T = censored_T,  event = event)
    colnames(data) <- c(xnames, "C", "censored_T", "event")
    data_calib <- data[1:n_calib,]
    data_test <- data[(n_calib+1) : (n_calib+n_test),]
    data <- rbind(data_fit,data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }


  if(setting == "hd_heterosc"){
    
    p <- length(xnames)

    ########################################
    ## Data generating models
    ########################################
    mu_x <- function(x) (beta * x[,1] + beta * sqrt(x[,3] * x[,5])) / 5 + 1
    sigma_x <- function(x) (x[,2] + 2) / 4
    gen_t <- function(x) exp(mu_x(x) + sigma_x(x) * rnorm(dim(x)[1]))
    gen_c <- function(x) rexp(rate = exp_rate * (x[,10] + 0.5), n = dim(x)[1])
    
    ## Generate training data
    X <- matrix(runif(n_train * p, min = xmin, max = xmax), n_train)
    T <- gen_t(X)
    C <- gen_c(X) 
    event <- (T<C)
    censored_T <- pmin(T,C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    colnames(data_fit) <- c(xnames, "C", "censored_T", "event")
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X <- matrix(runif((n_calib + n_test) * p, min = xmin, max = xmax), n_calib + n_test)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T<C)
    censored_T <- pmin(T,C)
    data <- data.frame(X, C = C, censored_T = censored_T,  event = event)
    colnames(data) <- c(xnames, "C", "censored_T", "event")
    data_calib <- data[1:n_calib,]
    data_test <- data[(n_calib+1) : (n_calib+n_test),]
    data <- rbind(data_fit,data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  

  ## Collect results 
  obj <- list(data_fit = data_fit, 
              data_calib = data_calib, 
              data_test = data_test, 
              data = data, 
              T_test = T_test)

  return(obj)

}
