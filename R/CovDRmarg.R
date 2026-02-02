# assume C is unknown, only (W, tilda_T and Delta) are known
CovDRSurvMarg<-function(
    covariates,
    follow.up.time,
    cal.covariates,
    cal.follow.up.time,
    candidate.taus=NULL,
    truncation.index=1,
    time.var,
    event.var,
    event.formula=NULL,
    censor.formula=NULL,
    # Q.formula=~.,
    event.method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas","survival_forest"),
    censor.method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas","survival_forest"),
    event.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"))),
    censor.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"))),
    Q.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    denom.survival.trunc=1e-3,
    conf.level=.95,
    p=1
){
    event.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))
    censor.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph"),
                    cens.SL.library=c("survSL.coxph")))

    event.surv.data <- cbind(covariates, follow.up.time)
    # cal.event.follow.up.time<-follow.up.time
    event.surv.cal.data <- cbind(cal.covariates, cal.follow.up.time)
    
    fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                           paste(as.character(event.formula),collapse=""),
                           collapse=""))

    fit_surv_arg<-c(
        list(method=event.method,formula=fmla_event,data=event.surv.data,newdata=event.surv.cal.data,time.var=time.var,event.var=event.var),
        event.control
    )
    pred_event_obj<-do.call(fit_surv, fit_surv_arg)

    # fit G = P(C | X) with superlearner 
    censor.follow.up.time<-follow.up.time%>%
                mutate("{event.var}":=1-.data[[event.var]],
                       "{time.var}":=left.shift.censoring(.data[[time.var]],.data[[event.var]]))
    censor.surv.data <- cbind(covariates, censor.follow.up.time)
    cal.censor.follow.up.time<-cal.follow.up.time%>%
                mutate("{event.var}":=1-.data[[event.var]],
                       "{time.var}":=left.shift.censoring(.data[[time.var]],.data[[event.var]]))
    censor.surv.cal.data <- cbind(cal.covariates, cal.censor.follow.up.time)

    fmla_censor <- as.formula(paste("Surv(",time.var,",",event.var,")",
                           paste(as.character(censor.formula),collapse=""),
                           collapse=""))
    fit_surv_arg<-c(
        list(method=censor.method,formula=fmla_censor,data=censor.surv.data,newdata=censor.surv.cal.data,time.var=time.var,event.var=event.var),
        censor.control
    )
    pred_censor_obj<-do.call(fit_surv, fit_surv_arg)    

    # P(T>=t∣X) event probability ( or Shat) at grid points
    pred_event_censor_obj<-list(event=pred_event_obj,censor=pred_censor_obj)

    # initP(C>=t∣X) event probability at grid points
    Ghat.minus <- pred_event_censor_obj$event$surv

    # find P(C>=t∣X) at each event grid points t probability
    # row is still pid
    for (i in seq_along(pred_event_censor_obj$event$time)) {
        # Find index j where censoring time >= current event time
        j <- find.first.TRUE.index(pred_event_censor_obj$censor$time >= pred_event_censor_obj$event$time[i],
                                   noTRUE = length(pred_event_censor_obj$censor$time) + 1)
        if (j == 1) {
            Ghat.minus[, i] <- 1
        } else {
            Ghat.minus[, i] <- pred_event_censor_obj$censor$surv[, j - 1]
        }
    }
    
    integrand <- pred_event_censor_obj$event$surv
    for (i in seq_along(pred_event_censor_obj$event$time)) {
        if (i == 1) {
            integrand[, i] <- (1 - pred_event_censor_obj$event$surv[, i]) /
                              pred_event_censor_obj$event$surv[, i] /
                              pmax(Ghat.minus[, i], denom.survival.trunc)
        }
        else if (i == length(pred_event_censor_obj$event$time)) {
            integrand[, i] <- integrand[, i-1]
        } else {
            integrand[, i] <- (pred_event_censor_obj$event$surv[, i-1] - pred_event_censor_obj$event$surv[, i]) /
                                pred_event_censor_obj$event$surv[, i-1] /
                              pred_event_censor_obj$event$surv[, i-1] /
                              pmax(Ghat.minus[, i-1], denom.survival.trunc)
        }
    }
    # For each element in a row, the value is the sum of all previous elements in that row (including the current one)
    integral <- rowCumsums(integrand)

    ## start calibration
    n<-nrow(cal.covariates)
    good_taus <- list()
    max_tau <- 0
    ############################################################################
    # doubly robust transformation and regression
    ############################################################################
    psi<-SE<-CI.lower<-numeric(length(candidate.taus))
    for(tau.index in 1:length(candidate.taus)){
        tau<-candidate.taus[tau.index]
        # Estimate the median survival time for each observation
        L_tau <- apply(pred_event_censor_obj$event$surv, 1, find_quantile_time, times = pred_event_censor_obj$event$time, quantile = 1-tau)
        
        results <- DRtransform(Ghat.minus, integral, cal.follow.up.time, pred_event_censor_obj, L_tau,denom.survival.trunc)

        S_hat_L_tau_wi = results$S_hat_L_tau_wi
        psi_n_tau_wi = results$psi_n_tau_wi

        naive.psi <- mean(S_hat_L_tau_wi)  # Compute as the average over the calibration set
        # Step 4: Compute the pathwise derivative values
        D_tau_values <- psi_n_tau_wi - naive.psi
        psi_n_tau <- naive.psi + mean(D_tau_values)

        if (psi_n_tau > 1 - alpha) {
            max_tau <- tau.index
        } else {
            break
        }
    }

    if(max_tau>0){
        if(max_tau==length(candidate.taus)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-candidate.taus[max_tau]
        cat(tau, psi_n_tau, "\n")
        return(tau)
    }
    else{
        cat("Try smaller tau.")
        return(0.0)
    }

}

# assume C is known, the tuple is (W, C, tilda_T)
CovDRSurvMargC<-function(
    covariates,
    follow.up.time,
    cal.covariates,
    cal.follow.up.time,
    candidate.taus=NULL,
    truncation.index=1,
    time.var,
    event.var,
    event.formula=NULL,
    censor.formula=NULL,
    # Q.formula=~.,
    event.method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas","survival_forest"),
    censor.method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas","survival_forest"),
    event.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"))),
    censor.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"))),
    Q.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    denom.survival.trunc=1e-3,
    conf.level=.95,
    p=1
){
    xnames <- paste0("X",1:p) 
    event.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))

    censor.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph"),
                    cens.SL.library=c("survSL.coxph")))

    event.surv.data <- cbind(covariates, follow.up.time)
    event.surv.cal.data <- cbind(cal.covariates, cal.follow.up.time)
    
    fmla_event <- as.formula(paste("Surv(",time.var,",",event.var,")",
                           paste(as.character(event.formula),collapse=""),
                           collapse=""))

    fit_surv_arg<-c(
        list(method=event.method,formula=fmla_event,data=event.surv.data,newdata=event.surv.cal.data,time.var=time.var,event.var=event.var),
        event.control
    )
    pred_event_obj<-do.call(fit_surv, fit_surv_arg)

    # fit G = P(C | X) with coxph 
    censor.follow.up.time<-follow.up.time%>%
                mutate("{event.var}":=1,
                       "{time.var}":=.data[["C"]])
    censor.surv.data <- cbind(covariates, censor.follow.up.time)
    cal.censor.follow.up.time<-cal.follow.up.time%>%
                mutate("{event.var}":=1,
                       "{time.var}":=.data[["C"]])
    censor.surv.cal.data <- cbind(cal.covariates, cal.censor.follow.up.time)

    fmla_censor <- as.formula(paste("Surv(",time.var,",",event.var,")",
                           paste(as.character(censor.formula),collapse=""),
                           collapse=""))
    fit_surv_arg<-c(
        list(method=censor.method,formula=fmla_censor,data=censor.surv.data,newdata=censor.surv.cal.data,time.var=time.var,event.var=event.var)
    )
    pred_censor_obj<-do.call(fit_surv, fit_surv_arg)  
    
    # P(T>=t∣X) event probability ( or Shat) at grid points
    pred_event_censor_obj<-list(event=pred_event_obj,censor=pred_censor_obj)

    # initP(C>=t∣X) event probability at grid points
    Ghat.minus <- pred_event_censor_obj$event$surv

    # find P(C>=t∣X) at each event grid points t probability
    # row is still pid
    for (i in seq_along(pred_event_censor_obj$event$time)) {
        # Find index j where censoring time >= current event time
        j <- find.first.TRUE.index(pred_event_censor_obj$censor$time >= pred_event_censor_obj$event$time[i],
                                   noTRUE = length(pred_event_censor_obj$censor$time) + 1)
        if (j == 1) {
            Ghat.minus[, i] <- 1
        } else {
            Ghat.minus[, i] <- pred_event_censor_obj$censor$surv[, j - 1]
        }
    }
    
    integrand <- pred_event_censor_obj$event$surv
    for (i in seq_along(pred_event_censor_obj$event$time)) {
        if (i == 1) {
            integrand[, i] <- (1 - pred_event_censor_obj$event$surv[, i]) /
                              pred_event_censor_obj$event$surv[, i] /
                              pmax(Ghat.minus[, i], denom.survival.trunc)
        }
        else if (i == length(pred_event_censor_obj$event$time)) {
            integrand[, i] <- integrand[, i-1]
        } else {
            integrand[, i] <- (pred_event_censor_obj$event$surv[, i-1] - pred_event_censor_obj$event$surv[, i]) /
                                pred_event_censor_obj$event$surv[, i-1] /
                              pred_event_censor_obj$event$surv[, i-1] /
                              pmax(Ghat.minus[, i-1], denom.survival.trunc)
        }
    }
    # For each element in a row, the value is the sum of all previous elements in that row (including the current one)
    integral <- rowCumsums(integrand)

    ## start calibration
    n<-nrow(cal.covariates)
    good_taus <- list()
    max_tau <- 0
    ############################################################################
    # doubly robust transformation and regression
    ############################################################################
    psi<-SE<-CI.lower<-numeric(length(candidate.taus))
    for(tau.index in 1:length(candidate.taus)){
        tau<-candidate.taus[tau.index]
        # Estimate the median survival time for each observation
        L_tau <- apply(pred_event_censor_obj$event$surv, 1, find_quantile_time, times = pred_event_censor_obj$event$time, quantile = 1-tau)
        
        results <- DRtransform(Ghat.minus, integral, cal.follow.up.time, pred_event_censor_obj, L_tau,denom.survival.trunc)

        S_hat_L_tau_wi = results$S_hat_L_tau_wi
        psi_n_tau_wi = results$psi_n_tau_wi

        naive.psi <- mean(S_hat_L_tau_wi)  # Compute as the average over the calibration set
        # Step 4: Compute the pathwise derivative values
        D_tau_values <- psi_n_tau_wi - naive.psi
        psi_n_tau <- naive.psi + mean(D_tau_values)

        if (psi_n_tau > 1 - alpha) {
            max_tau <- tau.index
        } else {
            break
        }
    }

    if(max_tau>0){
        if(max_tau==length(candidate.taus)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-candidate.taus[max_tau]
        cat(tau, psi_n_tau, "\n")
        return(tau)
    }
    else{
        cat("Try smaller tau.")
        return(0.0)
    }

}

