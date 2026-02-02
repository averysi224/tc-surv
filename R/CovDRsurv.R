#' @title Estimate (training-set conditional) survival probabilities with doubly robust transformation
#' @name CovDRsurv
#' @description
#' Estimate P(T > t | T > truncation time, covariates available at truncation time) for given t, where T is the time to event, using ly doubly robust transformation. Use a user-specified flexible method to fit survival curves of time to event/censoring at each stage and then use \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to regress pseudo-outcome on covariates in order to estimate P(T > t | T > truncation time, covariates available at truncation time).
#'
#' @param covariates a list of data frames of covarates in the order of check in times. Each data frame contains the covariates collected at a visit time. Data frames may have different numbers of variables (may collect different variables at different visit times) and different numbers of individuals (some individuals may have an event or is censored before a later visit time). All data frames must have a common character variable (see `id.var`) that identifies each individual but no other variables with common names. No missing data is allowed.
#' @param follow.up.time data frame of follow up times, i.e., times to event/censoring. Contains the follow up times and an indicator of event/(right-)censoring. Follow up times must be numeric. Indicator of event/censoring should be binary with 0=censored, 1=event.
#' @param candidate.taus times t for which P(T > t) given covariates are computed (T is the time to event). Default is all unique event times in `follow.up.time`. Will be sorted in ascending order.
#' @param truncation.index index of the visit time to which left-truncation is applied. The truncation time is `visit.times[truncation.index]`. Covariates available up to (inclusive) `visit.times[truncation.index]` are of interest. Default is 1, corresponding to no truncation.
#' @param time.var (character) name of the variable containing follow up times in the data frame `follow.up.time`.
#' @param event.var (character) name of the variable containing indicator of event/censoring in the data frame `follow.up.time`.
#' @param event.formula a list of formulas to specify covariates being used when estimating the conditional survival probabilities of time to event at each visit time. The length should be the number of check in times after `truncation.index` (inclusive). Default is `~ .` for all visit times, which includes main effects of all covariates available at each visit time.
#' @param censor.formula a list of formulas to specify covariates being used when estimating the conditional survival probabilities of time to censoring at each visit time. Similar to `event.formula`. If `event.method` is chosen as `"survSuperLearner"`, `censor.formula` is not used and `event.formula` is used instead.
#' @param Q.formula formula to specify covariates being used for estimating P(T > t | T > `visit.times[truncation.index]`, covariates available at `visit.times[truncation.index]`). Set to include intercept only (`~ 0` or `~ -1`) for marginal survival probability. Default is `~ .`, which includes main effects of all available covariates up to (inclusive) the `visit.times[truncation.index]`.
#' @param event.method one of `"survSuperLearner"`, `"rfsrc"`, `"ctree"`, `"rpart"`, `"cforest"`, `"coxph"`, `"coxtime"`, `"deepsurv"`, `"dnnsurv"`, `"akritas"`, `"survival_forest"`. The machine learning method to fit  survival curves of time to event in each time window. See the underlying wrappers \code{\link{fit_survSuperLearner}}, \code{\link{fit_rfsrc}}, \code{\link{fit_ctree}}, \code{\link{fit_rpart}}, \code{\link{fit_cforest}}, \code{\link{fit_coxph}}, \code{\link{fit_coxtime}}, \code{\link{fit_deepsurv}}, \code{\link{fit_dnnsurv}}, \code{\link{fit_akritas}}, \code{\link{fit_survival_forest}} for more details and the available options. Default is `"survSuperLearner"`, which may perform well with a decent amount of events and censoring but may fail if too few events or too little censoring in one time window. If at least one of `event.method` and `censor.method` is `"survSuperLearner"`, then survival curves of both time to event and time to censoring are fitted by \code{\link{fit_survSuperLearner}}.
#' @param censor.method one of `"survSuperLearner"`, `"rfsrc"`, `"ctree"`, `"rpart"`, `"cforest"`, `"coxph"`, `"coxtime"`, `"deepsurv"`, `"dnnsurv"`, `"akritas"`, `"survival_forest"`. The machine learning method to fit survival survival curves of time to censoring in each time window. Similar to `event.method`. Default is `"survSuperLearner"`, which may perform well with a decent amount of events and censoring but may fail if too few events or too little censoring in one time window. If at least one of `event.method` and `censor.method` is `"survSuperLearner"`, then survival curves of both time to event and time to censoring are fitted by \code{\link{fit_survSuperLearner}}.
#' @param event.control a returned value from \code{\link{fit_surv_option}} to control fitting survival curves of time to event. Ignored if `event.method` is chosen as `"survSuperLearner"`.
#' @param censor.control a returned value from \code{\link{fit_surv_option}} to control fitting survival curves of time to censoring. Ignored if `event.method` is chosen as `"survSuperLearner"`.
#' @param denom.survival.trunc the numeric truncation value for the survival function in the denominator. All denominators below `denom.survival.trunc` will be set to `denom.survival.trunc` for numerical stability.
#' @param conf.level 
#' @return a list of fitted `SuperLearner` models corresponding to each t in `tvals`.
#' @section Formula arguments:
#' All formulas should have covariates on the right-hand side and no terms on the left-hand side, e.g., `~ V1 + V2 + V3`. At each visit time, the corresponding formulas may (and usually should) contain covariates at previous visit times, and must only include available covariates up to (inclusive) that visit time. Interactions, polynomials and splines may be treated differently by different machine learning methods to estimate conditional survival curves.
#' @export
CovDRsurv<-function(
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
    conf.level=.95
){
    event.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))
    censor.control=fit_surv_option(
        option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names),
                    cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam",rfsrc.learners$names)))

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
        # browser()
        if (i == 1) {
            integrand[, i] <- (1 - pred_event_censor_obj$event$surv[, i]) /
                              pred_event_censor_obj$event$surv[, i] /
                              pmax(Ghat.minus[, i], denom.survival.trunc)
        }
        else if (i == length(pred_event_censor_obj$event$time)) {
            integrand[, i] <- integrand[, i-1]
        } else {
            # gfilter <- pmax(Ghat.minus[, i], denom.survival.trunc)
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

        # Step 5: Compute the centered influence function values
        IF <- D_tau_values 

        # Step 6: Compute the variance and standard error
        sigma2 <- sum(IF^2) / n
        SE[tau.index] <- sqrt(sigma2 / n)

        CI.lower[tau.index] <- psi_n_tau - SE[tau.index] * qnorm(conf.level)
        if (CI.lower[tau.index] > 1 - alpha) {
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
        cat(tau, CI.lower[max_tau], "\n")
        return(tau)
    }
    else{
        cat("Try smaller tau.")
        return(0.0)
    }

}

rfsrc.learners<-SuperLearner::create.Learner("survSL.rfsrc",tune=list(nodesize=c(10,15,20),ntime=c(150,200,250)),params=list(perf.type="none"))