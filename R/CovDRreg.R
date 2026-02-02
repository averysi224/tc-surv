#' @title Doubly Robust Transformation
#' @name DRtransform
#' @description
#' Given a `pred_event_censor` object in a time window, calculates the doubly robust transformation in the time window. The transformation is used as the outcome when estimating the conditional survival probability at the next visit time.
#' @param follow.up.time A data frame containing the observed follow-up times and event indicators for each individual. Must contain columns `time` and `event`, and rows must correspond to the individuals in `pred_event_censor_obj`.
#' @param pred_event_censor_obj A `pred_event_censor` object in the time window of interest.
#' @param L_tau A numeric vector of time points. Must be sorted in ascending order and all greater than the smallest time in `pred_event_censor_obj`.
#' @param next.visit.time The next visit time. Default is `Inf`, corresponding to the last time window.
#' @param denom.survival.trunc A numeric value for truncating the denominator survival probabilities to avoid division by very small numbers.
#' @return A matrix of transformations used for regression. Each row corresponds to an individual; each column corresponds to a value of `tvals`. Row names are elements in `pred_event_censor_obj$event$id`; column names are values of `tvals`.
#' @details
#' This function is designed to be called by other functions such as \code{\link{CovDRsurv}}, so inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages.
#' Should be applied on calibration set.
DRtransform <- function(Ghat.minus, integral, follow.up.time, pred_event_censor_obj, L_tau, denom.survival.trunc) {
    n_individuals <- nrow(pred_event_censor_obj$event$surv)
    psi_n_tau_wi <- matrix(nrow = n_individuals, ncol = 1)
    S_hat_L_tau_wi <- matrix(nrow = n_individuals, ncol = 1)
    
    # Ensure follow.up.time has the correct number of rows
    if (nrow(follow.up.time) != n_individuals) {
        stop("The number of rows in follow.up.time must match the number of individuals in pred_event_censor_obj.")
    }
    
    # calibration censored time
    Y <- follow.up.time$censored_T
    Delta <- follow.up.time$event
    
    for (i in seq_len(n_individuals)) {
        # Find Shat at Y
        # todo3: double check last first
        k.SY <- find.last.TRUE.index(pred_event_censor_obj$event$time <= Y[i], noTRUE = 0)
        Shat.Y <- if (k.SY == 0) 1 else pred_event_censor_obj$event$surv[i, k.SY]

        # Find Shat at L_tau(w)
        k.t <- find.last.TRUE.index(pred_event_censor_obj$event$time <= L_tau[i], noTRUE = 0)
        Shat.t <- if (k.t == 0) 1 else pred_event_censor_obj$event$surv[i, k.t]
        
        # Compute first IPW term in the bracket
        if (Delta[i] == 1 && Y[i] <= L_tau[i]) {
            # Find Ghat(Y-)
            # todo4: double check last first
            k.GY <- find.first.TRUE.index(pred_event_censor_obj$censor$time >= Y[i],
                                          noTRUE = length(pred_event_censor_obj$censor$time) + 1)
            Ghat.Yminus <- if (k.GY == 1) 1 else pred_event_censor_obj$censor$surv[i, k.GY - 1]
            
            IPW.term <- 1 / Shat.Y / pmax(Ghat.Yminus, denom.survival.trunc)
        } else {
            IPW.term <- 0
        }
        
        # Find integral at min(Y, )
        k.int <- min(k.t, k.SY)
        integral.Yt <- if (k.int == 0) 0 else integral[i, k.int]
        
        if (Shat.t != 0) {
            psi_n_tau_wi[i] <- Shat.t - Shat.t * (IPW.term - integral.Yt)
        } else {
            psi_n_tau_wi[i] <- Shat.t
        }  
        S_hat_L_tau_wi[i] <- Shat.t
    }    
    return(list(psi_n_tau_wi=psi_n_tau_wi, S_hat_L_tau_wi=S_hat_L_tau_wi))
}