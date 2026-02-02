#shift censoring time to the left a tiny bit
#used for survival analysis with event observed iff T<C (rather than traditionally T<=C)
#status: 1 if observed event; 0 if censoring
left.shift.censoring<-function(time,status){
    epsilon<-time%>%unique%>%sort%>%diff%>%min
    epsilon<-min(epsilon*.5,1e-5)
    time[status==0]<-time[status==0]-epsilon
    time
}

#find the first index of x that is TRUE
#noTRUE: return value when none of x is TRUE
find.first.TRUE.index<-function(x,noTRUE=length(x)+1){
    match(TRUE,x,nomatch=noTRUE)
}

find_quantile_time <- function(surv_probs, times, quantile) {
    index <- which.max(surv_probs <= quantile)
    return(times[index])
}

#find the last index of x that is TRUE (loop from the last entry back to the first to and find the first encountered TRUE)
#noTRUE: return value when none of x is TRUE
find.last.TRUE.index<-function(x,noTRUE=0){
    length(x)+1-match(TRUE,rev(x),nomatch=length(x)+1-noTRUE)
}

cens_prob_mat <- function(gpr_mdl, calib, t_vec, xnames){
    ## Computing the censoring scores for the calibration data
    predictions <- gpr_mdl$predict(as.matrix(calib[,names(calib) %in% xnames]),
                              se.fit = TRUE)

    mean_calib <- predictions$mean
    sd_calib <- predictions$se

    prob_matrix <- matrix(NA, nrow = length(mean_calib), ncol = length(t_vec))
    colnames(prob_matrix) <- paste0("t=", t_vec)

    # Calculate P(C > t_i | W = w_i) for each time point t_i
    for (j in seq_along(t_vec)) {
        t_i <- t_vec[j]
        # pr_calib <- 1 - pnorm((t_i - mean_calib) / sd_calib)
        pr_calib <- 1 - pnorm(t_i, mean = mean_calib, sd = sd_calib)
        prob_matrix[, j] <- pr_calib
    }

    return(prob_matrix)
}