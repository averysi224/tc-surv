#' @title Options of machine learning methods' wrappers for fitting conditional survival curves
#' @name fit_surv_option
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to the wrapped machine learning function. Will be used in a command like `do.call(machine.learning, option)` where `machine.learning` is the machine learning function being called. `formula` and `data` should not be specified. For \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, if `tune=TRUE`, then `mtry` and `nodesize` should not be specified either.
#' @param oob whether to use out-of-bag (OOB) fitted values from random forests (\code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, \code{\link[party:cforest]{party::cforest}}) and \code{\link[grf:survival_forest]{grf::survival_forest}}) when sample splitting is not used (`nfold=1`). Ignored otherwise.
#' @param tune whether to tune `mtry` and `nodesize` for \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. Ignored for other methods.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} is used and `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @param lambda bandwidth parameter for uniform smoothing kernel in nearest neighbours estimation for method `"akritas"`. The default value of 0.5 is arbitrary and should be chosen by the user
#' @export
fit_surv_option<-function(nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list(),lambda=0.5){
    assert_that(is.count(nfold))
    assert_that(is.flag(oob))
    assert_that(is.flag(tune))
    assert_that(is.number(lambda),lambda>0)
    out<-list(nfold=nfold,option=option,oob=oob,tune=tune,tune.option=tune.option,lambda=lambda)
    class(out)<-"fit_surv_option"
    out
}

fit_surv<-function(method=c("survSuperLearner","SuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","survival_forest","no_event"),...){
    method<-match.arg(method)
    if(method=="survSuperLearner"){
        fit_survSuperLearner(...)
    } else if(method=="rfsrc"){
        fit_rfsrc(...)
    } else if(method=="coxph"){
        fit_coxph(...)
    }
}

fit_survSuperLearner<-function(formula,data,newdata,time.var,event.var,nfold=1,option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc")),...){
    # .requireNamespace("survSuperLearner")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("time","event","X","newX") %in% names(option))){
        stop("option specifies time, event, X, or newX")
    }
    
    # formula<-as.formula(paste(as.character(formula)[-2],collapse=" "))
    if(nfold==1){
        time<-data%>%pull(time.var)
        event<-data%>%pull(event.var)
        
        # X<-model.frame(formula,data=data%>%select(!c(.data[[time.var]],.data[[event.var]])))
        X<-data%>%select(-all_of(c(time.var, event.var, "C")))
        # newX<-model.frame(formula,data=newdata%>%select(!c(.data[[time.var]],.data[[event.var]])))
        newX<-newdata%>%select(-all_of(c(time.var, event.var, "C")))

        new.time<-newdata%>%pull(time.var)
        new.times<-seq(0,max(new.time),length.out=1000) #t grid
        
        arg<-c(list(time=time,event=event,X=X,newX=newX,new.times=new.times),option)
        model<-do.call(survSuperLearner::survSuperLearner,arg)
        
        event.pred<-model$event.SL.predict
        return(list(time=new.times,surv=event.pred))
    }
}


#' @title Wrapper of `randomForestSRC::rfsrc`
#' @name fit_rfsrc
#' @param formula formula used by \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. We encourage using a named list. Will be passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} by running a command like `do.call(rfsrc, option)`. The user should not specify `formula` and `data`.
#' @param oob whether to use out-of-bag (OOB) fitted values from \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} when sample splitting is not used (`nfold=1`)
#' @param tune whether to tune `mtry` and `nodesize`.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_rfsrc<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list(),...){
    .requireNamespace("randomForestSRC")
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        if(tune){
            tune.arg<-c(
                list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
                tune.option
            )
            tune.output<-do.call(randomForestSRC::tune.rfsrc,tune.arg)
        }
        
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        if(tune){
            arg<-c(arg,list(mtry=tune.output$optimal["mtry"],nodesize=tune.output$optimal["nodesize"]))
        }
        model<-do.call(randomForestSRC::rfsrc,arg)
        if(oob){
            surv<-model$survival.oob
        }else{
            surv<-model$survival
        }
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=model$time.interest,surv=surv)
    }
}


#' @title Wrapper of `survival::coxph`
#' @name fit_coxph
#' @param formula formula used by \code{\link[survival:coxph]{survival::coxph}}. Currently \code{\link[survival]{strata}} is not supported.
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survival:coxph]{survival::coxph}}. We encourage using a named list. Will be passed to \code{\link[survival:coxph]{survival::coxph}} by running a command like `do.call(coxph, option)`. The user should not specify `formula` and `data`.
#' @param ... ignored
#' @param option a list containing optional arguments passed to \code{\link[survival:coxph]{survival::coxph}}. We encourage using a named list. Will be passed to \code{\link[survival:coxph]{survival::coxph}} by running a command like `do.call(coxph, option)`. The user should not specify `formula` and `data`.
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_coxph<-function(formula,data,newdata,time.var,event.var,nfold=1,option=list(),...){
    if(nfold==1){
        # time<-data%>%pull(time.var)
        # event<-data%>%pull(event.var)
        
        XY<-data%>%select(-all_of(c("C")))
        newX<-newdata%>%select(-all_of(c(time.var, event.var, "C")))

        new.time<-newdata%>%pull(time.var)
        new.times<-seq(0,max(new.time),length.out=1000) #t grid
        arg<-c(
            list(formula=formula,data=XY), #remove id.var to allow for . in formula
            option
        )
        cox.model<-do.call(survival::coxph,arg)
        surv<-lapply(new.times,function(t){
            predict(cox.model,newdata=newX%>%mutate("{time.var}":=t, event = 1),type="survival")
        })%>%do.call(what=cbind)
        # survival_probs <- predict(cox_model, newdata = data_cal, type = "survival")
        return(list(time=new.times,surv=surv))
    }
}
