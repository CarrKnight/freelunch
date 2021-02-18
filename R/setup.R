# run abc examples
# library(abc)
# library(abctools)
# library(R.utils)
# library(quantregForest)
# library(caret)
# library(knitrProgressBar)
# library(ranger)
# library(e1071) #caret likes this for some unidentified reason!
# library(rsample)
library(tidyverse)

DEGREE<-1 # regression degree
BOOTSTRAP_N = 200 # number of bootstraps
LOW_QUANTILE<-0.025
HIGH_QUANTILE<-0.975




.formula_maker<-function(yname,xnames,degree=DEGREE)
{
  formula(
    paste(
      paste0("`",yname,"`"),"~","(",
      paste(paste0("`",xnames,"`"),collapse = "+"),")",
      ifelse(degree>1,paste("^",degree),"")

    )
  )


}


## there are two use cases for target runs: you are testing yourself (so that the data contains the real parameters)
## or you are genuinely fitting real data at which point you obviously don't have parameter columns
## in the second case this function adds NAs everywhere
.add_NAs_to_target_runs_ifneeded<-function(target_runs,
                                           parameter_colnames){
  for(parameter in parameter_colnames){
    if(!any(parameter==colnames(target_runs))){
      target_runs[[parameter]]<-NA
    }
  }
  return(target_runs)

}


.check_that_summary_statistics_are_into_dataframe<-function(target_runs,summary_statistics_colnames){
  if(is.vector(target_runs) & !is.list(target_runs)) {
    stopifnot(length(target_runs)==length(summary_statistics_colnames))
    target_runs<-data.frame(t(target_runs))
    colnames(target_runs)<-summary_statistics_colnames
  }

  return(target_runs)

}

DEFAULT_TOL<-.1


#library(gam)
# this is in-loop cv parameters for caret for RF; not the outer-loop CV parameters
control<-caret::trainControl(method = "cv",
                             number = 5)



#   _____ ______ _____ _____ _____        _   _  ___   _     ___________  ___ _____ _____ _____ _   _
#  /  __ \| ___ \  _  /  ___/  ___|      | | | |/ _ \ | |   |_   _|  _  \/ _ \_   _|_   _|  _  | \ | |
#  | /  \/| |_/ / | | \ `--.\ `--. ______| | | / /_\ \| |     | | | | | / /_\ \| |   | | | | | |  \| |
#  | |    |    /| | | |`--. \`--. \______| | | |  _  || |     | | | | | |  _  || |   | | | | | | . ` |
#  | \__/\| |\ \\ \_/ /\__/ /\__/ /      \ \_/ / | | || |_____| |_| |/ /| | | || |  _| |_\ \_/ / |\  |
#   \____/\_| \_|\___/\____/\____/        \___/\_| |_/\_____/\___/|___/ \_| |_/\_/  \___/ \___/\_| \_/
#
#

## internal cross-validation function
.cross_validate<-function(total_data,ngroup=5,fitting_method,
                          cv_seed=0,
                          parameter_colnames,
                          summary_statistics_colnames,
                          compute_predictivity=TRUE,...){


  n<-total_data %>% nrow()
  leave.out <- trunc(n/ngroup)
  o <- R.utils::withSeed({
    sample(1:n)
  }, seed=cv_seed)
  groups <- vector("list", ngroup)
  for (j in 1:(ngroup - 1)) {
    jj <- (1 + (j - 1) * leave.out)
    groups[[j]] <- (o[jj:(jj + leave.out - 1)])
  }
  groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
  results<-list()


  for(group in 1:ngroup){
    message("*************************************")
    message(paste("cross validating group",group))

    test_data<- total_data %>%
      ungroup() %>%
      filter(row_number() %in% groups[[group]] )
    training_data<-
      suppressMessages(anti_join(total_data,test_data))

    results[[group]]<-
      fitting_method(training_runs = training_data,target_runs = test_data,
                     parameter_colnames,
                     summary_statistics_colnames,
                     ...)
  }

  predictivity<-NA
  performance<-NA
  ## now compute the same for average, so we can produce predictivity
  if(compute_predictivity){
    suppressMessages(averages<-
                       .cross_validate(total_data,ngroup=5,fit_average,cv_seed,
                                       parameter_colnames,
                                       summary_statistics_colnames,
                                       compute_predictivity=FALSE))
    errorone <- results %>% map_dfr(~bind_rows(.$rmse)) %>% colMeans()
    errortwo <- averages$rmse
    performance = 1 - errorone/errortwo

    errorone <- results %>% map_dfr(~bind_rows(.$prediction_residuals))
    errortwo <- averages$results %>% map_dfr(~bind_rows(.$prediction_residuals))
    predictivity = 1 - (errorone%>% map_dbl(~sum(.^2)))/(errortwo%>% map_dbl(~sum(.^2)))
  }


  out_of_sample_predictions<-
    map_dfr(results,~.$predictions) %>% mutate(id=row_number()) %>% gather("variable","prediction",-id)
  out_of_sample_errors<-
    map_dfr(results,~.$prediction_residuals) %>% mutate(id=row_number()) %>% gather("variable","error",-id)

  out_of_sample_errors<-inner_join(out_of_sample_predictions,out_of_sample_errors) %>% mutate(real=prediction+error)

  return(
    list(
      results=results,
      out_of_sample_errors = out_of_sample_errors,
      rmse= results %>% map_dfr(~bind_rows(.$rmse)) %>% colMeans,
      contained =  results %>% map_dfr(~bind_rows(.$contained)) %>% colMeans,
      interval_size =  results %>% map_dfr(~bind_rows(.$interval_size)) %>% colMeans,
      performance = performance,
      predictivity = predictivity
    )
  )

}


#   _     _                        ______                             _
#  | |   (_)                       | ___ \                           (_)
#  | |    _ _ __   ___  __ _ _ __  | |_/ /___  __ _ _ __ ___  ___ ___ _  ___  _ __
#  | |   | | '_ \ / _ \/ _` | '__| |    // _ \/ _` | '__/ _ \/ __/ __| |/ _ \| '_ \
#  | |___| | | | |  __/ (_| | |    | |\ \  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
#  \_____/_|_| |_|\___|\__,_|_|    \_| \_\___|\__, |_|  \___||___/___/_|\___/|_| |_|
#                                              __/ |
#

#' Tries to estimate the parameters through a linear regression of specified degree and uses bootstrap to generate
#' prediction intervals that make sense
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param degree the degree of the linear regression we are using to map back from summary statistics to parameters
#' @param bootstrap_n the number of boostrap resamples we run to generate prediction intervals (default is 200)
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ##notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'################# LINEAR REGRESSION
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#'fit_linear_regression(training_runs = training_data,
#'                      target_runs =  testing_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'fit_linear_regression(training_runs = training_data,
#'                      target_runs =  c(2,1),
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame, especially with parameters included, you get error and coverage stuff
#'fit_linear_regression(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
fit_linear_regression<-function(
  training_runs,target_runs,parameter_colnames,
  summary_statistics_colnames, # name of the X variables in the regression (summary statistics)
  degree,bootstrap_n = BOOTSTRAP_N){


  ##you are allowed to send summary statistics into a vector if you are just fitting once
  target_runs<-.check_that_summary_statistics_are_into_dataframe(target_runs,summary_statistics_colnames)
  ## check if the target runs need NAs added (happens when we aren't testing but properly fitting!)
  target_runs<-.add_NAs_to_target_runs_ifneeded(target_runs,parameter_colnames)


  regressions<-list()
  rmse<-rep.int(NA,length(parameter_colnames))
  interval_size<-rep.int(NA,length(parameter_colnames))
  contained<-rep.int(NA,length(parameter_colnames))
  predictions<-list()
  lows<-list()
  highs<-list()
  prediction_residuals<-list()
  for(i in 1:length(parameter_colnames)){
    y<-parameter_colnames[i]
    message(paste("fitting",y))
    formula1<-.formula_maker(yname = y,
                             xnames = summary_statistics_colnames,
                             degree = 1)

    regression<-lm(formula1,data=training_runs)
    residuals<-residuals(regression)/sqrt(1-hatvalues(regression))

    current_prediction<-predict(regression,  newdata=target_runs)
    predictions[[y]]<-  current_prediction

    rmse[i]<-  sqrt(mean((predictions[[y]]-target_runs[[y]])^2))

    # progress_bar<-NULL
    #  if((target_runs %>% nrow())>1)
    progress_bar<- knitrProgressBar::progress_estimated(BOOTSTRAP_N)
    #fit a regression function
    .boostrap_fit<-function(bootdata)
    {
      if(!is.null(progress_bar))
        knitrProgressBar::update_progress(progress_bar)
      return(
        lm(formula = formula1,data=rsample::analysis(bootdata))
      )

    }

    #fit many times
    booted<-rsample::bootstraps(training_runs,times=BOOTSTRAP_N) %>%
      mutate(model = purrr::map(splits,.boostrap_fit)) %>%
      mutate(predictions = purrr::map(model,predict,newdata=target_runs,se=F)) %>%
      mutate(residuals = purrr::map(model,
                                    function(x) sample(residuals(x)/sqrt(1-hatvalues(x)),
                                                       size=length(current_prediction))

      )) %>%
      filter(!is.null(predictions))
    #
    # residuals <-
    #   map(booted$residuals, ~ sample(.,size=10000/length(booted$residuals)))  %>% unlist()
    #
    #observe prediction errors
    ses<-
      booted %>% mutate(original=list(current_prediction)) %>%
      unnest(c(predictions,original,residuals)) %>% ungroup() %>%
      filter(is.finite(predictions)) %>% filter(!is.na(predictions)) %>%
      mutate(se = original-predictions) %>%
      mutate(full = se + residuals) %>%
      group_by(id) %>% mutate(index=dplyr::row_number())
    # bootstrap a big list of them
    test<-ses  %>% group_by(index) %>%
      summarise(error=list(sample(full,replace=TRUE,size = 10000)),.groups="keep")

    # bootstrap also residuals
    # test$residuals<-list(residuals)
    #now combine these with resampled errors and grab quantiles

    #should be done!
    lows[[y]]<- predictions[[y]]  + (map(test$error,
                                         function(x) quantile(x,LOW_QUANTILE)) %>% unlist())
    highs[[y]]<- predictions[[y]] + (map(test$error,
                                         function(x) quantile(x,HIGH_QUANTILE)) %>% unlist())
    prediction_residuals[[y]]<- target_runs[[y]] - predictions[[y]]

    contained[i]<- sum( (target_runs[[y]]>=lows[[y]]) & (target_runs[[y]]<=highs[[y]]))/length(predictions[[y]])
    interval_size[i]<-mean(highs[[y]]-lows[[y]])
  }



  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames


  return(
    list(
      predictions= as_tibble(predictions),
      lows = as_tibble(lows),
      highs = as_tibble(highs),
      prediction_residuals = as.data.frame(prediction_residuals),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )




}

#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_linear_regression} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param degree the degree of the linear regression we are using to map back from summary statistics to parameters
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_linear_regression} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_linear_regression} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_linear_regression} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ##notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#' cv_results<-cross_validate_lm(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#'
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_lm<-function(total_data,ngroup=5,
                            parameter_colnames,summary_statistics_colnames,
                            degree=DEGREE,
                            cv_seed=0){
  .cross_validate(total_data=total_data,cv_seed=cv_seed,
                  ngroup=ngroup,
                  fitting_method=fit_linear_regression,
                  parameter_colnames=parameter_colnames,
                  summary_statistics_colnames=summary_statistics_colnames,
                  compute_predictivity = TRUE,
                  degree=degree
  )
}


#    ___  ______  _____
#   / _ \ | ___ \/  __ \
#  / /_\ \| |_/ /| /  \/
#  |  _  || ___ \| |
#  | | | || |_/ /| \__/\
#  \_| |_/\____/  \____/
#
#


# performs ABC on each single row of the training set independently.
# this to give the full benefit to ABC
abc_testing<-function(training_runs,target_runs,tol=DEFAULT_TOL,parameter_colnames,
                      summary_statistics_colnames = NULL,
                      method="rejection",semiauto=FALSE,
                      #only used in semi-auto
                      satr=list(function(x){outer(x,Y=1:4,"^")}),# this comes from the COAL example; regress up to degree 4
                      ...)
{

  ##you are allowed to send summary statistics into a vector if you are just fitting once
  target_runs<-.check_that_summary_statistics_are_into_dataframe(target_runs,summary_statistics_colnames)
  ## check if the target runs need NAs added (happens when we aren't testing but properly fitting!)
  target_runs<-.add_NAs_to_target_runs_ifneeded(target_runs,parameter_colnames)



  #abctools works on matrices, not data.frames so we need to convert
  obsparam = target_runs %>%
    dplyr::select(parameter_colnames) %>% as.matrix()

  # print(summary_statistics_colnames)
  if(is.null(summary_statistics_colnames)){
    obs<- target_runs %>% dplyr::select(-parameter_colnames) %>%
      data.matrix()
  }else{
    obs<- target_runs %>% dplyr::select(all_of(summary_statistics_colnames)) %>%
      data.matrix()
  }
  progress_bar<- knitrProgressBar::progress_estimated(max(obs %>%nrow(),2))


  #store results here
  results<-data.frame(
    matrix(ncol=length(parameter_colnames),nrow=0)
  )
  colnames(results)<-parameter_colnames

  prediction_residuals<-list()

  #store intervals here
  highs<-data.frame(
    matrix(ncol=length(parameter_colnames),nrow=0)
  )
  lows<-data.frame(
    matrix(ncol=length(parameter_colnames),nrow=0)
  )
  colnames(highs)<-parameter_colnames
  colnames(lows)<-parameter_colnames

  # the x and the y of the training set
  param <-
    training_runs %>%
    dplyr::select(parameter_colnames) %>% data.matrix()

  if(is.null(summary_statistics_colnames)){
    sumstats<- training_runs %>% dplyr::select(-parameter_colnames) %>%
      data.matrix()
  }else{
    sumstats<- training_runs %>% dplyr::select(all_of(summary_statistics_colnames)) %>%
      data.matrix()
  }
  for(i in 1:nrow(target_runs))
  {
    rezult<-NULL

    if(semiauto)
    {
      rezult<-
        abctools::semiauto.abc(obs = obs[i,],
                               param = param ,
                               sumstats = sumstats,
                               method=method,
                               # this comes from the COAL example; regress up to degree 4
                               satr=satr,
                               verbose=FALSE,
                               plot=FALSE,
                               do.err=TRUE,
                               final.dens = TRUE,
                               tol=tol,
                               obspar = obsparam[i,],...)
      low<-rezult$post.sample %>% as.data.frame() %>% summarise_all(quantile,LOW_QUANTILE)
      medians<-rezult$post.sample %>% as.data.frame() %>% summarise_all(median)
      high<-rezult$post.sample %>% as.data.frame() %>% summarise_all(quantile,HIGH_QUANTILE)

    }
    else{
      rezult<-abc::abc(target = obs[i,],
                       param = param,
                       sumstat = sumstats
                       ,tol=tol,method=method,...)
      low<- summary(rezult, print = F)[2, ]
      medians<-summary(rezult, print = F)[3, ] #grabbed from the ABC code (the cv method)
      high<-summary(rezult, print = F)[6,]
    }

    names(medians) <- names(low) <- names(high) <-parameter_colnames
    medians<-data.frame(
      as.list(medians)
    )
    # 3 is median; 4 is mean; 5 is mode.
    results<-bind_rows(results,
                       medians)

    highs<-bind_rows(highs,
                     data.frame(as.list(high)))
    lows<-bind_rows(lows,
                    data.frame(as.list(low)))

    progress_bar$tick()
    progress_bar$print()

  }

  #compute prediction errors
  rmse<-vector(mode="numeric",length=length(parameter_colnames))
  contained<-vector(mode="numeric",length=length(parameter_colnames))
  interval_size<-vector(mode="numeric",length=length(parameter_colnames))

  for(i in 1:length(rmse))
  {
    rmse[i] <- sqrt(mean((obsparam[,i]-results[,i])^2))
    contained[i] <- sum( obsparam[,i] >= lows[,i] & obsparam[,i] <= highs[,i])/length(obsparam[,i])
    interval_size[i] <- mean(highs[,i]-lows[,i])
    prediction_residuals[[parameter_colnames[i]]]<-obsparam[,i]-results[,i]
  }
  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames

  return(
    list(
      medians= results,
      predictions = results, ##alias
      lows = lows,
      highs = highs,
      prediction_residuals = as_tibble(prediction_residuals),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )

}


########################## REJECTION

#' Tries to estimate the parameters through the classic rejection-sampling approximate bayesian computation method. The prediction is the median of the posterior
#' and the interval of prediction is the posterior (cut at 2.5th and 97.5th quantiles)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param \dots other arguments passed to the \code{abc::abc} method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#'fit_rejection_abc(training_runs = training_data,
#'                  target_runs =  testing_data,
#'                  parameter_colnames = c("paramone","paramtwo"),
#'                  summary_statistics_colnames = c("ssone","sstwo")
#')
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'
#'fit_rejection_abc(training_runs = training_data,
#'                  target_runs =  c(2,1),
#'                  parameter_colnames = c("paramone","paramtwo"),
#'                  summary_statistics_colnames = c("ssone","sstwo")
#')
#'## if you feed a full data.frame as target, especially with parameters included, you get error and coverage stuff
#'
#'fit_rejection_abc(training_runs = training_data,
#'                  target_runs =  training_data,
#'                  parameter_colnames = c("paramone","paramtwo"),
#'                  summary_statistics_colnames = c("ssone","sstwo")
#')
fit_rejection_abc<-function(training_runs,target_runs,
                            parameter_colnames,
                            summary_statistics_colnames,tol=DEFAULT_TOL,...){

  abc_testing(training_runs = training_runs,
              target_runs = target_runs,
              method="rejection",
              parameter_colnames = parameter_colnames,
              summary_statistics_colnames = summary_statistics_colnames,
              semiauto = FALSE,
              satr = NULL,
              tol=tol,
              ...)

}

#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_rejection_abc} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param \dots passed to the abc fitting method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_rejection_abc} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_rejection_abc} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_rejection_abc} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
# cv_results<-
# cross_validate_rejection_abc(training_data,ngroup = 5,
# parameter_colnames = c("paramone","paramtwo"),
# summary_statistics_colnames = c("ssone","sstwo"))
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_rejection_abc<-function(total_data,ngroup=5,
                                       parameter_colnames,
                                       summary_statistics_colnames,
                                       cv_seed=0,
                                       tol=DEFAULT_TOL,...){
  .cross_validate(total_data=total_data,
                  ngroup=ngroup,
                  fitting_method=fit_rejection_abc,
                  parameter_colnames = parameter_colnames,
                  summary_statistics_colnames = summary_statistics_colnames,
                  compute_predictivity = TRUE,
                  cv_seed=cv_seed,
                  ...)

}


############################# SEMIAUTO


#' Tries to estimate the parameters through the semi-automatic approximate bayesian computation method: rejection sampling
#' but on a projected space given by a set of linear regressions rather than the original summary statistics space.
#' The prediction is the median of the posterior
#' and the interval of prediction is the posterior (cut at 2.5th and 97.5th quantiles)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param degree the degree of the linear regression that projects summary statistics to the space where ABC then takes place
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param \dots other arguments passed to the \code{abc::abc} method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#' fit_semiauto_abc(training_runs = training_data,
#'                  target_runs =  testing_data,
#'                  parameter_colnames = c("paramone","paramtwo"),
#'                  summary_statistics_colnames = c("ssone","sstwo"),
#'                  degree = 1)
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'
#' fit_semiauto_abc(training_runs = training_data,
#'                  target_runs =  c(2,1),
#'                  parameter_colnames = c("paramone","paramtwo"),
#'                  summary_statistics_colnames = c("ssone","sstwo"),
#'                  degree = 2)
#'## if you feed a full data.frame as target, especially with parameters included, you get error and coverage stuff
#'
#'fit_semiauto_abc(training_runs = training_data,
#'                 target_runs =  training_data,
#'                 parameter_colnames = c("paramone","paramtwo"),
#'                 summary_statistics_colnames = c("ssone","sstwo"),
#'                 degree = 4)
#')
fit_semiauto_abc<-function(training_runs,target_runs,
                           parameter_colnames,
                           summary_statistics_colnames,
                           tol=DEFAULT_TOL,
                           degree = 1,
                           ...){

  abc_testing(training_runs = training_runs,
              target_runs = target_runs,
              method="rejection",
              parameter_colnames = parameter_colnames,
              summary_statistics_colnames = summary_statistics_colnames,
              semiauto = TRUE,
              satr =list(function(x){outer(x,Y=1:degree,"^")}),
              tol=tol,
              ...)

}

#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_semiauto_abc} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param degree the degree of the linear regression that projects summary statistics to the space where ABC then takes place
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param \dots passed to the abc fitting method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_semiauto_abc} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_semiauto_abc} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_semiauto_abc} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#'
#' cv_results<-
#'   cross_validate_semiauto_abc(training_data,ngroup = 5,degree=2,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_semiauto_abc<-function(total_data,ngroup=5,
                                      parameter_colnames,
                                      summary_statistics_colnames,
                                      cv_seed=0,
                                      tol=DEFAULT_TOL,
                                      degree = 1,
                                      ...){
  .cross_validate(total_data=total_data,
                  ngroup=ngroup,
                  fitting_method=fit_semiauto_abc,
                  parameter_colnames = parameter_colnames,
                  summary_statistics_colnames = summary_statistics_colnames,
                  compute_predictivity = TRUE,
                  cv_seed=cv_seed,
                  degree=degree,
                  ...)

}


############################## NEURAL NETWORK

#' Tries to estimate the parameters through the regression-adjusted approximate bayesian computation via neural network
#' The prediction is the median of the posterior
#' and the interval of prediction is the posterior (cut at 2.5th and 97.5th quantiles)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param \dots other arguments passed to the \code{abc::abc} method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#'fit_neural_network_abc(training_runs = training_data,
#'                       target_runs =  testing_data,
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'
#'fit_neural_network_abc(training_runs = training_data,
#'                       target_runs =  c(2,1),
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame as target, especially with parameters included, you get error and coverage stuff
#'
#'fit_neural_network_abc(training_runs = training_data,
#'                       target_runs =  training_data,
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'
fit_neural_network_abc<-function(training_runs,target_runs,
                                 parameter_colnames,
                                 summary_statistics_colnames,
                                 tol=DEFAULT_TOL,
                                 ...){

  abc_testing(training_runs = training_runs,
              target_runs = target_runs,
              method="neuralnet",
              parameter_colnames = parameter_colnames,
              summary_statistics_colnames = summary_statistics_colnames,
              semiauto = FALSE, hcor=TRUE,
              satr =NULL,
              tol=tol,
              ...)

}



#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_neural_network_abc} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param \dots passed to the abc fitting method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_neural_network_abc} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_neural_network_abc} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_neural_network_abc} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#'
#'cv_results<-
#'  cross_validate_neural_network_abc(training_data,ngroup = 5,
#'                                    parameter_colnames = c("paramone","paramtwo"),
#'                                    summary_statistics_colnames = c("ssone","sstwo"))
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_neural_network_abc<-function(total_data,ngroup=5,
                                            parameter_colnames,
                                            summary_statistics_colnames,
                                            cv_seed=0,
                                            tol=DEFAULT_TOL,...){
  .cross_validate(total_data=total_data,
                  ngroup=ngroup,
                  fitting_method=fit_neural_network_abc,
                  parameter_colnames = parameter_colnames,
                  summary_statistics_colnames = summary_statistics_colnames,
                  compute_predictivity = TRUE,
                  cv_seed=cv_seed,
                  ...)

}

############################## loclinear


#' Tries to estimate the parameters through the regression-adjusted approximate bayesian computation via loclinear regression
#' The prediction is the median of the posterior
#' and the interval of prediction is the posterior (cut at 2.5th and 97.5th quantiles)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param \dots other arguments passed to the \code{abc::abc} method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#'fit_loclinear_abc(training_runs = training_data,
#'                       target_runs =  testing_data,
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'
#'fit_loclinear_abc(training_runs = training_data,
#'                       target_runs =  c(2,1),
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame as target, especially with parameters included, you get error and coverage stuff
#'
#'fit_loclinear_abc(training_runs = training_data,
#'                       target_runs =  training_data,
#'                       parameter_colnames = c("paramone","paramtwo"),
#'                       summary_statistics_colnames = c("ssone","sstwo"))
#'
fit_loclinear_abc<-function(training_runs,target_runs,
                            parameter_colnames,
                            summary_statistics_colnames,
                            tol=DEFAULT_TOL,
                            ...){

  abc_testing(training_runs = training_runs,
              target_runs = target_runs,
              method="loclinear",
              parameter_colnames = parameter_colnames,
              summary_statistics_colnames = summary_statistics_colnames,
              semiauto = FALSE, hcor=TRUE,
              satr =NULL,
              tol=tol,
              transf = "none",
              ...)

}

#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_loclinear_abc} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param tol tolerance error for ABC: in practice what % of the original runs end up forming the posterior (default is 0.1)
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param \dots passed to the abc fitting method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_loclinear_abc} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_loclinear_abc} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_loclinear_abc} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#'
#'cv_results<-
#'  cross_validate_loclinear_abc(training_data,ngroup = 5,
#'                                    parameter_colnames = c("paramone","paramtwo"),
#'                                    summary_statistics_colnames = c("ssone","sstwo"))
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_loclinear_abc<-function(total_data,ngroup=5,
                                       parameter_colnames,
                                       summary_statistics_colnames,
                                       cv_seed=0,
                                       tol=DEFAULT_TOL,...){
  .cross_validate(total_data=total_data,
                  ngroup=ngroup,
                  fitting_method=fit_loclinear_abc,
                  parameter_colnames = parameter_colnames,
                  summary_statistics_colnames = summary_statistics_colnames,
                  compute_predictivity = TRUE,
                  cv_seed=cv_seed,
                  ...)

}


#loclinear

# ### FACADE that calls cross_validate with linear regression fit
# cross_validate_abc<-function(total_data,ngroup,tol=DEFAULT_TOL,parameter_colnames,
#                              method,semiauto=FALSE,cv_seed=0,
#                              satr=list(function(x){outer(x,Y=1:4,"^")}),
#                              ...){
#   .cross_validate(total_data=total_data,cv_seed=cv_seed,
#                  ngroup=ngroup,
#                  fitting_method=abc_testing,
#                  tol=tol,
#                  parameter_colnames=parameter_colnames,
#                  method=method,
#                  semiauto=semiauto,
#                  satr=satr)
# }

#   _____ _   _  ___   _   _ _____ _____ _      _____  ____________
#  |  _  | | | |/ _ \ | \ | |_   _|_   _| |    |  ___| | ___ \  ___|
#  | | | | | | / /_\ \|  \| | | |   | | | |    | |__   | |_/ / |_
#  | | | | | | |  _  || . ` | | |   | | | |    |  __|  |    /|  _|
#  \ \/' / |_| | | | || |\  | | |  _| |_| |____| |___  | |\ \| |
#   \_/\_\\___/\_| |_/\_| \_/ \_/  \___/\_____/\____/  \_| \_\_|
#
#




#' Tries to estimate the parameters through a quantile random forest and using the median as prediction and the 2.5th-97.5th
#' quantile as prediction intervals
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
# fit_quantile_random_forest(training_runs = training_data,
#                      target_runs =  testing_data,
#                      parameter_colnames = c("paramone","paramtwo"),
#                      summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'fit_quantile_random_forest(training_runs = training_data,
#'                      target_runs =  c(2,1),
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame, especially with parameters included, you get error and coverage stuff
#'fit_quantile_random_forest(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
fit_quantile_random_forest<-function(
  training_runs,target_runs,parameter_colnames,
  summary_statistics_colnames # name of the X variables in the regression (summary statistics)

){


  ##you are allowed to send summary statistics into a vector if you are just fitting once
  target_runs<-.check_that_summary_statistics_are_into_dataframe(target_runs,summary_statistics_colnames)
  ## check if the target runs need NAs added (happens when we aren't testing but properly fitting!)
  target_runs<-.add_NAs_to_target_runs_ifneeded(target_runs,parameter_colnames)

  #abctools works on matrices, not data.frames so we need to convert
  regressions<-list()
  rmse<-rep.int(NA,length(parameter_colnames))
  interval_size<-rep.int(NA,length(parameter_colnames))
  contained<-rep.int(NA,length(parameter_colnames))
  predictions<-list()
  lows<-list()
  highs<-list()
  prediction_residuals<-list()
  for(i in 1:length(parameter_colnames)){
    y<-parameter_colnames[i]
    formula1<-.formula_maker(yname = y,
                             xnames = summary_statistics_colnames,
                             degree = 1)

    framed<-model.frame(formula1,data=training_runs)
    test_frame<-model.frame(formula1,data=target_runs,na.action =stats::na.pass )

    #run the regressions
    train_X<-
      framed[,-1]
    train_Y<-
      framed[,1]

    test_X<-
      test_frame[,-1]
    test_Y<-
      test_frame[,1]

    if(is.null(ncol(train_X)))
      train_X <- train_X %>% as.data.frame()

    regression<- quantregForest::quantregForest(x =train_X,y=train_Y )



    if(is.null(ncol(test_X)))
      test_X <- test_X %>% as.data.frame()

    predictions[[y]]<- predict(regression,  test_X %>% as.data.frame(), what=0.5)
    rmse[i]<-  sqrt(mean((predictions[[y]]-test_Y)^2))
    prediction_residuals[[y]]<-test_Y - predictions[[y]]

    lows[[y]]<- predict(regression,  test_X, what=LOW_QUANTILE)
    highs[[y]]<- predict(regression,  test_X, what=HIGH_QUANTILE)


    contained[i]<- sum( (test_Y>=lows[[y]]) & (test_Y<=highs[[y]]))/length(predictions[[y]])
    interval_size[i]<-mean(highs[[y]]-lows[[y]])
  }



  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames


  return(
    list(
      predictions= as_tibble(predictions),
      lows = as_tibble(lows),
      highs = as_tibble(highs),
      prediction_residuals= as_tibble(prediction_residuals),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )



}


#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_quantile_random_forest} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_quantile_random_forest} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_quantile_random_forest} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_quantile_random_forest} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#' cv_results<-cross_validate_quantile_random_forest(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#'
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_quantile_random_forest<-function(total_data,ngroup=5,
                                                parameter_colnames,summary_statistics_colnames,
                                                cv_seed=0){
  .cross_validate(total_data=total_data,cv_seed=cv_seed,
                  ngroup=ngroup,
                  fitting_method=fit_quantile_random_forest,
                  parameter_colnames=parameter_colnames,
                  summary_statistics_colnames=summary_statistics_colnames)
}

#  ______                _                  ______                  _
#  | ___ \              | |                 |  ___|                | |
#  | |_/ /__ _ _ __   __| | ___  _ __ ___   | |_ ___  _ __ ___  ___| |_
#  |    // _` | '_ \ / _` |/ _ \| '_ ` _ \  |  _/ _ \| '__/ _ \/ __| __|
#  | |\ \ (_| | | | | (_| | (_) | | | | | | | || (_) | | |  __/\__ \ |_
#  \_| \_\__,_|_| |_|\__,_|\___/|_| |_| |_| \_| \___/|_|  \___||___/\__|
#
#



#' Tries to estimate the parameters through a classic random forest and using bootstrapped residuals + infinitesimal jacknife standard error
#' to build prediction intervals. Internally it uses caret and ranger packages for some self-calibration.
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param fast when set to true the fitting will not try to find the best hyper-parameters through the \code{caret} package and will
#' instead use the default \code{ranger} parameters. This is FALSE by default to replicate the paper's result, but in my experience
#' setting it to TRUE speeds up the process tenfold for a minor hit in performance.
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
# fit_random_forest(training_runs = training_data,
#                      target_runs =  testing_data,
#                      parameter_colnames = c("paramone","paramtwo"),
#                      summary_statistics_colnames = c("ssone","sstwo"), fast=TRUE)
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'fit_random_forest(training_runs = training_data,
#'                      target_runs =  c(2,1),
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"),fast=FALSE)
#'## if you feed a full data.frame, especially with parameters included, you get error and coverage stuff
#'fit_random_forest(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"),fast=FALSE)
fit_random_forest<-function(
  training_runs,target_runs,parameter_colnames,
  summary_statistics_colnames, # name of the X variables in the regression (summary statistics)
  fast=FALSE
){
  ##you are allowed to send summary statistics into a vector if you are just fitting once
  target_runs<-.check_that_summary_statistics_are_into_dataframe(target_runs,summary_statistics_colnames)
  ## check if the target runs need NAs added (happens when we aren't testing but properly fitting!)
  target_runs<-.add_NAs_to_target_runs_ifneeded(target_runs,parameter_colnames)


  #abctools works on matrices, not data.frames so we need to convert
  regressions<-list()
  rmse<-rep.int(NA,length(parameter_colnames))
  interval_size<-rep.int(NA,length(parameter_colnames))
  contained<-rep.int(NA,length(parameter_colnames))
  predictions<-list()
  lows<-list()
  highs<-list()

  prediction_residuals<-list()


  if(length(parameter_colnames)>1)
  {
    progress_bar<- knitrProgressBar::progress_estimated(length(parameter_colnames))
    #knitrProgressBar::update_progress()
    #progress_bar$print()
  }

  for(i in 1:length(parameter_colnames)){
    y<-parameter_colnames[i]
    formula1<-.formula_maker(yname = y,
                             xnames = summary_statistics_colnames,
                             degree = 1)


    if(fast == FALSE & length(summary_statistics_colnames)>1)
    {
      #use caret to find out the best hyper-parameters
      regression<- caret::train(formula1,data=training_runs,method="ranger")

      #ugh, run it once more (with best hyper-parameters!)
      regression<- ranger::ranger(formula1,data=training_runs,
                                  keep.inbag = TRUE,
                                  mtry = regression$bestTune$mtry,
                                  splitrule = regression$bestTune$splitrule,
                                  min.node.size = regression$bestTune$min.node.size)
    }
    else{
      regression<- ranger::ranger(formula1,data=training_runs,
                                  keep.inbag = TRUE)
    }

    #matrix with rows--> X and columns--> tree prediction at that X
    ## may complain about lack of calibration even with calibrate set to false
    suppressWarnings(
      best_prediction<-predict(regression,  data=target_runs,
                               type="se",se.method="infjack",calibrate=FALSE)
    )
    #if infinitesimal jacknife fails, we have to switch to plain jackkniffe (but pay the computational cost!)
    if(sum(is.na(best_prediction$se))>0)
      best_prediction<-predict(regression,  data=target_runs,
                               type="se",se.method="jack")
    #get the residuals --> OUT OF BAG!

    residuals<-training_runs[[y]]-regression$predictions


    test <- best_prediction$se %>%
      as.list() %>%
      map(~rnorm(mean=0,sd=.,n=1000)+ sample(residuals,size = 1000,replace=TRUE))

    low<- best_prediction$predictions +
      test %>% map_dbl(~quantile(.,probs=c(0.025)))
    high<- best_prediction$predictions +
      test %>% map_dbl(~quantile(.,probs=c(0.975)))


    test_Y<-
      target_runs[[y]]


    predictions[[y]]<- best_prediction$predictions
    rmse[i]<-  sqrt(mean((predictions[[y]]-test_Y)^2))
    prediction_residuals[[y]]<-test_Y - predictions[[y]]

    lows[[y]]<- low
    highs[[y]]<- high


    contained[i]<- sum( (test_Y>=lows[[y]]) & (test_Y<=highs[[y]]))/length(
      predictions[[y]])
    interval_size[i]<-mean(highs[[y]]-lows[[y]])

    if(length(parameter_colnames)>1)
    {
      progress_bar$tick()
      progress_bar$print()
    }

  }



  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames


  return(
    list(
      predictions= as_tibble(predictions),
      lows = as_tibble(lows),
      highs = as_tibble(highs),
      prediction_residuals = as_tibble(prediction_residuals),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )



}





#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_random_forest} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param fast when set to true the fitting will not try to find the best hyper-parameters through the \code{caret} package and will
#' instead use the default \code{ranger} parameters. This is FALSE by default to replicate the paper's result, but in my experience
#' setting it to TRUE speeds up the process tenfold for a minor hit in performance.
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_random_forest} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_random_forest} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_random_forest} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#' cv_results<-cross_validate_random_forest(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#'
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_random_forest<-function(total_data,ngroup=5,
                                       parameter_colnames,summary_statistics_colnames,
                                       cv_seed=0,fast=FALSE){
  .cross_validate(total_data=total_data,cv_seed=cv_seed,
                  ngroup=ngroup,
                  fitting_method=fit_random_forest,
                  parameter_colnames=parameter_colnames,
                  summary_statistics_colnames=summary_statistics_colnames,
                  fast=fast)
}




#   _____    ___   ___  ___
#  |  __ \  / _ \  |  \/  |
#  | |  \/ / /_\ \ | .  . |
#  | | __  |  _  | | |\/| |
#  | |_\ \ | | | | | |  | |
#   \____/ \_| |_/ \_|  |_/
#
#

#' Tries to estimate the parameters through a generalized additive model and using bootstrapped residuals + estimated standard error
#' to build prediction intervals. Internally it uses mgcv package.
#'
#' It may run into numerical issues with a very large number of summary statistics; if so setting "bam=TRUE" helps.
#' It may also automatically ignore summary statistics for which you don't get enough variation (>10 unique values)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param bam defauls to FALSE, when it's true it switches optimizer mgcv uses (akin to using mgcv:bam AND setting discrete to True) to
#' something more robust.
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#' fit_gam(training_runs = training_data,
#'                      target_runs =  testing_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'fit_gam(training_runs = training_data,
#'                      target_runs =  c(2,1), bam=TRUE,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame, especially with parameters included, you get error and coverage stuff
#'fit_gam(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
fit_gam<-function(
  training_runs,target_runs,parameter_colnames,
  summary_statistics_colnames, # name of the X variables in the regression (summary statistics)
  bam=FALSE # should we use BAM instead of GAM? Useful for big data
){
  # print(bam)

  ##you are allowed to send summary statistics into a vector if you are just fitting once
  target_runs<-.check_that_summary_statistics_are_into_dataframe(target_runs,summary_statistics_colnames)
  ## check if the target runs need NAs added (happens when we aren't testing but properly fitting!)
  target_runs<-.add_NAs_to_target_runs_ifneeded(target_runs,parameter_colnames)


  regressions<-list()
  rmse<-rep.int(NA,length(parameter_colnames))
  interval_size<-rep.int(NA,length(parameter_colnames))
  contained<-rep.int(NA,length(parameter_colnames))
  predictions<-list()
  lows<-list()
  highs<-list()
  prediction_residuals<-list()
  if(length(parameter_colnames)>1)
  {
    progress_bar<- knitrProgressBar::progress_estimated(length(parameter_colnames))
    progress_bar$print()
  }

  for(i in 1:length(parameter_colnames)){
    y<-parameter_colnames[i]
    #  message(paste("fitting parameter",y))

    #we need to avoid putting too high k on variables whose unique values are only a few
    unique_values<-
      training_runs %>% select(-all_of(parameter_colnames)) %>% gather(variable,value) %>% group_by(variable) %>%
      summarise(uniques=length(unique(value)),.groups="keep")  %>% filter(uniques<10)

    #with more than 10 unique observations, they go in the formula directly
    valids<-setdiff(summary_statistics_colnames,unique_values$variable)
    #adding s makes it non-parametric; caret does it for us when we call train
    #but in the bootstrap phase we just call gam directly
    gam_formula<-
      paste(
        y,"~","s(",
        paste(valids,collapse = ")+s("),")"

      )
    #with 1 unique observation, they don't go in the formula at all
    unique_values<-unique_values %>% filter(uniques>1)
    if(nrow(unique_values)>0)
      gam_formula<-paste(gam_formula,"+",
                         paste("s(",unique_values$variable,",k=",pmin(10,unique_values$uniques),")",collapse = "+")
      )

    gam_formula<-formula(gam_formula)


    if(!bam)
    {
      regression<- mgcv::gam(gam_formula,data=training_runs,select=TRUE)
    } else{
      regression<- mgcv::bam(gam_formula,data=training_runs,select=TRUE,discrete = TRUE)

    }
    residuals<-residuals(regression)

    current_prediction<-predict(regression,  newdata=target_runs,se=T)
    predictions[[y]]<-  current_prediction$fit

    rmse[i]<-  sqrt(mean((predictions[[y]]-target_runs[[y]])^2))
    prediction_residuals[[y]]<-target_runs[[y]] - predictions[[y]]
    residuals<-residuals(regression)


    test <- current_prediction$se %>%
      as.list() %>%
      map(~rnorm(mean=0,sd=.,n=10000)+ sample(residuals,size = 10000,replace=TRUE))

    low<- current_prediction$fit +
      test %>% map_dbl(~quantile(.,probs=c(0.025)))
    high<- current_prediction$fit +
      test %>% map_dbl(~quantile(.,probs=c(0.975)))


    #should be done!
    lows[[y]]<- low
    highs[[y]]<- high

    contained[i]<- sum( (target_runs[[y]]>=lows[[y]]) & (target_runs[[y]]<=highs[[y]]))/length(predictions[[y]])
    interval_size[i]<-mean(highs[[y]]-lows[[y]])

    if(length(parameter_colnames)>1)
    {
      knitrProgressBar::update_progress(progress_bar)
      # progress_bar$print()
    }
  }



  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames
  contained

  return(
    list(
      predictions= as_tibble(predictions),
      lows = as_tibble(lows),
      highs = as_tibble(highs),
      prediction_residuals = as_tibble(prediction_residuals),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )



}

#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_gam} and the last group as \code{target_runs}. Then repeats for all groups.
#' Finally it collects performance and predictivity to get a good handle on identification and estimation quality
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param bam defauls to FALSE, when it's true it switches optimizer mgcv uses (akin to using mgcv:bam AND setting discrete to True) to
#' something more robust.
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @param \dots passed to the fitting method
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item performance: 1 - the ratio of out-of-sample RMSE of the \code{fit_gam} and RMSE of just using the average (named vector of length equal to the number of parameters)
#'   \item predictivity: 1 - sum of squared residuals of the \code{fit_gam} and sum of squared residuals when using just the average value (named vector of length equal to the number of parameters)
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_gam} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#' cv_results<-cross_validate_gam(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#'
#' ## if we look at performance we can see that paramone has been identified and paramtwo has not:
#' ## (in the paper I use either 0.1 or 0.3 as minimum performance below which you just didn't identify the parameter)
#' cv_results$performance
#' ## different numbers but same story looking at predictivity
#' cv_results$predictivity
#' ## however it looks like coverage is almost perfect (the 95% interval contains out-of-sample the real parameters 94.X% of the time)
#' cv_results$contained
cross_validate_gam<-function(total_data,ngroup=5,
                             parameter_colnames,summary_statistics_colnames,
                             bam=FALSE,
                             cv_seed=0,...){
  .cross_validate(total_data=total_data,cv_seed=cv_seed,
                  ngroup=ngroup,
                  fitting_method=fit_gam,
                  parameter_colnames=parameter_colnames, bam=bam,
                  summary_statistics_colnames=summary_statistics_colnames,...)
}



#' This method is useful for other cross-validation methods to compute predictivity/performance. On its own it's pretty pointless.
#'
#' Performs a full-cross validation where it splits all the runs in \code{total_data}
#' into a set of groups, uses all but one group as \code{training_runs} in
#' \code{fit_average} and the last group as \code{target_runs}. Then repeats for all groups.
#'
#' It doesn't collect any performance indicator since by construction here they all end up being 0.
#'
#' @param total_data A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param ngroup Number of groups to split the total data in when cross-validating. Defaults to 5
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param cv_seed random seed controlling how CV groups are formed: keep it constant to compare cross-validations across methods
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item rmse - Average out of sample RMSE for all the groups in the CV  (named vector of length equal to the number of parameters)
#'   \item contained - The out-of-sample percentage of times the real parameters are contained in the estimated interval;
#'   \item interval_size - The average range between lower and upper bound of the prediction intervals for each parameter, i.e. how wide our confidence bounds are on average (named vector of length equal to the number of parameters)
#'   \item results list of \code{fit_gam} results, one for each CV group
#'   \item out_of_sample_errors a tidy data-frame containing for each row an out-of-sample prediction and error made for one parameter; useful for debugging
#' }
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#' ### or simply do a full cross-validation
#' cross_validate_average<-cross_validate_gam(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
#'
cross_validate_average<-function(total_data,ngroup=5,
                                 parameter_colnames,summary_statistics_colnames,
                                 cv_seed=0){
  .cross_validate(total_data=total_data,cv_seed=cv_seed,
                  ngroup=ngroup,
                  fitting_method=fit_average,
                  parameter_colnames=parameter_colnames,
                  summary_statistics_colnames=summary_statistics_colnames)
}




#' "Estimates" parameters by just always returning the average observed in training. Not useful at all except as a baseline to compare
#' other policies (if they do not do better than this, they haven't identified anything)
#'
#' @param training_runs A data.frame where each row represents a separate simulation run. The data.frame ought to
#' contain both the summary statistics (output of the model) and the parameters that generated them.
#' @param target_runs The "real" summary statistics for which we are trying to estimate which parameters generated them.
#' Target runs can be a data.frame with one or more rows of summary statistics and optionally also the generating parameters if we
#' are testing the estimation method. Alternatively a vector of length equal to the number of summary statistics
#' @param parameter_colnames The name of the columns in the training_runs data.frame that represent the parameters (inputs) of the
#' simulation model
#' @param summary_statistics_colnames The name of the columns in the training_runs data.frame that represent the summary statistics
#' (outputs) of the simulation model
#' @param bam defauls to FALSE, when it's true it switches optimizer mgcv uses (akin to using mgcv:bam AND setting discrete to True) to
#' something more robust.
#' @return A list containing the basic estimation results
#' #' \itemize{
#'   \item predictions - The estimated values for each parameter (data.frame with as many rows as in \code{target_runs})
#'   \item lows - Estimated lower bound (2.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item highs - Estimated upper bound (97.5th quantile) for each estimated parameter (data.frame with as many rows as in \code{target_runs})
#'   \item prediction_residuals - Difference between "real" parameters and estimated ones for each row in the \code{target_runs};
#'   if target_runs didn't come with parameters, then only NAs are shown (data.frame with as many rows as in \code{target_runs})
#'   \item rmse - The root mean square error of the estimated parameter compared to the real one; only computed if target_runs also
#'   provided the real parameters otherwise NAs (named vector of length equal to the number of parameters)
#'   \item contained - The percentage of times the real parameters are contained in the estimated interval; NA if \code{target_runs}
#'   doesn't provide the real parameters (named vector of length equal to the number of parameters)
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### classic arguments: a data.frame for training and one for testing plus the column names
#' fit_average(training_runs = training_data,
#'                      target_runs =  testing_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## you can fit just a vector of summary statistics (still need to provide the names though)
#'fit_average(training_runs = training_data,
#'                      target_runs =  c(2,1), bam=TRUE,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## if you feed a full data.frame, especially with parameters included, you get error and coverage stuff
#'fit_average(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
fit_average<-function(
  training_runs,target_runs,parameter_colnames,
  summary_statistics_colnames # name of the X variables in the regression (summary statistics)
){


  #abctools works on matrices, not data.frames so we need to convert
  rmse<-rep.int(NA,length(parameter_colnames))
  interval_size<-rep.int(NA,length(parameter_colnames))
  contained<-rep.int(NA,length(parameter_colnames))
  predictions<-list()
  lows<-list()
  highs<-list()
  prediction_residuals<-list()

  for(i in 1:length(parameter_colnames)){
    y<-parameter_colnames[i]


    predictions[[y]]<-  mean(training_runs[[y]])

    rmse[i]<-  sqrt(mean((predictions[[y]]-target_runs[[y]])^2))
    prediction_residuals[[y]]<- target_runs[[y]] - predictions[[y]]


    #should be done!
    lows[[y]]<- predictions[[y]]  - 1.96 * sd(training_runs[[y]])
    highs[[y]]<- predictions[[y]]  + 1.96 * sd(training_runs[[y]])


    contained[i]<- sum( (target_runs[[y]]>=lows[[y]]) & (target_runs[[y]]<=highs[[y]]))/length(target_runs[[y]])
    interval_size[i]<-mean(highs[[y]]-lows[[y]])
  }



  names(rmse)<- names(contained) <- names(interval_size) <-parameter_colnames


  return(
    list(
      predictions= as_tibble(predictions),
      prediction_residuals = as.data.frame(prediction_residuals),
      lows = as_tibble(lows),
      highs = as_tibble(highs),
      rmse = rmse,
      contained = contained,
      interval_size = interval_size
    )
  )
}


#' Helper function to produce a tidy data frame containing upper,lower,estimate and real value for each parameter
#' in each simulated run. Useful for plotting
#'
#' @param estimation the output of one of the \code{fit_} or \code{cross_validate_} calls in this package
#' @return A data-frame containing the following columns
#' #' \itemize{
#'   \item type whether the row is the estimated parameter, its bounds or the real value
#'   \item run the simulation run index
#'   \item variable which parameter are we estimating
#'   \item value the estimated value (NA if we don't know the REAL value of a parameter when type is "real")
#' }
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### fit a gam
#'estimation<- fit_gam(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'
#'## we can clean it to be more presentable
#'tidy_up_estimation(estimation)
#'
#'## we can do the same with a cv result
#'cv_results<-
#'cross_validate_rejection_abc(training_data,ngroup = 5,
#'                             parameter_colnames = c("paramone","paramtwo"),
#'                             summary_statistics_colnames = c("ssone","sstwo"))
#' tidy_up_estimation(cv_results)
tidy_up_estimation<-function(estimation){

  ## check if it's a cross-validation error (in which case you need to pre-process)
  if("results" %in% names(estimation)){
    ## it's cross validation!
    fit_result<-
      list(
        lows= map_dfr(estimation$results,~.$lows),
        highs= map_dfr(estimation$results,~.$highs),
        predictions=map_dfr(estimation$results,~.$predictions),
        prediction_residuals=map_dfr(estimation$results,~.$prediction_residuals)
      )

  } else{
    ## no need to pre-process
    fit_result<-estimation
  }
  return(
    bind_rows(fit_result$lows %>%
                mutate(type="lower_bound",id=row_number()),
              fit_result$predictions %>%
                mutate(type="estimate",id=row_number()),
              (fit_result$predictions + fit_result$prediction_residuals) %>%
                mutate(type="real",id=row_number()),
              fit_result$highs %>% mutate(type="higher_bound",id=row_number())) %>%
      gather("variable","value",-type,-id)  %>% spread(type,value) %>%
      rename(run=id) %>%
      arrange(run,variable)
  )


}

#' Plotting function that samples at random a few estimated simulation runs
#' and shows an error plot highlighting
#' the size of the 95th confidence interval (black error bars) and whether the real value (red dots) are in
#' them.
#'
#' Useful to get a quick overview about the quality and size of the confidence intervals
#'
#' @param estimation the output of one of the \code{fit_} or \code{cross_validate_} calls in this package
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### fit a gam
#'estimation<- fit_gam(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## we can check how the confidence intervals look like on a subset of runs
#'## because target_runs were just the training_data again this is showing IN-SAMPLE CI
#'plot_confidence_intervals(estimation)
#'## look in particular at how non-committal  paramtwo CIs are! Clearly GAM knows
#'## that it is NOT identifying that parameter
#'
#'## we can do the same plot for cross-validations results; which
#'## would make it OUT-OF-SAMPLE
#'cv_results<-
#'  cross_validate_rejection_abc(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
# plot_confidence_intervals(cv_results,
#                          limit_number_of_runs_to_plot_to = 95,
#                          random_seed=1) + theme_minimal()
# wider CIs with ABC, but not much difference besides that!
plot_confidence_intervals<-function(estimation,
                                    limit_number_of_runs_to_plot_to=100,
                                    random_seed=0){



  tidied<-tidy_up_estimation(estimation)

  ##remove NAs
  tidied<-tidied %>% na.omit()
  if(nrow(tidied)==0)stop("There are no real values provided; can't plot")


  runs_to_plot<- o <- R.utils::withSeed({
    sample(tidied$run,size=
             min(limit_number_of_runs_to_plot_to,
                 max(tidied$run)))
  }, seed=random_seed)



  return(
    tidied %>% group_by(variable) %>%
      filter(run %in% runs_to_plot) %>% mutate(run=rank(real,ties.method = "first")) %>%
      ggplot(aes(x=run)) +
      geom_errorbar(aes(ymin=lower_bound,ymax=higher_bound),alpha=0.8) +
      geom_point(aes(y=real),col="red",size=2) +
      facet_wrap(variable~.,scales = "free") +
      xlab("") + ylab("Parameter")
  )

}

#' Produces a scatterplot of real parameters (y axis) and estimated parameters (x axis) to quickly
#' look at quality of point prediction. A perfect method would have all these dots on the 45 degree line
#' which is highlighted here with a red dashed line.
#'
#' I find this plot particularly useful to spot partial identifications: when estimation is possible for
#' a sub-interval of the parameter range
#' #'
#' @param estimation the output of one of the \code{fit_} or \code{cross_validate_} calls in this package
#' @return A list of ggplots
#' @export
#'
#' @examples
#' ##generate some fake data where paramone,paramtwo ---> ssone,sswto;
#' ## notice that paramtwo is basically unidentifiable!
#' paramone<-rnorm(n=5000)
#' paramtwo<-runif(n=5000,min=2,max=5)
#' ssone<-2*paramone + rnorm(n=5000)
#' sstwo<- paramone/paramtwo  + rnorm(n=5000)
#' training_data<-
#'   data.frame(
#'     paramone,
#'     paramtwo,
#'     ssone,
#'     sstwo
#'   )
#'   ## this would be the "real" data, what we want to estimate our model with!
#' testing_data<-data.frame(
#'   ssone=2,
#'  sstwo=0.25
#')
#'
#'
#'### fit a gam
#'estimation<- fit_gam(training_runs = training_data,
#'                      target_runs =  training_data,
#'                      parameter_colnames = c("paramone","paramtwo"),
#'                      summary_statistics_colnames = c("ssone","sstwo"))
#'## we can check how the prediction match real parameters
#'## because target_runs were just the training_data again this is showing IN-SAMPLE errors
#'plot_point_prediction_quality(estimation)
#'## notice that basically GAM for paramtwo just returns the average (i.e. "I don't know!")
#'
#'## we can theme it, but we need to map (since this is patchwork output)
# to_theme<-plot_point_prediction_quality(estimation)
# to_theme[[1]] <-to_theme[[1]] + theme_light()
# to_theme[[2]] <-to_theme[[2]] + theme_dark()
# to_theme
#'## we can do the same plot for cross-validations results; which
#'## would make it OUT-OF-SAMPLE
#'cv_results<-
#'  cross_validate_rejection_abc(training_data,ngroup = 5,
#'                               parameter_colnames = c("paramone","paramtwo"),
#'                               summary_statistics_colnames = c("ssone","sstwo"))
# plot_point_prediction_quality(cv_results)
# very similar as above
plot_point_prediction_quality<-function(estimation){
  tidied<-tidy_up_estimation(estimation)
  parameters<- tidied %>% pull(variable) %>% unique()
  plot_list<-map(parameters,
                 function(x) {
                   filtered<- tidied %>% filter(variable==x)

                   small_x<- min(filtered$real)
                   high_x<-max(filtered$real)

                   filtered %>%
                     ggplot(aes(x=estimate,y=real)) + ylab("Real value") + xlab("Estimate") +
                     geom_point() +
                     geom_abline(slope=1,lwd=2,lty=2,col="red",alpha=0.7) +
                     xlim(small_x,high_x) + ylim(small_x,high_x) +ggtitle(x)
                 }
  )
  patchwork::wrap_plots(plot_list)
}

## todo: would like to add mxnet fits, but seems impossible to install
## would like to add BACCO fits, but no idea on how to choose GP hyper-pameters
