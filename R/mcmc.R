#' importFrom dplyr "%>%"


#' @title Adaptive MCMC algorithm
#'
#' @description MCMC which adapts its proposal distribution for faster convergence following:
#' Sherlock, C., Fearnhead, P. and Roberts, G.O. The Random Walk Metrolopois: Linking Theory and Practice Through a Case Study. Statistical Science 25, no.2 (2010): 172-190.
#'
#' @param lprior A function returning the log prior probability of the parameters 
#' @param llikelihood A function returning the log likelihood of the parameters given the data
#' @param nburn Number of iterations of burn in
#' @param initial Vector with starting parameter values
#' @param nbatch Number of batches to run (number of samples to return)
#' @param blen Length of each batch
#' @param outfun A function that is called for each batch. Can be useful to log certain values. 
#' @param acceptfun A function that is called whenever a sample is accepted. 
#' @param verbose Output debugging information
#' @param ... Extra parameters passed to the log likelihood function
#' 
#' @return Returns a list with the accepted samples and the corresponding llikelihood values and the return of the optional outfun
#'
#' @seealso \code{\link{adaptive.mcmc.cpp}} Used internally by this function.
adaptive.mcmc <- function(lprior, llikelihood, nburn, 
                          initial, nbatch, blen = 1, outfun = NULL, 
                          acceptfun = NULL, verbose = FALSE, ...)
{
  if (is.null(outfun))
    outfun <- function() { NULL }
  if (is.null(acceptfun))
    acceptfun <- function() { NULL }
  
  adaptive.mcmc.cpp(lprior, function(pars) llikelihood(pars, ...), outfun,
                    acceptfun, nburn, initial, nbatch, blen, verbose)
}


#' MCMC based inference of the parameter values given the different data sets
#'
#' @param demography A vector with the population size by each age {0,1,..}
#' @param ili The number of Influenza-like illness cases per week
#' @param mon_pop The number of people monitored for ili
#' @param n_pos The number of positive samples for the given strain (per week)
#' @param n_samples The total number of samples tested 
#' @param vaccine_calendar A vaccine calendar valid for that year
#' @param polymod_data Contact data for different age groups
#' @param initial Vector with starting parameter values
#' @param mapping Optional age group mapping from model age groups to data age groups
#' @param nburn Number of iterations of burn in
#' @param nbatch Number of batches to run (number of samples to return)
#' @param blen Length of each batch
#' 
#' @return Returns a list with the accepted samples and the corresponding llikelihood values and a matrix (contact.ids) containing the ids (row number) of the contacts data used to build the contact matrix.
#'
#' @export
inference <- function(demography, ili, mon_pop, n_pos, n_samples, 
        vaccine_calendar, polymod_data, initial, mapping,
        risk_ratios, nburn = 0, nbatch = 1000, blen = 1 )
{
  warning("Should accept age_group_mapping, risk_group_mapping and combine it (and sort it). Then pass it to cpp as group_mapping matrix, with each row a to, from and weight.")
  if (missing(mapping))
    mapping <- age_group_mapping(c(1,5,15,25,45,65), c(5,15,45,65))
  if (missing(risk_ratios)) {
    risk_ratios <- matrix(c(
      0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
      0, 0, 0, 0, 0, 0, 0                          
    ), ncol = 7, byrow = T)
  }
  if (class(risk_ratios) == "matrix") {
    no_risk_groups <- nrow(risk_ratios) + 1
    rv <- c(rep(1, ncol(risk_ratios)) - colSums(risk_ratios), t(risk_ratios))
    risk_ratios <- data.frame(
      AgeGroup = rep(unique(mapping$from), no_risk_groups),
      value = rv
    ) %>% dplyr::group_by(AgeGroup) %>% dplyr::mutate(RiskGroup = factor(row_number())) %>%
      dplyr::ungroup()
  } else if (class(risk_ratios) != "data.frame") {
    risk_ratios <- data.frame(value = risk_ratios) 
  }
  if (missing(no_risk_groups))
    no_risk_groups <- nrow(risk_ratios)/unique(mapping$from)
  
  .inference_cpp(demography, ili, mon_pop, n_pos, n_samples, vaccine_calendar, polymod_data, initial, mapping,
               risk_ratios$value, no_risk_groups, nburn, nbatch, blen)
}

#' Aggregate model results at different time points
#' 
#' This function is useful to convert the mcmc results into aggregated results, such as mean, variance etc.
#' 
#' @param func The function that gets called for each set of parameters (e.g. infectionODEs)
#' @param batch Posterior parameters samples resulting from mcmc. Each row is a set of parameters
#' @param aggregate The aggragation function (i.e. mean, var etc). If the function returns a list, then the names are used in the return value
#' @param ... Extra parameters passed to func.
#' 
#' @return Returns a data.frame in long format, with for each time point, interval, column.ID the value
#' 
aggregateModel <- function( func, batch, aggregate, ... )
{
  values <- list()
  for( k in 1:nrow(batch) )
  {
    df <- func(batch[k,], ...)
    values[[k]] <- df
  }
  
  column.ID <- c()
  row.ID <- c()
  interval <- c()
  value <- c()
  for (i in rownames(values[[1]]))
    for (j in colnames(values[[1]]))
    {
      agg <- aggregate(sapply( values, function(v) v[i,j] ))
      vs <- as.vector(unlist(agg))
      column.ID <- c( column.ID, rep(j,length(vs)) )
      row.ID <- c( row.ID, rep(i,length(vs)) )
      
      if (class(agg)=="list")
      {
        nms <- names(agg)
        for (l in 1:length(agg))
             interval <- c(interval, rep( nms[[l]], length( agg[[l]] ) ) )
      }
      else
        interval <- c( interval, seq(1,length(vs)) )
      value <- c( value, vs )
    }
  return(
      data.frame(list(
      "column.ID"=column.ID,
      "row.ID"=row.ID,
      "aggregate"=interval,
      "value"=value
    ) ) )
}

#' @title Calculate the (equal tailed) credible interval at different time points
#' 
#' @description This function is useful to convert the mcmc results into credible intervals, which is needed for plotting
#' your results. Calls aggregateModel to aggregate the results by the passed intervals.
#' 
#' @param func The function that gets called for each set of parameters (e.g. infectionODEs)
#' @param batch Posterior parameters samples resulting from mcmc. Each row is a set of parameters
#' @param intervals The intervals to calculate (by default 0 (median) and 0.98)
#' @param ... Extra parameters passed to func
#' 
#' @return Returns a data.frame in long format, with for each time point, interval, column.ID the value
#' 
credible.interval.model <- function( func, batch, intervals = c( 0.0, 0.98 ), ... )
{
  agg.f <- function( v )
  {
    xs <- sort(v)
    int <- lapply( intervals, function(i) 
    { 
      perc <- i/2
      list(0.5-perc,0.5+perc)
    })
    ls <- lapply( int, function(i) xs[unlist(i)*length(xs)] )
    names(ls) <- intervals
    return(ls)
  }
  aggregateModel( func, batch, agg.f, ... )
}