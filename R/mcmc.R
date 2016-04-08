#' Adaptive MCMC algorithm
#'
#' MCMC which adapts its proposal distribution for faster convergence following:
#' Sherlock, C., Fearnhead, P. and Roberts, G.O. The Random Walk Metrolopois: Linking Theory and Practice Through a Case Study. Statistical Science 25, no.2 (2010): 172-190.
#'
#' @param lprior A function returning the log prior probability of the parameters 
#' @param llikelihood A function returning the log likelihood of the parameters given the data
#' @param nburn Number of iterations of burn in
#' @param initial Vector with starting parameter values
#' @param nbatch Number of batches to run (number of samples to return)
#' @param blen Length of each batch
#' @param ... Extra parameters passed to the log likelihood function
#' 
#' @return Returns a list with the accepted samples and the corresponding llikelihood values
#'
#' @seealso \code{\link{adaptive.mcmc.cpp}} Used internally by this function.
adaptive.mcmc <- function( lprior, llikelihood, nburn, initial, nbatch, blen = 1, ...  )
{
    adaptive.mcmc.cpp( lprior, function(pars) llikelihood( pars, ... ),
                  nburn, initial, nbatch, blen )
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
#' @return Returns a data.table in long format, with for each time point, interval, column.ID the value
#' 
aggregate.model <- function( func, batch, aggregate, ... )
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
      column.ID <- c( column.ID, rep(j,length(agg)) )
      row.ID <- c( row.ID, rep(i,length(agg)) )
      
      if (class(agg)=="list")
        interval <- c( interval, names(agg) )
      else
        interval <- c( interval, seq(1,length(agg)) )
      value <- c( value, unlist(agg) )
    }
  return(
      data.frame(list(
      "column.ID"=column.ID,
      "row.ID"=row.ID,
      "aggregate"=interval,
      "value"=value
    ) ) )
}

#' Calculate the credible interval at different time points
#' 
#' This function is useful to convert the mcmc results into credible intervals, which is needed for plotting
#' your results. Calls aggregate.model to aggregate the results by the passed intervals.
#' 
#' @param func The function that gets called for each set of parameters (e.g. infectionODEs)
#' @param batch Posterior parameters samples resulting from mcmc. Each row is a set of parameters
#' @param intervals The cut offs at which to calculate the intervals
#' @param ... Extra parameters passed to func.
#' 
#' @return Returns a data.table in long format, with for each time point, interval, column.ID the value
#' 
credible.interval.model <- function( func, batch, intervals = c( 0.01, 0.5, 0.99 ), ... )
{
  agg.f <- function( v )
  {
    xs <- sort(v)
    ls <- lapply( intervals, function(i) xs[i*length(xs)] )
    names(ls) <- intervals
    return(ls)
  }
  aggregate.model( func, batch, agg.f, ... )
}