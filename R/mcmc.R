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
#' @details
#' The method we use here combines data from numerous sources that are then used to compute the likelihood of the predicted 
#' number of influenza cases in a given week. Given the data and the likelihood function we use MCMC to obtain the posterior 
#' distribution of the parameters of an underlying epidemiological model (see also: \code{\link{infectionODEs}}).
#' 
#' When running inference there are four main steps needed are 1) prepare the data, 2) load a vaccination calendar (\code{\link{as_vaccination_calendar}})
#' 3) decide on parameterisation of the model (\url{https://blackedder.github.io/flu-evidence-synthesis/modelling.html}) and 4) run the inference using
#' this function.
#' 
#' The initial parameters vector should contain values for the parameters (in order):
#' 
#' * Ascertainment probabilty for each age group (epsilon)
#' * Outside infection (psi)
#' * Transmissibility
#' * Susceptibility for each age group
#' * Initial number of infections (log transformed)
#' 
#' If your model is more complex and the number of age groups and risk groups are different between the epidemiological model (vaccination calendar) and the influenza data then you need to 
#' provide (one or more of) the following extra variables to the function: \code{parameter_map} (see also: \code{\link{parameter_mapping}}), \code{age_group_map} (see also: 
#' \code{\link{age_group_mapping}}) and \code{risk_group_map} (see also: \code{\link{risk_group_mapping}}). 
#' See \url{https://blackedder.github.io/flu-evidence-synthesis/inference.html} for more details.
#' 
#' @md
#'
#' @param demography A vector with the population size by each age {0,1,..}
#' @param ili The number of Influenza-like illness cases per week
#' @param mon_pop The number of people monitored for ili
#' @param n_pos The number of positive samples for the given strain (per week)
#' @param n_samples The total number of samples tested 
#' @param vaccine_calendar A vaccine calendar valid for that year
#' @param polymod_data Contact data for different age groups
#' @param initial Vector with starting parameter values
#' @param parameter_map Optional mapping parameter (by description and age group) to the relevant index
#' in the initial vector. Needed parameters are: epsilon (ascertainment) with a separate value per data
#' age group, transmissibility, psi (infection from outside sources), susceptibility (with a value per age group)
#' and log of initial_infected population \code{\link{parameter_mapping}}.
#' @param age_groups Optional age groups upper limits used in your model and data. If you use different age groups for the model and the data you need
#' to provide a age_group_map instead.
#' @param age_group_map Optional age group mapping from model age groups to data age groups (\code{\link{age_group_mapping}})
#' @param risk_group_map Optional risk group mapping from model risk groups to data risk groups (\code{\link{risk_group_mapping}}).
#' This parameter is not needed if only one risk group is modelled
#' @param risk_ratios A matrix with the fraction in the risk groups. The leftover fraction is assumed to be low risk. (\code{\link{stratify_by_risk}})
#' @param lprior Optional function returning the log prior probability of the parameters. If no function is passed then a flat prior is used.
#' @param nburn Number of iterations of burn in
#' @param nbatch Number of batches to run (number of samples to return)
#' @param blen Length of each batch
#' 
#' @return Returns a list with the accepted samples and the corresponding llikelihood values and a matrix (contact.ids) containing the ids (row number) of the contacts data used to build the contact matrix.
#'
#' @seealso \code{\link{infectionODEs}}; \code{\link{age_group_mapping}}; \code{\link{risk_group_mapping}}; \code{\link{parameter_mapping}}; \url{https://blackedder.github.io/flu-evidence-synthesis/inference.html}
#'
#' @export
inference <- function(demography, ili, mon_pop, n_pos, n_samples, 
        vaccine_calendar, polymod_data, initial, parameter_map, age_groups, age_group_map,
        risk_group_map, risk_ratios, lprior, nburn = 0, nbatch = 1000, blen = 1 )
{
  uk_defaults <- F
  if (any(n_samples>ili))
    stop("The model assumes that the virological samples are a subsample of patients diagnosed as ILI cases. The ili counts should always be larger than or equal to n_samples") 
  no_risk_groups <- vaccine_calendar$no_risk_groups
  no_age_groups <- vaccine_calendar$no_age_groups
  if (no_risk_groups >= 2 && no_age_groups == 7 && ncol(ili) == 5)
    uk_defaults <- T
  if (missing(age_group_map))
    if (missing(age_groups))
      if (uk_defaults) { # Using UK defaults
        age_group_map <- age_group_mapping(c(1,5,15,25,45,65), c(5,15,45,65))
      } else
        stop("Need to provide either age_groups or age_group_map")
    else 
      age_group_map <- age_group_mapping(age_groups, age_groups)
  if (missing(risk_ratios)) {
    if (uk_defaults) {
        risk_ratios <- matrix(c(0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, rep(0,no_age_groups*(no_risk_groups-2))), ncol = 7, byrow = T)
    } else {
      if (no_risk_groups > 1)
        stop("No risk ratios supplied.")
      risk_ratios <- rep(1, no_age_groups)
    }
  }
  if (class(risk_ratios) == "matrix") {
    rv <- c(rep(1, ncol(risk_ratios)) - colSums(risk_ratios), t(risk_ratios))
    risk_ratios <- data.frame(
      AgeGroup = rep(unique(age_group_map$from), no_risk_groups),
      value = rv
    ) %>% dplyr::group_by(AgeGroup) %>% dplyr::mutate(RiskGroup = factor(row_number())) %>%
      dplyr::ungroup()
  } else if (class(risk_ratios) != "data.frame") {
    risk_ratios <- data.frame(value = risk_ratios) 
  }
  if (missing(risk_group_map)) {
    risk_group_map <- risk_group_mapping(from = factor(1:no_risk_groups),
                                 to = factor(1:(ncol(ili)/length(unique(age_group_map$to)))))
  }
  
  mapping <- data.frame()
  ## Add risk group
  for (i in 1:nrow(risk_group_map)) {
    mapping <- rbind(mapping, age_group_map %>% dplyr::mutate(join = i))
  }
  mapping <- mapping %>% dplyr::left_join(
    risk_group_map %>% dplyr::mutate(join = row_number()) %>%
      dplyr::rename(from_risk = from, to_risk = to, weight_risk = weight), by = "join") %>%
    dplyr::mutate(weight = weight*weight_risk) %>% dplyr::select(-weight_risk, -join)
  
  from_i <- mapping %>% dplyr::select(from, from_risk) %>% dplyr::distinct() %>% dplyr::mutate(from_i = row_number() - 1)
  to_j <- mapping %>% dplyr::select(to, to_risk) %>% dplyr::distinct() %>% dplyr::mutate(to_j = row_number() - 1)
  mapping <- mapping %>% dplyr::left_join(from_i, by = c("from", "from_risk")) %>% 
    dplyr::left_join(to_j, by = c("to", "to_risk")) %>% dplyr::select(from_i, to_j, weight)
  
  if (missing(parameter_map)) {
    if (uk_defaults) {
      parameter_map <-
        parameter_mapping(
          epsilon = c(1,1,2,2,3),
          psi = 4,
          transmissibility = 5,
          susceptibility = c(6,6,6,7,7,7,8),
          initial_infected = 9)
    } else if (length(initial) == 2*no_age_groups + 3) {
      parameter_map <- parameter_mapping(parameters = initial)
    } else {
      stop("Missing parameter map")
    }
  }
  
  # Go over parameter_map. Shorten list and also make sure min(index) = 0
  m <- min(unlist(parameter_map))
  batch_cols <- rep("NULL", length(initial)) 
  for(n in names(parameter_map)) {
    parameter_map[[substr(n, 1, 1)]] <- parameter_map[[n]] - m
    for (i in parameter_map[[n]]) {
      if (!is.null(names(initial))) 
        batch_cols[i-m+1] <- names(initial)[i]
      else
        batch_cols[i-m+1] <- n
    }
  }
  
  # Add numbering for duplicate parameter names
  b_cols <- data.frame(value = batch_cols) %>%
    dplyr::group_by(value) %>%
    dplyr::mutate(ext = row_number(), mx = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(value = paste0(value, ifelse(mx > 1, paste0("_", ext),"")))
  
  pass_prior = T
  if (missing(lprior)) {
    lprior <- function(pars) {} # Dummy function
    pass_prior = F
  }
  
  results <- .inference_cpp(demography, sort(unique(age_group_limits(as.character(age_group_map$from)))),
                 as.matrix(ili), as.matrix(mon_pop), as.matrix(n_pos), as.matrix(n_samples), vaccine_calendar, polymod_data, initial, 
                 as.matrix(mapping), risk_ratios$value, 
                 parameter_map$e, parameter_map$p, parameter_map$t, parameter_map$s, parameter_map$i, lprior, pass_prior,
                 no_age_groups, no_risk_groups, uk_defaults, nburn, nbatch, blen)
  if (is.null(names(initial))) {
    colnames(results$batch) <- b_cols$value
  } else
    colnames(results$batch) <- batch_cols
  results
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