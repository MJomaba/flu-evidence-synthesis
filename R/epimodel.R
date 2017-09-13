#' Run the SEIR model for the given parameters
#'
#' @param population The population size of the different age groups, subdivided into risk groups 
#' @param initial_infected The corresponding number of initially infected
#' @param vaccine_calendar A vaccine calendar valid for that year
#' @param contact_matrix Contact rates between different age groups
#' @param susceptibility Vector with susceptibilities of each age group
#' @param transmissibility The transmissibility of the strain
#' @param infection_delays Vector with the time of latent infection and time infectious
#' @param interval Optional: interval (in days) between data points (used if dates are not provided)
#' @param dates Optional: dates to return values for.
#' @return A data frame with number of new cases after each interval during the year
#' 
#' @seealso \code{\link{infectionODEs.cpp}} Used internally by this function.
infectionODEs <- function(population, initial_infected, vaccine_calendar, contact_matrix,
                          susceptibility, transmissibility, infection_delays, interval = 7,
                          dates = NULL )
{
  if (is.null(dates))
  {
    yr <- 1970
    if (length(vaccine_calendar$dates)>0)
      yr <- data.table::year(vaccine_calendar$dates[1])
    start.date <- as.Date(getTimeFromWeekYear(35,yr))
    dates <- c(start.date)
    #latest.date <- start.date
    latest.date <- start.date + interval
    while(!(data.table::year(latest.date) > data.table::year(start.date) & 
              data.table::yday(latest.date) >= data.table::yday(start.date)))
    {
      dates <- c(dates, latest.date)
      latest.date <- latest.date + interval
    }
  }
  #print(dates)
  #if (class(dates[1]!=Date))
  #  stop( "Dates must be of class Date" );
  
  infectionODEs.cpp(population, initial_infected, vaccine_calendar, contact_matrix,
                    susceptibility, transmissibility, infection_delays, dates )
}

#' Mapping parameters to the model
#' 
#' @description To reduce complexity it is common to map certain parameters to multiple age groups. For example
#' in the UK model we use the same susceptibility for age groups 1,2 and 3. Similarly for fitting purposes we
#' need a vector of parameters. That vector will hold different parameters: ascertainment (epsilon; by age group), 
#' psi, transmissibility, susceptibility (by age group) and initial infected. This method makes this easier.
#' 
#' In the UK model we use 9 parameters, the first three are the epsilon for age groups 1,2 and 3,4 and 5. This is then followed by the 
#' psi parameter (4th parameter in the parameter vector). Next is the transmissibility parameter. Then three susceptibility parameters, used in
#' age groups: 1,2,3 and 4,5,6 and 7. Finally followed by initial_infected as the 9th parameter in the parameters vector. 
#' \code{parameter_map(epsilon = c(1,1,2,2,3), psi = 4, transmissibility = 5, susceptibility = c(6,6,6,7,7,7,8), initial_infected = 9}
#' 
#' @param epsilon A vector holding the indeces of the ascertainment parameters for each age group in the parameters vector
#' @param psi Index of the psi parameter in the parameters vector
#' @param transmissibility Index of the transmissibility parameter in the parameters vector
#' @param susceptibility A vector holding the indeces of the susceptibility parameter for each age group in the parameters vector
#' @param initial_infected Index of the initial_infected parameter in the parameters vector
#' @param parameters Optional parameters vector. For simple cases it is possible to infer the parameter map just from the (named) parameters vector. In that case
#' no other variables need to be passed to this function.
#' 
#' @return A list specifying the parameter mapping as used in \code{\link{inference}} and 
#' \code{\link{vaccination_scenario}}
parameter_mapping <- function(epsilon, psi, transmissibility, susceptibility, 
                              initial_infected, parameters) {
  if (missing(parameters)) {
    if (!missing(psi) && !missing(transmissibility) && !missing(susceptibility) &&
        !missing(initial_infected)) {
      parameter_map <- list(
        psi = psi,
        transmissibility = transmissibility,
        susceptibility = susceptibility,
        initial_infected = initial_infected
      )
      if (!missing(epsilon))
        parameter_map$epsilon = epsilon
      return(parameter_map)
    } else {
      stop("Missing values for one of the following: psi, transmissibility, susceptibility or initial_infected")
    }
  }
  if (!missing(parameters) && (!is.null(names(parameters)) || is.character(parameters))) {
    parameter_map <- list()
    if (is.character(parameters))
      ns <- parameters
    else
      ns <- names(parameters)
    browser()
    for (i in 1:length(parameters)) {
      if (substr(ns[i],1,1) == "e")
        parameter_map[["epsilon"]] <- c(parameter_map[["epsilon"]], i)
      if (substr(ns[i],1,1) == "p")
        parameter_map[["psi"]] <- c(parameter_map[["psi"]], i)
      if (substr(ns[i],1,1) == "t")
        parameter_map[["transmissibility"]] <- c(parameter_map[["transmissibility"]], i)
      if (substr(ns[i],1,1) == "s")
        parameter_map[["susceptibility"]] <- c(parameter_map[["susceptibility"]], i)
      if (substr(ns[i],1,1) == "i")
        parameter_map[["initial.infected"]] <- c(parameter_map[["initial.infected"]], i)
    }
  } else {
    parameter_map <- list()
    no_age_groups <- (length(parameters) - 3)/2
    parameter_map <- list(
      epsilon = seq(1,no_age_groups),
      psi = no_age_groups + 1,
      transmissibility = no_age_groups + 2,
      susceptibility = seq(no_age_groups + 3, 2*no_age_groups + 2),
      initial.infected = 2*no_age_groups + 3
    )
  }
  parameter_map
}
