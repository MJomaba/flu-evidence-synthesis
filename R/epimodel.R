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
      yr <- year(vaccine_calendar$dates[1])
    start.date <- as.Date(getTimeFromWeekYear(35,yr))
    dates <- c(start.date)
    #latest.date <- start.date
    latest.date <- start.date + interval
    while( !( year(latest.date)>year(start.date)&
              #month(latest.date)>=month(start.date)&
              yday(latest.date)>=yday(start.date) ) )
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
