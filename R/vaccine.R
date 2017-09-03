#' Function to read legacy file format for vaccine calendar
#'
#' @param file The file to load
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
#'
read.legacy.vaccine.file <- function( file )
{
  results <- list()
  # fill results with vaccine_calendars, which are of type:
  # list( "efficacy"=c(), "calendar" = matrix())
  count <- 0
  started <- FALSE
  for (line in readLines(file))
  {
    if (count==2 & started==FALSE)
    {
      count <- 0
      started <- TRUE
    }
    if (started)
    {
      if (count==2)
      {
        efficacy <- as.numeric(utils::read.csv(text=line,sep=" ",header=FALSE)[1,])
        str.calendar <- ""
      } else if (count>3 & count<127) {
        str.calendar <- paste(str.calendar,line,sep="\n")
      }
      if (count==126)
      {
        calendar <- as.matrix(utils::read.csv(text=str.calendar,sep=" ",header=FALSE))
        results[[length(results)+1]] <- list("efficacy"=efficacy,
                                         "calendar"=calendar)
        count <- -1
      }
    }
    count <- count + 1
  }
  results
}

#' @title
#' Create a vaccination calendar object
#' 
#' @description 
#' Helper function to create (valid) vaccination calendars. It will take the given data, do some
#' sanity checking and convert it into a list object that can be passed to other functions.
#'
#' @param efficacy Vector holding the vaccine efficacy for each age group and risk group
#' @param dates Vector containing the dates at which the coverage was measured.  
#' @param coverage Data frame with the fraction of coverage for each age and risk group in each column. Different rows are the different dates of measurements
#' @param legacy An optional legacy object to convert. When passed a legacy object this function will sanity check the given object
#' @param no_risk_groups The total number of risk groups (optional)
#' @param no_age_groups The total number of age groups (optional)
#' @param starting_year The year of the start of the season. Only used when passed a legacy file
#' @seealso \url{http://blackedder.github.io/flu-evidence-synthesis/vaccination.html}
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
#'
as_vaccination_calendar <- function(efficacy = NULL, dates = NULL, coverage = NULL, legacy = NULL, no_risk_groups = NULL, no_age_groups = NULL, starting_year = NULL) {
  if (!is.null(legacy))
  {
    # We are dealing with a legacy object
    vc <- legacy
    if (is.null(no_risk_groups))
    {
      no_risk_groups <- 1
      if (is.null(no_age_groups))
      {
        if (length(vc$efficacy) == ncol(vc$calendar))
        {
          if (length(vc$efficacy) %% 3 == 0)
            no_risk_groups <- 3
          if (length(vc$efficacy) %% 2 == 0)
            no_risk_groups <- 2
        }
        else
          no_risk_groups <- ncol(vc$calendar)/length(vc$efficacy)
      }
      else
        no_risk_groups <- ncol(vc$calendar)/no_age_groups
    }
    if (is.null(no_age_groups))
      no_age_groups <- ncol(vc$calendar)/no_risk_groups
    
    if (no_age_groups != round(no_age_groups) 
        || no_risk_groups != round(no_risk_groups))
      stop("no_age_groups: ", no_age_groups," and no_risk_groups: ", no_risk_groups, 
           " need to be integer values")
    
    if (length(vc$efficacy) != ncol(vc$calendar))
      vc$efficacy <- rep(vc$efficacy, no_risk_groups)
    if (no_risk_groups < 3) # append zeros
      vc$efficacy <- c(vc$efficacy, rep(0, no_age_groups*(3 - no_risk_groups)))
    if (no_risk_groups > 3)
      stop("Only three risk groups supported currently")
    
    # Dates (if not given, and 123 rows are there -> use default uk values)
    if (is.null(starting_year))
      starting_year <- 1970
    if (nrow(vc$calendar) == 123 & is.null(vc$dates)) 
    {
      vc[["dates"]] <- c(as.Date(paste0(starting_year,"-10-01")),
                         as.Date(paste0(starting_year,"-11-01")),
                         as.Date(paste0(starting_year,"-12-01")),
                         as.Date(paste0(starting_year + 1,"-01-01")),
                         as.Date(paste0(starting_year + 1,"-02-01")))
      vc[["calendar"]] <- matrix(c(vc[["calendar"]][1,],
                                   vc[["calendar"]][32,],
                                   vc[["calendar"]][62,],
                                   vc[["calendar"]][93,]),ncol = 21, byrow = TRUE)
    }
    # Also make sure they are dates
    
    # Calendar (convert it to a matrix)
    # If dates 1 longer -> add row of zeros to calendar
    if (length(vc$dates) == nrow(vc$calendar) + 1)
    {
      vc$calendar <- rbind(vc$calendar, rep(0, ncol(vc$calendar)))
    }
    else if (!all(vc$calendar[nrow(vc$calendar),] == 0))
    {
      warning("No date given when vaccination is stopped.")
    }
    # If no_risk_groups<3, fill with extra zero columns
    if (no_risk_groups < 3)
    {
      zeros <- matrix(rep(0,(3 - no_risk_groups)*no_age_groups*nrow(vc$calendar)), nrow = nrow(vc$calendar))
      vc$calendar <- cbind( vc$calendar, zeros )
    }
    
    if (!all(vc$calendar >= 0))
        stop("Vaccination rates should be positive")
    
    
    if (length(vc$dates) != nrow(vc$calendar))
      stop("Incompatible number of dates and rows in the calendar")
    
    rates <-  vc$calendar[1:(nrow(vc$calendar) - 1),]*
      as.numeric((vc$dates[2:(length(vc$dates))] - vc$dates[1:(length(vc$dates) - 1)]))
    if (class(rates) == "matrix")
      rates <- colSums(rates)
    
    # Check that summed rates (by column,multiplied by number of days) are below 1
    if (!all(rates <= 1))
      stop( "Total fraction of people vaccinated greater than 1" )
    vc$no_risk_groups <- no_risk_groups
    vc$no_age_groups <- no_age_groups
    return(vc)
  } else {
    vc <- list("efficacy" = efficacy, "dates" = dates)
    if (!all(coverage[1,] == 0)) {
      stop("The first row of the coverage table should be 0. This is required to no the beginning of the vaccination program (using the first element of the dates vector)")
    }
    
    diff.dates <- dates[2:length(dates)] - dates[1:(length(dates) - 1)]
    diff.coverage <- coverage[2:nrow(coverage),] - coverage[1:(nrow(coverage) - 1),]
    vc$calendar <- as.matrix(diff.coverage)/as.double(diff.dates)
    # Note that if passed data we can recursively call it with the legacy argument to 
    # normalize the result
    return(as_vaccination_calendar(legacy = vc, no_risk_groups = no_risk_groups, no_age_groups = no_age_groups))
  }
}

#' Calculate number of influenza cases given a vaccination strategy
#'
#' @param vaccine_calendar A vaccine calendar valid for that year
#' @param parameters The parameters to use. Both a vector or a data frame with each row a set of parameters 
#' (e.g. a batch of inferred parameters by adaptive_mcmc$batch) are accepted.
#' @param contact_ids Optional: The contact_ids used to infer the contact matrix. Similar to the \code{parameters} this can be a
#' vector or a data frame.
#' @param incidence_function An optional function that takes a \code{vaccine_calendar}, \code{parameters} and optionally
#' \code{contact_ids} and returns the incidence over time. If none is provided then the default
#' \code{infectionODEs} is used. In that case both \code{polymod_data} and \code{demography} data need to be
#' provided as well.
#' @param time_column An optional time column that will be removed from the data returned by the incidence_function
#' @param ... Further parameters that will be passed to the \code{incidence_function}
#' @param verbose Whether to display warnings. Default is TRUE.
#' 
#' @seealso \code{\link{infectionODEs}}
#' 
#' @return The total incidence in a year per age and risk group
vaccination_scenario <- function(vaccine_calendar, parameters,  
                                 contact_ids, incidence_function,
                                 time_column, ..., verbose = T) {
  if (missing(incidence_function)) {
    var_names <- names(sys.call())
    if (!"polymod_data" %in% var_names) {
      stop("No polymod_data set provided")
    } else {
      polymod_data <- eval(match.call()[["polymod_data"]])
    }
    if (missing(contact_ids)) {
      stop("No contact_ids set provided")
    }
    if (!"demography" %in% var_names) {
      stop("No demography provided, i.e. a vector with population size by age (starting at age is zero)")
    } else {
      demography <- eval(match.call()[["demography"]])
    }
    time_column = "Time"
    incidence_function <- function(vaccine_calendar, parameters, contact_ids, ...) {
      if (!"age_group_limits" %in% var_names) {
        if (verbose)
          warning("Missing age_group_limits, using default: c(1,5,15,25,45,65)")
        age_group_limits <- c(1,5,15,25,45,65)
      } else {
        age_group_limits <- eval(match.call()[["age_group_limits"]])
      }
      contacts <- contact_matrix(as.matrix(polymod_data[contact_ids,]),
                                 demography, age_group_limits )
      
      age.groups <- stratify_by_age( demography, 
                                     age_group_limits )
      
      # Fraction of each age group classified as high risk
      # We can classify a third risk group, but we are not doing
      # that here (the second row is 0 in our risk.ratios matrix)
      if (!"risk_ratios" %in% var_names) {
        if (verbose)
          warning("Missing risk_ratios, using default UK based values")
        risk_ratios <- matrix(c(
          0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
          0, 0, 0, 0, 0, 0, 0                          
        ), ncol = 7, byrow = T)
      } else {
        risk_ratios <- eval(match.call()[["risk_ratios"]])
      }
      
      # Only print warnings on first use of the function
      verbose <<- F
      
      # Population sizes in each age and risk group
      popv <- stratify_by_risk(
        age.groups, risk_ratios );
      
      # Population size initially infected by age and risk group
      initial.infected <- rep( 10^parameters[9], 7 ) 
      initial.infected <- stratify_by_risk(
        initial.infected, risk_ratios );
      
      # Run simulation
      # Note that to reduce complexity 
      # we are using the same susceptibility parameter for multiple age groups
      infectionODEs(popv, initial.infected,
                    vaccine_calendar,
                    contacts,
                    c(parameters[6], parameters[6], parameters[6],
                      parameters[7], parameters[7], parameters[7], parameters[8]),
                    transmissibility = parameters[5],
                    c(0.8,1.8), 7)
    }
  }
 
  # One set of parameters
  if (is.null(nrow(parameters))) {
    if (missing(contact_ids)) {
      result <- incidence_function(vaccine_calendar, parameters, ...)
    } else {
      result <- incidence_function(vaccine_calendar, parameters, contact_ids, ...)
    }
    if (!missing(time_column) && !is.null(time_column))
      result[[time_column]] <- NULL
    return(colSums(result))
  } else {
    # Table of parameters and optionally contact_ids
    if (missing(time_column))
      time_column <- NULL
    if (missing(contact_ids)) {
      return(apply(parameters, 1, function(pars) 
        vaccination_scenario(parameters = pars, vaccine_calendar = vaccine_calendar,
                             incidence_function = incidence_function, 
                             time_column = time_column, ...)
        ))
    } else {
      pc <- cbind(parameters, contact_ids)
      return(t(apply(pc, 1, function(pars_contacts) 
        vaccination_scenario(parameters = pars_contacts[1:ncol(parameters)], 
                             contact_ids = pars_contacts[(ncol(parameters)+1):length(pars_contacts)],
                             vaccine_calendar = vaccine_calendar,
                             incidence_function = incidence_function, 
                             time_column = time_column, ...))
        ))
    }
  }
}