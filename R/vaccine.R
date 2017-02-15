#' Number of Influenza-like illness cases per week
#'
#' @format A list with both the number of cases and the total number of monitored people
#' \describe{
#'   \item{ili}{Matrix with total number of cases per week (row) and age group (column)}
#'   \item{total.monitored}{The number of people monitored over the same time/age group}
#' }
"ili"

#' Number of confirmed positive samples per week
#'
#' Generally only for a fraction of the Influenza-like illness cases an actual (blood) sample is taken to confirm whether/which strain the subject is infected with.
#'
#' @format A list with both the number of positive samples and the total number samples taken
#' \describe{
#'   \item{positive}{Matrix with total number of positive samples per week (row) and age group (column)}
#'   \item{total.samples}{The total number of samples tested over the same time/age group}
#' }
"confirmed.samples"

#' Number of people of each age for the UK in 1999
#'
"age_sizes"

#' MCMC Result produced by the inference function
#'
"inference.results"

#' DEPRECATED: MCMC Result produced by the inference function
#' 
#' @keywords internal
"mcmcsample"

#' Polymod contact data for the UK
#'
#' @format Data frame containing all polymod data for the UK
#' \describe{
#'    \item{V1}{Age of the subject}
#'    \item{V2}{Data taken in the weekend (1) or during the week (0)}
#'    \item{V3-V9}{Number of contacts with persons from the different age groups}
#' }
"polymod_uk"

#' Data on the uptake rates of vaccination in the UK
#' @keywords internal
"uptake"

#' Vaccination uptake rate for the UK in 1999
"vaccine_calendar"

#' Data on vaccine coverage during the 2007-08 season
"coverage"


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