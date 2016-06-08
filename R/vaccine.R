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
"uptake"

#' Vaccination uptake rate for the UK in 1999
"vaccine_calendar"

#' Function to create a vaccine calendar object
#'
#' @param coverage A vector containing the total coverage for the vaccine for each age/risk group
#' @param efficacy The efficacy of the vaccine for that year (subdivided by each age group)
#' @param uptake The uptake of the vaccine over different time periods
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
#' 
vaccine.calendar <- function( coverage, efficacy, uptake )
{
  .Deprecated("as.vaccination.calendar")
    new.vacc.cal<-matrix(rep(0,123*21),ncol=21)

    LR.cov<-as.numeric(coverage[1:7])
    HR.cov<-as.numeric(coverage[8:14])

    #October
    for(j in 1:31)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[2],6),uptake[3],rep(uptake[2],6),uptake[3],rep(uptake[2],6),uptake[3])/(31*100)  
    #November
    for(j in 32:61)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3],rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3],rep(uptake[4]-uptake[2],6),uptake[5]-uptake[3])/(30*100)
    #December
    for(j in 62:92)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5],rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5],rep(uptake[6]-uptake[4],6),uptake[7]-uptake[5])/(31*100)
    #Januari
    for(j in 93:123)
        new.vacc.cal[j,]<-c(LR.cov,HR.cov,LR.cov)*c(rep(100-uptake[6],6),100-uptake[7],rep(100-uptake[6],6),100-uptake[7],rep(100-uptake[6],6),100-uptake[7])/(31*100)

    cal <- list( "efficacy"=efficacy, "calendar"=new.vacc.cal )
    cal
}

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
#' @param legacy An optional legacy object to convert. Will sanity check the given object
#' @param no.risk.groups The total number of risk groups (optional)
#' @param no.age.groups The total number of age groups (optional)
#' @param starting.year The year of the start of the season. Only used when passed a legacy file
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
#'
as.vaccination.calendar <- function( legacy = NULL, no.risk.groups = NULL, no.age.groups = NULL, starting.year = NULL )
{
  if (!is.null(legacy))
  {
    # We are dealing with a legacy object
    vc <- legacy
    if (is.null(no.risk.groups))
    {
      no.risk.groups <- 1
      if (is.null(no.age.groups))
      {
        if (length(vc$efficacy)==ncol(vc$calendar))
        {
          if (length(vc$efficacy)%%3==0)
            no.risk.groups <- 3
          if (length(vc$efficacy)%%2==0)
            no.risk.groups <- 2
        }
        else
          no.risk.groups <- ncol(vc$calendar)/length(vc$efficacy)
      }
      else
        no.risk.groups <- ncol(vc$calendar)/no.age.groups
    }
    if (is.null(no.age.groups))
      no.age.groups <- ncol(vc$calendar)/no.risk.groups
    
    if (no.age.groups != round(no.age.groups) 
        || no.risk.groups != round(no.risk.groups))
      stop("no.age.groups: ", no.age.groups," and no.risk.groups: ", no.risk.groups, 
           " need to be integer values")
    
    if (length(vc$efficacy)!=ncol(vc$calendar))
      vc$efficacy <- rep(vc$efficacy, no.risk.groups)
    if (no.risk.groups < 3) # append zeros
      vc$efficacy <- c(vc$efficacy, rep(0, no.age.groups*(3-no.risk.groups)))
    if (no.risk.groups > 3)
      stop("Only three risk groups supported currently")
    
    # Dates (if not given, and 123? rows are there -> use default uk values)
    if (is.null(starting.year))
      starting.year <- 1970
    if (nrow(vc$calendar)==123 & is.null(vc$dates)) 
    {
      vc[["dates"]] <- c(as.Date(paste0(starting.year,"-10-01")),
                         as.Date(paste0(starting.year,"-11-01")),
                         as.Date(paste0(starting.year,"-12-01")),
                         as.Date(paste0(starting.year+1,"-01-01")),
                         as.Date(paste0(starting.year+1,"-02-01")))
      vc[["calendar"]] <- matrix(c(vc[["calendar"]][1,],
                                   vc[["calendar"]][32,],
                                   vc[["calendar"]][62,],
                                   vc[["calendar"]][93,]),ncol=21,byrow=TRUE)
    }
    # Also make sure they are dates
    
    # Calendar (convert it to a matrix)
    # If dates 1 longer -> add row of zeros to calendar
    if (length(vc$dates) == nrow(vc$calendar)+1)
    {
      vc$calendar <- rbind(vc$calendar, rep(0, ncol(vc$calendar)))
    }
    else if (!equal(vc$calendar[nrow(vc$calendar),], rep(0,ncol(vc$calendar))))
    {
      warning("No date given when vaccination is stopped.")
    }
    # If no.risk.groups<3, fill with extra zero columns
    if (no.risk.groups < 3)
    {
      zeros <- matrix(rep(0,(3-no.risk.groups)*no.age.groups*nrow(vc$calendar)), nrow = nrow(vc$calendar))
      vc$calendar <- cbind( vc$calendar, zeros )
    }
    
    
    if (length(vc$dates) != nrow(vc$calendar))
      stop("Incompatible number of dates and rows in the calendar")
    
    rates <-  vc$calendar[1:(nrow(vc$calendar)-1),]*as.numeric((vc$dates[2:(length(vc$dates))]-vc$dates[1:(length(vc$dates)-1)]))
    if(class(rates)=="matrix")
      rates <- colSums(rates)
    
    # Check that summed rates (by column,multiplied by number of days) are below 1
    if (!all(rates<=1))
      stop( "Total fraction of people vaccinated greater than 1" )
    
    return(vc)
  }
  # Note that if passed data we can recursively call it with the legacy argument to 
  # normalize the result
}