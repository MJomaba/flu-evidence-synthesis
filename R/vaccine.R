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

#' Function to read legacy file format for vaccines
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
        efficacy <- as.numeric(read.csv(text=line,sep=" ",header=FALSE)[1,])
        str.calendar <- ""
      } else if (count>3 & count<127) {
        str.calendar <- paste(str.calendar,line,sep="\n")
      }
      if (count==126)
      {
        calendar <- as.matrix(read.csv(text=str.calendar,sep=" ",header=FALSE))
        results[[length(results)+1]] <- list("efficacy"=efficacy,
                                         "calendar"=calendar)
        count <- -1
      }
    }
    count <- count + 1
  }
  results
}

