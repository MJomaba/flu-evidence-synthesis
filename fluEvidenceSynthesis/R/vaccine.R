
#' Number of Influenza-like illness cases per week
#'
#' @format A list with both the number of cases and the total number of monitored people
#' \describe{
#'   \item{ili}{Matrix with total number of cases per week (row) and age group (column)}
#'   \item{total.monitored}{The number of people monitored over the same time/age group}
#' }
"ili"

#' Number of people of each age for the UK in 1999
#'
"age_sizes"

#' MCMC Result produced by the inference function
#'
"mcmcsample"

#' Function to create a vaccine calendar object
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

