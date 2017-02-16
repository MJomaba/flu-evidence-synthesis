#' Separate the population into age groups
#'
#' @description Superseded by stratify_by_age function
#'
#' @param age_sizes A vector containing the population size by age (first element is number of people of age 1 and below)
#' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
#'
#' @keywords internal
#'
#' @return A vector with the population in each age group.
separate.into.age.groups <- function(age_sizes, limits = c(1, 5, 15, 25, 45, 65)) {
  .Deprecated("stratify_by_age")
  stratify_by_age(age_sizes, limits)
}

#' Separate the population into risk groups
#'
#' @description Superseded by stratify_by_risk function
#' 
#' @param age_groups A vector containing the population size of each age group
#' @param risk A matrix with the fraction in the risk groups. The leftover fraction is assumed to be low risk
#'
#' @keywords internal
#'
#' @return A vector with the population in the low risk groups, followed by the other risk groups. The length is equal to the number of age groups times the number of risk groups (including the low risk group).
separate.into.risk.groups <- function(age_groups, risk) {
  .Deprecated("stratify_by_risk")
  stratify_by_risk(age_groups, risk)
}


#' Age as age group
#'
#' @description Renamed to as_age_group
#'
#' @param age The relevant age 
#' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
#' 
#' @keywords internal
#'
#' @return An integer representing the age group
as.age.group <- function(age, limits = as.numeric(c(1, 5, 15, 25, 45, 65))) {
  .Deprecated("as_age_group")
  as_age_group(age, limits)
}

#' Function to create a vaccine calendar object
#'
#' @description Deprecated
#' @param coverage A vector containing the total coverage for the vaccine for each age/risk group
#' @param efficacy The efficacy of the vaccine for that year (subdivided by each age group)
#' @param uptake The uptake of the vaccine over different time periods
#' 
#' @keywords internal
#' 
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
vaccine.calendar <- function( coverage, efficacy, uptake )
{
  .Deprecated("as_vaccination_calendar")
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

#' @title
#' Create a vaccination calendar object
#' 
#' @description Renamed to as_vaccination_calendar
#' 
#' @param efficacy Vector holding the vaccine efficacy for each age group and risk group
#' @param dates Vector containing the dates at which the coverage was measured.  
#' @param coverage Data frame with the fraction of coverage for each age and risk group in each column. Different rows are the different dates of measurements
#' @param legacy An optional legacy object to convert. When passed a legacy object this function will sanity check the given object
#' @param no_risk_groups The total number of risk groups (optional)
#' @param no_age_groups The total number of age groups (optional)
#' @param starting_year The year of the start of the season. Only used when passed a legacy file
#' @seealso  \code{\link{as_vaccination_calendar}}
#' @keywords internal
#' @return A list that contains the \code{calendar} and \code{efficacy} of the vaccine for that year
#'
as.vaccination.calendar <- function(efficacy = NULL, dates = NULL, coverage = NULL, legacy = NULL, no_risk_groups = NULL, no_age_groups = NULL, starting_year = NULL) {
  .Deprecated("as_age_group")
  as_vaccination_calendar(efficacy, dates, coverage, legacy, no_risk_groups, 
                          no_age_groups, starting_year)
}

