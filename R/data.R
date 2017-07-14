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

#' Demography of the UK in 1999
#' 
#' @description A vector containing the number of people of a certain age in the UK. The first element
#' is people of age 0, second is age 1 etc.
"demography"


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
#'    \item{Age}{Age of the subject}
#'    \item{Weekend}{Data taken in the weekend (1) or during the week (0)}
#'    \item{...}{Number of contacts with individuals from the age groups}
#' }
"polymod_uk"

#' Data on the uptake rates of vaccination in the UK
#' @keywords internal
"uptake"

#' Vaccination uptake rate for the UK in 1999
"vaccine_calendar"

#' Data on vaccine coverage during the 2007-08 season
"coverage"
