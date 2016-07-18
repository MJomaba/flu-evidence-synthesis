#' Separate the population into age groups
#'
#' @description Superseded by stratify_by_age function
#'
#' @param age_sizes A vector containing the population size by age (first element is number of people of age 1 and below)
#' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
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
#' @return A vector with the population in the low risk groups, followed by the other risk groups. The length is equal to the number of age groups times the number of risk groups (including the low risk group).
separate.into.risk.groups <- function(age_groups, risk) {
  .Deprecated("stratify_by_risk")
  stratify_by_risk(age_groups, risk)
}