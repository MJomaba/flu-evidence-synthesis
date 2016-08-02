.proportion_by_group <- function(proportion, incidence, no_risk_groups = NULL, no_age_groups = NULL) {
  if (length(proportion) == 1 || length(proportion) == length(incidence))
    return(proportion*incidence)
  if (is.null(no_risk_groups) && is.null(no_age_groups)) {
    stop("Either no_risk_groups or no_age_groups needs to be specified")
  }
  if (is.null(no_risk_groups)) {
    no_risk_groups <- length(incidence)/no_age_groups
  } else {
    no_age_groups <- length(incidence)/no_risk_groups
  }
  if (no_age_groups == no_risk_groups && length(proportion) == no_age_groups) {
    stop(paste("Unable to infer whether proportion is by age or by risk groups if no_age_groups is equal to no_risk_groups.",
               "Please provide the proportion for all age and risk groups (proportion vector should be the same length as", 
               "the incidence vector)."))
  }
  if (length(proportion) == no_age_groups) {
    return(rep(proportion, no_risk_groups)*incidence)
  }
  if (length(proportion) == no_risk_groups) {
    return(as.vector(sapply(proportion, function(p) rep(p, no_age_groups)))*incidence)
  }
  stop("Length of proportion and incidence vectors are incompatible")
}

#' @title Calculate number of hospitalisations from incidence data
#' 
#' @description Uses the provided proportion to calculate the number of hosptitalisations from the incidence data. The expected
#' proportion of hospitalisations needs to be provided by the user and is likely to be country/disease specific.
#' 
#' @param proportion The expected proportion of cases resulting in hospitalisations. 
#' This can be a constant value for each age/risk groups, or a vector with the proportion for each age/risk group. Finally if
#' either no_risk_groups or no_age_groups is specified it can also be a vector with the proportion by risk groups or by age group.
#' @param incidence Vector containing the expected incidence for each age/risk group as returned by for example \code{\link{vaccinationScenario}}.
#' @param no_risk_groups The total number of risk groups (optional)
#' @param no_age_groups The total number of age groups (optional)
#' @return A vector with the number of hospitalisations for each age/risk group
#' @export
hospitalisations <- function(proportion, incidence, no_risk_groups = NULL, no_age_groups = NULL) {
  .proportion_by_group(proportion, incidence, no_risk_groups, no_age_groups)
}

#' @title Calculate number of deaths from incidence data
#' 
#' @description Uses the provided proportion to calculate the mortality from the incidence data. The expected
#' proportion of hospitalisations needs to be provided by the user and is likely to be country/disease specific.
#' 
#' @param proportion The expected proportion of cases resulting in death. 
#' This can be a constant value for each age/risk groups, or a vector with the proportion for each age/risk group. Finally if
#' either no_risk_groups or no_age_groups is specified it can also be a vector with the proportion by risk groups or by age group.
#' @param incidence Vector containing the expected incidence for each age/risk group as returned by for example \code{\link{vaccinationScenario}}.
#' @param no_risk_groups The total number of risk groups (optional)
#' @param no_age_groups The total number of age groups (optional)
#' @return A vector with the number of deaths for each age/risk group
#' @export
mortality <- function(proportion, incidence, no_risk_groups = NULL, no_age_groups = NULL) {
  .proportion_by_group(proportion, incidence, no_risk_groups, no_age_groups)
}

# Calculate final coverage from a passed vaccination_calendar (used in vaccine_doses)
.final_coverage <- function(vaccination_calendar) {
  dim <- length(vaccination_calendar$dates)
  durations <- vaccination_calendar$dates[2:dim] - vaccination_calendar$dates[1:(dim - 1)]
  return(as.vector(colSums(durations*vaccination_calendar$calendar[1:(dim - 1),])))
}

#' @title Calculate number of vaccine doses needed based on vaccine calendar
#' 
#' @description Calculates the number of doses required to vaccinate the population under the provided vaccination calendar.
#' 
#' @param vaccination_calendar The vaccination calendar. This object should have the same layout as a 
#' vaccination_calendar returned by \code{\link{as.vaccination.calendar}}.
#' @param population A vector with the size of the population in each age and risk group
#' @return The number of doses by age and risk group.
vaccine_doses <- function(vaccination_calendar, population) {
  cov <- .final_coverage(vaccination_calendar)
  if (length(population) == length(cov)) {
    return(cov*population)
  } else {
    dim <- length(population)
    if (all(cov[(dim + 1):length(cov)] == 0)) {
      # Maintain compatibility and just ignore zeros
      return(cov[1:dim]*population)
    } else {
      stop("Vaccination calendar number of groups incompatible with population vector")
    }
  }
}