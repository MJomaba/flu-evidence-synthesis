#' @title Calculate a mapping from one sets of age groups to another set of age groups
#' 
#' @description Sometimes data/results are stratified in different ways. For example in the UK data set vaccination
#' is divided into age groups: <1, <5, <15, <25, <45, <65 and 65+, while virological data is 
#' provided as <5, <15, <45, <65 and 65+. This function calculates a mapping from one set of age groups
#' to another. If the age groups we want to map to is a subset of the original age groups (as in the UK)
#' then this is straightforward. If instead one or more of the from age group needs to be split over multiple to age groups
#' then the function weights the age group by the size of the population (demography).
#' 
#' @param from The upper limit to each age groups you want to map from
#' @param to The upper limit to each age groups you want to map to
#' @param demography Optional vector with the demography of the population. This value is needed if the to age groups are
#' not a subset of the from age groups. First value is the population of 0 year olds.
#' Second of 1 year olds etc...
#' @return A data frame containing which age groups map to which other age group. And what proportion (weight)
#' of the age groups maps to the other age grouping.
#' @export
age_group_mapping <- function(from, to, demography) {
  i <- 1
  j <- 1
  is <- c()
  js <- c()
  weights <- c()
  while(TRUE) {
    if (j == length(to) + 1 && i == length(from) + 1) {
      is <- c(is, i)
      js <- c(js, j)
      weights <- c(weights, 1)
      break
    } else if (j == length(to) + 1) {
      is <- c(is, i)
      js <- c(js, j)
      weights <- c(weights, 1)
      i <- i + 1
    } else if (from[i] <= to[j]) {
      is <- c(is, i)
      js <- c(js, j)
      weights <- c(weights, 1)
      if (from[i] == to[j])
        j <- j + 1
      i <- i + 1
    } else if (from[i] > to[j]) {
      is <- c(is, i)
      js <- c(js, j)
      if (missing(demography))
        stop("Mapping incompatible age groups requires demography data")
      fa <- function(ag) if (ag < 0 || ag >= length(demography)) 0 else demography[ag+1]
      low_w <- sum(sapply(seq(max(-1, from[i - 1], na.rm = T), to[j] - 1), fa))
      high_w <- sum(sapply(seq(to[j], from[i] - 1), fa))
      j <- j + 1
      is <- c(is, i)
      js <- c(js, j)
      i <- i + 1
      weights <- c(weights, low_w/(low_w+high_w), high_w/(low_w+high_w))
      #weights <- c(weights, low_w, high_w)
    } else {
      warning("Invalid state in age_groups_mapping()")
    }
  }
  data.frame(from = factor(is, labels = age_group_levels(from)), 
             to = factor(js, labels = age_group_levels(to)), 
             weight = weights)
}

#' @title Calculate a mapping from one sets of risk groups to another set of risk groups
#' 
#' @description Sometimes data/results are stratified in different ways. For example in the UK data set vaccination
#' is divided into risk groups: low risk and high risk, while for the virological data is not stratified by risk. 
#' This function provides a standardised way of mapping one to another. 
#' 
#' @param from The risk groups you want to map from
#' @param to The risk groups you want to map to
#' @param weights Optional vector with the weight/fraction of how each from group should be grouped into a to group. For example if your epidemiological model
#' groups everything together, but the data is divided into two equal risk groups you would use: 
#' \code{risk_group_mapping(from = "All", to = c("R1", "R2"), weights = c(0.5, 0.5))}. The UK example uses the mapping: 
#' \code{risk_group_mapping(from = c("High risk", "Low risk", "Pregnant"), to = "All")}
#' @return A data frame containing which risk groups map to which other risk group. And what proportion (weight)
#' of the risk groups maps to the other risk grouping.
#' @export
risk_group_mapping <- function(from, to, weights) {
  if (length(from) == 1 && !missing(to))
    from <- rep(from, length(to))
  if (missing(to))
    to <- rep("All", length(from))
  if (length(to) == 1)
    to <- rep(to, length(from))
  if (missing(weights))
    weights <- rep(1, length(from))
  # Do some checking that weights are valid
  if (sum(weights) != length(unique(from)))
    warning(paste0("Total weight should sum up to: ", length(unique(from))))
  data.frame(from = from,
             to = to,
             weight = weights)
}

#' @title Stratify age groups into different risk groups
#' 
#' @description Stratifies the age groups and returns the population size of each age group and risk group.
#'
#' @param age_groups A vector containing the population size of each age group
#' @param risk A matrix with the fraction in the risk groups. The leftover fraction is assumed to be low risk
#'
#' @return A vector with the population in the low risk groups, followed by the other risk groups. The length is equal to the number of age groups times the number of risk groups (including the low risk group).
#'
#' @export
stratify_by_risk <- function(age_groups, risk_ratios, no_risk_groups)
{
  if (class(age_groups) != "data.frame")
    age_groups <- data.frame(value = age_groups) %>%
      dplyr::mutate(AgeGroup = row_number())
  if (class(risk_ratios) == "matrix") {
    rv <- c(rep(1, ncol(risk_ratios)) - colSums(risk_ratios), t(risk_ratios))
    if (missing(no_risk_groups))
      no_risk_groups <- length(rv)/nrow(age_groups)
    risk_ratios <- data.frame(
      AgeGroup = rep(age_groups$AgeGroup, no_risk_groups),
      value = rv
    ) %>% dplyr::group_by(AgeGroup) %>% dplyr::mutate(RiskGroup = factor(row_number())) %>%
      dplyr::ungroup()
  }
  if (missing(no_risk_groups))
    no_risk_groups <- nrow(risk_ratios)/length(age_groups)
  
  .stratify_by_risk(age_groups$value, risk_ratios$value, no_risk_groups)
}