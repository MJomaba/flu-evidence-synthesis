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
  uni <- sort(union(from, to))
  i <- 1
  j <- 1
  is <- c()
  js <- c()
  if (length(uni) == length(from)) {
    wi <- rep(1, length(from) + 1)
    wj <- wi
  } else {
    # to are not strict subgroups, so need to have demography
    if (missing(demography))
      stop("Mapping incompatible age groups requires demography data")
    wi <- stratify_by_age(demography, from)
    wj <- stratify_by_age(demography, uni)
  }
  for (k in 1:length(uni)) {
    ag <- uni[k]
    if (i <= length(from) && from[i] < ag) {
      i <- i + 1
      cw_i <- 0
    }
    is <- c(is, i)
    if (j <= length(to) && to[j] < ag) {
      j <- j + 1
      cw_j <- 0
    }
    js <- c(js, j)
  }
  is <- c(is, length(from) + 1)
  js <- c(js, length(to) + 1)
  weights <- wj/wi[is]
  data.frame(from = factor(is, labels = age_group_levels(from)), 
             to = factor(js, labels = age_group_levels(to)), 
             weight = as.numeric(weights))
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
#' @param no_risk_groups Optional number of risk groups.
#' @param labels Optional names for the risk groups.
#'
#' @return A named vector with the population in the low risk groups, followed by the other risk groups. The length is equal to the number of age groups times the number of risk groups (including the low risk group).
#'
#' @export
stratify_by_risk <- function(age_groups, risk_ratios, no_risk_groups, labels)
{

  if (class(age_groups) != "data.frame") {
    age_labels <- paste0("AG", seq(1, length(age_groups)))
    if (!is.null(names(age_groups)) && 
        length(unique(names(age_groups))) == length(age_groups))
      age_labels <- names(age_groups) 
    age_groups <- data.frame(value = age_groups) %>%
      dplyr::mutate(AgeGroup = age_labels)
  }
     
  if (class(risk_ratios) == "matrix") {
    rv <- c(rep(1, ncol(risk_ratios)) - colSums(risk_ratios), t(risk_ratios))
    if (missing(no_risk_groups))
      no_risk_groups <- length(rv)/nrow(age_groups)
    risk_ratios <- data.frame(
      value = rv
    )
  } else if (is.null(dim(risk_ratios))) {
    risk_ratios <- data.frame(
      value = risk_ratios
    )
  }
  if (missing(no_risk_groups)) {
    if (missing(labels))
      no_risk_groups <- nrow(risk_ratios)/length(age_groups)
    else
      no_risk_groups <- length(labels)
  }
  if (missing(labels)) {
    labels <- paste0("RG", seq(1, no_risk_groups))
  }
  
  v <- .stratify_by_risk(age_groups$value, risk_ratios$value, no_risk_groups)
  
  if (is.null(age_groups[["AgeGroup"]]))
    age_groups <- age_groups %>%
      dplyr::mutate(AgeGroup = paste0("AG", row_number()))
  if (is.null(risk_ratios[["AgeGroup"]]))
    risk_ratios <- risk_ratios %>% 
      dplyr::mutate(AgeGroup = rep(age_groups$AgeGroup, no_risk_groups))
  if (is.null(risk_ratios[["RiskGroup"]])) {
    risk_ratios <- risk_ratios %>% 
      dplyr::group_by(AgeGroup) %>%
      dplyr::mutate(RiskGroup = labels) %>%
      dplyr::ungroup()
  }
  
  names(v) <- paste(risk_ratios$RiskGroup, risk_ratios$AgeGroup, sep = " ")
  v
}