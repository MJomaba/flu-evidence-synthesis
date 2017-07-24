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