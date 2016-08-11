context("CostEffectiveness")

test_that("We can calculate the number of certain health outcome", 
  {
    expect_equal(public_health_outcome(2, c(3,4)), c(6,8))
    expect_equal(public_health_outcome(c(2,3), c(3,4)), c(6,12))
    
    # Require no of groups if unclear
    expect_error(public_health_outcome(c(2, 3), c(1, 4, 9, 13, 5, 6)))
    
    # Override no of risk groups
    expect_equal(public_health_outcome(c(2,3,4), c(1,4,9,13,5,6), no_risk_groups = 3), c(2*1,2*4,3*9,3*13,4*5,4*6))
    
    # Override no of age groups
    expect_equal(public_health_outcome(c(2,3,4), c(1,4,9,13,5,6), no_age_groups = 3), c(2*1,3*4,4*9,2*13,3*5,4*6))
  }
)

test_that("We can pass a list to public_health_outcome", 
  {
    res <- public_health_outcome(list("gp" = 2, "hosp" = c(2, 3)), c(3, 4))
    expect_equal(names(res), c("gp", "hosp"))
    expect_equal(res$gp, c(6, 8))
    expect_equal(res$hosp, c(6, 12))
  }
)


test_that("We can calculate final coverage from a vaccination calendar", 
  {
    if (!exists(".final_coverage"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".final_coverage not available (hidden)")
    }
    data(coverage)
    # Coverage rates for respectively low risk <65, low risk 65+,
    # high risk <65 and 65+. Original is in percentages. Here converted to fraction
    cov <- coverage[,c("Under.65","X65","at.risk.under.65","X65")]/100.0
    
    vaccine_calendar <- as.vaccination.calendar(efficacy = c(0.7, 0.4, 0.7, 0.4), 
                                                dates = coverage$Date,
                                                coverage = cov, 
                                                no_age_groups = 2, no_risk_groups = 2)
    
    expect_equal(.final_coverage(vaccine_calendar), c(0.068, 0.736, 0.466, 0.736, 0, 0))
    expect_equal(vaccine_doses(vaccine_calendar, c(2,3,4,5,0,0)), c(2*0.068, 3*0.736, 4*0.466, 5*0.736, 0, 0))
    expect_equal(vaccine_doses(vaccine_calendar, c(2,3,4,5)), c(2*0.068, 3*0.736, 4*0.466, 5*0.736))
    expect_error(vaccine_doses(vaccine_calendar, c(2,3,4)))
  }
)
 