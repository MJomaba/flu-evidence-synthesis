context("CostEffectiveness")

test_that("We can calculate the number of hospitalisations", 
  {
    expect_equal(hospitalisations(2, c(3,4)), c(6,8))
    expect_equal(hospitalisations(c(2,3), c(3,4)), c(6,12))
    
    # Require no of groups if unclear
    expect_error(hospitalisations(c(2, 3), c(1, 4, 9, 13, 5, 6)))
    
    # Override no of risk groups
    expect_equal(hospitalisations(c(2,3,4), c(1,4,9,13,5,6), no_risk_groups = 3), c(2*1,2*4,3*9,3*13,4*5,4*6))
    
    # Override no of age groups
    expect_equal(hospitalisations(c(2,3,4), c(1,4,9,13,5,6), no_age_groups = 3), c(2*1,3*4,4*9,2*13,3*5,4*6))
  }
)

test_that("We can calculate the mortality", 
  {
    expect_equal(mortality(2, c(3,4)), c(6,8))
    expect_equal(mortality(c(2,3), c(3,4)), c(6,12))
    
    # Require no of groups if unclear
    expect_error(mortality(c(2, 3), c(1, 4, 9, 13, 5, 6)))
    
    # Override no of risk groups
    expect_equal(mortality(c(2,3,4), c(1,4,9,13,5,6), no_risk_groups = 3), c(2*1,2*4,3*9,3*13,4*5,4*6))
    
    # Override no of age groups
    expect_equal(mortality(c(2,3,4), c(1,4,9,13,5,6), no_age_groups = 3), c(2*1,3*4,4*9,2*13,3*5,4*6))
  }
)