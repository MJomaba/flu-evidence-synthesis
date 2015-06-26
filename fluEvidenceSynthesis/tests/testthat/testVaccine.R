context("Vaccine")

test_that("We can load age data", 
          {
              data("age_sizes")
              expect_that( nrow(age_sizes), equals(86) )
          }
)

test_that("We can load contact data", 
          {
              data("polymod_uk")
              expect_that( nrow(polymod_uk), equals(597) )
              expect_that( ncol(polymod_uk), equals(9) )
              expect_true( polymod_uk$weekend[2] )
          }
)

test_that("We can call scenario",
          {
              data("age_sizes")
              expect_true( vaccination_scenario( age_sizes=age_sizes[,1] ) )
          }
)
