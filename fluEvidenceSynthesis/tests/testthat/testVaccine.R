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

test_that("We can load vaccine calendar",
          {
              data("vaccine_calendar")
              expect_that( vaccine_calendar$efficacy[1], equals( 0.7 ) )
              expect_that( length(vaccine_calendar$efficacy), equals( 7 ) )
              expect_that( ncol(vaccine_calendar$calendar), equals( 21 ) )
              expect_that( nrow(vaccine_calendar$calendar), equals( 123 ) )
          }
)

test_that("We can load mcmc sample",
          {
              data("mcmcsample")
              expect_that( length(mcmcsample[["contact_ids"]]), equals( 597 ) )

          }
)


test_that("We can call scenario",
          {
              skip("Not yet implemented")
              data("age_sizes")
              data("vaccine_calendar")
              data("mcmcsample")
              expect_true( 
                vaccinationScenario( age_sizes=age_sizes[,1], 
                    vaccine_calendar=vaccine_calendar,
                    sample=mcmcsample ) 
                )
          }
)
