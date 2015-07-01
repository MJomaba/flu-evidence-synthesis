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
              expect_that( polymod_uk[2,2], equals(1) )
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
              data("age_sizes")
              data("vaccine_calendar")
              data("mcmcsample")
              data("polymod_uk")

              total_size <- vaccinationScenario( age_sizes=age_sizes[,1], 
                    vaccine_calendar=vaccine_calendar,
                    sample=mcmcsample,
                    polymod_uk=as.matrix(polymod_uk) )

              expect_equal( total_size,
                    c( 103841.212, 560157.389, 1284767.061, 2607491.160, 3884830.785, 2428974.589, 328566.100, 2181.114, 28518.728, 118383.065, 168843.285, 264320.534, 380895.545, 268826.809, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 )
                               )
          }
)
