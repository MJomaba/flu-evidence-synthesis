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
              data("inference.results")
              data("polymod_uk")

              total_size <- vaccinationScenario( age_sizes=age_sizes[,1], 
                    vaccine_calendar=vaccine_calendar,
                    polymod_data=as.matrix(polymod_uk),
                    contact_ids = inference.results$contact.ids[1000,],
                    parameters = inference.results$batch[1000,],
                    )

              exp.total.size <- c( 119474.807, 636440.213, 1469251.696, 2881440.126, 4408817.625, 2796622.448, 498176.149, 2541.504, 36427.619, 155706.821, 264664.169, 416935.424, 530054.119, 359250.707, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 )
              for( i in 1:length(total_size) )
              {
                  if (total_size[i] == 0)
                      expect_that( exp.total.size[i], equals(total_size[i]) )
                  else
                  {
                      ratio <- exp.total.size[i]/total_size[i]
                      expect_less_than( ratio, 1.03 )
                      expect_more_than( ratio, 0.97 )
                  }
              }
          })
