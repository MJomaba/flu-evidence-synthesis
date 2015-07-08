context("Inference")

test_that("We can run inference", 
          {
              skip( "Currently often hangs. Needs to be reenabled" )
              data("age_sizes")
              data("vaccine_calendar")
              data("polymod_uk")
              data("mcmcsample")
              data("ili")
              data("confirmation")

              set.seed(100)
              results <- inference( age_sizes=age_sizes[,1],
                                   vaccine_calendar=vaccine_calendar,
                                   polymod_data=as.matrix(polymod_uk),
                                   init_state=mcmcsample, 
                                   ili=as.matrix(ili),
                                   mon_pop=as.matrix(mon_pop),
                                   n_pos=as.matrix(n_pos),
                                   n_samples=as.matrix(n_samples),
                                   mcmc_chain_length=1000,
                                   burn_in=1000, thinning=1 )

              expect_that( length(results), equals( 1000 ) )
          }
)
 
