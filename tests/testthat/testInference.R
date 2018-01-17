context("Inference")

test_that("Likelihood function returns the correct value",
  {
      if (!exists(".infection.model"))
      {
          # In certain situation this function is hidden 
          # (i.e. for devtools::check(), but not for devtools::test()

          skip(".infection.model not available (hidden)")
      }
      data("demography")
      data("vaccine_calendar")
      data("polymod_uk")
      data("mcmcsample")
      data("ili")
      data("confirmed.samples")

      week <- 18
      age.group <- 3
      population <- sum(demography[16:45]) # age group 3 is from 16 to 45

      inf.model <- .infection.model(demography,vaccine_calendar,
                                   as.matrix(polymod_uk[mcmcsample$contact_ids+1,]),
                                   mcmcsample$parameters$susceptibility,
                                   mcmcsample$parameters$transmissibility,
                                   mcmcsample$parameters$init_pop,
                                   c(0.8,1.8),
                                   7)

      # Convert results that are 7 groups to age group 3 (out of 5),
      # i.e group 4 and 5 is equal to group 3
      predicted <- rowSums(inf.model[week,5:6])+rowSums(inf.model[week,12:13])

      ll <- log_likelihood_cases(
                              c(mcmcsample$parameters$epsilon[age.group]),
                              mcmcsample$parameters$psi,
                              as.matrix(predicted), c(population),
                              as.matrix(ili$ili[week,age.group]),
                              as.matrix(ili$total.monitored[week,age.group]),
                              as.matrix(confirmed.samples$positive[week,age.group]),
                              as.matrix(confirmed.samples$total.samples[week,age.group]))
      expect_lt(ll, 121.9881)
      expect_gt(ll, 120.9881)
  }
)

test_that("We can run inference", 
  {
      library(moments)
      data("demography")
      data("vaccine_calendar")
      data("polymod_uk")
      data("ili")
      data("confirmed.samples")

      expect_equal( length(vaccine_calendar$efficacy), 7 )

      set.seed(100)
      results <- inference(demography = demography,
                           vaccine_calendar=vaccine_calendar,
                           polymod_data=as.matrix(polymod_uk),
                           initial=c(0.01188150,0.01831852,0.05434378,
                             1.049317e-05,0.1657944,
                             0.3855279,0.9269811,0.5710709,
                             -0.1543508), 
                           ili=ili$ili,
                           mon_pop=ili$total.monitored,
                           n_pos=confirmed.samples$positive,
                           n_samples=confirmed.samples$total.samples,
                           nbatch=1000,
                           nburn=1000, blen=1 )

      expect_that( nrow(results$batch), equals( 1000 ) )
      expect_that( nrow(results$contact.ids), equals( 1000 ) )
      expect_true(all(results$contact.ids > 0))
      expect_true(all(results$contact.ids <= nrow(polymod_uk)))

      # Contact ids are mixing
      expect_false(identical(results$contact.ids[1,], results$contact.ids[1000,]))
      #mean
      m1 <- moment(results$llikelihoods,central=FALSE)
      expect_lt(m1, 2268 )
      expect_gt(m1, 2242 )
      m2 <- moment(results$llikelihoods,central=TRUE,2)
      expect_lt(m2, 5.5 )
      expect_gt(m2, 0.9 )
      m3 <- moment(results$llikelihoods,central=TRUE,3)
      expect_lt(m3, 1.0 )
      expect_gt(m3, -7.0 )
  }
)

test_that("dmultinom and dmultinom.cpp return same value", 
    {
        dp <- dmultinom( c(5,4,3), 12, c(0.4, 0.5, 0.1) )
        expect_identical( dmultinom.cpp( c(5,4,3), 12, c(0.4,0.5,0.1) ), dp )
    }
)

test_that("Likelihood works for large values", {
  expect_gt(log_likelihood_cases(c(0.00622462018361167),0.0414822583510223,matrix(170776.505911481),31742426,matrix(589),matrix(383614),matrix(36),matrix(100), depth = 8),
        426)    
  expect_lt(log_likelihood_cases(c(0.00622462018361167),0.0414822583510223,matrix(170776.505911481),31742426,matrix(589),matrix(383614),matrix(36),matrix(100), depth = 8),
        428)    
  expect_gt(log_likelihood_cases(c(0.00622462018361167),0.0414822583510223,matrix(170776.505911481),31742426,matrix(589),matrix(383614),matrix(136),matrix(210), depth = 8)
         , -10e10)
})

test_that("Group mapping works correctly", {
  # Create scenario with two age groups and two risk groups (same size)
})