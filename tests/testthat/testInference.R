context("Inference")

test_that("Likelihood function returns the correct value",
  {
      data("age_sizes")
      data("vaccine_calendar")
      data("polymod_uk")
      data("mcmcsample")
      data("ili")
      data("confirmed.samples")

      week <- 18
      age.group <- 3
      population <- sum(age_sizes$V1[16:45]) # age group 3 is from 16 to 45

      inf.model <- infection.model(age_sizes$V1,vaccine_calendar,
                                   as.matrix(polymod_uk[mcmcsample$contact_ids+1,]),
                                   mcmcsample$parameters$susceptibility,
                                   mcmcsample$parameters$transmissibility,
                                   mcmcsample$parameters$init_pop,
                                   c(0.8,1.8),
                                   7)

      # Convert results that are 7 groups to age group 3 (out of 5),
      # i.e group 4 and 5 is equal to group 3
      predicted <- rowSums(inf.model[week,5:6])+rowSums(inf.model[week,12:13])

      ll <- llikelihood.cases(mcmcsample$parameters$epsilon[age.group],
                              mcmcsample$parameters$psi,
                              predicted, population,
                              ili$ili[week,age.group],
                              ili$total.monitored[week,age.group],
                              confirmed.samples$positive[week,age.group],
                              confirmed.samples$total.samples[week,age.group])
      expect_less_than(ll, 121.9881)
      expect_more_than(ll, 120.9881)
  }
)

test_that("We can run inference", 
  {
      library(moments)
      data("age_sizes")
      data("vaccine_calendar")
      data("polymod_uk")
      data("ili")
      data("confirmed.samples")

      expect_equal( length(vaccine_calendar$efficacy), 7 )

      set.seed(100)
      results <- inference( age_sizes=age_sizes$V1,
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
      # Contact ids are mixing
      expect_false(identical(results$contact.ids[1,], results$contact.ids[1000,]))
      #mean
      m1 <- moment(results$llikelihoods,central=FALSE)
      expect_less_than(m1, 2268 )
      expect_more_than(m1, 2242 )
      m2 <- moment(results$llikelihoods,central=TRUE,2)
      expect_less_than(m2, 3.5 )
      expect_more_than(m2, 0.9 )
      m3 <- moment(results$llikelihoods,central=TRUE,3)
      expect_less_than(m3, -1.0 )
      expect_more_than(m3, -5.0 )
  }
)
