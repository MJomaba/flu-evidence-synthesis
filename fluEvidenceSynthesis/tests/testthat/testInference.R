context("Inference")

test_that("We can run inference", 
  {
      library(moments)
      data("age_sizes")
      data("vaccine_calendar")
      data("polymod_uk")
      data("mcmcsample")
      data("ili")
      data("confirmed.samples")

      set.seed(100)
      results <- inference( age_sizes=age_sizes$V1,
                           vaccine_calendar=vaccine_calendar,
                           polymod_data=as.matrix(polymod_uk),
                           init_sample=mcmcsample, 
                           ili=ili$ili,
                           mon_pop=ili$total.monitored,
                           n_pos=confirmed.samples$positive,
                           n_samples=confirmed.samples$total.samples,
                           mcmc_chain_length=1000,
                           burn_in=1000, thinning=1 )

      expect_that( length(results), equals( 1000 ) )
      # Contact ids are mixing
      expect_false(identical(results[[1]]$contact_ids, results[[1000]]$contact_ids))
      all.moments(sapply(results, function(x) x$likelihood))
      #mean
      m1 <- moment(sapply(results, function(x) x$likelihood),central=FALSE)
      expect_less_than(m1, 2268 )
      expect_more_than(m1, 2262 )
      m2 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,2)
      expect_less_than(m2, 3.5 )
      expect_more_than(m2, 0.9 )
      m3 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,3)
      expect_less_than(m3, -1.0 )
      expect_more_than(m3, -5.0 )
  }
)

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

      inf.model <- infection.model(age_sizes$V1,vaccine_calendar,as.matrix(polymod_uk),mcmcsample,7)

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
 
