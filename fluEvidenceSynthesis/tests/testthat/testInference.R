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
      expect_less_than(m1, 2265 )
      expect_more_than(m1, 2264 )
      m2 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,2)
      expect_less_than(m2, 1 )
      expect_more_than(m2, 0.9 )
      m3 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,3)
      expect_less_than(m3, -1.2 )
      expect_more_than(m3, -1.5 )
  }
)
 
