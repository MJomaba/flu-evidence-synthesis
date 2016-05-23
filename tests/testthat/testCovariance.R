context("Covariance")

test_that("We calculate (biased) covariance correctly", 
  {
    if (!exists(".updateMeans"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".updateMeans not available (hidden)")
    }
      library(MASS)
      set.seed(100)
      realCov <- matrix(rnorm(5*5,0,0.3),ncol=5)
      for (i in 1:5)
          realCov[i,i] <- runif(1,1,5)
      realMeans <- runif(5,-5,5)
      samples <- mvrnorm(1000,realMeans,realCov)

      means <- samples[1,]
      cov <- matrix(rep(0,5*5),ncol=5)
      prev <- 0
      for (i in 2:nrow(samples))
      {
          means <- .updateMeans( 
            means, samples[i,], i )
          cov <- .updateCovariance( 
            cov, samples[i,],means,i)
          if (i == 2)
          {
              # cov() returns unbiased, we use biased
              expect_that( max(abs(2*cov-cov(samples[1:i,]))),
                        equals(0) )
              prev <- max(abs(cov-cov(samples[1:i,])))
          }
          else if (i %% 20 == 0)
          {
              # In general (except for numerical errors) the difference 
              # between biased and unbiased should decrease
              diff <- max(abs(cov-cov(samples[1:i,])))
              expect_lt( diff, prev )
              prev <- diff
          }
      }
  }
)


