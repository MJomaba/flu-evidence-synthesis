context("Demography")

test_that("We can separate into risk groups", 
  {
      popv <- separate.into.risk.groups( 
                c(1000,1000,1000),
                matrix( c(0.4,0.4,0.4,0.2,0.2,0.2), ncol=3 ) )
      expect_identical( popv,
            c( 400, 400, 400, 400, 400, 400, 200, 200, 200 ) )

  })
