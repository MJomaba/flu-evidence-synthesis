context("Demography")

test_that("We can separate into risk groups", 
  {
      popv <- separate.into.risk.groups( 
                c(3000,2000,1000),
                matrix( c(0.4,0.4,0.3,0.1,0.1,0.1), ncol=3, byrow=T ) )
      expect_identical( popv,
            c( 1500, 1000, 600, 1200, 800, 300, 300, 200, 100 ) )

  })
