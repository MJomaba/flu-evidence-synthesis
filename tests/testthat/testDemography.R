context("Demography")

test_that("We can separate into risk groups", 
  {
      popv <- stratify_by_risk( 
                c(3000,2000,1000),
                matrix( c(0.4,0.4,0.3,0.1,0.1,0.1), ncol=3, byrow=T ) )
      expect_identical( popv,
            c( 1500, 1000, 600, 1200, 800, 300, 300, 200, 100 ) )

  })

test_that("We can convert age to age group", {
  expect_equal(as.numeric(as_age_group( 6, c(1,5,10) )),3)
  expect_equal(as.numeric(as_age_group( 5, c(1,5,10) )),3)
  expect_equal(as.numeric(as_age_group( 4, c(1,5,10) )),2)
  expect_equal(as.numeric(as_age_group( 0, c(1,5,10) )),1)
  expect_equal(as.numeric(as_age_group( 10, c(1,5,10) )),4)
})

test_that("We can convert age vector to age group", {
  v <- as_age_group(c(4,5), c(1,5,10))
  expect_equal(as.numeric(v), c(2, 3))
  expect_is(v, "factor")
  expect_equal(levels(v), c("[0,1)", "[1,5)", "[5,10)", "[10,+)"))
  
  v <- as_age_group(c(4,11), c(1,5,10))
  expect_equal(as.numeric(v), c(2, 4))
  
  v <- as_age_group(11, c(1,5,10))
  expect_equal(as.numeric(v), 4)
})
 
test_that("We can map age groups", {
  mp <- age_group_mapping(c(1,5), c(5))
  expect_identical(as.numeric(mp$from), c(1, 2, 3))
  expect_identical(as.numeric(mp$to), c(1, 1, 2))
  expect_identical(mp$weight, c(1,1,1))
  
  mp <- age_group_mapping(c(1,5,8), c(5))
  expect_identical(as.numeric(mp$from), c(1, 2, 3, 4))
  expect_identical(as.numeric(mp$to), c(1, 1, 2, 2))
  expect_identical(mp$weight, c(1, 1, 1, 1))
  
  mp <- age_group_mapping(c(1,3,4), c(2), demography = c(100, 200, 400, 800, 1600))
  expect_identical(as.numeric(mp$from), c(1, 2, 2, 3, 4))
  expect_identical(as.numeric(mp$to), c(1, 1, 2, 2, 2))
  expect_identical(mp$weight, c(1, 200/600, 400/600, 1, 1))
  
  mp <- age_group_mapping(c(1,5), c(3), demography = c(100, 200, 400, 800, 1600, 3200))
  expect_identical(as.numeric(mp$from), c(1, 2, 2, 3))
  expect_identical(as.numeric(mp$to), c(1, 1, 2, 2))
  expect_identical(mp$weight, c(1, 600/(600+2400), 2400/(600+2400), 1))
})