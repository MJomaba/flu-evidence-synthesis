context("Contacts")

test_that("We can create a contact matrix", 
  {
      data(age_sizes)
      data(polymod_uk)
      age.group.limits <- c( 1,5,15,25,45,65 )
      cm <- contract.matrix( as.matrix(polymod_uk), age_sizes[1,], 
                            age.group.limits );
      expect_equal( cm[1,1], 1 )

      cm <- contract.matrix( as.matrix(polymod_uk), age_sizes[1,] );
      expect_equal( cm[1,1], 1 )
  }
)


test_that("We can use custom age limits", 
  {
      expect_equal( 0, 1 )
  }
)

