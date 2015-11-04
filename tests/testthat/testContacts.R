context("Contacts")

test_that("We can create a contact matrix", 
  {
      data(age_sizes)
      data(polymod_uk)
      age.group.limits <- c( 1,5,15,25,45,65 )
      cm <- contact.matrix( as.matrix(polymod_uk), age_sizes[,1], 
                            age.group.limits );
      expect_false( cm[2,2] == Inf )

      #cm <- contract.matrix( as.matrix(polymod_uk), age_sizes[1,] );
      #expect_equal( cm[1,1], 1 )
  }
)


test_that("We can use custom age limits", 
  {
      skip("Not implemented yet")
      expect_equal( 0, 1 )
  }
)

