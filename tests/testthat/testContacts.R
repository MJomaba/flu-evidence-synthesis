context("Contacts")

test_that("We can create a contact matrix", 
  {
      data(age_sizes)
      data(polymod_uk)
      age.group.limits <- c( 1,5,15,25,45,65 )
      cm <- contact.matrix( as.matrix(polymod_uk), age_sizes[,1], 
                           age.group.limits )
      expect_lt( sum(cm), 7.06338e-06 + 1e-8 )
      expect_gt( sum(cm), 7.06338e-06 - 1e-8 )
      for (i in 1:nrow(cm))
        for (j in 1:nrow(cm))
          expect_equal( cm[i,j], cm[j,i])

      cm2 <- contact.matrix( as.matrix(polymod_uk), age_sizes[,1] )
      expect_identical( cm, cm2 )
  }
)


test_that("We can use custom age limits", 
  {
      data(age_sizes)
      poly <- matrix( c( 5, 0, 4, 2, 
                        5, 1, 3, 1,
                        25, 0, 1, 5, 
                        25, 1, 2, 6 ), nrow = 4, byrow= TRUE)

      age.group.limits <- c( 15 )
      cm <- contact.matrix( poly, age_sizes[,1], 
                           age.group.limits )
      expect_equal( nrow(cm), 2 )
      expect_equal( ncol(cm), 2 )

      expect_equal( sum(cm>0), 4 ) # All greater than 0
      expect_lt( sum(cm), 6.682946e-07 + 1e-8 )
      expect_gt( sum(cm), 6.682946e-07 - 1e-8 )
  }
)

