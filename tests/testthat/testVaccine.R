context("Vaccine")

test_that("We can load age data", 
          {
              data("age_sizes")
              expect_that( nrow(age_sizes), equals(86) )
          }
)

test_that("We can load contact data", 
          {
              data("polymod_uk")
              expect_that( nrow(polymod_uk), equals(597) )
              expect_that( ncol(polymod_uk), equals(9) )
              expect_that( polymod_uk[2,2], equals(1) )
          }
)

test_that("We can load vaccine calendar",
          {
              data("vaccine_calendar")
              expect_that( vaccine_calendar$efficacy[1], equals( 0.7 ) )
              expect_that( length(vaccine_calendar$efficacy), equals( 7 ) )
              expect_that( ncol(vaccine_calendar$calendar), equals( 21 ) )
              expect_that( nrow(vaccine_calendar$calendar), equals( 123 ) )
          }
)

test_that("We can load mcmc sample",
          {
              data("mcmcsample")
              expect_that( length(mcmcsample[["contact_ids"]]), equals( 597 ) )

          }
)

test_that("We can call scenario",
          {
              data("age_sizes")
              data("vaccine_calendar")
              data("inference.results")
              data("polymod_uk")

              total_size <- vaccinationScenario( age_sizes=age_sizes[,1], 
                    vaccine_calendar=vaccine_calendar,
                    polymod_data=as.matrix(polymod_uk),
                    contact_ids = inference.results$contact.ids[1000,],
                    parameters = inference.results$batch[1000,]
                    )

              exp.total.size <- c(46500.6746557518, 303225.043210381, 649632.503513013, 2471269.51343945, 4089296.63404964, 2176800.79303793, 421066.463090916, 989.173167588841, 17354.8114395601, 68842.3418775971, 226975.890171658, 386713.740264845, 412465.102974455, 303654.162949902, 0, 0, 0, 0, 0, 0, 0)
              #exp.total.size <- c( 119474.807, 636440.213, 1469251.696, 2881440.126, 4408817.625, 2796622.448, 498176.149, 2541.504, 36427.619, 155706.821, 264664.169, 416935.424, 530054.119, 359250.707, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 )
              for( i in 1:length(total_size) )
              {
                  if (total_size[i] == 0)
                      expect_that( exp.total.size[i], equals(total_size[i]) )
                  else
                  {
                      ratio <- exp.total.size[i]/total_size[i]
                      expect_lt( ratio, 1.03 )
                      expect_gt( ratio, 0.97 )
                  }
              }
          })

test_that("Vaccination scenario and infectionODEs give similar results", 
{
  data("age_sizes")
  data("vaccine_calendar")
  data("inference.results")
  data("polymod_uk")

  total_size <- vaccinationScenario( age_sizes=age_sizes[,1], 
        vaccine_calendar=vaccine_calendar,
        polymod_data=as.matrix(polymod_uk),
        contact_ids = inference.results$contact.ids[1000,],
        parameters = inference.results$batch[1000,]
        )

  test.vac <- vaccine_calendar
  # We would take first of the month, but because we want to directly
  # compare with the old way we use their dates
  #test.vac[["dates"]] <- c(as.Date("1999-10-01"), as.Date("1999-11-01"),
  #                         as.Date("1999-12-01"), as.Date("2000-01-01"),
  #                         as.Date("2000-02-01"))
  test.vac[["dates"]] <- c(as.Date("1970-10-07"), as.Date("1970-11-07"),
                           as.Date("1970-12-07"), as.Date("1971-01-07"),
                           as.Date("1971-02-07"))
  test.vac[["calendar"]] <- matrix(c(test.vac[["calendar"]][1,],
                                     test.vac[["calendar"]][32,],
                                     test.vac[["calendar"]][62,],
                                     test.vac[["calendar"]][93,]),ncol=21,byrow=TRUE)
  
  age.groups <- separate.into.age.groups( age_sizes[,1], 
                                          c(1,5,15,25,45,65) )
  
  risk.ratios <- matrix( c(
    0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
    0, 0, 0, 0, 0, 0, 0                          
  ), ncol=7, byrow=T )
  
  popv <- separate.into.risk.groups(
    age.groups, risk.ratios );
  
  initial.infected <- rep( 10^inference.results$batch[1000,9], 7 )
  
  initial.infected <- separate.into.risk.groups(
    initial.infected, risk.ratios );
  
  # Need to separate into age groups... 
  odes <- infectionODEs( popv, initial.infected,
                         test.vac,
                         contact.matrix( as.matrix(polymod_uk[inference.results$contact.ids[1000,]+1,]),
                                         age_sizes[,1], c(1,5,15,25,45,65) ),
                         c(rep(inference.results$batch[1000,6],3),rep(inference.results$batch[1000,7],3),inference.results$batch[1000,8]),
                         inference.results$batch[1000,5],
                         c(0.8,1.8), 1 )
  expect_lt(abs(sum(total_size-colSums(odes[,2:ncol(odes)]))),1e-3)
})

test_that("We can specify efficacy by risk group", {
    data(age_sizes) 
    data(polymod_uk)

    ag <- separate.into.age.groups(age_sizes$V1, limits=c(65)) # c( 43670500, 8262600 )
    population <- separate.into.risk.groups( ag, matrix(c(0.01,0.4),nrow=1) ) # c( 43233795, 4957560, 436705, 3305040 )
    ag <- c(1000,1000)
    initial.infected <- separate.into.risk.groups( ag, matrix(c(0.01,0.4),nrow=1) )
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )

    # Polymod data is subdivided in seven age groups
    poly <- polymod_uk[,c(1,2,3,9)]
    poly[,3] <- rowSums(polymod_uk[,3:8])

    contacts <- contact.matrix(as.matrix(poly), age_sizes$V1, c(65))
    susceptibility <- c( 0.7, 0.3 ) # Different for different ages
    transmissibility <- 0.17 # Same for all ages
    infection_delays <- c( 0.8, 1.8 ) # 0.8 and 1.8 day.

    odes1 <- infectionODEs( population, initial.infected, 
                vaccine_calendar,  contacts, 
                susceptibility, transmissibility, infection_delays, 7 )
    expect_equal(ncol(odes1), 5)
    expect_equal(format(odes1$Time[1],format="%Y"),"2010" );
    
    comp <- mapply( function(x,y) difftime(y,x)-7, odes1$Time[1:(nrow(odes1)-1)], odes1$Time[2:nrow(odes1)] )

    # Up to two time diffs will be slightly more/less than 7 days due to summer time
    expect_lt(sum(comp!=0),3)
    expect_lt( sum(abs(comp)), 2.1/24 )
    
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3,0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    odes2 <- infectionODEs( population, initial.infected, 
                            vaccine_calendar,  contacts, 
                            susceptibility, transmissibility, infection_delays, 7 )
    
    expect_equal(odes1,odes2)
    
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3,0.0,0.0),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    odes3 <- infectionODEs( population, initial.infected, 
                            vaccine_calendar,  contacts, 
                            susceptibility, transmissibility, infection_delays, 7 )
    
    expect_gt(sum(odes3[,2:ncol(odes2)]),sum(odes2[,2:ncol(odes2)]))
    
    vaccine_calendar <- list(
      "efficacy" = c(0.0,0.0,0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    odes4 <- infectionODEs( population, initial.infected, 
                            vaccine_calendar,  contacts, 
                            susceptibility, transmissibility, infection_delays, 7 )
    
    expect_gt(sum(odes4[,2:ncol(odes2)]),sum(odes2[,2:ncol(odes2)]))
    expect_gt(sum(odes4[,2:ncol(odes2)]),sum(odes3[,2:ncol(odes2)]))
    
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3,0.0,0.3),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    odes5 <- infectionODEs( population, initial.infected, 
                            vaccine_calendar,  contacts, 
                            susceptibility, transmissibility, infection_delays, 7 )
    
    expect_equal( as.numeric(which.max(
      colSums(odes5[,2:ncol(odes5)])/colSums(odes1[,2:ncol(odes1)]))), 3 )
})


test_that("as.vaccination.calendar should normalize given vaccination calendars", {
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.001,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    
    vc <- as.vaccination.calendar(legacy = vaccine_calendar)
    expect_equal( vc$efficacy, c(0.7, 0.3, 0.7, 0.3, 0, 0) )
    expect_equal( ncol(vc$calendar), 6 )
    expect_equal( vc$calendar[nrow(vc$calendar),], rep(0,6) )
})

test_that("as.vaccination.calendar fails when more than 100% gets vaccinated (if sum rates > 1)", {
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.02,0.007),nrow=1),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2011-02-01")) # begin and end date
    )
    
    expect_error(as.vaccination.calendar(legacy = vaccine_calendar))
    
    vaccine_calendar <- list(
      "efficacy" = c(0.7,0.3),
      "calendar" = matrix(c(0,0.007,0.02,0.007,0,0.007,0.02,0.007),nrow=2,byrow=T),
      "dates" =  c(as.Date("2010-10-01"), as.Date("2010-12-30"), as.Date("2011-02-01")) # begin and end date
    )
    
    expect_error(as.vaccination.calendar(legacy = vaccine_calendar))
})
 
    # Add test with legacy vaccine.calendar
test_that("as.vaccination.calendar correctly converts legacy vaccine.calendar ", {
  data( vaccine_calendar )
  expect_equal( nrow(vaccine_calendar$calendar), 123 )
  expect_null( vaccine_calendar$dates )
    
  vc <- as.vaccination.calendar(legacy = vaccine_calendar)
  expect_equal( length(vc$dates), 5 )
  expect_equal(vc$efficacy, rep(vaccine_calendar$efficacy,3))
})
