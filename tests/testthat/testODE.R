context("ODE")

test_that("We correctly convert week and year into a date", {
    expect_that(as.character(getTimeFromWeekYear( 35, 1974 )), equals("1974-08-26 12:00:00"))
    expect_that(as.character(getTimeFromWeekYear( 3, 1973 )), equals("1973-01-15 12:00:00"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals("1973-12-31 12:00:00"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals(as.character(getTimeFromWeekYear( 1, 1974 ))))
    expect_that(as.character(getTimeFromWeekYear( 3, 1972 )), equals("1972-01-17 12:00:00"))
})

test_that("Runge Kutta ode solver works correctly",
{
    if (!exists(".runRKF"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".runRKF not available (hidden)")
    }
    #library("deSolve")
    #predprey <- function(t,y,prs) 
    #{
    #  r <- 1.5
    #  a <- 1
    #  g <- 1
    #  m <- 3
    #  list(c(r*y[1]-a*y[1]*y[2],g*y[1]*y[2]-m*y[2]))
    #}
    #yini <- c(10,5)
    #sol <- deSolve::ode(y=yini, times=seq(0,20,0.1),func=predprey,parms=c())

    sol <- .runStep()
    rkf <- .runRKF(step_size=0.1)
    expect_less_than( sum(abs(sol[,2:3]-rkf[,2:3]))/nrow(sol),0.05 )
})

test_that("Runge Kutta ode solver works with a large initial integration dt",
{
    if (!exists(".runRKF"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".runRKF not available (hidden)")
    }
    sol <- .runStep(step_size=0.1)
    rkf <- .runRKF(step_size=0.1,h_step=1.0)
    expect_less_than( sum(abs(sol[,2:3]-rkf[,2:3]))/nrow(sol),0.05 )
})

test_that("Return the correct ODE results", {
    data("age_sizes")
    data("vaccine_calendar")
    data("polymod_uk")
    data("mcmcsample")
    odes <- infection.model( age_sizes[,1],
                            vaccine_calendar,
                            as.matrix(polymod_uk[mcmcsample$contact_ids+1,]),
                            mcmcsample$parameters$susceptibility,
                            mcmcsample$parameters$transmissibility,
                            mcmcsample$parameters$init_pop,
                            c(0.8,1.8) )
    # Check all sums
    sums <- c(119002.481089555, 633877.330662713, 1463334.37622813, 2870258.66403954, 4391396.3601774, 2786119.98841863, 496487.103476314, 2531.54735113481, 36283.6707667877, 155096.824666694, 263690.71146717, 415426.823720658, 528526.506836539, 358207.414064678)
    odeSums <- colSums(odes[,2:15])
    for( i in 1:length(sums) )
    {
        ratio <- odeSums[i]/sums[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }

    # Check max value
    maxs <- c(3132.20022124067, 16720.2726361168, 38489.4726884184, 77065.0627636948, 115970.502880153, 73496.0794524567, 13004.8663674856, 66.6225150332303, 956.806903893294, 4077.66057862465, 7073.81188424645, 10955.7871500801, 13891.6328086064, 9360.94718096195)
    odeMaxs <- sapply( seq(2,15), function(x) max(odes[,x]) ) 
    for( i in 1:length(maxs) )
    {
        ratio <- odeMaxs[i]/maxs[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }
    # Check time of max
    maxts <- c(118, 117, 118, 114, 117, 117, 118, 118, 117, 118, 114, 117, 117, 118)
    odeMaxts <- sapply( seq(2,15), function(x) which.max(odes[,x]) )
    for( i in 1:length(maxts) )
    {
        dt <- abs(odeMaxts[i]-maxts[i])
        expect_less_than( dt, 5 )
    }
})

test_that("Return the correct ODE results, when using a different interval", {
    data("age_sizes")
    data("vaccine_calendar")
    data("polymod_uk")
    data("mcmcsample")
    odes <- infection.model( age_sizes[,1],
                            vaccine_calendar,
                            as.matrix(polymod_uk[mcmcsample$contact_ids+1,]),
                            mcmcsample$parameters$susceptibility,
                            mcmcsample$parameters$transmissibility,
                            mcmcsample$parameters$init_pop,
                            c(0.8,1.8),
                  7 )

    # Check all sums
    sums <- c(119002.481089555, 633877.330662713, 1463334.37622813, 2870258.66403954, 4391396.3601774, 2786119.98841863, 496487.103476314, 2531.54735113481, 36283.6707667877, 155096.824666694, 263690.71146717, 415426.823720658, 528526.506836539, 358207.414064678)
    odeSums <- colSums(odes[,2:15])
    for( i in 1:length(sums) )
    {
        ratio <- odeSums[i]/sums[i]
        expect_less_than( ratio, 1.05 )
        expect_more_than( ratio, 0.97 )
    }

    # Check max value
    maxs <- c(3132.20022124067, 16720.2726361168, 38489.4726884184, 77065.0627636948, 115970.502880153, 73496.0794524567, 13004.8663674856, 66.6225150332303, 956.806903893294, 4077.66057862465, 7073.81188424645, 10955.7871500801, 13891.6328086064, 9360.94718096195)
    odeMaxs <- sapply( seq(2,15), function(x) max(odes[,x]) ) 
    for( i in 1:length(maxs) )
    {
        ratio <- odeMaxs[i]/(7*maxs[i])
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }
    # Check time of max
    maxts <- c(118, 117, 118, 114, 117, 117, 118, 118, 117, 118, 114, 117, 117, 118)
    odeMaxts <- sapply( seq(2,15), function(x) which.max(odes[,x]) )
    for( i in 1:length(maxts) )
    {
        dt <- abs(odeMaxts[i]*7-maxts[i]) # Take into account different time
        expect_less_than( dt, 6 )
    }
})

test_that("ODE works correctly with the new vaccine date vector", {
    data("age_sizes")
    data("polymod_uk")
    data("mcmcsample")
    data("vaccine_calendar")
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
    odes <- infection.model( age_sizes[,1],
                            test.vac,
                            as.matrix(polymod_uk[mcmcsample$contact_ids+1,]),
                            mcmcsample$parameters$susceptibility,
                            mcmcsample$parameters$transmissibility,
                            mcmcsample$parameters$init_pop,
                            c(0.8,1.8) )

    # Check all sums
    sums <- c(119002.481089555, 633877.330662713, 1463334.37622813, 2870258.66403954, 4391396.3601774, 2786119.98841863, 496487.103476314, 2531.54735113481, 36283.6707667877, 155096.824666694, 263690.71146717, 415426.823720658, 528526.506836539, 358207.414064678)
    odeSums <- colSums(odes[,2:15])
    for( i in 1:length(sums) )
    {
        ratio <- odeSums[i]/sums[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }

    # Check max value
    maxs <- c(3132.20022124067, 16720.2726361168, 38489.4726884184, 77065.0627636948, 115970.502880153, 73496.0794524567, 13004.8663674856, 66.6225150332303, 956.806903893294, 4077.66057862465, 7073.81188424645, 10955.7871500801, 13891.6328086064, 9360.94718096195)
    odeMaxs <- sapply( seq(2,15), function(x) max(odes[,x]) ) 
    for( i in 1:length(maxs) )
    {
        ratio <- odeMaxs[i]/maxs[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }
    # Check time of max
    maxts <- c(118, 117, 118, 114, 117, 117, 118, 118, 117, 118, 114, 117, 117, 118)
    odeMaxts <- sapply( seq(2,15), function(x) which.max(odes[,x]) )
    for( i in 1:length(maxts) )
    {
        dt <- abs(odeMaxts[i]-maxts[i])
        expect_less_than( dt, 5 )
    }
})

test_that("We can set population sizes etc with infectionODEs", {
    data("age_sizes")
    data("polymod_uk")
    data("mcmcsample")
    data("vaccine_calendar")
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

    initial.infected <- rep( 10^mcmcsample$parameters$init_pop, 7 )
    initial.infected <- separate.into.risk.groups(
              initial.infected, risk.ratios );

    # Need to separate into age groups... 
    odes <- infectionODEs( popv, initial.infected,
                            test.vac,
                            contact.matrix( as.matrix(polymod_uk[mcmcsample$contact_ids+1,]), age_sizes[,1], c(1,5,15,25,45,65) ),
                            mcmcsample$parameters$susceptibility,
                            mcmcsample$parameters$transmissibility,
                            c(0.8,1.8), 1 )

    # Check all sums
    sums <- c(119002.481089555, 633877.330662713, 1463334.37622813, 2870258.66403954, 4391396.3601774, 2786119.98841863, 496487.103476314, 2531.54735113481, 36283.6707667877, 155096.824666694, 263690.71146717, 415426.823720658, 528526.506836539, 358207.414064678)
    odeSums <- colSums(odes[,2:15])
    for( i in 1:length(sums) )
    {
        ratio <- odeSums[i]/sums[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }

    # Check max value
    maxs <- c(3132.20022124067, 16720.2726361168, 38489.4726884184, 77065.0627636948, 115970.502880153, 73496.0794524567, 13004.8663674856, 66.6225150332303, 956.806903893294, 4077.66057862465, 7073.81188424645, 10955.7871500801, 13891.6328086064, 9360.94718096195)
    odeMaxs <- sapply( seq(2,15), function(x) max(odes[,x]) ) 
    for( i in 1:length(maxs) )
    {
        ratio <- odeMaxs[i]/maxs[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
    }
    # Check time of max
    maxts <- c(118, 117, 118, 114, 117, 117, 118, 118, 117, 118, 114, 117, 117, 118)
    odeMaxts <- sapply( seq(2,15), function(x) which.max(odes[,x]) )
    for( i in 1:length(maxts) )
    {
        dt <- abs(odeMaxts[i]-maxts[i])
        expect_less_than( dt, 5 )
    }
})

test_that("Second risk group works as expected", {
     data("age_sizes")
    data("polymod_uk")
    data("mcmcsample")
    data("vaccine_calendar")
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
    test.vac[["calendar"]][,15:21] <- test.vac[["calendar"]][,8:14] # Use same calendar for Preg as for high risk


    age.groups <- separate.into.age.groups( age_sizes[,1], 
                                           c(1,5,15,25,45,65) )

    risk.ratios <- matrix( c(
        0, 0, 0, 0, 0, 0, 0,                          
        0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45
                          ), ncol=7, byrow=T )

    popv <- separate.into.risk.groups(
              age.groups, risk.ratios );

    initial.infected <- rep( 10^mcmcsample$parameters$init_pop, 7 )
    initial.infected <- separate.into.risk.groups(
              initial.infected, risk.ratios );

    # Need to separate into age groups... 
    odes <- infectionODEs( popv, initial.infected,
                            test.vac,
                            contact.matrix( as.matrix(polymod_uk[mcmcsample$contact_ids+1,]), age_sizes[,1], c(1,5,15,25,45,65) ),
                            mcmcsample$parameters$susceptibility,
                            mcmcsample$parameters$transmissibility,
                            c(0.8,1.8), 1 )

    # Check all sums
    sums <- c(119002.481089555, 633877.330662713, 1463334.37622813, 2870258.66403954, 4391396.3601774, 2786119.98841863, 496487.103476314, 0,0,0,0,0,0,0, 2531.54735113481, 36283.6707667877, 155096.824666694, 263690.71146717, 415426.823720658, 528526.506836539, 358207.414064678)
    odeSums <- colSums(odes[,2:22])
    for( i in 1:length(sums) )
    {
        if( odeSums[i] != 0 & sums[i] != 0 )
        {
            ratio <- odeSums[i]/sums[i]
            expect_less_than( ratio, 1.03 )
            expect_more_than( ratio, 0.97 )
        }
    }

    # Check max value
    maxs <- c(3132.20022124067, 16720.2726361168, 38489.4726884184, 77065.0627636948, 115970.502880153, 73496.0794524567, 13004.8663674856, 0,0,0,0,0,0,0, 66.6225150332303, 956.806903893294, 4077.66057862465, 7073.81188424645, 10955.7871500801, 13891.6328086064, 9360.94718096195)
    odeMaxs <- sapply( seq(2,22), function(x) max(odes[,x]) ) 
    for( i in 1:length(maxs) )
    {
      if( odeMaxs[i] != 0 & maxs[i] != 0 )
      {
        ratio <- odeMaxs[i]/maxs[i]
        expect_less_than( ratio, 1.03 )
        expect_more_than( ratio, 0.97 )
      }
    }
    # Check time of max
    maxts <- c(118, 117, 118, 114, 117, 117, 118,
               0,0,0,0,0,0,0,
               118, 117, 118, 114, 117, 117, 118)
    odeMaxts <- sapply( seq(2,22), function(x) which.max(odes[,x]) )
    for( i in 1:length(maxts) )
    {
        dt <- abs(odeMaxts[i]-maxts[i])
        expect_less_than( dt, 5 )
    }
})


test_that("infectionODEs works with less than 3 risk groups", {
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

    odes <- infectionODEs( population, initial.infected, 
                vaccine_calendar,  contacts, 
                susceptibility, transmissibility, infection_delays, 7 )
    expect_equal(ncol(odes), 5)
    expect_equal(format(odes$Time[1],format="%Y"),"2010" );
    
    comp <- mapply( function(x,y) difftime(y,x)-7, odes$Time[1:(nrow(odes)-1)], odes$Time[2:nrow(odes)] )

    # Up to two time diffs will be slightly more/less than 7 days due to summer time
    expect_less_than(sum(comp!=0),3)
    expect_less_than( sum(abs(comp)), 2.1/24 )
})

