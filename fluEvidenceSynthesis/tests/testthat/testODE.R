context("ODE")

test_that("We correctly convert week and year into a date", {
    expect_that(as.character(getTimeFromWeekYear( 35, 1974 )), equals("1974-08-26"))
    expect_that(as.character(getTimeFromWeekYear( 3, 1973 )), equals("1973-01-15"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals("1973-12-31"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals(as.character(getTimeFromWeekYear( 1, 1974 ))))
    expect_that(as.character(getTimeFromWeekYear( 3, 1972 )), equals("1972-01-17"))
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
    if (!exists(".runSEIRModel"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".runSEIRModel not available (hidden)")
    }
    data("age_sizes")
    data("vaccine_calendar")
    data("polymod_uk")
    data("mcmcsample")
    odes <- .runSEIRModel( age_sizes[,1],
                  vaccine_calendar,
                  as.matrix(polymod_uk),
                  mcmcsample )

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

test_that("ODE works correctly with the new vaccine date vector", {
    if (!exists(".runSEIRModel"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".runSEIRModel not available (hidden)")
    }
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
    odes <- .runSEIRModel( age_sizes[,1],
                  test.vac,
                  as.matrix(polymod_uk),
                  mcmcsample )

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
