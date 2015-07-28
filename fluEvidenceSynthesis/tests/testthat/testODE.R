context("ODE")

test_that("We correctly convert week and year into a date", {
    expect_that(as.character(getTimeFromWeekYear( 35, 1974 )), equals("1974-08-26"))
    expect_that(as.character(getTimeFromWeekYear( 3, 1973 )), equals("1973-01-15"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals("1973-12-31"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals(as.character(getTimeFromWeekYear( 1, 1974 ))))
    expect_that(as.character(getTimeFromWeekYear( 3, 1972 )), equals("1972-01-17"))
})

test_that("Runge Kutta ode solver works correctly")
{
    if (!exists(".runRKF"))
    {
        # In certain situation this function is hidden 
        # (i.e. for devtools::check(), but not for devtools::test()

        skip(".runRKF not available (hidden)")
    }
    library(deSolve)
    predprey <- function(t,y,prs) 
    {
      r <- 1.5
      a <- 1
      g <- 1
      m <- 3
      list(c(r*y[1]-a*y[1]*y[2],g*y[1]*y[2]-m*y[2]))
    }

    yini <- c(10,5)
    sol <- ode(y=yini, times=seq(0,20,0.1),func=predprey,parms=c(),method="ode45")
    rkf <- runRKF()
    expect_less_than( sum(abs(sol[,2:3]-rkf[,2:3]))/nrow(sol),0.05 )
}

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
    # Highest densities around 125 days, so we better check a couple around then
    #paste("c(", paste(as.character(odes[100,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[100,2:22]), equals(c(1185.42601715722, 6792.58863150519, 15078.5929324777, 39389.685725592, 49211.7694979212, 31336.4720297433, 5144.74181590203, 25.2311639178007, 389.215084386994, 1600.69354707852, 3625.56049288265, 4675.38259229678, 6011.25250068733, 3738.37677692342, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[125,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[125,2:22]), equals(c(2953.69086684003, 15357.8323022363, 35905.4817919461, 61904.8742125427, 104597.134598715, 66229.3838862123, 12131.8289518787, 62.8190393915206, 878.627906990143, 3802.70675021195, 5677.8173773949, 9870.96897884856, 12483.5243800358, 8720.80061692352, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[150,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[150,2:22]), equals(c(514.859952355205, 2518.85353671921, 6113.29524858995, 8125.35395798149, 16518.2262227089, 10442.5487538174, 2071.02875613771, 10.9467538291165, 144.017536691814, 646.866431468683, 744.224569009413, 1554.79711275651, 1954.76439001486, 1483.78466453535, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[175,]),collapse=", "), ")", sep="")
    expect_that( as.numeric(odes[175,2:22]), equals(c(42.4306250847742, 206.145425771812, 502.78743373508, 647.817533358771, 1345.61592861698, 850.046745985084, 170.25923829655, 0.902011284344771, 11.7830443139393, 53.1778897621223, 59.2956887296562, 126.496047217001, 158.580963174792, 121.85333202487, 0, 0, 0, 0, 0, 0, 0)))
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

    # Highest densities around 125 days, so we better check a couple around then
    #paste("c(", paste(as.character(odes[100,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[100,2:22]), equals(c(1185.42601715722, 6792.58863150519, 15078.5929324777, 39389.685725592, 49211.7694979212, 31336.4720297433, 5144.74181590203, 25.2311639178007, 389.215084386994, 1600.69354707852, 3625.56049288265, 4675.38259229678, 6011.25250068733, 3738.37677692342, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[125,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[125,2:22]), equals(c(2953.69086684003, 15357.8323022363, 35905.4817919461, 61904.8742125427, 104597.134598715, 66229.3838862123, 12131.8289518787, 62.8190393915206, 878.627906990143, 3802.70675021195, 5677.8173773949, 9870.96897884856, 12483.5243800358, 8720.80061692352, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[150,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[150,2:22]), equals(c(514.859952355205, 2518.85353671921, 6113.29524858995, 8125.35395798149, 16518.2262227089, 10442.5487538174, 2071.02875613771, 10.9467538291165, 144.017536691814, 646.866431468683, 744.224569009413, 1554.79711275651, 1954.76439001486, 1483.78466453535, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[175,]),collapse=", "), ")", sep="")
    expect_that( as.numeric(odes[175,2:22]), equals(c(42.4306250847742, 206.145425771812, 502.78743373508, 647.817533358771, 1345.61592861698, 850.046745985084, 170.25923829655, 0.902011284344771, 11.7830443139393, 53.1778897621223, 59.2956887296562, 126.496047217001, 158.580963174792, 121.85333202487, 0, 0, 0, 0, 0, 0, 0)))
})
