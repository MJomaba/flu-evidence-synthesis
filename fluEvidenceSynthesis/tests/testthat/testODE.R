context("ODE")

test_that("We correctly convert week and year into a date", {
    expect_that(as.character(getTimeFromWeekYear( 35, 1974 )), equals("1974-08-26"))
    expect_that(as.character(getTimeFromWeekYear( 3, 1973 )), equals("1973-01-15"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals("1973-12-31"))
    expect_that(as.character(getTimeFromWeekYear( 53, 1973 )), equals(as.character(getTimeFromWeekYear( 1, 1974 ))))
    expect_that(as.character(getTimeFromWeekYear( 3, 1972 )), equals("1972-01-17"))
})

test_that("Return the correct ODE results", {
    data("age_sizes")
    data("vaccine_calendar")
    data("polymod_uk")
    data("mcmcsample")
    odes <- runSEIRModel( age_sizes[,1],
                  vaccine_calendar,
                  as.matrix(polymod_uk),
                  mcmcsample )
    # Highest densities around 125 days, so we better check a couple around then
    #paste("c(", paste(as.character(odes[100,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[100,]), equals(c(1160.29815611962, 6651.54303696748, 14762.0587218382, 38640.5640137718, 48205.7611098923, 30699.171008947, 5038.38829144536, 24.6968035807556, 381.147892544313, 1567.18171893041, 3556.91858512515, 4580.56072507346, 5891.52233256964, 3662.10109127876, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[125,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[125,]), equals(c(2972.889441366, 15475.3586008999, 36156.6628793019, 62641.3500626147, 105472.377345004, 66785.4625576766, 12216.3713066405, 63.2276085737547, 885.358827022898, 3829.35555657518, 5745.47153673136, 9953.91418158275, 12589.5012810659, 8782.03691110781, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[150,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[150,]), equals(c(527.13439460105, 2579.40274757502, 6259.43155341605, 8326.56723514911, 16917.3546437955, 10694.9703021943, 2120.51713280515, 11.2077508736787, 147.480088198045, 662.333594129373, 762.661300126587, 1592.39345085831, 2002.10924865605, 1519.26249860887, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[175,]),collapse=", "), ")", sep="")
    expect_that( as.numeric(odes[175,]), equals(c(43.5300237643213, 211.490118566854, 515.817292452149, 664.651978715041, 1380.51679689959, 872.094216392817, 174.671202722151, 0.925382847053637, 12.0885410370055, 54.5560077181041, 60.8365701976313, 129.77694021797, 162.694041789986, 125.010943743409, 0, 0, 0, 0, 0, 0, 0)))
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
    odes <- runSEIRModel( age_sizes[,1],
                  test.vac,
                  as.matrix(polymod_uk),
                  mcmcsample )

    # Highest densities around 125 days, so we better check a couple around then
    #paste("c(", paste(as.character(odes[100,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[100,]), equals(c(1160.29815611962, 6651.54303696748, 14762.0587218382, 38640.5640137718, 48205.7611098923, 30699.171008947, 5038.38829144536, 24.6968035807556, 381.147892544313, 1567.18171893041, 3556.91858512515, 4580.56072507346, 5891.52233256964, 3662.10109127876, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[125,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[125,]), equals(c(2972.889441366, 15475.3586008999, 36156.6628793019, 62641.3500626147, 105472.377345004, 66785.4625576766, 12216.3713066405, 63.2276085737547, 885.358827022898, 3829.35555657518, 5745.47153673136, 9953.91418158275, 12589.5012810659, 8782.03691110781, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[150,]),collapse=", "), ")", sep="") 
    expect_that( as.numeric(odes[150,]), equals(c(527.13439460105, 2579.40274757502, 6259.43155341605, 8326.56723514911, 16917.3546437955, 10694.9703021943, 2120.51713280515, 11.2077508736787, 147.480088198045, 662.333594129373, 762.661300126587, 1592.39345085831, 2002.10924865605, 1519.26249860887, 0, 0, 0, 0, 0, 0, 0)))
    #paste("c(", paste(as.character(odes[175,]),collapse=", "), ")", sep="")
    expect_that( as.numeric(odes[175,]), equals(c(43.5300237643213, 211.490118566854, 515.817292452149, 664.651978715041, 1380.51679689959, 872.094216392817, 174.671202722151, 0.925382847053637, 12.0885410370055, 54.5560077181041, 60.8365701976313, 129.77694021797, 162.694041789986, 125.010943743409, 0, 0, 0, 0, 0, 0, 0)))

})
