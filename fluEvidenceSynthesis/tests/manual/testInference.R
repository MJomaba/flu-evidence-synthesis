context("Inference")

test_that("We can run inference", 
    {
        library(moments)
        data("age_sizes")
        data("vaccine_calendar")
        data("polymod_uk")
        data("mcmcsample")
        data("ili")
        data("confirmed.samples")

        set.seed(100)
        results <- inference( age_sizes=age_sizes[,1],
                              vaccine_calendar=vaccine_calendar,
                              polymod_data=as.matrix(polymod_uk),
                              init_state=mcmcsample, 
                              ili=ili$ili,
                              mon_pop=ili$total.monitored,
                              n_pos=confirmed.samples$positive,
                              n_samples=confirmed.samples$total.samples,
                              mcmc_chain_length=10000000,
                              burn_in=1000000, thinning=10000 )

        m1 <- moment(sapply(results, function(x) x$likelihood),central=FALSE)
        expect_less_than(m1, 2266 )
        expect_more_than(m1, 2264 )
        m2 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,2)
        expect_less_than(m2, 12 )
        expect_more_than(m2, 9 )
        m3 <- moment(sapply(results, function(x) x$likelihood),central=TRUE,3)
        expect_less_than(m3, -4.5 )
        expect_more_than(m2, -6 )
    }
)
 
