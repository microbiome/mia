context("estimateDominance")

test_that("estimateDominance", {

    data(esophagus)

    #Calculates index that is not accepted, should get an error.
    expect_error((estimateDominance(esophagus, index = "shannon")))

    #Calculates DBP in two different ways, should get the same result.
    expect_true(all(estimateDominance(esophagus)@colData$DBP == estimateDominance(esophagus, index = "DBP")@colData$DBP))

    #Calculates absolute in two different ways, should get the same result.
    expect_true(all(estimateDominance(esophagus, index = NULL, relative = FALSE, rank = 1)@colData$absolute ==
                        estimateDominance(esophagus, index = "absolute")@colData$absolute))

    #Calculates relative in two different ways, should get the same result.
    expect_true(all(estimateDominance(esophagus, index = NULL, relative = TRUE, rank = 1)@colData$relative ==
                        estimateDominance(esophagus, index = "relative")@colData$relative))

    #Calculates DBP in two different ways, should get the same result.
    expect_true(all(estimateDominance(esophagus, index = NULL, relative = TRUE, rank = 1)@colData$DBP ==
                        estimateDominance(esophagus, index = "DBP")@colData$DBP))

    #Calculates DMN in two different ways, should get the same result.
    expect_true(all(estimateDominance(esophagus, index = NULL, relative = TRUE, rank = 2)@colData$DMN ==
                        estimateDominance(esophagus, index = "DMN")@colData$DMN))

})
