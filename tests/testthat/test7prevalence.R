context("prevalence")
test_that("prevalence", {

    data(GlobalPatterns)
    
    # Output is always frequencies between 0 to 1	
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, sort=TRUE, as_relative=TRUE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)

    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, sort=TRUE, as_relative=FALSE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)    

    # If we look at 
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include_lowest=TRUE, sort=TRUE, as_relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0/100, include_lowest=FALSE, sort=TRUE, as_relative=TRUE)        
    expect_true(all(pr1 == pr2))

})
