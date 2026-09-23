
test_that("applyByModule", {
    
    set.seed(123)
    tse <- makeTSE()
    
    assayNames(tse) <- "counts"
    
    expect_error(
        applyByModule(
            tse,
            "rows",
            "var1",
            getAlpha,
            index = c("observed", "shannon")
        )
    )
    
    expect_error(
        applyByModule(
            tse,
            "rows",
            "var1",
            addAlpha,
            index = c("observed", "shannon"),
            min.group.size = 2
        ),
        "No 'group' has size above 'min.group.size'.",
        fixed = TRUE
    )
    
    rowData(tse)$var1 <- rep(letters[1:2], 5)
    
    tse <- applyByModule(
        tse,
        "rows",
        "var1",
        getAlpha,
        index = c("observed", "shannon"),
        min.group.size = 2
    )
    
    namess <- interaction(c("observed", "shannon"), c("a", "b"))
    expect_contains(names(colData(tse)), namess)
    
    tse <- applyByModule(
        tse,
        "rows",
        "var1",
        addAlpha,
        index = c("observed", "shannon"),
        min.group.size = 2
    )
    
    # expect_in(altExpNames)
    
})