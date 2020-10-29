context("merge")
test_that("merge", {
    # .check_f
    expect_error(SEtup:::.norm_f(),
                 'argument "f" is missing')
    expect_error(SEtup:::.norm_f(6),
                 'argument "f" is missing')
    expect_error(SEtup:::.norm_f(6,5),
                 "'f' must be a factor or character vector")
    f <- factor(c(rep("a",3),rep("b",3)))
    expect_true(is.factor(SEtup:::.norm_f(6,f)))
    # .check_archetype
    expect_error(SEtup:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(SEtup:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_error(SEtup:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_equal(SEtup:::.norm_archetype(f, 1),c(1,1))
    expect_equal(SEtup:::.norm_archetype(f, c(1,2)),c(1,2))
    expect_error(SEtup:::.norm_archetype(f, c(1,2,3)),
                 "length of 'archetype' must have the same length as levels")
    expect_error(SEtup:::.norm_archetype(f, c(5)),
                 "'archetype' out of bounds for some levels of 'f'")
    # .norm_archetype
    expect_error(SEtup:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(SEtup:::.norm_archetype(f),
                 'argument "archetype" is missing')
    actual <- SEtup:::.norm_archetype(f, c(1,2))
    expect_equal(actual, c(1,2))
    actual <- SEtup:::.norm_archetype(f, c(1))
    expect_equal(actual, c(1,1))
    # .get_element_pos
    expect_error(SEtup:::.get_element_pos(),
                 'argument "archetype" is missing')
    expect_error(SEtup:::.get_element_pos(f),
                 'argument "archetype" is missing')
    actual <- SEtup:::.get_element_pos(f, archetype = SEtup:::.norm_archetype(f, 1))
    expect_equal(actual,c(a = 1, b = 4))
    actual <- SEtup:::.get_element_pos(f, archetype = SEtup:::.norm_archetype(f, 2))
    expect_equal(actual,c(a = 2, b = 5))
    actual <- SEtup:::.get_element_pos(f, archetype = c(2,1))
    expect_equal(actual,c(a = 2, b = 4))
    # .merge_rows
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    expect_error(SEtup:::.merge_rows(),
                 'argument "f" is missing')
    x <- SummarizedExperiment(assays = list(mat = mat))
    xr <- SummarizedExperiment(assays = list(mat = mat),
                               rowRanges = gr)
    xrl <- SummarizedExperiment(assays = list(mat = mat),
                                rowRanges = unname(grl))
    expect_error(SEtup:::.merge_rows(x),
                 'argument "f" is missing')
    FUN_check_x <- function(x,archetype=1){
        actual <- mergeRows(x, f, archetype)
        expect_s4_class(actual,class(x))
        expect_equal(dim(actual),c(2,10))
    }
    lapply(list(x,xr,xrl),FUN_check_x)
    lapply(list(x,xr,xrl),FUN_check_x,archetype=2)
    #
    f <- factor(c(rep("a",3),rep("b",3)))
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    me <- MicrobiomeExperiment(assays = list(mat = mat),
                               rowRanges = unname(grl))
    FUN_check_x <- function(x,archetype=1){
      actual <- mergeRows(x, f, archetype)
      expect_s4_class(actual,class(x))
      expect_equal(dim(actual),c(2,10))
    }
    lapply(list(xtse,me),FUN_check_x)
    lapply(list(xtse,me),FUN_check_x,archetype=2)
})
