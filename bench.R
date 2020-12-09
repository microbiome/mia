library(phyloseq)
library(MicrobiomeExperiment)

d <- data("GlobalPatterns", package = "phyloseq")
GPp <- get(d)
d <- data("GlobalPatterns", package = "MicrobiomeExperiment")
GPm <- get(d)

# accessors
bench::mark(tax_table(GPp),rowData(GPm), check = FALSE, min_iterations = 100)
bench::mark(sample_data(GPp),colData(GPm), check = FALSE, min_iterations = 100)
bench::mark(otu_table(GPp),assay(GPm,"counts"), check = FALSE, min_iterations = 100)
bench::mark(phy_tree(GPp),rowTree(GPm), check = FALSE, min_iterations = 100)

# object creation
data <- assay(GPm,"counts")
rowData <- as.data.frame(rowData(GPm))
colData <- as.data.frame(colData(GPm))
rowTree <- rowTree(GPm)
fun_p1 <- function(){
    phyloseq(sample_data(colData),
             otu_table(data, taxa_are_rows = TRUE),
             tax_table(as.matrix(rowData)),
             phy_tree(rowTree))
}
fun_m1 <- function(){
    MicrobiomeExperiment(assays = list(data),
                         rowData = rowData,
                         colData = colData,
                         rowTree = rowTree)
}
bench::mark(fun_p1(), fun_m1(), check = FALSE, min_iterations = 10)

# subsetting
# rowwise
# this shouldn't work
bench::mark(subset_taxa(GPp, c(TRUE,TRUE,rep(FALSE,24))),
            GPm[1:2,],
            check = FALSE, min_iterations = 10)
# this should
bench::mark(subset_taxa(GPp, c(TRUE,TRUE,rep(FALSE,19214))),
            GPm[1:2,],
            check = FALSE, min_iterations = 10)
# colwise
bench::mark(subset_samples(GPp, c(TRUE,TRUE,rep(FALSE,24))),
            GPm[,1:2],
            check = FALSE, min_iterations = 10)


# cooercion
bench::mark(as.data.frame(sample_data(GPp)),
            as.data.frame(colData(GPm)),
            check = FALSE, min_iterations = 100)

bench::mark(as.data.frame(tax_table(GPp)),
            as.data.frame(rowData(GPm)),
            check = FALSE, min_iterations = 100)

bench::mark(as.data.frame(otu_table(GPp)),
            as.data.frame(assay(GPm,"counts")),
            check = FALSE, min_iterations = 100)

bench::mark(as.matrix(sample_data(GPp)),
            as.matrix(colData(GPm)),
            check = FALSE, min_iterations = 100)

bench::mark(as.matrix(tax_table(GPp)),
            as.matrix(rowData(GPm)),
            check = FALSE, min_iterations = 100)

bench::mark(as.matrix(otu_table(GPp)),
            as.matrix(assay(GPm,"counts")),
            check = FALSE, min_iterations = 100)

################################################################################
## bigger objects
colData <- do.call(rbind,lapply(seq_len(75),function(i){colData}))
data <- do.call(cbind,lapply(seq_len(75),function(i){data}))
rowData <- do.call(rbind,lapply(seq_len(3),function(i){rowData}))
data <- do.call(rbind,lapply(seq_len(3),function(i){data}))
rownames(colData) <- make.unique(rownames(colData))
colData$add <- seq_len(nrow(colData))
colData[["X.SampleID"]] <- make.unique(as.character(colData[["X.SampleID"]]))
rownames(data) <- make.unique(rownames(data))
rownames(rowData) <- rownames(data)

# this object should not be valid since it reports only 26 samples for
# sample_data(p)
p <- phyloseq(sample_data(colData),
              otu_table(data, taxa_are_rows = TRUE),
              tax_table(as.matrix(rowData)))
m <- MicrobiomeExperiment(assays = list(data),
                          rowData = rowData,
                          colData = colData)

fun_p2 <- function(){
    phyloseq(sample_data(colData),
             otu_table(data, taxa_are_rows = TRUE),
             tax_table(as.matrix(rowData)))
}
fun_m2 <- function(){
    MicrobiomeExperiment(assays = list(data),
                         rowData = rowData,
                         colData = colData)
}
fun_p2()
fun_m2()
# creation of big object
bench::mark(fun_p2(), fun_m2(), check = FALSE, min_iterations = 10)

# subsetting bigger objects
# this shouldn't work
bench::mark(subset_taxa(p, c(TRUE,TRUE,rep(FALSE,24))),
            m[1:2,],
            check = FALSE, min_iterations = 10)
# this should
bench::mark(subset_taxa(p, c(TRUE,TRUE,rep(FALSE,57646))),
            m[1:2,],
            check = FALSE, min_iterations = 10)
# colwise
bench::mark(subset_samples(p, c(TRUE,TRUE,rep(FALSE,1948))),
            m[,1:2],
            check = FALSE, min_iterations = 10)
