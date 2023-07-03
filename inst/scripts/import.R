# Load mia R package
library(mia)

# Set file paths for input data
biom_file_path <- "PATH_TO_BIOM_FILE"
sample_meta_file_path <- "PATH_TO_SAMPLE_METADATA_FILE"
tree_file_path <- "PATH_TO_PHYLOGENETIC_TREE_FILE"

# Load data from the biom file into a SummarizedExperiment container
se <- loadFromBiom(biom_file_path, removeTaxaPrefixes = TRUE, rankFromPrefix = TRUE)

# Read sample metadata from file and add column names if necessary
sample_meta <-
    read.table(
        sample_meta_file_path,
        sep = ",",
        header = FALSE,
        row.names = 1
    )
# Add headers for the columns (if they seem to be missing)
colnames(sample_meta) <- c("patient_status", "cohort",
                           "patient_status_vs_cohort", "sample_name")

# Add sample metadata to colData slot of the TSE object
# Note that the data must be given in a DataFrame format (required for our purposes)
colData(tse) <- cbind(colData(tse), DataFrame(sample_meta))

# Read the phylogenetic tree from file and assign it to the rowTree slot of the TSE object
tree <- ape::read.tree(tree_file_path)
rowTree(se) <- tree

# Save the final TSE object as an R data file
saveRDS(se, file = "OUTPUT_FILE_NAME.rds")
