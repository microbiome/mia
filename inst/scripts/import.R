# Load mia R package
library(mia)

# Set file paths for input data
biom_file_path <- "PATH_TO_BIOM_FILE"
sample_meta_file_path <- "PATH_TO_SAMPLE_METADATA_FILE"
tree_file_path <- "PATH_TO_PHYLOGENETIC_TREE_FILE"

# Load data from the biom file into a SummarizedExperiment container
x1 <- biomformat::read_biom(biom_file_path)
se <- makeTreeSEFromBiom(x1)

# Convert the SummarizedExperiment to a TreeSummarizedExperiment container
tse <- as(se, "TreeSummarizedExperiment")

# We notice that the rowData fields do not have descriptibve names.
# Hence, let us rename the columns in rowData
names(rowData(tse)) <- c("Kingdom", "Phylum", "Class", "Order", 
                         "Family", "Genus")


# We also notice that the taxa names are of form "c__Bacteroidia" etc.
# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(tse), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so convert this back to DataFrame format. 
# and assign the cleaned data back to the TSE rowData
rowData(tse) <- DataFrame(rowdata_modified)

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
colData(tse) <- DataFrame(sample_meta)

# Read the phylogenetic tree from file and assign it to the rowTree slot of the TSE object
tree <- ape::read.tree(tree_file_path)
rowTree(tse) <- tree

# Save the final TSE object as an R data file
saveRDS(tse, file = "OUTPUT_FILE_NAME.rds")
