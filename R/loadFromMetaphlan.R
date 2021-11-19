

file <- "~/Downloads/merged_abundance_table.txt"

readLines(file)

# Read the table
table <- read.table(file, header = TRUE)

# Get the lowest level of each row
levels <- sapply(table[["clade_name"]], FUN = .get_level)

.get_level <- function(string){
    levels <- gregexpr("([kpcofgs]+)__", string)[[1]]
    lowest_level_ind <- levels[length(levels)]
    lowest_level <- substr(string, start = lowest_level_ind, stop = lowest_level_ind)
    return(lowest_level)
}

order <- c("s", "g", "f", "o", "c", "p", "k")
# Order the data and get the lowest level of data
lowest_level_found <- levels[order(match(levels, order))][1]
lowest_level_found 

# Get those rows that include information at lowest level
table <- table[grepl(paste0(lowest_level_found, "__"), table[["clade_name"]]), ]
table

# Get the data that belongs to rowData
rowdata_columns <- c("clade_name", "NCBI_tax_id")
rowdata <- table[, colnames(table) %in% rowdata_columns, drop = FALSE]

# Get those columns that belong to assay
assay_columns <- colnames(table)[!colnames(table) %in% rowdata_columns]
assay <- table[, assay_columns, drop = FALSE]

# Store taxonomic ids to add them later
tax_id <- rowdata$NCBI_tax_id

# Psrse taxonomic levels
rowdata <- .parse_taxonomy(rowdata, sep = "\\|", column_name = "clade_name")

# Add taxonomic ids
rowdata$NCBI_tax_id <- tax_id

assay

# Create SE
se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay), 
                                           rowData = rowdata)
# Add rownames from the lowest taxonomic level
rownames(se) <- .get_taxonomic_label(se)

se
rowData(se)
