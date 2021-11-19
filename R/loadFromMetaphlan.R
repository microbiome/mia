

file <- "~/Downloads/merged_abundance_table.txt"

readLines(file)

table <- read.table(file, header = TRUE)

levels <- sapply(table[["clade_name"]], FUN = .get_level)

.get_level <- function(string){
    levels <- gregexpr("([kpcofgs]+)__", string)[[1]]
    lowest_level_ind <- levels[length(levels)]
    lowest_level <- substr(string, start = lowest_level_ind, stop = lowest_level_ind)
    return(lowest_level)
}

order <- c("s", "g", "f", "o", "c", "p", "k")

lowest_level_found <- levels[order(match(levels, order))][1]
lowest_level_found 

table <- table[grepl(paste0(lowest_level_found, "__"), table[["clade_name"]]), ]
table
rowdata <- table[, c("clade_name", "NCBI_tax_id")]

assay_columns <- colnames(table)[!colnames(table) %in% c("clade_name", "NCBI_tax_id")]

assay <- table[, assay_columns]

taxonomy <- rowdata[, "clade_name", drop = FALSE]
taxonomy
rowdata[, "clade_name"] <- NULL
rowdata
colnames(taxonomy)[colnames(taxonomy) == "clade_name"] <- "Taxon"

taxonomy <- .parse_q2taxonomy(taxonomy, sep = "\\|", removeTaxaPrefixes = FALSE)
rowdata <- cbind(taxonomy, rowdata)

assay

se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay), 
                                           rowData = rowdata)
rownames(se) <- .get_taxonomic_label(se)

se
