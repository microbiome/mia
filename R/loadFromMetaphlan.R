

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
# 
# taxonomy <- .parse_taxonomy(taxonomy, sep = "\\|", column_name = "clade_name", removeTaxaPrefixes = FALSE)
# 
# # .parse_taxonomy <- function(taxa_table, column_name = "Taxon", sep = "\\|",
#                             removeTaxaPrefixes = FALSE, ...){
#     
#     #  work with any combination of taxonomic ranks available
#     all_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#     all_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
#     
#     # split the taxa strings
#     taxa_split <- CharacterList(strsplit(taxa_tab[,"Taxon"],sep))
#     # extract present prefixes
#     taxa_prefixes <- lapply(taxa_split, substr, 1L, 3L)
#     # match them to the order given by present_prefixes
#     taxa_prefixes_match <- lapply(taxa_prefixes, match, x = all_prefixes)
#     taxa_prefixes_match <- IntegerList(taxa_prefixes_match)
#     # get the taxa values
#     if(removeTaxaPrefixes){
#         taxa_split <- lapply(taxa_split,
#                              gsub,
#                              pattern = "([kpcofgs]+)__",
#                              replacement = "")
#         taxa_split <- CharacterList(taxa_split)
#     }
#     # extract by order matches
#     taxa_split <- taxa_split[taxa_prefixes_match]
#     #
#     if(length(unique(lengths(taxa_split))) != 1L){
#         stop("Internal error. Something went wrong")
#     }
#     taxa_tab <- DataFrame(as.matrix(taxa_split))
#     colnames(taxa_tab) <- all_ranks
# }


rowdata <- cbind(taxonomy, rowdata)

assay

se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay), 
                                           rowData = rowdata)
#rownames(se) <- .get_taxonomic_label(se)

se
