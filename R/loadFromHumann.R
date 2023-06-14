

file <- "data/pathcoverage.tsv"
pathabundance <- read.delim(path, check.names = FALSE)


loadFromHumann(file)

loadFromHumann <- function(file, colData = NULL, ...){
    ################################ Input check ################################
    if(!mia:::.is_non_empty_string(file)){
        stop("'file' must be a single character value.",
             call. = FALSE)
    }
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if(!is.null(colData) && !mia:::.is_non_empty_string(colData)){
        stop("'colData' must be a single character value or NULL.",
             call. = FALSE)
    }
    ############################## Input check end #############################
    # Read metaphlan data
    data <- .read_humann(file)
    # Create TreeSE from the data
    tse <- .create_tse_from_humann(data, colData = colData, ...)
    return(tse)
}

################################ HELP FUNCTIONS ################################

# Read Humann file, catch error if it occurs
.read_humann <- function(file){
    # Read the table. Catch error and give more informative message
    table <- tryCatch(
        {
            read.delim(file, check.names = FALSE)
        },
        error = function(condition){
            stop("Error while reading ", file,
                 "\nPlease check that the file is in merged Metaphlan file format.",
                 call. = FALSE)
        }
    )
    # IN the first column name, there is "# " prefix. Remove it
    colnames(table)[1] <- gsub("# ", "", colnames(table)[1])
    # Add rownames
    rownames(table) <- table[, 1] 
    # Check that file is in right format
    if( .check_humann(table) ){
        stop("Error while reading ", file,
             "\nPlease check that the file is in merged humann file format.",
             call. = FALSE)
    }
    return(table)
}

# Check that metaphlan file contains correct information
.check_humann <- function(data){
    # Get rowdata columns
    rowdata_col <- c("Pathway", "Gene Family")
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata_columns <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay_columns <- data[ , -rowdata_id, drop = FALSE]
    # Initialize result 
    result <- TRUE
    
    # Check rowdata column names that they contain right information, and check that 
    # rest of the columns represents abundances in samples.
    # If these requirements are met, give FALSE. Otherwise, give TRUE.
    if( any(colnames(rowdata_columns) %in% c("Pathway", "Gene Family")) && 
        is.numeric(unlist(assay_columns)) ){
        result <- FALSE
    }
    return(result)
}


.create_tse_from_humann <- function(data, colData, assay.type = "counts", ...){
    # Get rowdata columns
    rowdata_col <- c("Pathway", "Gene Family")
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay <- data[ , -rowdata_id, drop = FALSE]
    
    # Parse rowdata. The data includes gene/pathway info along with taxonomy
    # info. Separate the information to own columns.
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(rowdata)))
    # Get the column and split gene/pathway and taxonomy info
    rowdata_temp <- strsplit(rowdata[[1]], "\\|",)
    # Are some rows missing taxonomy info? Get their indices.
    missing_taxa <- lengths(rowdata_temp) == 1
    # Create a df from the list of splitted data.
    rowdata_temp <- data.frame(t(data.frame(rowdata_temp)))
    # Replace missing taxonomy info with NA
    rowdata_temp[missing_taxa, 2] <- NA
    # Replace column names
    colnames(rowdata_temp) <- c(colnames(rowdata), "Taxonomy")
    
    # Now we have rowdata that includes gene/pathway info in one column and
    # taxonomy info in other. Let's parse the taxonomy info so that species
    # genus etc levels are in unique columns.
    taxonomy <- mia:::.parse_taxonomy(
        rowdata_temp, column_name = "Taxonomy", sep = "\\.", ...)
    
    # Convert all data to DataFrame
    rowdata <- DataFrame(rowdata)
    rowdata_temp <- DataFrame(rowdata_temp)
    
    # Create rowdata from the information that is parsed
    rownames(rowdata) <- rowdata[ , 1]
    colnames(rowdata) <- paste0(colnames(rowdata), "_long")
    rowdata[["Pathway"]] <- rowdata_temp[["Pathway"]]
    rowdata <- cbind(rowdata, taxonomy)
    
    # Create assays
    assay <- as.matrix(assay)
    assays <- SimpleList(counts = assay)
    names(assays) <- assay.type
    
    # If coldata was provided, get it
    if( is.character(colData) ){
        colData <- read.table(file = colData, header = TRUE, sep = "\t")
        colData(tse) <- colData
    } else if( is.null(colData) ) {
        # If NULL, initialize DF
        colData <- DataFrame(row.names = colnames(assay))
    }
    # Ensure that the class is DF
    colData <- DataFrame(colData)
    
    # Create TreeSE
    tse <- TreeSummarizedExperiment(
        assays = assays, rowData = rowdata, colData = colData)
    return(tse)
}
