# Load necessary packages
library(mia)
library(stringr)
library(S4Vectors)
 
### Defining FILE PATHS ###
path_to_assay <- "Tito2024QMP_assay.csv"
path_to_coldata <- "Tito2024QMP_coldata.csv"
path_to_rowdata <- "Tito2024QMP_rowdata.csv"
 
### Importing FILES ###
Tito2024QMP_assay <- read.csv(path_to_assay, row.names = 1, check.names = FALSE)
colnames(Tito2024QMP_assay) <- gsub("\\.", "_", colnames(Tito2024QMP_assay))
 
Tito2024QMP_coldata <- read.csv(path_to_coldata, row.names = 1)
Tito2024QMP_rowdata <- read.csv(path_to_rowdata, row.names = 1)

 
# Replace periods with underscores in the Species column values
if ("species" %in% colnames(Tito2024QMP_rowdata)) {
  Tito2024QMP_rowdata$species <- str_replace_all(Tito2024QMP_rowdata$species, "\\.", "_")
}

# Convert empty strings to NA in coldata and rowdata
Tito2024QMP_coldata[Tito2024QMP_coldata == ""] <- NA
Tito2024QMP_rowdata[Tito2024QMP_rowdata == ""] <- NA

# Add row names as a new column in CRC_coldata
Tito2024QMP_coldata$sampleID <- as.character(rownames(Tito2024QMP_coldata))
            
# Convert CRC_coldata and CRC_rowdata to DataFrame objects
Tito2024QMP_coldata <- DataFrame(Tito2024QMP_coldata)
Tito2024QMP_rowdata <- DataFrame(Tito2024QMP_rowdata)

# Convert colData fields to factors except for sampleID
Tito2024QMP_coldata[] <- lapply(seq_along(Tito2024QMP_coldata), function(i) {
    if (is.character(Tito2024QMP_coldata[[i]]) && colnames(Tito2024QMP_coldata)[i] != "sampleID") {
       return(as.factor(Tito2024QMP_coldata[[i]]))
       } else {
          return(Tito2024QMP_coldata[[i]])
       }
})
                
# Reassign the list back to the DataFrame
Tito2024QMP_coldata <- DataFrame(Tito2024QMP_coldata)
                  
### Constructing TSE ###
Tito2024QMP <- TreeSummarizedExperiment(
    assays = SimpleList(counts = as.matrix(Tito2024QMP_assay)),
    colData = Tito2024QMP_coldata,
    rowData = Tito2024QMP_rowdata
)
                    
### SAVE into RDA ###
save(Tito2024QMP, file = "Tito2024QMP.rda")