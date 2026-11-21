# Fix the sample metadata of the Tengeler2020 dataset (mia issue #904).
#
# The original import (inst/scripts/import.R, 2023) attached the sample
# metadata with cbind() in the row order of the metadata file, not by sample
# name. The counts kept the column order of the BIOM file, so the sample
# names and colData (patient_status, cohort, ...) were misaligned with the
# counts. This script matches the metadata to the counts by sample name and
# keeps everything else (features, rowData, tree) unchanged.
#
# Source files (from the OMA book repository, inst/extdata):
#   Aggregated_humanization2.biom          counts and taxonomy
#   Mapping_file_ADHD_aggregated.csv       sample metadata
#   Data_humanization_phylo_aggregation.tre phylogenetic tree

library(mia)

extdata <- "PATH_TO/OMA/inst/extdata"

data("Tengeler2020", package = "mia")
tse <- Tengeler2020

# Sample order of the counts, from the BIOM file
biom <- importBIOM(
    file.path(extdata, "Aggregated_humanization2.biom"),
    prefix.rm = TRUE, rank.from.prefix = TRUE, artifact.rm = TRUE)
# The counts of the packaged dataset are the BIOM counts in BIOM order
stopifnot(identical(unname(assay(tse, "counts")), unname(assay(biom))))

# Sample metadata, matched by sample name
meta <- read.csv(
    file.path(extdata, "Mapping_file_ADHD_aggregated.csv"),
    comment.char = "", check.names = FALSE)
meta <- meta[match(colnames(biom), meta[["#SampleID"]]), ]
stopifnot(identical(meta[["#SampleID"]], colnames(biom)))

colnames(tse) <- colnames(biom)
colData(tse) <- DataFrame(
    patient_status = meta[["Treatment"]],
    cohort = meta[["Cohort"]],
    patient_status_vs_cohort = meta[["TreatmentxCohort"]],
    sample_name = meta[["Description"]],
    row.names = colnames(biom))

Tengeler2020 <- tse
save(Tengeler2020, file = "data/Tengeler2020.rda", compress = "xz")
