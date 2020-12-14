#' Estimate Dominance
#'
#' This function calculates the community dominance index.
#'
#' This includes the \sQuote{DBP}, \sQuote{DMN}, \sQuote{Absolute},
#' \sQuote{Relative}, \sQuote{Core Abundance}, \sQuote{Gini}
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param abund_values A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to use for prevalence calculation.
#'
#' @param index Specifies the indices which are calculated.
#'
#' @param rank,... additional arguments
#' \itemize{
#'   \item{If \code{!is.null(rank)} arguments are passed on to
#'   \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'   \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#' }
#'
#' @param as_relative logical scalar: Should the detection threshold be applied
#'   on compositional (relative) abundances? (default: \code{TRUE})
#'
#' @param aggregate (Optional, default = TRUE) Aggregate the top members or not.
#' If aggregate=TRUE, then the sum of relative abundances is returned.
#' Otherwise the relative abundance is returned for the single taxa with
#' the indicated rank.
#'
#' @param name A name for the column of the colData where the calculated
#' Dominance indices should be stored in.
#'
#' @details
#' \code{estimateDominance} calculates the community dominance indices.
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[=estimateDiversity]{estimateDiversity}}}
#' }
#'
#' @name estimateDominance
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus, package = "MicrobiomeExperiment")
#' esophagus <- as(esophagus, "MicrobiomeExperiment")
#'
#' #Calculates simpson dominance index
#' esophagus <- estimateDominance(esophagus, index="simpson")
#' #Shows all indices
#' colData(esophagus)
#'
#' #Tries to calculate index that is not accepted. Gets an error.
#' esophagus <- estimateDominance(esophagus, index="does_not_exist")
#' #Tries to calculate index that is not accepted
#' # and index that is accepted. Gets an error.
#' esophagus <- estimateDominance(esophagus, index=c("DMN", "does_not_exist"))
#' #Shows all indices
#' colData(esophagus)
#'
#' #Indices must be written correctly (e.g. DBP, not dbp), unless gets an error.
#' esophagus <- estimateDominance(esophagus, index="dbp")
#' #Calculates DBP and Core Abundance indices
#' esophagus <- estimateDominance(esophagus, index=c("DBP", "core_abundance"))
#' #Shows all indices
#' colData(esophagus)
#' #Shows DBP index
#' colData(esophagus)$DBP
#'
#' #Deletes DBP index
#' colData(esophagus)$DBP <- NULL
#' #Shows all indices, DBP is deleted
#' colData(esophagus)
#' #Deletes all indices
#' colData(esophagus) <- NULL
#'
#' #Names of columns can be chosen, but the length of arguments must match.
#' esophagus <- estimateDominance(esophagus,
#'     index=c("DBP", "core_abundance"),
#'     name = c("index1", "index2"))
#' #Shows all indices
#' colData(esophagus)
#' #If they do not match, gets an error.
#' esophagus <- estimateDominance(esophagus,
#'     index="simpson",
#'     name = c("index3", "index4"))
#' #Shows all indices
#' colData(esophagus)
#' #Deletes all indices
#' colData(esophagus) <- NULL
#'
#' #Calculates all indices
#' esophagus <- estimateDominance(esophagus)
#' #Shows all indices
#' colData(esophagus)
#'
NULL

#' @rdname estimateDominance
#' @export
setGeneric("estimateDominance",signature = c("x"),
           function(x, abund_values = "counts", index = c("DBP", "DMN", "absolute", "relative", "simpson", "core_abundance", "gini"), rank=1, as_relative=TRUE, aggregate=TRUE,
                    name = index, ...)
               standardGeneric("estimateDominance"))


#' @rdname estimateDominance
#' @export
setMethod("estimateDominance", signature = c(x = "MicrobiomeExperiment"),
          function(x, abund_values = "counts", index = c("DBP", "DMN", "absolute", "relative", "simpson", "core_abundance", "gini"), rank=1, as_relative=TRUE, aggregate=TRUE,
                   name = index, ...){

              # input check
              index <- match.arg(index, several.ok = TRUE)

              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }

              #Initialize table that is used to store indices
              tab <- NULL

              #If the indices vector is not null
              if (!is.null(index)) {
                  #loop through index list
                  for (idx in index) {
                      #Get the specific index, and save it to table
                      tab <- cbind(tab,
                          .dominance_help(x,
                              abund_values=abund_values, index=idx, rank=rank,
                              as_relative=as_relative,
                              aggregate=aggregate))
                  }

                  # Add indices names to table's columns' names
                  colnames(tab) <- index

                  # Convert table to data frame
                  tab <- as.data.frame(tab)

                  # Adds index data to original object
                  if (length(name) == 1 && name == "all") {
                      name <- colnames(tab)
                  }

                  #colData(x) <- cbind(colData(x), tab)
                  x <- .add_indices_to_coldata(x, tab, name)
              }

              # Return object
              return(x)


          }
)




#---------------------------Help functions----------------------------------------------------------------

.dominance_help <- function(x, abund_values = "counts", index="all", rank=1, as_relative=TRUE,
                           aggregate=TRUE) {
    #Stores the absolute abundances to "otu" variable
    otu <- assay(x,abund_values)

    #if index does not have any values
    if (is.null(index)) {
        rank <- rank
    } else if (index == "absolute") {
        #Rank=1 by default but can be tuned
        as_relative <- FALSE
    } else if (index %in% c("relative")) {
        #Rank=1 by default but can be tuned
        as_relative <- TRUE
    } else if (index %in% c("dbp")) {
        #Berger-Parker
        rank <- 1
        as_relative <- TRUE
    } else if (index %in% c("dmn")) {
        #McNaughton's dominance
        rank <- 2
        as_relative <- TRUE
        aggregate <- TRUE
        #If index is "Simpson", calculates the Simpson to all the samples
    } else if (index %in% c("simpson")) {
        tmp <- apply(otu, 2, function(x) {
            .simpson_dominance(x)})
        return(tmp)
    } else if (index %in% c("core_abundance")) {
        prevalence <- getPrevalentAbundance(x, detection=0, as_relative=TRUE)
        return(prevalence)
        #If index is "Gini" calculates the gini index to all the samples
    } else if (index == "gini") {

        tmp <- apply(otu, 2, function(x) {
            .get_gini(x)})

        return(tmp)
    }

    if (rank == 1 && as_relative) {
        index <- "dbp"
    } else if (rank == 2 && as_relative && aggregate) {
        index <- "dmn"
    }

    if (as_relative) {
        #Calculates the relative abundance to all the columns
        otu <- apply(otu, 2, function(x) {
            x/sum(x, na.rm=TRUE)
        })
    }

    #Aggregate or not
    if (!aggregate) {
        do <- apply(otu, 2, function(x) {
            sort(x, decreasing = TRUE)[[rank]]
        })
    } else {
        do <- apply(otu, 2, function(x) {
            sum(sort(x, decreasing = TRUE)[seq_len(rank)])
        })
    }

    #Add otus' names to the table
    names(do) <- colnames(otu)

    #Creates a matrix that has one column. Values of index is saved to that column.
    if (is.vector(do)) {
        do <- as.matrix(do, ncol=1)
        colnames(do) <- index
    }

    #Returns a matrix that has one column. Rows are otus, column includes values of index.
    return(do)

}

# x: Species count vector
.simpson_dominance <- function(x, zeroes=TRUE) {

    if (!zeroes) {
        x[x > 0]
    }

    # Relative abundances
    p <- x/sum(x)

    # Simpson index (has interpretation as dominance)
    lambda <- sum(p^2)

    # More advanced Simpson dominance (Simpson 1949) However let us not use
    # this as it is not in [0,1] and it is very highly correlated with the
    # simpler variant lambda Species richness (number of species)
    # S <- length(x) sum(p * (p - 1)) / (S * (S - 1))

    lambda

}

.get_gini <- function(x) {

    # Gini index for each sample
    gini_for_sample <- .calculate_gini(x)

    return(gini_for_sample)

}



.calculate_gini <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

#Gets ME object and index data frame as parameters
.add_indices_to_coldata <- function(x, indicesDF, names){
    #Stores sample names
    sampleNames <- rownames(indicesDF)

    #Checks if the length of names and indices match. If not, gives an error.
    if(length(names) != length(colnames(indicesDF))){
        stop("Number of names and indices do not match.")

    }

    #Loops through columns of data frame. One by one adds index data
    #to the ME object.
    order_of_name = 1
    for(indexName in colnames(indicesDF)){
        #Stores coldata to colData variable
        colData <- colData(x)

        #Transpose and converts index data to numeric values
        indexData <- as.numeric(t(indicesDF[[indexName]]))

        #Adds sample names to index values
        names(indexData) <- sampleNames

        #Add index data to the colData variable. The name is from "names" variable
        colData[[names[order_of_name]]] <- indexData
        order_of_name = order_of_name + 1

        #Adds the new colData to the ME object
        colData(x) <- colData
    }

    return(x)

}
