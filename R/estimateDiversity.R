#' Estimate (alpha) diversity measures
#'
#' Several functions for calculating (alpha) diversity indices, including 
#' the \code{vegan} package options and some others.
#'
#' The available indices include the \sQuote{Coverage}, 
#' \sQuote{Faith's phylogenetic diversity}, \sQuote{Fisher alpha},
#' \sQuote{Gini-Simpson}, 
#' \sQuote{Inverse Simpson}, \sQuote{log-modulo skewness}, and \sQuote{Shannon} 
#' indices. See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object or \code{\link{TreeSummarizedExperiment}}.
#' The latter is recommended for microbiome data sets and tree-based alpha diversity indices.
#' 
#' @param tree A phylogenetic tree that is used to calculate 'faith' index.
#'   If \code{x} is a \code{TreeSummarizedExperiment}, \code{rowTree(x)} is 
#'   used by default.
#'
#' @param assay_name the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param abund_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay_name} instead. At some point \code{abund_values}
#'   will be disabled.)
#'
#' @param index a \code{character} vector, specifying the diversity measures
#'   to be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices.
#'   
#' @param tree_name a single \code{character} value for specifying which
#'   rowTree will be used to calculate faith index. 
#'   (By default: \code{tree_name = "phylo"})
#'   
#' @param node_lab NULL or a character vector specifying the links between rows and 
#'   node labels of \code{tree}. If a certain row is not linked with the tree, missing 
#'   instance should be noted as NA. When NULL, all the rownames should be found from
#'   the tree. (By default: \code{node_lab = NULL})
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{threshold}{ A numeric value in the unit interval,
#'   determining the threshold for coverage index. By default,
#'   \code{threshold} is 0.9.}
#'   \item{quantile}{ Arithmetic abundance classes are evenly cut up to to
#'   this quantile of the data. The assumption is that abundances higher than
#'   this are not common, and they are classified in their own group.
#'   By default, \code{quantile} is 0.5.}
#'   \item{num_of_classes}{ The number of arithmetic abundance classes
#'   from zero to the quantile cutoff indicated by \code{quantile}. 
#'   By default, \code{num_of_classes} is 50.}
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#'
#' @details
#'
#' Alpha diversity is a joint quantity that combines elements or community richness
#' and evenness. Diversity increases, in general, when species richness or
#' evenness increase.
#'
#' By default, this function returns all indices.
#'
#' \itemize{
#' 
#' \item{'coverage' }{Number of species needed to cover a given fraction of
#' the ecosystem (50 percent by default). Tune this with the threshold
#' argument.}
#' 
#' \item{'faith' }{Faith's phylogenetic alpha diversity index measures how
#' long the taxonomic distance is between taxa that are present in the sample.
#' Larger values represent higher diversity. Using this index requires
#' rowTree. (Faith 1992)}
#' 
#' \item{'fisher' }{Fisher's alpha; as implemented in
#' \code{\link[vegan:diversity]{vegan::fisher.alpha}}. (Fisher et al. 1943)}
#' 
#' \item{'gini_simpson' }{Gini-Simpson diversity i.e. \eqn{1 - lambda},
#' where \eqn{lambda} is the
#' Simpson index, calculated as the sum of squared relative abundances.
#' This corresponds to the diversity index
#' 'simpson' in \code{\link[vegan:diversity]{vegan::diversity}}.
#' This is also called Gibbs–Martin, or Blau index in sociology,
#' psychology and management studies. The Gini-Simpson index (1-lambda)
#' should not be
#' confused with Simpson's dominance (lambda), Gini index, or
#' inverse Simpson index (1/lambda).}
#' 
#' \item{'inverse_simpson' }{Inverse Simpson diversity:
#' \eqn{1/lambda} where \eqn{lambda=sum(p^2)} and p refers to relative
#' abundances.
#' This corresponds to the diversity index
#' 'invsimpson' in vegan::diversity. Don't confuse this with the
#' closely related Gini-Simpson index}
#'
#' \item{'log_modulo_skewness' }{The rarity index characterizes the
#' concentration of species at low abundance. Here, we use the skewness of
#' the frequency 
#' distribution of arithmetic abundance classes (see Magurran & McGill 2011).
#' These are typically right-skewed; to avoid taking log of occasional
#' negative skews, we follow Locey & Lennon (2016) and use the log-modulo
#' transformation that adds a value of one to each measure of skewness to
#' allow logarithmization.}
#'
#' \item{'shannon' }{Shannon diversity (entropy).}
#' 
#' }
#'
#' @references
#'
#' Beisel J-N. et al. (2003)
#' A Comparative Analysis of Diversity Index Sensitivity.
#' _Internal Rev. Hydrobiol._ 88(1):3-15.
#' \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#'
#' Bulla L. (1994)
#' An  index of diversity and its associated diversity measure.
#' _Oikos_ 70:167--171
#'
#' Faith D.P. (1992)
#' Conservation evaluation and phylogenetic diversity.
#' _Biological Conservation_ 61(1):1-10.
#'
#' Fisher R.A., Corbet, A.S. & Williams, C.B. (1943)
#' The relation between the number of species and the number of individuals in
#' a random sample of animal population.
#' _Journal of Animal Ecology_ *12*, 42-58.
#' 
#' Locey K.J. & Lennon J.T. (2016)
#' Scaling laws predict global microbial diversity.
#' _PNAS_ 113(21):5970-5975.
#'
#' Magurran A.E., McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment.
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Smith B. & Wilson JB. (1996)
#' A Consumer's Guide to Diversity Indices.
#' _Oikos_ 76(1):70-82.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateEvenness]{estimateEvenness}}}
#'   \item{\code{\link[mia:estimateDominance]{estimateDominance}}}
#'   \item{\code{\link[vegan:diversity]{diversity}}}
#'   \item{\code{\link[vegan:specpool]{estimateR}}}
#' }
#'
#' @name estimateDiversity
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # All index names as known by the function
#' index <- c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher", 
#' "faith",  "log_modulo_skewness")
#'
#' # Corresponding polished names
#' name <- c("Shannon","GiniSimpson","InverseSimpson", "Coverage", "Fisher", 
#' "Faith",  "LogModSkewness")
#'
#' # Calculate diversities
#' tse <- estimateDiversity(tse, index = index)
#'
#' # The colData contains the indices with their code names by default
#' colData(tse)[, index]
#'
#' # Removing indices
#' colData(tse)[, index] <- NULL
#' 
#' # 'threshold' can be used to determine threshold for 'coverage' index
#' tse <- estimateDiversity(tse, index = "coverage", threshold = 0.75)
#' # 'quantile' and 'num_of_classes' can be used when
#' # 'log_modulo_skewness' is calculated
#' tse <- estimateDiversity(tse, index = "log_modulo_skewness",
#'        quantile = 0.75, num_of_classes = 100)
#'
#' # It is recommended to specify also the final names used in the output.
#' tse <- estimateDiversity(tse,
#'   index = c("shannon", "gini_simpson", "inverse_simpson", "coverage",
#'                "fisher", "faith", "log_modulo_skewness"),
#'    name = c("Shannon", "GiniSimpson",  "InverseSimpson",  "Coverage",
#'                "Fisher", "Faith", "LogModSkewness"))
#'
#' # The colData contains the indices by their new names provided by the user
#' colData(tse)[, name]
#'
#' # Compare the indices visually
#' pairs(colData(tse)[, name])
#'
#' # Plotting the diversities - use the selected names
#' library(scater)
#' plotColData(tse, "Shannon")
#' # ... by sample type
#' plotColData(tse, "Shannon", "SampleType")
#' \dontrun{
#' # combining different plots
#' library(patchwork)
#' plot_index <- c("Shannon","GiniSimpson")
#' plots <- lapply(plot_index,
#'                plotColData,
#'                object = tse,
#'                x = "SampleType",
#'                colour_by = "SampleType")
#' plots <- lapply(plots,"+",
#'    theme(axis.text.x = element_text(angle=45,hjust=1)))
#' names(plots) <- plot_index
#' plots$Shannon + plots$GiniSimpson + plot_layout(guides = "collect")
#' }
NULL

#' @rdname estimateDiversity
#' @export
setGeneric("estimateDiversity",signature = c("x"),
        function(x, assay_name = abund_values, abund_values = "counts",
                index = c("coverage", "fisher", "gini_simpson", 
                        "inverse_simpson", "log_modulo_skewness", "shannon"),
                    name = index, ...)
                    standardGeneric("estimateDiversity"))

#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x="SummarizedExperiment"),
    function(x, assay_name = abund_values, abund_values = "counts",
            index = c("coverage", "fisher", "gini_simpson", 
                    "inverse_simpson", "log_modulo_skewness", "shannon"),
                    name = index, ..., BPPARAM = SerialParam()){

        # input check
        index<- match.arg(index, several.ok = TRUE)
        
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                "same length than 'index'.",
                call. = FALSE)
        }
        .check_assay_present(assay_name, x)
        .require_package("vegan")

        dvrsts <- BiocParallel::bplapply(index,
                                        .get_diversity_values,
                                        x = x,
                                        mat = assay(x, assay_name),
                                        BPPARAM = BPPARAM,
                                        ...)
        .add_values_to_colData(x, dvrsts, name)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x="TreeSummarizedExperiment"),
    function(x, assay_name = abund_values, abund_values = "counts",
            index = c("coverage", "faith", "fisher", "gini_simpson", 
                    "inverse_simpson", "log_modulo_skewness", "shannon"),
            name = index, tree_name = "phylo", 
            ..., BPPARAM = SerialParam()){
        # input check
        # Check tree_name
        if( !.is_non_empty_character(tree_name) ){
            stop("'tree_name' must be a character specifying a rowTree of 'x'.",
                 call. = FALSE)
        }
        # Check indices
        index <- match.arg(index, several.ok = TRUE)
        if(!.is_non_empty_character(name) || length(name) != length(index)){
            stop("'name' must be a non-empty character value and have the ",
                "same length than 'index'.",
                call. = FALSE)
        }
        
        # If 'faith' is one of the indices
        if( "faith" %in% index ){
            # Get the name of "faith" index
            faith_name <- name[index %in% "faith"]
            # Store original names
            name_original <- name
            # And delete it from name
            name <- name[!index %in% "faith"]

            # Delete "faith" from indices
            index <- index[!index %in% "faith"]
            
            # Faith will be calculated
            calc_faith <- TRUE
        } else{
            # Faith will not be calculated
            calc_faith <- FALSE
        }
        
        # If index list contained other than 'faith' index, the length of the
        # list is over 0
        if( length(index)>0){
            # Calculates all indices but not 'faith'
            x <- callNextMethod()
        }
        # If 'faith' was one of the indices, 'calc_faith' is TRUE
        if( calc_faith ){
            # Get tree to check whether faith can be calculated
            tree <- rowTree(x, tree_name)
            # Check if faith can be calculated. Give warning and do not run estimateFaith
            # if there is no rowTree and other indices were also calculated. Otherwise, 
            # run estimateFaith. (If there is no rowTree --> error)
            if( (is.null(tree) || is.null(tree$edge.length)) &&
                length(index) >= 1 ){
                warning("Object does not have a tree called 'tree_name' or the tree does not ",
                        "have any branches. \nThe 'faith' alpha diversity index. ",
                        "cannot be calculated without rowTree. Therefore it is excluded ",
                        "from the results. \nYou can consider adding rowTree to include this index.",            
                        call. = FALSE)
            } else{
                x <- estimateFaith(x, name = faith_name, tree_name = tree_name, ...)
                # Ensure that indices are in correct order
                colnames <- colnames(colData(x))
                colnames <- c(colnames[ !colnames %in% name_original ], name_original)
                colData(x) <- colData(x)[ , colnames]
            }
        }
        return(x)
    }
)

#' @rdname estimateDiversity
#' @export
setGeneric("estimateFaith",signature = c("x", "tree"),
            function(x, tree = "missing", 
                    assay_name = abund_values, abund_values = "counts",
                    name = "faith", ...)
            standardGeneric("estimateFaith"))

#' @rdname estimateDiversity
#' @export
setMethod("estimateFaith", signature = c(x="SummarizedExperiment", tree="phylo"),
    function(x, tree, assay_name = abund_values, abund_values = "counts",
            name = "faith", node_lab = NULL, ...){
        # Input check
        # Check 'tree'
        # IF there is no rowTree gives an error
        if( is.null(tree) || is.null(tree$edge.length) ){
            stop("'tree' is NULL or it does not have any branches.",
                "'faith' is not possible to calculate.",
                call. = FALSE)
        }
        # Check 'assay_name'
        .check_assay_present(assay_name, x)
        # Check that it is numeric
        if( !is.numeric(assay(x, assay_name)) ){
            stop("The abundance matrix specificied by 'assay_name' must be numeric.",
                 call. = FALSE)
        }
        # Check 'name'
        if(!.is_non_empty_character(name)){
            stop("'name' must be a non-empty character value.",
                call. = FALSE)
        }
        # Check that node_lab is NULL or it specifies links between rownames and 
        # node labs
        if( !( is.null(node_lab) || 
               is.character(node_lab) && length(node_lab) == nrow(x) ) ){
            stop("'node_lab' must be NULL or a vector specifying links between ",
                 "rownames and node labs of 'tree'.",
                 call. = FALSE)
        }
        # Get the abundance matrix
        mat <- assay(x, assay_name)
        # Check that it is numeric
        if( !is.numeric(mat) ){
            stop("The abundance matrix specificied by 'assay_name' must be numeric.",
                 call. = FALSE)
        }
        # Subset and rename rows of the assay to correspond node_labs
        if( !is.null(node_lab) ){
            # Subset 
            mat <- mat[ !is.na(node_lab), ]
            node_lab <- node_lab[ !is.na(node_lab) ]
            # Rename
            rownames(mat) <- node_lab
        }
        # Calculates Faith index
        faith <- list(.calc_faith(mat, tree))
        # Adds calculated Faith index to colData
        .add_values_to_colData(x, faith, name)
    }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateFaith", signature = c(x="TreeSummarizedExperiment", tree="missing"),
    function(x, assay_name = abund_values, abund_values = "counts",
            name = "faith", tree_name = "phylo", ...){
        # Check tree_name
        if( !.is_non_empty_character(tree_name) ){
            stop("'tree_name' must be a character specifying a rowTree of 'x'.",
                 call. = FALSE)
        }
        # Gets the tree
        tree <- rowTree(x, tree_name)
        if( is.null(tree) || is.null(tree$edge.length)){
            stop("rowTree(x, tree_name) is NULL or the tree does not have any branches. ",
            "Faith's diversity cannot be calculated.",
                call. = FALSE)
        }
        # Get node labs
        node_lab <- rowLinks(x)[ , "nodeLab" ]
        node_lab[ rowLinks(x)[, "whichTree"] != tree_name ] <- NA
        # Give a warning, data will be subsetted
        if( any(is.na(node_lab)) ){
            warning("The rowTree named 'tree_name' does not include all the ",
                    "rows which is why 'x' is subsetted when 'faith' alpha ",
                    "diversity index is calculated.",
                    call. = FALSE)
        }
        # Calculates the Faith index
        estimateFaith(x, tree, name = name, node_lab = node_lab, ...)
    }
)


################################################################################

.calc_shannon <- function(mat, ...){
    vegan::diversity(t(mat), index="shannon")
}

# NOTE: vegan::diversity(x, index = "simpson")
# gives Simpson diversity, also called Gini-Simpson
# index: 1-lambda, where lambda is the Simpson index
# (lambda). This may cause confusion if your familiarity
# with diversity indices is limited.
# Moreover, Simpson's lambda is simply the
# squared sum of relative abundances so we can
# just use that for clarity and simplicity.
#.get_simpson <- function(x, ...){
.simpson_lambda <- function(mat, ...){

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Squared sum of relative abundances
    colSums2(rel^2)
}

.calc_gini_simpson <- function(mat, ...){
    1 - .simpson_lambda(mat, ...)
}

.calc_inverse_simpson <- function(mat, ...){
    1 / .simpson_lambda(mat, ...)
}

.calc_coverage <- function(mat, threshold = 0.9, ...){

    # Threshold must be a numeric value between 0-1
    if( !( is.numeric(threshold) && (threshold >= 0 && threshold <= 1) ) ){
        stop("'threshold' must be a numeric value between 0-1.",
            call. = FALSE)
    }

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Number of groups needed to have threshold (e.g. 50 %) of the
    # ecosystem occupied
    coverage <- apply(rel, 2, function(x) {
        min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))
    })
    names(coverage) <- colnames(rel)
    coverage
}

.calc_fisher <- function(mat, ...){
    vegan::fisher.alpha(t(mat))
}

.calc_faith <- function(mat, tree, ...){

    # Gets vector where number represent nth sample
    samples <- seq_len(ncol(mat))

    # Repeats taxa as many times there are samples, i.e. get all the
    # taxa that are
    
    # analyzed in each sample.
    taxa <- rep(rownames(mat), length(samples))

    # Gets those taxa that are present/absent in each sample.
    # Gets one big list that combines
    # taxa from all the samples.
    present_combined <- taxa[ mat[, samples] > 0 ]
    absent_combined <- taxa[ mat[, samples] == 0 ]
    
    # Gets how many taxa there are in each sample. 
    # After that, determines indices of samples' first taxa with cumsum.
    split_present <- as.vector(cumsum(colSums(mat > 0)))
    split_absent <- as.vector(cumsum(colSums(mat == 0)))
    
    # Determines which taxa belongs to which sample by first determining
    # the splitting points,
    # and after that giving every taxa number which tells their sample.
    split_present <- as.factor(cumsum((seq_along(present_combined)-1) %in%
                        split_present))
    split_absent <- as.factor(cumsum((seq_along(absent_combined)-1) %in%
                        split_absent))
    
    # Assigns taxa to right samples based on their number that they got from
    # previous step, and deletes unnecessary names.
    present <- unname(split(present_combined, split_present))
    absent <- unname(split(absent_combined, split_absent))

    # Assign NA to all samples
    faiths <- rep(NA,length(samples))

    # If there is one taxon present, then faith is the age of taxon
    f <- lengths(present) == 1
    faiths[f] <- tree$ages[which(tree$edge[, 2] == which(tree$tip.label == present[f]))]

    # If all the taxa are present, then faith is the sum of all edges of taxa
    faiths[lengths(absent) == 0] <- sum(tree$edge.length)

    # If there are taxa that are not present,
    f <- lengths(absent) > 0
    # absent taxa are dropped
    trees <- lapply(absent[f], ape::drop.tip, phy = tree)
    # and faith is calculated based on the subset tree
    faiths[f] <- vapply(trees, function(t){sum(t$edge.length)},numeric(1))

    # IF there are no taxa present, then faith is 0
    faiths[lengths(present) == 0] <- 0

    return(faiths)
}

.calc_log_modulo_skewness <- function(mat, quantile = 0.5, num_of_classes = 50, ...){
    # quantile must be a numeric value between 0-1
    if( !( is.numeric(quantile) && (quantile >= 0 && quantile <= 1) ) ){
        stop("'quantile' must be a numeric value between 0-1.",
            call. = FALSE)
    }
    # num_of_classes must be a positive numeric value
    if( !( is.numeric(num_of_classes) && num_of_classes > 0 ) ){
        stop("'num_of_classes' must be a positive numeric value.",
            call. = FALSE)
    }
    # Determine the quantile point.
    quantile_point <- quantile(max(mat), quantile)
    # Tabulate the arithmetic abundance classes. Use the same classes
    # for all samples for consistency
    cutpoints <- c(seq(0, quantile_point, length=num_of_classes), Inf)
    # Calculates sample-wise frequencies. How many taxa in each interval?
    freq_table <- table(cut(mat, cutpoints), col(mat))
    # Calculates the skewness of frequency table. Returns skewness for each
    # sample
    r <- .calc_skewness(freq_table)
    # Return log-modulo
    log(1 + r)
}

#' @importFrom DelayedMatrixStats rowSums2 rowMeans2
.calc_skewness <- function(x) {
    # Transposes the table
    x <- t(x)
    # Each value is substracted by sample-wise mean, which is raised to the
    # power of 3.
    # Then the sample-wise sum is taken from these values. 
    numerator <- rowSums2((x - rowMeans2(x))^3)
    # Sample-wise sum is divided by number of taxa that are not NA.
    numerator <- numerator/rowSums2(!is.na(x))
    # Each value is substracted by sample-wise mean, which is raises to the
    # power of 2.
    # Then the sample-wise sum is taken from these values. 
    denominator <- rowSums2((x - rowMeans2(x))^2)
    # Sample-wise sum is divided by number of taxa that are not NA. Then
    # these values
    # are raised to the power of 3/2.
    denominator <- (denominator/rowSums2(!is.na(x)))^(3/2)
    # Result
    result <- numerator/denominator
    return(result)
}

#' @importFrom SummarizedExperiment assay assays
.get_diversity_values <- function(index, x, mat, tree, ...){
    FUN <- switch(index,
                        shannon = .calc_shannon,
                        gini_simpson = .calc_gini_simpson,
                        inverse_simpson = .calc_inverse_simpson,
                        coverage = .calc_coverage,
                        fisher = .calc_fisher,
                        faith = .calc_faith,
                        log_modulo_skewness = .calc_log_modulo_skewness
                        )

    FUN(x = x, mat = mat, tree = tree, ...)
}
