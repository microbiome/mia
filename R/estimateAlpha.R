#' Estimate alpha indices using rarefaction
#' 
#' The function estimates alpha diversity measures optionally using n rounds of rarefaction,
#' given the rarefaction depth, then stores results at \code{\link{colData}}.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' 
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param index a \code{character} vector, specifying the alpha diversity measures
#'   to be calculated
#'   
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices specifying the alpha diversity measures used.
#'   
#' @param ... optional arguments.
#' 
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#' 
#' @param rarify logical scalar: Should the alpha diversity measures be estimated
#'   using rarefaction? (default: \code{FALSE})
#' 
#' @param seed a single \code{integer} value as the seed used for the nround
#'  rarefaction.
#' 
#' @param nrounds a single \code{integer} value for the number of rarefaction
#' rounds.
#' 
#' @param rarefaction_depth a \code{double} value as for the minimim size or 
#' rarefaction_depth. (default: \code{min(colSums(assay(x, "counts")), na.rm = TRUE)})
#' 
#' @return \code{x} with additional \code{\link{colData}} named after the index 
#' used.
#' 
#' @examples
#' 
#' data("GlobalPatterns")
#' tse <- GlobalPatterns
#' 
#' # Calculate the default Shannon index with no rarefaction
#' tse <- estimateAlpha(tse, assay.type = "counts", index = "shannon")
#' 
#' # Shows the estimated Shannon index
#' colData(tse)$shannon_diversity
#'
#'# Calculate observed richness with 10 rarefaction rounds
#' tse <- estimateAlpha(tse, assay.type = "counts", index = "observed_richness",
#' rarify=TRUE, nrounds=10)
#' 
#' # Shows the estimated observed richness
#' colData(tse)$observed_richness
#' 
#' @importFrom dplyr %>% 
#' 
#' @rdname estimateAlpha
#' @export
estimateAlpha <- function(x, assay.type = "counts", assay_name = NULL,
                          index = c("coverage_diversity", "fisher_diversity",
                                    "faith_diversity", "faith",
                                    "gini_simpson_diversity", "inverse_simpson_diversity",
                                    "log_modulo_skewness_diversity", "shannon_diversity",
                                    "absolute_dominance", "dbp_dominance",
                                    "core_abundance_dominance", "gini_dominance",
                                    "dmn_dominance", "relative_dominance",
                                    "simpson_lambda_dominance",
                                    "camargo_evenness", "pielou_evenness",
                                    "simpson_evenness", "evar_evenness",
                                    "bulla_evenness",
                                    "ace_richness", "chao1_richness", "hill_richness",
                                    "observed_richness"),
                          name = index,
                          ...,
                          BPPARAM = SerialParam(),
                          rarify=FALSE,
                          seed = 123,
                          nrounds=10,
                          rarefaction_depth=min(colSums(assay(x, "counts")), na.rm = TRUE)){
    if(!.is_non_empty_string(index)) {
        stop("'index' should be a non empty string.",
             call. = FALSE)
    }
    if(!.is_a_bool(rarify)){
        stop("'rarify' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_an_integer(seed)) {
        stop("'seed' must be an interger.",
             call. = FALSE)
    }
    if(!.is_an_integer(nrounds)) {
        stop("'nrounds' must be an integer.",
             call. = FALSE)
    }
    if(!(is.double(rarefaction_depth) & rarefaction_depth > 0)) {
        stop("'rarefaction_depth' must be a non-zero positive double.",
             call. = FALSE)
    }
    diversity_indices <- c("coverage_diversity", "coverage",
                           "faith_diversity", "faith",
                           "fisher_diversity", "fisher",
                           "gini_simpson_diversity", "gini_simpson",
                           "inverse_simpson_diversity", "inverse_simpson",
                           "log_modulo_skewness_diversity", "log_modulo_skewness",
                           "shannon_diversity", "shannon")
    dominance_indices <- c("absolute_dominance", "absolute",
                           "dbp_dominance", "dbp",
                           "core_abundance_dominance", "core_abundance",
                           "gini_dominance", "gini", 
                           "dmn_dominance", "dmn",
                           "relative_dominance", "relative",
                           "simpson_lambda_dominance", "simpson_lambda")
    evenness_indices <- c("camargo_evenness", "camargo",
                         "pielou_evenness", "pielou",
                         "simpson_evenness",
                         "evar_evenness", "evar",
                         "bulla_evenness", "bulla")
    richness_indices <- c("ace_richness", "ace",
                          "chao1_richness", "chao1",
                          "hill_richness", "hill",
                          "observed_richness", "observed")
    FUN <- NULL
    if(index %in% diversity_indices) {
        name <- .parse_name(index, name, "diversity")
        index <- gsub("_diversity", "", index)
        FUN <- estimateDiversity
    } else if(index %in% dominance_indices) {
        name <-  .parse_name(index, name, "dominance")
        index <- gsub("_dominance", "", index)
        FUN <- estimateDominance
    } else if (index %in% evenness_indices) {
        name <- .parse_name(index, name, "evenness")
        if (index!="simpson_evenness") {
            index <- gsub("_evenness", "", index)
        }
        FUN <- estimateEvenness
    } else if (index %in% richness_indices) {
        name <- .parse_name(index, name, "richness")
        index <- gsub("_richness", "", index)
        FUN <- estimateRichness
    } else {
        stop("'index' is coresponding to none of the alpha diversity measures.",
             call. = FALSE)
    }
    
    if (rarify) {
        .alpha_rarefaction(x, nrounds = nrounds, seed = seed,
                           args.sub = list(assay.type=assay.type,
                                           min_size=rarefaction_depth,
                                           verbose=FALSE),
                           FUN=FUN,
                           args.fun=list(index=index,
                                         assay.type="subsampled",
                                         ...,
                                         BPPARAM=BPPARAM),
                           name=name)
    } else {
        suppressWarnings(do.call(FUN, list(x, assay.type=assay.type, assay_name=assay_name,
                          index=index, name=name, ..., BPPARAM=BPPARAM)))
    }

}

.alpha_rarefaction <- function(x,
                               nrounds=1L,
                               seed=123,
                               args.sub=list(assay.type="counts",
                                             min_size=min(colSums(assay(x, "counts")),
                                                          na.rm = TRUE),
                                             verbose=FALSE),
                               FUN=estimateDiversity,
                               args.fun=list(index="shannon",
                                             assay.type="subsampled",
                                             ...,
                                             BPPARAM=BPPARAM),
                               name = args.fun$index) {
    set.seed(seed)
    colData(x)[, name] <- lapply(seq(nrounds), function(i){
        x_sub <- do.call(subsampleCounts, append(list(x), args.sub))
        suppressWarnings(x_sub <- do.call(FUN, append(list(x_sub), args.fun)))
        colData(x_sub)[, args.fun$index, drop=FALSE]
    }) %>% as.data.frame() %>% rowMeans() %>% as.data.frame()
    return(x)
}

.parse_name <- function(index, name, measure) {
    # don't change name if defined by user
    if (name==index) {
        if (measure %in% unlist(strsplit(index, "\\_"))) {
            name = index
        } else {
            name = paste0(index, "_", measure)
        }
    } else {
        return(name)
    }
}

################################# Alpha Functions ##############################

#' @rdname estimateDiversity
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
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
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
#'   \item{only.tips}{ A boolean value specifying whether to remove internal
#'   nodes when Faith's inex is calculated. When \code{only.tips=TRUE}, those
#'   rows that are not tips of tree are removed.
#'   (By default: \code{only.tips=FALSE})}
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
#' rowTree. (Faith 1992)
#' 
#' If the data includes features that are not in tree's tips but in
#' internal nodes, there are two options. First, you can keep those features,
#' and prune the tree to match features so that each tip can be found from
#' the features. Other option is to remove all features that are not tips.
#' (See \code{only.tips} parameter)}
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
#' @export
setGeneric("estimateDiversity",signature = c("x"),
           function(x, assay.type = "counts", assay_name = NULL,
                    index = c("coverage_diversity", "coverage",
                              "faith_diversity", "faith",
                              "fisher_diversity", "fisher",
                              "gini_simpson_diversity", "gini_simpson",
                              "inverse_simpson_diversity", "inverse_simpson",
                              "log_modulo_skewness_diversity", "log_modulo_skewness",
                              "shannon_diversity", "shannon"),
                    name = index, ...)
               standardGeneric("estimateDiversity"))

#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x="SummarizedExperiment"),
          function(x, assay.type = "counts", assay_name = NULL,
                   index = c("coverage_diversity", "coverage",
                             "faith_diversity", "faith",
                             "fisher_diversity", "fisher",
                             "gini_simpson_diversity", "gini_simpson",
                             "inverse_simpson_diversity", "inverse_simpson",
                             "log_modulo_skewness_diversity", "log_modulo_skewness",
                             "shannon_diversity", "shannon"),
                   name = index, ..., BPPARAM = SerialParam()){
              .Deprecated(old="estimateDiversity", new="estimateAlpha",
                          "Now estimateDiversity is deprecated. Use estimateAlpha instead.")
              if (!is.null(assay_name)) {
                  .Deprecated(old="assay_name", new="assay.type",
                              "Now assay_name is deprecated. Use assay.type instead.")
              }
              
              # input check
              index<- match.arg(index, several.ok = TRUE)
              
              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }
              .check_assay_present(assay.type, x)
              .require_package("vegan")
              
              dvrsts <- BiocParallel::bplapply(index,
                                               .get_diversity_values,
                                               x = x,
                                               mat = assay(x, assay.type),
                                               BPPARAM = BPPARAM,
                                               ...)
              .add_values_to_colData(x, dvrsts, name)
          }
)

#' @rdname estimateDiversity
#' @export
setMethod("estimateDiversity", signature = c(x="TreeSummarizedExperiment"),
          function(x, assay.type = "counts", assay_name = NULL,
                   index = c("coverage_diversity", "coverage",
                             "faith_diversity", "faith",
                             "fisher_diversity", "fisher",
                             "gini_simpson_diversity", "gini_simpson",
                             "inverse_simpson_diversity", "inverse_simpson",
                             "log_modulo_skewness_diversity", "log_modulo_skewness",
                             "shannon_diversity", "shannon"),
                   name = index, tree_name = "phylo", 
                   ..., BPPARAM = SerialParam()){
              .Deprecated(old="estimateDiversity", new="estimateAlpha",
                          "Now estimateDiversity is deprecated. Use estimateAlpha instead.")
              # input check
              # Check tree_name
              if( !.is_non_empty_string(tree_name) ){
                  stop("'tree_name' must be a character specifying a rowTree of 'x'.",
                       call. = FALSE)
              }
              if (!is.null(assay_name)) {
                  .Deprecated(old="assay_name", new="assay.type",
                              "Now assay_name is deprecated. Use assay.type instead.")
              }	
              # Check indices
              index <- match.arg(index, several.ok = TRUE)
              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }
              
              # If 'faith' is one of the indices
              if( "faith" %in% unlist(strsplit(index, "\\_")) ){
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
                      warning("Faith diversity has been excluded from the results ",
                              "since it cannot be calculated without rowTree. ",
                              "This requires a rowTree in the input argument x. ",
                              "Make sure that 'rowTree(x)' is not empty, or ",
                              "make sure to specify 'tree_name' in the input ",
                              "arguments. Warning is also provided if the tree does ",
                              "not have any branches. You can consider adding ",
                              "rowTree to include this index.", 
                              call. = FALSE)
                  } else {
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

#' @rdname estimateFaith
#' @export
setGeneric("estimateFaith",signature = c("x", "tree"),
           function(x, tree = "missing", 
                    assay.type = "counts", assay_name = NULL,
                    name = "faith", ...)
               standardGeneric("estimateFaith"))

#' @rdname estimateFaith
#' @export
setMethod("estimateFaith", signature = c(x="SummarizedExperiment", tree="phylo"),
          function(x, tree, assay.type = "counts", assay_name = NULL,
                   name = "faith", node_lab = NULL, ...){
              .Deprecated(old="estimateFaith", new="estimateAlpha",
                          "Now estimateFaith is deprecated. Use estimateAlpha instead.")
              # Input check
              # Check 'tree'
              # IF there is no rowTree gives an error
              if( is.null(tree) || is.null(tree$edge.length) ){
                  stop("'tree' is NULL or it does not have any branches.",
                       "The Faith's alpha diversity index is not possible to calculate.",
                       call. = FALSE)
              }
              # Check 'assay.type'
              .check_assay_present(assay.type, x)
              # Check that it is numeric
              if( !is.numeric(assay(x, assay.type)) ){
                  stop("The abundance matrix specificied by 'assay.type' must be numeric.",
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
              mat <- assay(x, assay.type)
              # Check that it is numeric
              if( !is.numeric(mat) ){
                  stop("The abundance matrix specificied by 'assay.type' must be numeric.",
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
              faith <- list(.calc_faith(mat, tree, ...))
              # Adds calculated Faith index to colData
              .add_values_to_colData(x, faith, name)
          }
)

#' @rdname estimateFaith
#' @export
setMethod("estimateFaith", signature = c(x="TreeSummarizedExperiment", tree="missing"),
          function(x, assay.type = "counts", assay_name = NULL,
                   name = "faith", tree_name = "phylo", ...){
              .Deprecated(old="estimateFaith", new="estimateAlpha",
                          "Now estimateFaith is deprecated. Use estimateAlpha instead.")
              # Check tree_name
              if( !.is_non_empty_character(tree_name) ){
                  stop("'tree_name' must be a character specifying a rowTree of 'x'.",
                       call. = FALSE)
              }
              # Gets the tree
              tree <- rowTree(x, tree_name)
              if( is.null(tree) || is.null(tree$edge.length)){
                  stop("rowTree(x, tree_name) is NULL or the tree does not have any branches. ",
                       "The Faith's alpha diversity index cannot be calculated.",
                       call. = FALSE)
              }
              # Get node labs
              node_lab <- rowLinks(x)[ , "nodeLab" ]
              node_lab[ rowLinks(x)[, "whichTree"] != tree_name ] <- NA
              # Give a warning, data will be subsetted
              if( any(is.na(node_lab)) ){
                  warning("The rowTree named 'tree_name' does not include all the ",
                          "rows which is why 'x' is subsetted when the Faith's alpha ",
                          "diversity index is calculated.",
                          call. = FALSE)
              }
              # Calculates the Faith index
              estimateFaith(x, tree, name = name, node_lab = node_lab, ...)
          }
)

#' @rdname estimateDominance
#' Estimate dominance measures
#'
#' This function calculates community dominance indices.
#' This includes the \sQuote{Absolute}, \sQuote{Berger-Parker},
#' \sQuote{Core abundance},
#' \sQuote{Gini}, \sQuote{McNaughton’s}, \sQuote{Relative}, and
#' \sQuote{Simpson's} indices.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the sample-wise estimates.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param index a \code{character} vector, specifying the indices to be
#'   calculated.
#'
#' @param ntaxa Optional and only used for the \code{Absolute} and
#'   \code{Relative} dominance indices: The n-th position of the dominant taxa
#'   to consider (default: \code{ntaxa = 1}). Disregarded for the indices
#'   \dQuote{dbp},
#'   \dQuote{core_abundance}, \dQuote{Gini}, \dQuote{dmn}, and \dQuote{Simpson}.
#'
#' @param aggregate Optional and only used for the \code{Absolute}, \code{dbp},
#'   \code{Relative}, and \code{dmn} dominance indices:
#'   Aggregate the values for top members selected by \code{ntaxa} or not. If
#'   \code{TRUE}, then the sum of relative abundances is returned. Otherwise the
#'   relative abundance is returned for the single taxa with the indicated rank
#'   (default: \code{aggregate = TRUE}). Disregarded for the indices
#'   \dQuote{core_abundance}, \dQuote{gini}, \dQuote{dmn}, and \dQuote{simpson}.
#'
#' @param name A name for the column(s) of the colData where the calculated
#'   Dominance indices should be stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'   (Currently not used)
#'
#' @param ... additional arguments currently not used.
#'
#' @details
#'
#' A dominance index quantifies the dominance of one or few species in a
#' community. Greater values indicate higher dominance.
#'
#' Dominance indices are in general negatively correlated with alpha diversity
#' indices (species richness, evenness, diversity, rarity). More dominant
#' communities are less diverse.
#'
#' \code{estimateDominance} calculates the following community dominance
#' indices:
#'
#' \itemize{
#' 
#' \item{'absolute' }{Absolute index equals to the absolute abundance of the
#' most dominant n species of the sample (specify the number with the argument
#' \code{ntaxa}). Index gives positive integer values.}
#' 
#' \item{'dbp' }{Berger-Parker index (See Berger & Parker 1970) calculation
#' is a special case of the 'relative' index. dbp is the relative abundance of
#' the most
#' abundant species of the sample. Index gives values in interval 0 to 1,
#' where bigger value represent greater dominance.
#'
#' \deqn{dbp = \frac{N_1}{N_{tot}}}{%
#' dbp = N_1/N_tot} where \eqn{N_1} is the absolute abundance of the most
#' dominant species and \eqn{N_{tot}} is the sum of absolute abundances of all
#' species.}
#' 
#' \item{'core_abundance' }{ Core abundance index is related to core species.
#' Core species are species that are most abundant in all samples, i.e., in
#' whole data set. Core species are defined as those species that have
#' prevalence over 50\%. It means that in order to belong to core species,
#' species must be prevalent in 50\% of samples. Core species are used to
#' calculate the core abundance index. Core abundance index is sum of relative
#' abundances of core species in the sample. Index gives values in interval
#' 0 to 1, where bigger value represent greater dominance.
#'
#' \deqn{core_abundance = \frac{N_{core}}{N_{tot}}}{%
#' core_abundance = N_core/N_tot} where \eqn{N_{core}} is the sum of absolute
#' abundance of the core species and \eqn{N_{tot}} is the sum of absolute
#' abundances of all species.}
#' 
#' \item{'gini' }{ Gini index is probably best-known from socio-economic
#' contexts (Gini 1921). In economics, it is used to measure, for example, how
#' unevenly income is distributed among population. Here, Gini index is used
#' similarly, but income is replaced with abundance. 
#' 
#' If there is small group of species
#' that represent large portion of total abundance of microbes, the inequality
#' is large and Gini index closer to 1. If all species has equally large
#' abundances, the equality is perfect and Gini index equals 0. This index
#' should not be confused with Gini-Simpson index, which quantifies diversity.}
#'
#' \item{'dmn' }{McNaughton’s index is the sum of relative abundances of the two
#' most abundant species of the sample (McNaughton & Wolf, 1970). Index gives
#' values in the unit interval:
#'
#' \deqn{dmn = (N_1 + N_2)/N_tot}
#'
#' where \eqn{N_1} and \eqn{N_2} are the absolute
#' abundances of the two most dominant species and \eqn{N_{tot}} is the sum of
#' absolute abundances of all species.}
#'
#' \item{'relative' }{ Relative index equals to the relative abundance of the
#' most dominant n species of the sample (specify the number with the
#' argument \code{ntaxa}).
#' This index gives values in interval 0 to 1.
#'
#' \deqn{relative = N_1/N_tot}
#'
#' where \eqn{N_1} is the absolute abundance of the most
#' dominant species and \eqn{N_{tot}} is the sum of absolute abundances of all
#' species.}
#'
#' \item{'simpson_lambda' }{ Simpson's (dominance) index or Simpson's lambda is
#' the sum of squared relative abundances. This index gives values in the unit interval.
#' This value equals the probability that two randomly chosen individuals
#' belongs to the
#' same species. The higher the probability, the greater the dominance (See
#' e.g. Simpson 1949).
#'
#' \deqn{lambda = \sum(p^2)}
#'
#' where p refers to relative abundances.
#'
#' There is also a more advanced Simpson dominance index (Simpson 1949).
#' However, this is not provided and the simpler squared sum of relative
#' abundances is used instead as the alternative index is not in the unit
#' interval and it is highly
#' correlated with the simpler variant implemented here.}
#' 
#' }
#'
#' @references
#'
#' Berger WH & Parker FL (1970)
#' Diversity of Planktonic Foraminifera in Deep-Sea Sediments.
#' _Science_ 168(3937):1345-1347. doi: 10.1126/science.168.3937.1345
#'
#' Gini C (1921)
#' Measurement of Inequality of Incomes.
#' _The Economic Journal_ 31(121): 124-126. doi: 10.2307/2223319
#'
#' McNaughton, SJ and Wolf LL. (1970).
#' Dominance and the niche in ecological systems.
#' _Science_ 167:13, 1--139
#'
#' Simpson EH (1949)
#' Measurement of Diversity.
#' _Nature_ 163(688). doi: 10.1038/163688a0
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateEvenness]{estimateEvenness}}}
#'   \item{\code{\link[mia:estimateDiversity]{estimateDiversity}}}
#' }
#'
#' @name estimateDominance
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#'
#' # Calculates Simpson's lambda (can be used as a dominance index)
#' esophagus <- estimateDominance(esophagus, index="simpson_lambda")
#'
#' # Shows all indices
#' colData(esophagus)
#'
#' # Indices must be written correctly (e.g. dbp, not dbp), otherwise an error
#' # gets thrown
#' \dontrun{esophagus <- estimateDominance(esophagus, index="dbp")}
#' # Calculates dbp and Core Abundance indices
#' esophagus <- estimateDominance(esophagus, index=c("dbp", "core_abundance"))
#' # Shows all indices
#' colData(esophagus)
#' # Shows dbp index
#' colData(esophagus)$dbp
#' # Deletes dbp index
#' colData(esophagus)$dbp <- NULL
#' # Shows all indices, dbp is deleted
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices
#' esophagus <- estimateDominance(esophagus)
#' # Shows all indices
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices with explicitly specified names
#' esophagus <- estimateDominance(esophagus,
#'     index = c("dbp", "dmn", "absolute", "relative",
#'               "simpson_lambda", "core_abundance", "gini"),
#'     name  = c("BergerParker", "McNaughton", "Absolute", "Relative",
#'               "SimpsonLambda", "CoreAbundance", "Gini")
#' )
#' # Shows all indices
#' colData(esophagus)
#' @export
setGeneric("estimateDominance",signature = c("x"),
           function(x,
                    assay.type = assay_name, assay_name = "counts",
                    index = c("absolute_dominance", "absolute",
                              "dbp_dominance", "dbp",
                              "core_abundance_dominance", "core_abundance",
                              "gini_dominance", "gini", 
                              "dmn_dominance", "dmn",
                              "relative_dominance", "relative",
                              "simpson_lambda_dominance", "simpson_lambda"),
                    ntaxa = 1,
                    aggregate = TRUE,
                    name = index,
                    ...,
                    BPPARAM = SerialParam())
               standardGeneric("estimateDominance"))
#' @rdname estimateDominance
#' @export
setMethod("estimateDominance", signature = c(x = "SummarizedExperiment"),
          function(x,
                   assay.type = assay_name, assay_name = "counts",
                   index = c("absolute_dominance", "absolute",
                             "dbp_dominance", "dbp",
                             "core_abundance_dominance", "core_abundance",
                             "gini_dominance", "gini", 
                             "dmn_dominance", "dmn",
                             "relative_dominance", "relative",
                             "simpson_lambda_dominance", "simpson_lambda"),
                   ntaxa = 1,
                   aggregate = TRUE,
                   name = index,
                   ...,
                   BPPARAM = SerialParam()){
              .Deprecated(old="estimateDominance", new="estimateAlpha",
                          "Now estimateDominance is deprecated. Use estimateAlpha instead.")
              # Input check
              # Check assay.type
              .check_assay_present(assay.type, x)
              # Check indices
              index <- match.arg(index, several.ok = TRUE)
              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }
              
              # Check aggregate
              if(!.is_a_bool(aggregate)){
                  stop("'aggregate' must be TRUE or FALSE.", call. = FALSE)
              }
              
              # Calculates dominance indices
              dominances <- BiocParallel::bplapply(index,
                                                   FUN = .get_dominance_values,
                                                   mat = assay(x,assay.type),
                                                   ntaxa = ntaxa,
                                                   aggregate = aggregate,
                                                   BPPARAM = BPPARAM)
              
              # Add dominance indices to colData
              .add_values_to_colData(x, dominances, name)
          }
)

#' @rdname estimateEvenness
#' #' Estimate Evenness measures
#'
#' This function calculates community evenness indices.
#' These include the \sQuote{Camargo}, \sQuote{Pielou}, \sQuote{Simpson},
#' \sQuote{Evar} and \sQuote{Bulla} evenness measures.
#' See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} used for
#'   calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param index a \code{character} vector, specifying the evenness measures to be
#'   calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{threshold}{ a numeric threshold. assay values below or equal
#'     to this threshold will be set to zero.}
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#'
#' @details
#' Evenness is a standard index in community ecology, and it quantifies how evenly the abundances
#' of different species are distributed. The following evenness indices are provided:
#'
#' By default, this function returns all indices.
#'
#' The available evenness indices include the following (all in lowercase):
#' \itemize{
#'   \item{'camargo' }{Camargo's evenness (Camargo 1992)}
#'   \item{'simpson_evenness' }{Simpson’s evenness is calculated as inverse Simpson diversity (1/lambda) divided by
#'   observed species richness S: (1/lambda)/S.}
#'   \item{'pielou' }{Pielou's evenness (Pielou, 1966), also known as Shannon or Shannon-Weaver/Wiener/Weiner
#'     evenness; H/ln(S). The Shannon-Weaver is the preferred term; see Spellerberg and Fedor (2003).}
#'   \item{'evar' }{Smith and Wilson’s Evar index (Smith & Wilson 1996).}
#'   \item{'bulla' }{Bulla’s index (O) (Bulla 1994).}
#' }
#'   
#' Desirable statistical evenness metrics avoid strong bias towards very
#' large or very small abundances; are independent of richness; and range
#' within the unit interval with increasing evenness (Smith & Wilson 1996).
#' Evenness metrics that fulfill these criteria include at least camargo,
#' simpson, smith-wilson, and bulla. Also see Magurran & McGill (2011)
#' and Beisel et al. (2003) for further details.
#'
#' @references
#'
#' Beisel J-N. et al. (2003)
#' A Comparative Analysis of Evenness Index Sensitivity.
#' _Internal Rev. Hydrobiol._ 88(1):3-15.
#' URL: \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#'
#' Bulla L. (1994)
#' An  index  of  evenness  and  its  associated  diversity  measure.
#' _Oikos_ 70:167--171.
#'
#' Camargo, JA. (1992)
#' New diversity index for assessing structural alterations in aquatic communities.
#' _Bull. Environ. Contam. Toxicol._ 48:428--434.
#'
#' Locey KJ and Lennon JT. (2016)
#' Scaling laws predict global microbial diversity.
#' _PNAS_ 113(21):5970-5975; doi:10.1073/pnas.1521291113.
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Pielou, EC. (1966)
#' The measurement of diversity in different types of
#' biological collections. _J Theoretical Biology_ 13:131--144.
#'
#' Smith B and Wilson JB. (1996)
#' A Consumer's Guide to Evenness Indices.
#' _Oikos_ 76(1):70-82.
#'
#' Spellerberg and Fedor (2003).
#' A tribute to Claude Shannon (1916 –2001) and a plea for more rigorous use of species richness,
#' species diversity and the ‘Shannon–Wiener’ Index.
#' _Alpha Ecology & Biogeography_ 12, 177–197.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateDominance]{estimateDominance}}}
#'   \item{\code{\link[mia:estimateDiversity]{estimateDiversity}}}
#' }
#'
#' @name estimateEvenness
#'
#' @examples
#' data(esophagus)
#' tse <- esophagus
#'
#' # Specify index and their output names
#' index <- c("pielou", "camargo", "simpson_evenness", "evar", "bulla")
#' name  <- c("Pielou", "Camargo", "SimpsonEvenness",  "Evar", "Bulla")
#'
#' # Estimate evenness and give polished names to be used in the output
#' tse <- estimateEvenness(tse, index = index, name = name)
#'
#' # Check the output
#' head(colData(tse))
#'
#' @export
setGeneric("estimateEvenness",signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts",
                    index = c("camargo_evenness", "camargo",
                              "pielou_evenness", "pielou",
                              "simpson_evenness",
                              "evar_evenness", "evar",
                              "bulla_evenness", "bulla"),
                    name = index, ...)
               standardGeneric("estimateEvenness"))

#' @rdname estimateEvenness
#' @export
setMethod("estimateEvenness", signature = c(x = "SummarizedExperiment"),
          function(x, assay.type = assay_name, assay_name = "counts",
                   index = c("camargo_evenness", "camargo",
                             "pielou_evenness", "pielou",
                             "simpson_evenness",
                             "evar_evenness", "evar",
                             "bulla_evenness", "bulla"),
                   name = index, ..., BPPARAM = SerialParam()){
              .Deprecated(old="estimateEvenness", new="estimateAlpha",
                          "Now estimateEvenness is deprecated. Use estimateAlpha instead.")
              # input check
              index <- match.arg(index, several.ok = TRUE)
              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }
              .check_assay_present(assay.type, x)
              #
              vnss <- BiocParallel::bplapply(index,
                                             .get_evenness_values,
                                             mat = assay(x, assay.type),
                                             BPPARAM = BPPARAM, ...)
              .add_values_to_colData(x, vnss, name)
          }
)

#' @rdname estimateRichness
#' Estimate richness measures
#'
#' Several functions for calculation of community richness indices available via
#' wrapper functions. They are implemented via the \code{vegan} package.
#'
#' These include the \sQuote{ace}, \sQuote{Chao1}, \sQuote{Hill}, and 
#' \sQuote{Observed} richness measures.
#' See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object.
#'
#' @param assay.type the name of the assay used for calculation of the
#'   sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param index a \code{character} vector, specifying the richness measures
#'   to be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in.
#'
#' @param detection a numeric value for selecting detection threshold
#' for the abundances. The default detection threshold is 0.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... additional parameters passed to \code{estimateRichness}
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @details
#'
#' The richness is calculated per sample. This is a standard index in community
#' ecology, and it provides an estimate of the number of unique species in the
#' community. This is often not directly observed for the whole community but
#' only for a limited sample from the community. This has led to alternative
#' richness indices that provide different ways to estimate the species
#' richness.
#'
#' Richness index differs from the concept of species diversity or evenness in
#' that it ignores species abundance, and focuses on the binary presence/absence
#' values that indicate simply whether the species was detected.
#'
#' The function takes all index names in full lowercase. The user can provide
#' the desired spelling through the argument \code{\link{name}} (see examples).
#'
#' The following richness indices are provided.
#'
#' \itemize{
#'   
#'   \item{'ace' }{Abundance-based coverage estimator (ACE) is another
#'   nonparametric richness
#'   index that uses sample coverage, defined based on the sum of the
#'   probabilities
#'   of the observed species. This method divides the species into abundant
#'   (more than 10
#'   reads or observations) and rare groups
#'   in a sample and tends to underestimate the real number of species. The
#'   ACE index
#'   ignores the abundance information for the abundant species,
#'   based on the assumption that the abundant species are observed regardless
#'   of their
#'   exact abundance. We use here the bias-corrected version
#'   (O'Hara 2005, Chiu et al. 2014) implemented in
#'   \code{\link[vegan:specpool]{estimateR}}.
#'   For an exact formulation, see \code{\link[vegan:specpool]{estimateR}}.
#'   Note that this index comes with an additional column with standard
#'   error information.}
#'   
#'   \item{'chao1' }{This is a nonparametric estimator of species richness. It
#'   assumes that rare species carry information about the (unknown) number
#'   of unobserved species. We use here the bias-corrected version
#'   (O'Hara 2005, Chiu et al. 2014) implemented in
#'   \code{\link[vegan:specpool]{estimateR}}. This index implicitly
#'   assumes that every taxa has equal probability of being observed. Note
#'   that it gives a lower bound to species richness. The bias-corrected
#'   for an exact formulation, see \code{\link[vegan:specpool]{estimateR}}.
#'   This estimator uses only the singleton and doubleton counts, and
#'   hence it gives more weight to the low abundance species.
#'   Note that this index comes with an additional column with standard
#'   error information.}
#'   
#'   \item{'hill' }{Effective species richness aka Hill index
#'   (see e.g. Chao et al. 2016).
#'   Currently only the case 1D is implemented. This corresponds to the exponent
#'   of Shannon diversity. Intuitively, the effective richness indicates the
#'   number of
#'   species whose even distribution would lead to the same diversity than the
#'   observed
#'   community, where the species abundances are unevenly distributed.}
#'   
#'   \item{'observed' }{The _observed richness_ gives the number of species that
#'   is detected above a given \code{detection} threshold in the observed sample
#'   (default 0). This is conceptually the simplest richness index. The
#'   corresponding index in the \pkg{vegan} package is "richness".}
#'   
#' }
#'
#'
#' @references
#'
#' Chao A. (1984)
#' Non-parametric estimation of the number of classes in a population.
#' _Scand J Stat._ 11:265–270.
#'
#' Chao A, Chun-Huo C, Jost L (2016).
#' Phylogenetic Diversity Measures and Their Decomposition:
#' A Framework Based on Hill Numbers. Biodiversity Conservation and
#' Phylogenetic Systematics,
#' Springer International Publishing, pp. 141–172,
#' doi:10.1007/978-3-319-22461-9_8.
#'
#' Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014).
#' Improved nonparametric lower bound of species richness via a modified
#' Good-Turing frequency formula.
#' _Biometrics_ 70, 671-682.
#'
#' O'Hara, R.B. (2005).
#' Species richness estimators: how many species can dance on the head of a pin?
#' _J. Anim. Ecol._ 74, 375-386.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[vegan:specpool]{estimateR}}}
#' }
#'
#' @name estimateRichness
#'
#' @export
#'
#' @author Leo Lahti. Contact: \url{microbiome.github.io}
#'
#' @examples
#' data(esophagus)
#'
#' # Calculates all richness indices by default
#' esophagus <- estimateRichness(esophagus)
#'
#' # Shows all indices
#' colData(esophagus)
#'
#' # Shows Hill index
#' colData(esophagus)$hill
#'
#' # Deletes hill index
#' colData(esophagus)$hill <- NULL
#'
#' # Shows all indices, hill is deleted
#' colData(esophagus)
#'
#' # Delete the remaining indices
#' colData(esophagus)[, c("observed", "chao1", "ace")] <- NULL
#'
#' # Calculates observed richness index and saves them with specific names
#' esophagus <- estimateRichness(esophagus,
#'     index = c("observed", "chao1", "ace", "hill"),
#'      name = c("Observed", "Chao1", "ACE", "Hill"))
#'
#' # Show the new indices
#' colData(esophagus)
#'
#' # Deletes all colData (including the indices)
#' colData(esophagus) <- NULL
#'
#' # Calculate observed richness excluding singletons (detection limit 1)
#' esophagus <- estimateRichness(esophagus, index="observed", detection = 1)
#'
#' # Deletes all colData (including the indices)
#' colData(esophagus) <- NULL
#'
#' # Indices must be written correctly (all lowercase), otherwise an error
#' # gets thrown
#' \dontrun{esophagus <- estimateRichness(esophagus, index="ace")}
#'
#' # Calculates Chao1 and ACE indices only
#' esophagus <- estimateRichness(esophagus, index=c("chao1", "ace"),
#'                                           name=c("Chao1", "ACE"))
#'
#' # Deletes all colData (including the indices)
#' colData(esophagus) <- NULL
#'
#' # Names of columns can be chosen arbitrarily, but the length of arguments
#' # must match.
#' esophagus <- estimateRichness(esophagus,
#'                                index = c("ace", "chao1"),
#'                                name = c("index1", "index2"))
#' # Shows all indices
#' colData(esophagus)
#'
#' @export
setGeneric("estimateRichness",signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts",
                    index = c("ace_richness", "ace",
                              "chao1_richness", "chao1",
                              "hill_richness", "hill",
                              "observed_richness", "observed"),
                    name = index,
                    detection = 0,
                    ...,
                    BPPARAM = SerialParam())
               standardGeneric("estimateRichness"))

#' @rdname estimateRichness
#' @export
setMethod("estimateRichness", signature = c(x = "SummarizedExperiment"),
          function(x,
                   assay.type = assay_name, assay_name = "counts",
                   index = c("ace_richness", "ace",
                             "chao1_richness", "chao1",
                             "hill_richness", "hill",
                             "observed_richness", "observed"),
                   name = index,
                   detection = 0,
                   ...,
                   BPPARAM = SerialParam()){
              .Deprecated(old="estimateRichness", new="estimateAlpha",
                          "Now estimateRichness is deprecated. Use estimateAlpha instead.")
              # Input check
              # Check assay.type
              .check_assay_present(assay.type, x)
              # Check indices
              index <- match.arg(index, several.ok = TRUE)
              if(!.is_non_empty_character(name) || length(name) != length(index)){
                  stop("'name' must be a non-empty character value and have the ",
                       "same length than 'index'.",
                       call. = FALSE)
              }
              # Calculates richness indices
              richness <- BiocParallel::bplapply(index,
                                                 FUN = .get_richness_values,
                                                 mat = assay(x, assay.type),
                                                 detection = detection,
                                                 BPPARAM = BPPARAM)
              # Add richness indices to colData
              .add_values_to_colData(x, richness, name)
          }
)
   
################################# Utils  #######################################

## Diversity helper function

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

.calc_faith <- function(mat, tree, only.tips = FALSE, ...){
    # Input check
    if( !.is_a_bool(only.tips) ){
        stop("'only.tips' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Remove internal nodes if specified
    if( only.tips ){
        mat <- mat[ rownames(mat) %in% tree$tip.label, ]
    }
    # To ensure that the function works with NA also, convert NAs to 0.
    # Zero means that the taxon is not present --> same as NA (no information)
    mat[ is.na(mat) ] <- 0
    
    # Gets vector where number represent nth sample
    samples <- seq_len(ncol(mat))
    
    # Repeats taxa as many times there are samples, i.e. get all the
    # taxa that are analyzed in each sample.
    taxa <- rep(rownames(mat), length(samples))
    
    # Gets those taxa that are present/absent in each sample.
    # Gets one big list that combines
    # taxa from all the samples.
    present_combined <- taxa[ mat[, samples] > 0 ]
    
    # Gets how many taxa there are in each sample. 
    # After that, determines indices of samples' first taxa with cumsum.
    split_present <- as.vector(cumsum(colSums(mat > 0)))
    
    # Determines which taxa belongs to which sample by first determining
    # the splitting points,
    # and after that giving every taxa number which tells their sample.
    split_present <- as.factor(cumsum((seq_along(present_combined)-1) %in%
                                          split_present))
    
    # Assigns taxa to right samples based on their number that they got from
    # previous step, and deletes unnecessary names.
    present <- unname(split(present_combined, split_present))
    
    # If there were samples without any taxa present/absent, the length of the
    # list is not the number of samples since these empty samples are missing.
    # Add empty samples as NULL.
    names(present) <- names(which(colSums2(mat) > 0))
    present[names(which(colSums2(mat) == 0))] <- list(NULL)
    present <- present[colnames(mat)]
    
    # Assign NA to all samples
    faiths <- rep(NA,length(samples))
    
    # If there are no taxa present, then faith is 0
    ind <- lengths(present) == 0
    faiths[ind] <- 0
    
    # If there are taxa present
    ind <- lengths(present) > 0
    # Loop through taxa that were found from each sample
    faiths_for_taxa_present <- lapply(present[ind], function(x){
        # Trim the tree
        temp <- .prune_tree(tree, x)
        # Sum up all the lengths of edges
        temp <- sum(temp$edge.length)
        return(temp)
    })
    faiths_for_taxa_present <- unlist(faiths_for_taxa_present)
    faiths[ind] <- faiths_for_taxa_present
    return(faiths)
}

# This function trims tips until all tips can be found from provided set of nodes
#' @importFrom ape drop.tip
.prune_tree <- function(tree, nodes){
    # Get those tips that can not be found from provided nodes
    remove_tips <- tree$tip.label[!tree$tip.label %in% nodes]
    # As long as there are tips to be dropped, run the loop
    while( length(remove_tips) > 0 ){
        # Drop tips that cannot be found. Drop only one layer at the time. Some
        # dataset might have taxa that are not in tip layer but they are higher
        # higher rank. IF we delete more than one layer at the time, we might
        # loose the node for those taxa. --> The result of pruning is a tree
        # whose all tips can be found provided nodes i.e., rows of TreeSE. Some
        # taxa might be higher rank meaning that all rows might not be in tips
        # even after pruning; they have still child-nodes.
        tree <- drop.tip(tree, remove_tips, trim.internal = FALSE, collapse.singles = FALSE)
        # If all tips were dropped, the result is NULL --> stop loop
        if( is.null(tree) ){
            break
        }
        # Again, get those tips of updated tree that cannot be found from provided nodes
        remove_tips <- tree$tip.label[!tree$tip.label %in% nodes]
    }
    return(tree)
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


## Dominance helper function

.gini_dominance <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    x <- as.vector(x)
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

.calc_gini_dominance <- function(mat, ...){
    apply(mat, 2L, .gini_dominance)
}

.calc_core_dominance <- function(mat, ...){
    getPrevalentAbundance(mat, detection = 0, as_relative = TRUE)
}

.calc_dominance <- function(mat, ntaxa, aggregate, index){
    
    # Check ntaxa
    if(!(ntaxa>0 && ntaxa<3)){
        stop("'ntaxa' must be a numerical value 1 or 2.", call. = FALSE)
    }
    #
    if (index == "absolute") {
        # ntaxa=1 by default but can be tuned
        as_relative <- FALSE
    } else if (index == "relative") {
        # ntaxa=1 by default but can be tuned
        as_relative <- TRUE
    } else if (index == "dbp") {
        # Berger-Parker: if selected fix the following values
        ntaxa <- 1
        as_relative <- TRUE
    } else if (index == "dmn") {
        # McNaughton's dominance: if selected fix the following values
        ntaxa <- 2
        aggregate <- TRUE
        as_relative <- TRUE
    }
    
    if (as_relative) {
        # Calculates the relative abundance per sample
        mat <- .calc_rel_abund(mat)
    }
    
    # Aggregate or not
    if (!aggregate) {
        idx <- apply(mat, 2L,
                     function(mc) {
                         order(as.vector(mc), decreasing = TRUE)[[ntaxa]]
                     })
    } else {
        idx <- apply(mat, 2L,
                     function(mc) {
                         order(as.vector(mc), decreasing = TRUE)[seq_len(ntaxa)]
                     })
        idx <- split(as.vector(idx),
                     unlist(lapply(seq_len(length(idx) / ntaxa),rep.int,ntaxa)))
    }
    
    ans <- lapply(mapply(function(i,j,x){x[i,j]},
                         i = idx,
                         j = seq_len(ncol(mat)),
                         MoreArgs = list(x = mat),
                         SIMPLIFY = FALSE),
                  sum)
    ans <- unlist(ans)
    
    # Adds sample names to the table
    names(ans) <- colnames(mat)
    ans
}

.get_dominance_values <- function(index, mat, ntaxa = 1, aggregate = TRUE, ...) {
    
    FUN <- switch(index,
                  simpson_lambda = .simpson_lambda,
                  core_abundance = .calc_core_dominance,
                  gini = .calc_gini_dominance,
                  absolute = .calc_dominance,
                  relative = .calc_dominance,
                  dbp = .calc_dominance,
                  dmn = .calc_dominance
    )
    
    FUN(index, mat = mat, ntaxa = ntaxa, aggregate = aggregate, ...)
    
}

## evenness helper function

.calc_bulla_evenness <- function(mat) {
    # Species richness (number of species)
    S <- colSums2(mat > 0, na.rm = TRUE)
    
    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)
    
    i <- seq_len(nrow(p))
    O <- vapply(i,function(i){sum(pmin(p[i,], 1/S[i]))},numeric(1))
    
    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}

# Camargo's evenness x: species counts zeroes: include zeros Inspired
# by code from Pepijn de Vries and Zhou Xiang at
# researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here
.calc_camargo_evenness <- function(mat) {
    N <- colSums2(mat > 0, na.rm = TRUE)
    
    seq <- IntegerList(lapply(N - 1,seq_len))
    
    x <- mapply(
        function(i, n, s){
            xx <- 0
            for (j in s) {
                xx <- xx + sum(abs(mat[(j + 1):n,i] - mat[j,i]))
            }
            xx
        },
        seq_along(N),
        N,
        seq)
    # Return
    1 - x/(colSums2(mat, na.rm = TRUE) * N)
}

# x: Species count vector
.calc_simpson_evenness <- function(mat) {
    
    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)
    
    # Simpson evenness (Simpson diversity per richness)
    .calc_inverse_simpson(mat)/S
}

# x: Species count vector
.calc_pielou_evenness <- function(mat) {
    # Remove zeroes
    mat[mat == 0] <- NA
    
    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)
    
    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)
    
    # Shannon index
    H <- (-rowSums2(p * log(p), na.rm = TRUE))
    
    # Simpson evenness
    H/log(S)
}

# Smith and Wilson’s Evar index
.calc_evar_evenness <- function(mat) {
    N <- colSums2(mat, na.rm = TRUE)
    
    # Log abundance
    a <- log(mat)
    a[is.na(a) | is.infinite(a)] <- 0
    
    # Richness
    S <- colSums2(mat > 0, na.rm = TRUE)
    
    c <- colSums2(a, na.rm = TRUE)/S
    d <- t((t(a) - c)^2/S)
    d[mat == 0] <- 0
    
    f <- colSums2(d, na.rm = TRUE)
    
    (1 - 2/pi * atan(f))
}

.get_evenness_values <- function(index, mat, threshold = 0, ...){
    
    if(!is.numeric(threshold) || length(threshold) != 1L){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if(threshold > 0){
        mat[mat <= threshold] <- 0
    }
    
    FUN <- switch(index,
                  camargo = .calc_camargo_evenness,
                  pielou = .calc_pielou_evenness,
                  simpson_evenness = .calc_simpson_evenness,
                  evar = .calc_evar_evenness,
                  bulla = .calc_bulla_evenness)
    
    FUN(mat = mat, ...)
}

## Richness helper function

.calc_observed <- function(mat, detection, ...){
    # vegan::estimateR(t(mat))["S.obs",]
    colSums(mat > detection)
}

.calc_chao1 <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
        mat <- matrix(mat, nrow = nrow(mat))
    }
    
    ans <- t(vegan::estimateR(t(mat))[c("S.chao1","se.chao1"),])
    colnames(ans) <- c("","se")
    ans
}

.calc_ace <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
        mat <- matrix(mat, nrow = nrow(mat))
    }
    
    ans <- t(vegan::estimateR(t(mat))[c("S.ACE","se.ACE"),])
    colnames(ans) <- c("","se")
    ans
}

.calc_hill <- function(mat, ...){
    # Exponent of Shannon diversity
    exp(vegan::diversity(t(mat), index="shannon"))
}

.get_richness_values <- function(index, mat, detection, ...) {
    
    FUN <- switch(index,
                  observed = .calc_observed,
                  chao1 = .calc_chao1,
                  ace = .calc_ace,
                  hill = .calc_hill
    )
    
    FUN(mat = mat, detection = detection, ...)
    
}
