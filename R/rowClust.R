# REQUIRED PACKAGES: "fastcluster", "stats", "coop"
rowClust <- function(tse,
                     abund_values = "counts",
                     metric = "spearman", # c("pearson", "kendall", "spearman" || 
                                          #    "cosine" ||
                                          #    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
                     method = "ward.D2", # c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid","median")
                     name = "clust_id",
                     k = NULL, # number of clusters
                     h = NULL, # height of cut-off
                     ... # arguments passed to stats::dist (mainly for the p arg of "minkowski")
                     ) {
  # Getting Assay data
  if(abund_values %in% assayNames(tse))
    x <- assay(tse, abund_values)
  else
    stop(paste0("Assay ",abund_values," not Found!"))
  
  # Calculating distance
  if (metric %in% c("pearson", "kendall", "spearman"))
    x <- as.dist(cor(t(x), method = metric))
  else if (metric=="cosine")
    x <- as.dist(1 - coop::cosine(t(x)))
  else if (metric %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    x <- stats::dist(x, method = metric, ...)
  else
    stop("metric not recognized! (look at methods at help(cor) or help(dist))")
  
  # Clustering
  if (method %in% c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid","median"))
    hc <- fastcluster::hclust(x, method = method)
  else
    stop("Clustering method unknown! (look at help(hclust))")
  
  if (is.null(k) & is.null(h))
    stop("Specify a cut-off point; k(number of clusters) | h(height)")
  else
    clusters <- stats::cutree(hc, k=k, h=h)
    rowData(tse)[,name] <- as.factor(clusters)
    return(tse)
}
