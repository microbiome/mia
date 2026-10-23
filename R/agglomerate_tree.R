.agglomerate_tree_into_mrca <- function(df, node_num, groups){
    # For each group, find most recent common ancestor
    mrcas <- c()
    group_unique <- groups |> unique()
    for( group in group_unique ){
        nodes <- node_num[ groups %in% group ]
        mrcas <- c(mrcas, .find_mrca_for_group(df, nodes))
    }
    # Now we know for each group their MRCA. Now the task is to determine
    # which internal nodes are kept. That is why we have to find full paths
    # for these MRCAs.
    all_nodes <- c()
    for( mrca in mrcas ){
        full_path <- .find_full_path_for_node(df, mrca)
        all_nodes <- c(all_nodes, full_path)
    }
    # Agglomerate the tree, i.e., select only MRCAs and their internal nodes
    new_df <- df[df[["node"]] %in% all_nodes, ]

    # Add to data.frame information on the linkage between nodes and groups
    new_df[["link_to_row"]] <- group_unique[ match(new_df[["node"]], mrcas) ]

    return(new_df)
}

# This function gets group of nodes, and it returns their MRCA
.find_mrca_for_group <- function(df, nodes){
    # For each node, find their all parents
    paths <- list()
    for( i in nodes |> length() |> seq_len() ){
        paths[[i]] <- .find_full_path_for_node(df, nodes[[i]])
    }
    # Now we have full paths for each node. We can find the MRCA by selecting
    # the the most recent node that all paths share.
    mrca <- .mrca_from_paths(paths)
    return(mrca)
}

# This function gets single node as input, and it returns all parents of it.
.find_full_path_for_node <- function(df, node){
    all_nodes <- c(node)
    while( TRUE ){
        # Get info on current branch
        temp <- df[ df[["node"]] == node, ]
        # Add the parent to list
        all_nodes <- c(all_nodes, temp[["parent"]])
        # If parent node is same as current node or it is NA, it means that we
        # are now in root. End here.
        if( is.na(node) || temp[["parent"]] == node ){
            break
        }
        # Go one branch up
        node <- temp[["parent"]]
    }
    return(all_nodes)
}

# This function gets as input, multiple paths. It finds the MRCA from them.
.mrca_from_paths <- function(path_list) {
    # Determine which nodes are shared
    shared <- .shared_nodes(path_list)
    mrca <- NA_integer_
    # If shared nodes were found
    if( length(shared) > 0L ){
        # First element is the mrca
        mrca <- shared[[1L]]
    }
    return(mrca)
}

# This function returns which nodes are shared by paths.
.shared_nodes <- function(path_list) {
    # Get nodes from first path as "shared"
    shared <- path_list[[1]]
    # If the group contained more than one path
    if( length(path_list) > 1L ){
        # Go through one by one each path.
        for (i in 2:length(path_list)) {
            # Check which nodes are shared between current shared nodes and
            # this path.
            shared <- shared[shared %in% path_list[[i]]]
            # If the list of shared nods is empty, i.e., no shared nodes,
            # we can stop already here.
            if( length(shared) == 0L ) {
                break
            }
        }
    }
    return(shared)
}
