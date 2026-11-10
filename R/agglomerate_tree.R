.agglomerate_tree_into_mrca <- function(df, node_num, groups){
    # Looking up parents from a vector is much faster than searching them from
    # the table.
    parent <- .get_parents(df)
    # For each group, find most recent common ancestor
    mrcas <- c()
    group_unique <- groups |> unique()
    for( group in group_unique ){
        nodes <- node_num[ groups %in% group ]
        mrcas <- c(mrcas, .find_mrca_for_group(parent, nodes))
    }
    # Now we know for each group their MRCA. Now the task is to determine
    # which internal nodes are kept. That is why we have to find full paths
    # for these MRCAs.
    all_nodes <- c()
    for( mrca in mrcas ){
        full_path <- .find_full_path_for_node(parent, mrca)
        all_nodes <- c(all_nodes, full_path)
    }
    # Agglomerate the tree, i.e., select only MRCAs and their internal nodes
    new_df <- df[df[["node"]] %in% all_nodes, ]

    # Add to data.frame information on the linkage between nodes and groups
    new_df[["link_to_row"]] <- group_unique[ match(new_df[["node"]], mrcas) ]

    return(new_df)
}

# This function gets group of nodes, and it returns their MRCA
.find_mrca_for_group <- function(parent, nodes){
    # For each node, find their all parents
    paths <- vector("list", length(nodes))
    for( i in seq_along(nodes) ){
        paths[[i]] <- .find_full_path_for_node(parent, nodes[[i]])
    }
    # Now we have full paths for each node. We can find the MRCA by selecting
    # the most recent node that all paths share.
    mrca <- .mrca_from_paths(paths)
    return(mrca)
}

# This function gets single node as input, and it returns the node and all its
# parents up to the root. parent is a vector that gives the parent of each
# node; the parent of the root is the root itself.
.find_full_path_for_node <- function(parent, node){
    all_nodes <- node
    # Go one branch up until the root is reached
    while( parent[[node]] != node ){
        node <- parent[[node]]
        all_nodes <- c(all_nodes, node)
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
            # If the list of shared nodes is empty, i.e., no shared nodes,
            # we can stop already here.
            if( length(shared) == 0L ) {
                break
            }
        }
    }
    return(shared)
}

# This function checks that groups are monophyletic, i.e., each group forms
# its own subtree that does not include nodes from other groups. The method
# that agglomerates the tree into MRCAs assumes this.
.check_monophyly <- function(df, nodes, groups){
    not_monophyletic <- .get_non_monophyletic_groups(df, nodes, groups)
    is_monophyletic <- length(not_monophyletic) == 0L
    return(is_monophyletic)
}

# This function returns the groups that are not monophyletic. Nodes that are
# not linked to any group are ignored. The input is the tree in table format
# (.tree_to_table()).
#
# The tree is traversed once from tips to root. Each node gets the group of its
# descendants, or it is marked as mixed if its descendants belong to multiple
# groups. A group is monophyletic if it has exactly one topmost node, i.e., a
# node whose parent is mixed or that is the root.
.get_non_monophyletic_groups <- function(df, nodes, groups){
    parent <- .get_parents(df)
    n_nodes <- length(parent)
    is_root <- parent == seq_len(n_nodes)
    # Number of branches from the root to each node
    level <- .get_levels_of_nodes(parent)
    # Group of each node. NA means that the node is not linked to any group.
    # A node is mixed if its descendants (or the node itself) belong to
    # multiple groups.
    node_group <- rep(NA_character_, n_nodes)
    mixed <- rep(FALSE, n_nodes)
    for( i in seq_along(nodes) ){
        node <- nodes[[i]]
        # Features from different groups can be linked to the same node
        mixed[node] <- mixed[node] || isTRUE(node_group[node] != groups[[i]])
        node_group[node] <- groups[[i]]
    }
    # Go through the nodes from tips to root, and pass the group to the parent.
    # order() gives node numbers sorted by level, the deepest nodes first, so
    # children are always handled before their parents. The root has no
    # parent, so it is skipped.
    ord <- order(level, decreasing = TRUE)
    ord <- ord[ !is_root[ord] ]
    for( i in ord ){
        p <- parent[i]
        # The parent is mixed if the child is mixed or if they have different
        # groups
        mixed[p] <- mixed[p] || mixed[i] ||
            isTRUE(node_group[p] != node_group[i])
        # The parent gets the group of the child if it does not have one yet
        if( is.na(node_group[p]) ){
            node_group[p] <- node_group[i]
        }
    }
    # Find topmost node of each group. Each group must have exactly one. If it
    # has none, the group's nodes are all mixed with other groups.
    is_top <- !mixed & !is.na(node_group) & (is_root | mixed[parent])
    n_top <- factor(node_group[is_top], levels = unique(groups)) |> table()
    not_monophyletic <- names(n_top)[n_top != 1L]
    return(not_monophyletic)
}

# This function returns the parent of each node. Node numbers go from 1 to n,
# so they can be used as indices. The parent of the root is the root itself.
.get_parents <- function(df){
    parent <- seq_len(max(df[["node"]]))
    parent[ df[["node"]] ] <- df[["parent"]]
    return(parent)
}

# This function gets the parent of each node as input. The parent of the root
# is the root itself. It returns the level of each node, i.e., the number of
# branches from the root (root = 0). Unlike depth, the level does not depend on
# branch lengths, so children always have a higher level than their parents.
.get_levels_of_nodes <- function(parent){
    is_root <- parent == seq_along(parent)
    if( sum(is_root) != 1L ){
        stop("The tree must have a single root.", call. = FALSE)
    }
    level <- rep(NA_integer_, length(parent))
    level[is_root] <- 0L
    # Go one level down at a time until all nodes have a level
    while( anyNA(level) ){
        ind <- is.na(level) & !is.na(level[parent])
        if( !any(ind) ){
            stop("All nodes must be connected to the root.", call. = FALSE)
        }
        level[ind] <- level[ parent[ind] ] + 1L
    }
    return(level)
}
