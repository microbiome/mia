/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <iostream>
#include <vector>

#include <Rcpp.h>

#include "assay.h"
#include "tree.h"
#include "propmap.h"

//' Calculate Faith's PD 
//' 
//' This function calculates Faith's phylogenetic diversity for a given assay
//' and rowTree, using a C++ implementation of the Stacked Faith's Phylogenetic
//' Diversity (SFPhD) algorithm.
//' 
//' @details
//' This function makes several assumptions about the contents of
//' \code{assay} and \code{rowTree}, namely that:
//' \itemize{
//'  \item \code{assay} and \code{rowTree} are both non-empty.
//'  \item \code{assay} has row and column names.
//'  \item \code{rowTree}'s nodes are arranged in cladewise order.
//' }
//' These checks should all be handled in the surrounding R code.
//' 
//' The values returned by this function are equivalent to the values returned
//' by \code{picante::pd()} with the parameter \code{include.root=TRUE}.
//' 
//' The C++ code was adapted from an implementation by the Unifrac team (see
//' \url{https://genome.cshlp.org/content/31/11/2131} or
//' \url{https://github.com/biocore/unifrac}), which is licensed under the BSD
//' 3-Clause license. 
//' 
//' @param assay An R numeric matrix containing the assay of a \code{TreeSE}
//' object.
//' @param rowTree An \code{ape::phylo} object containing the rowTree of a
//' \code{TreeSE} object.
//' @return A vector containing Faith's PD values.
//' 
//' @keywords internal
// [[Rcpp::export(.faith_cpp)]]
Rcpp::NumericVector faith_cpp(const Rcpp::NumericMatrix & assay,
                              const Rcpp::List & rowTree){
    
    su::BPTree tree = su::BPTree(rowTree);      
    su::Assay table = su::Assay(assay);
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(),
                                            table.obs_ids.end());
    
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    
    su::PropMap propmap(table.n_samples);
    
    uint32_t node;
    std::vector<double> node_proportions;
    double length;
    
    std::vector<double> results = std::vector<double>(table.n_samples, 0.0);
    
    // for node in postorderselect
    const unsigned int max_k = (tree_sheared.nparens>1) ? ((tree_sheared.nparens / 2) - 1) : 0;
    for(unsigned int k = 0; k < max_k; k++) {
        node = tree_sheared.postorderselect(k);
        
        // get branch length
        length = tree_sheared.lengths[node];
        
        // get node proportions and set intermediate scores
        node_proportions = set_proportions(tree_sheared, node, table, propmap, false);
        
        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            results[sample] += (node_proportions[sample] > 0) * length;
        }
    }
    
    Rcpp::NumericVector faith = Rcpp::NumericVector(results.size());
    
    for(unsigned int i = 0; i < results.size(); i++){
        faith[i] = results[i];
    }
    
    return faith;
}
