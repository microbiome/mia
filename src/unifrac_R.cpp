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
#include "stripemap.h"

#include "unifrac.h"


// Calculate Unifrac
//
// @keywords internal
// [[Rcpp::export(.unifrac_cpp)]]
Rcpp::List unifrac_cpp(const Rcpp::NumericMatrix & assay,
                              const Rcpp::List & rowTree){
    // 
    // std::unordered_set<std::string> to_keep(table.obs_ids.begin(),
    //                                         table.obs_ids.end());
    // 
    // su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    // 
    // su::PropMap propmap(table.n_samples);
    // 
    // uint32_t node;
    // std::vector<double> node_proportions;
    // double length;
    // 
    // std::vector<double> results = std::vector<double>(table.n_samples, 0.0);
    // 
    // 
    // // For node in postorderselect
    // const unsigned int max_k = (tree_sheared.nparens>1) ?
    // ((tree_sheared.nparens / 2) - 1) : 0;
    // 
    // for( unsigned int k = 0; k < max_k; k++ ){
    //     node = tree_sheared.postorderselect(k);
    //     
    //     // Get branch length
    //     length = tree_sheared.lengths[node];
    //     
    //     // Get node proportions and set intermediate scores
    //     node_proportions = set_proportions(tree_sheared, node, table, propmap,
    //                                        false);
    //     
    //     for( unsigned int sample = 0; sample < table.n_samples; sample++ ){
    //         // Calculate contribution of node to score
    //         results[sample] += (node_proportions[sample] > 0) * length;
    //     }
    // } 
    
    su::BPTree tree = su::BPTree(rowTree);
    su::Assay table = su::Assay(assay);
    std::string method = "unweighted";
    
    Rcpp::Rcout << "Start\n";
    
    su::mat_t results = su::one_off(table, tree, method, 1.0, false, false);
    
    Rcpp::Rcout << "All done\n";
    
    //condensed_form is the main values, returned in result
    //Sample_ids can be handled with a map?
    //n_samples, cf_size, is_upper_triangle are single values that can be passed in some other way?
    
    
    //Rcpp::NumericVector cf = Rcpp::NumericVector(results.cf_size);
    
    return Rcpp::List::create(Rcpp::Named("n_samples") = results.n_samples,
                              Rcpp::Named("is_upper_triangle") = results.is_upper_triangle,
                              Rcpp::Named("cf_size") = results.cf_size,
                              Rcpp::Named("c_form") = results.condensed_form);
}

