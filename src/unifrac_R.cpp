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
#include "unifrac.hpp"
#include "unifrac_task.hpp"

// Calculate Unifrac
//
// @keywords internal
// [[Rcpp::export(.unifrac_cpp)]]
Rcpp::NumericVector unifrac_cpp(const Rcpp::NumericMatrix & assay,
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
    
    su::mat_t results = one_off(table, tree, method, false, 1.0, false);
    
    //condensed_form is the main values, returned in result
    //Sample_ids can be handled with a map?
    //n_samples, cf_size, is_upper_triangle are single values that can be passed in some other way?
    
    /*
    Rcpp::NumericVector unifrac = Rcpp::NumericVector(results.size());
    
    for( unsigned int i = 0; i < results.cf_size; i++ ){
        unifrac[i] = results.condensed_form[i];
    }
    */
    
    return Rcpp::List::create(Rcpp::Named("n_samples") = results.n_samples,
                              Rcpp::Named("is_upper_triangle") = results.is_upper_triangle,
                              Rcpp::Named("cf_size") = results.cf_size,
                              Rcpp::Named("c_form") = results.condensed_form);
}






/* Compute UniFrac - condensed form
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads to use.
 * result <mat_t**> the resulting distance matrix in condensed form, this is initialized within the method so using **
 *
 * one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */

su::mat_t one_off(const su::Assay & table,
                       const su::BPTree & tree,
                       std::string unifrac_method,
                       double alpha,
                       bool variance_adjust,
                       bool bypass_tips) {
    
    //Check that method is valid - pass it as something other than string?
    //SET_METHOD(unifrac_method, unknown_method)
    
    //Stripes relate to matrix calculations
    const unsigned int stripe_stop = (table.n_samples + 1) / 2;
    
    //Originally std::vector of double pointers - this is where the data travels?
    su::StripeMap dm_stripes(table.n_samples);
    su::StripeMap dm_stripes_total(table.n_samples);
    
    su::task_parameters task;
    
    task.tid = 0;
    task.start = 0;
    task.stop = stripe_stop;
    task.bypass_tips = bypass_tips;
    task.n_samples = n_samples;
    task.g_unifrac_alpha = alpha;
    
    //Main action
    //Calls either unifrac or _vaw depending on variance_adjust
    //makes use of std::ref?
    //Versions for accelerated and cpu - let's go with cpu for now
    //method is "unweighted" by default, let's start with that and see what else may be needed
    
    //Threads are passed here after being created with a vector
    //This does nothing except pass the number of threads, however
    su::process_stripes(table, tree_sheared, method, variance_adjust,
                        dm_stripes, dm_stripes_total, task);
    

    //Only use of threading in this version of code was for stripes to condensed form
    //Basically each thread calls stripes_to_condensed_form
    //Which is just a bunch of binomial calculations
    
    su::mat_t result;
    result->n_samples = table.n_samples;
    result->cf_size = su::comb_2(table.n_samples);
    result->sample_ids = table.sample_ids;
    result->condensed_form = std::vector<double>(su::comb_2(table.n_samples),
                                                 0.0);
    result->is_upper_triangle = true;
    
    return su::stripes_to_condensed_form(dm_stripes,
                                   table.n_samples,
                                   result,
                                   task.start,
                                   task.stop);
}
        
