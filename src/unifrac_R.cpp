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

#include <chrono>
/*
#include <chrono>
 auto start = std::chrono::high_resolution_clock::now();
 auto stop = std::chrono::high_resolution_clock::now();
 auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 Rcpp::Rcout << "Main thread: " << duration.count() << "\n";
 
 
 start = std::chrono::high_resolution_clock::now();
 stop = std::chrono::high_resolution_clock::now();
 duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
 Rcpp::Rcout << "Condensed form: " << duration.count() << "\n";
 */

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
                              const Rcpp::List & rowTree,
                              bool weighted,
                              bool normalized,
                              bool bypass_tips){
    
    // Normalized matches the results given by weighted
    
    auto start = std::chrono::high_resolution_clock::now();
    
    su::BPTree tree = su::BPTree(rowTree);
    su::Assay table = su::Assay(assay);
    std::string method = "unweighted";
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(),
                                            table.obs_ids.end());
    
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    
    su::mat_t results = su::one_off(table, tree_sheared, weighted, normalized, bypass_tips);
    
    //condensed_form is the main values, returned in result
    //Sample_ids can be handled with a map?
    //n_samples, cf_size, is_upper_triangle are single values that can be passed in some other way?
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    Rcpp::Rcout << "Main thread: " << duration.count() << "\n";
    
    return Rcpp::List::create(Rcpp::Named("n_samples") = results.n_samples,
                              Rcpp::Named("is_upper_triangle") = results.is_upper_triangle,
                              Rcpp::Named("cf_size") = results.cf_size,
                              Rcpp::Named("c_form") = results.condensed_form);
}

