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
// This function calculates Unifrac distances for a given assay and rowTree,
// using a C++ implementation of the Striped Unifrac algorithm.
//
// @details
// This function makes several assumptions about the contents of
// \code{assay} and \code{rowTree}, namely that:
// \itemize{
//  \item \code{assay} and \code{rowTree} are both non-empty.
//  \item \code{assay} has row and column names.
//  \item \code{rowTree}'s nodes are arranged in cladewise order.
// }
// These checks should all be handled in the surrounding R code.
//
// The C++ code was adapted from an implementation by the Unifrac team
// (Armstrong et al. 2021), which is licensed under the BSD 3-Clause license.
//
// @param assay An R numeric matrix containing the assay of a \code{TreeSE}
// object.
// @param rowTree An \code{ape::phylo} object containing the rowTree of a
// \code{TreeSE} object.
// @param weighted Boolean: Whether to calculate unweighted or weighted Unifrac.
// @param bypass_tips Boolean: Whether to bypass tips during calculations. This
// speeds up calculations considerably, and does not seem to have a noticeable
// effect on the results.
// @return A vector containing Unifrac distances.
//
// @keywords internal
// [[Rcpp::export(.unifrac_cpp)]]
Rcpp::NumericVector unifrac_cpp(const Rcpp::NumericMatrix & assay,
                              const Rcpp::List & rowTree,
                              bool weighted,
                              bool bypass_tips){
    
    su::BPTree tree = su::BPTree(rowTree);
    su::Assay table = su::Assay(assay);
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(),
                                            table.obs_ids.end());
    
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    
    su::mat_t results = su::one_off(table, tree_sheared, weighted, bypass_tips);
    
    unsigned int n = results.condensed_form.size();
    Rcpp::NumericVector unifrac = Rcpp::NumericVector(n);
    
    for( unsigned int i = 0; i < n; i++ ){
        unifrac[i] = results.condensed_form[i];
    }
    
    return unifrac;
}

