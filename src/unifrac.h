/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC
#define __UNIFRAC 1

#include <stack>
#include <string>
#include <vector>
#include <unordered_map>

#include "assay.h"
#include "tree.h"
#include "propmap.h"
#include "unifrac_task.h"

    namespace su {
        
        typedef struct mat {
            unsigned int n_samples;
            unsigned int cf_size;
            bool is_upper_triangle;
            std::vector<double> condensed_form;
            std::vector<std::string> sample_ids;
        } mat_t;
        
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
                              bool weighted,
                              bool bypass_tips);
        
        // Chooses the right task for the job and constructs a unifracTT
        void unifrac(const su::Assay &table,
                     const su::BPTree &tree,
                     su::StripeMap & dm_stripes,
                     su::StripeMap & dm_stripes_total,
                     bool weighted,
                     const su::task_parameters task_p);
        
        // Works the vectors
        template<class TaskT>
        inline void unifracTT(const su::Assay & table,
                              const su::BPTree & tree,
                              const bool want_total,
                              su::StripeMap & dm_stripes,
                              su::StripeMap & dm_stripes_total,
                              const su::task_parameters & task_p);
        
        inline uint64_t comb_2(uint64_t N) {
            // based off of _comb_int_long
            // https://github.com/scipy/scipy/blob/v0.19.1/scipy/special/_comb.pyx
            
            // Compute binom(N, k) for integers.
            //
            // we're disregarding overflow as that practically should not
            // happen unless the number of samples processed is in excess
            // of 4 billion 
            uint64_t val, j, M, nterms;
            uint64_t k = 2;
            
            M = N + 1;
            nterms = k < (N - k) ? k : N - k;
            
            val = 1;
            
            for(j = 1; j < nterms + 1; j++) {
                val *= M - j;
                val /= j;
            }
            return val;
        }
        
        // Stripes to condensed form for the results
        std::vector<double> stripes_to_condensed_form(su::StripeMap & stripes,
                                            uint32_t n,
                                            unsigned int start,
                                            unsigned int stop);
    }
    
#endif
