/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "unifrac.h"
#include "propmap.h"
#include "stripemap.h"
#include "tree.h"

/*
#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>
*/




su::Method su::set_method(std::string requested_method) {
    if(requested_method == "unweighted")                                                               
        return unweighted;                                                                                           
    else if(requested_method == "weighted_normalized")                                                 
        return weighted_normalized;                                                                                  
    else if(requested_method == "weighted_unnormalized")                                               
        return weighted_unnormalized;                                                                                
    else if(requested_method == "generalized")                                                         
        return generalized;                                                                                          
    /*else if(std::strcmp(requested_method, "unweighted_fp32") == 0)                                                     
     method = unweighted_fp32;                                                                                      
     else if(std::strcmp(requested_method, "weighted_normalized_fp32") == 0)                                            
     method = weighted_normalized_fp32;                                                                             
     else if(std::strcmp(requested_method, "weighted_unnormalized_fp32") == 0)                                          
     method = weighted_unnormalized_fp32;                                                                           
     else if(std::strcmp(requested_method, "generalized_fp32") == 0)                                                    
     method = generalized_fp32;   */                                                                                  
    else {                                                                                                             
        return unknown;                                                                                                
    }             
}                                    



su::mat_t su::one_off(const su::Assay & table,
                      const su::BPTree & tree,
                      std::string unifrac_method,
                      double alpha,
                      bool variance_adjust,
                      bool bypass_tips) {
    
    Rcpp::Rcout << "Start one_off\n";
    
    //Check that method is valid - pass it as something other than string?
    su::Method method = set_method(unifrac_method);
    
    //Number of stripes to be used, basically half of samples
    const unsigned int stripe_stop = (table.n_samples + 1) / 2;
    
    //Originally std::vector of double pointers - this is where the data travels?
    su::StripeMap dm_stripes(table.n_samples);
    su::StripeMap dm_stripes_total(table.n_samples);
    
    su::task_parameters task;
    
    //thread id - currently single-threaded
    task.tid = 0;
    //Stripes to start and stop on - single task, so the entire thing
    task.start = 0;
    task.stop = stripe_stop;
    
    task.bypass_tips = bypass_tips;
    task.n_samples = table.n_samples;
    task.g_unifrac_alpha = alpha;
    
    //Main action
    //Calls either unifrac or _vaw depending on variance_adjust
    //makes use of std::ref?
    //Versions for accelerated and cpu - let's go with cpu for now
    //method is "unweighted" by default, let's start with that and see what else may be needed
    
    
    //This could potentially be threaded
    //Wasn't in the code because doesn't work with openacc/openmp?
    
    su::unifrac(std::ref(table),
                std::ref(tree),
                method,
                std::ref(dm_stripes),
                std::ref(dm_stripes_total),
                task,
                variance_adjust);
    
    Rcpp::Rcout << "unifrac done\n";
    
    //Only use of std::thread in this version of code was for stripes to condensed form
    //Basically each thread calls stripes_to_condensed_form
    //Which is just a bunch of binomial calculations
    
    su::mat_t result;
    result.n_samples = table.n_samples;
    result.cf_size = su::comb_2(table.n_samples);
    result.sample_ids = table.sample_ids;
    result.is_upper_triangle = true;
    result.condensed_form =  su::stripes_to_condensed_form(dm_stripes,
                                                            table.n_samples,
                                                            task.start,
                                                            task.stop);
    
    Rcpp::Rcout << "one_off done\n";
    
    return result;
}




void su::unifrac(const su::Assay &table,
                 const su::BPTree &tree,
                 su::Method unifrac_method,
                 su::StripeMap & dm_stripes,
                 su::StripeMap & dm_stripes_total,
                 const su::task_parameters task_p,
                 bool variance_adjust)
{
    
    Rcpp::Rcout << "Start unifrac\n";
    
    if(variance_adjust)
    {
        /*
         switch(unifrac_method) {
    case su::unweighted:
         unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<double>,double>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    case su::weighted_normalized:
         unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<double>,double>(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    case su::weighted_unnormalized:
         unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<double>,double>( table, tree, false, dm_stripes,dm_stripes_total,task_p);
         break;
    case su::generalized:
         unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<double>,double>(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    case su::unweighted_fp32:
         unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<float >,float >(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    case su::weighted_normalized_fp32:
         unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<float >,float >(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    case su::weighted_unnormalized_fp32:
         unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<float >,float >( table, tree, false, dm_stripes,dm_stripes_total,task_p);
         break;
    case su::generalized_fp32:
         unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<float >,float >(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
         break;
    default:
         fprintf(stderr, "Unknown unifrac task\n");
         exit(1);
         break;
         }
         */
    }
    else
    {
        switch(unifrac_method)
        {
        case su::unweighted:
            unifracTT<su::UnifracUnweightedTask>(
                table, tree, true,  dm_stripes,dm_stripes_total,
                task_p );
            break;
            /*case su::weighted_normalized:
             unifracTT<su::UnifracNormalizedWeightedTask<double>,double>(
             table, tree, true,  dm_stripes,dm_stripes_total,
             task_p );
             break;
        case su::weighted_unnormalized:
             unifracTT<su::UnifracUnnormalizedWeightedTask<double>,
             double>(table, tree, false, dm_stripes,
             dm_stripes_total, task_p );
             break;
        case su::generalized:
             unifracTT<su::UnifracGeneralizedTask<double>,double>(
             table, tree, true,  dm_stripes,dm_stripes_total,
             task_p );
             break;
             */
        default:
            fprintf(stderr, "Unknown unifrac task\n");
        exit(1);
        break;
        }
    }
}





template<class TaskT>
inline void su::unifracTT(const su::Assay & table,
                      const su::BPTree & tree,
                      const bool want_total,
                      su::StripeMap & dm_stripes,
                      su::StripeMap & dm_stripes_total,
                      const su::task_parameters & task_p)
{
    
    Rcpp::Rcout << "Start unifracTT\n";
    
    if(table.n_samples != task_p.n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }
    
    const unsigned int n_samples = task_p.n_samples;
    const uint64_t  n_samples_r = ((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK; // round up
    
    
    //su::PropStackMulti<TFloat> propstack_multi(table.n_samples);
    su::PropMap propmap(table.n_samples);
    
    const unsigned int max_emb =  TaskT::RECOMMENDED_MAX_EMBS;
    
    
    
    Rcpp::Rcout << "Start taskObj\n";
    
    TaskT taskObj(dm_stripes, dm_stripes_total, max_emb, task_p);
    
    Rcpp::Rcout << "taskObj done\n";
    
    std::vector<double> lengths = std::vector<double>(max_emb);
    
    /*
     * The values in the example vectors correspond to index positions of an
     * element in the resulting distance matrix. So, in the example below,
     * the following can be interpreted:
     *
     * [0 1 2]
     * [1 2 3]
     *
     * As comparing the sample for row 0 against the sample for col 1, the
     * sample for row 1 against the sample for col 2, the sample for row 2
     * against the sample for col 3.
     *
     * In other words, we're computing stripes of a distance matrix. In the
     * following example, we're computing over 6 samples requiring 3
     * stripes.
     *
     * A; stripe == 0
     * [0 1 2 3 4 5]
     * [1 2 3 4 5 0]
     *
     * B; stripe == 1
     * [0 1 2 3 4 5]
     * [2 3 4 5 0 1]
     *
     * C; stripe == 2
     * [0 1 2 3 4 5]
     * [3 4 5 0 1 2]
     *
     * The stripes end up computing the following positions in the distance
     * matrix.
     *
     * x A B C x x
     * x x A B C x
     * x x x A B C
     * C x x x A B
     * B C x x x A
     * A B C x x x
     *
     * However, we store those stripes as vectors, ie
     * [ A A A A A A ]
     *
     * We end up performing N / 2 redundant calculations on the last stripe
     * (see C) but that is small over large N.
     */
    
    Rcpp::Rcout << "Start calcs\n";
    
    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens / 2) - 1;
    
    // num_prop_chunks = 1
    while (k<max_k)
    {
        const unsigned int k_start = k;
        unsigned int filled_emb = 0;
        
        // ck = 0
        // chunk the progress to maximize cache reuse
        const unsigned int tstart = 0;
        const unsigned int tend = 0; // end of propstack?
        unsigned int my_filled_emb = 0;
        unsigned int my_k=k_start;
        
        while ((my_filled_emb<max_emb) && (my_k<max_k)) {
            const uint32_t node = tree.postorderselect(my_k);
            my_k++;
            
            //TFloat *node_proportions = propstack.pop(node);
            //su::set_proportions_range(node_proportions, tree, node, table, tstart, tend, propstack);
            
            //calculate proportions range for given node
            su::set_proportions_range(tree, node, table, tstart, tend, propmap);
            
            //propstack pop ERASES any existing vector for node and gives a blank one
                //creates memory leaks if node isn't pushed before popping? 
            //get just returns the given vector
            //push removes the vector from use
            //Any time a propstack vector is modified, remember to do a propmap update
            
            if(task_p.bypass_tips && tree.isleaf(node))
                continue;
            
            lengths[filled_emb] = tree.lengths[node];
            filled_emb++;
            
            //store the proportions inside the taskobject's continuous buffer
            //Shouldn't modify node_proportions
            std::vector<double> node_proportions = propmap.get(node);
            
            //Rcpp::Rcout << "start embed_proportions_range\n";
            taskObj.embed_proportions_range(node_proportions, tstart, tend, my_filled_emb);
            //Rcpp::Rcout << "embed_proportions_range done\n";
            my_filled_emb++;
        }
         
        k=my_k;
        
        //This is used to keep track of filled embeds over different threads?
        //Does nothing without openacc
        //taskObj.sync_embedded_proportions(filled_emb);
        
        Rcpp::Rcout << "start taskObj._run\n";
        taskObj._run(filled_emb,lengths);
        Rcpp::Rcout << "taskObj._run done\n";
        
        filled_emb=0;
    }
    
    
    Rcpp::Rcout << "calcs done\n";
    
    
    //I suppose want_total is used if you want the results as a percentage of the total?
    if(want_total) {
        const uint64_t start_idx = task_p.start;
        const uint64_t stop_idx = task_p.stop;
        
        for(uint64_t i = start_idx; i < stop_idx; i++){
            /*
             std::vector<double> dm_stripes_buf = std::vector<double>  ;
             std::vector<double> dm_stripes_total_buf = taskObj.dm_stripes_total.get(idx);
            std::copy(std::begin(taskObj.dm_stripes.buf),
                      std::end(taskObj.dm_stripes.buf),
                      std::begin(dm_stripes_buf) + (emb8<<8));
            */
            
            std::vector<double> dm_stripes_buf = taskObj.dm_stripes.buf;
            std::vector<double> dm_stripes_total_buf = taskObj.dm_stripes_total.buf;
            
            for(uint64_t j = 0; j < n_samples; j++) {
                uint64_t idx = (i-start_idx)*n_samples_r+j;
                dm_stripes_buf[idx] = dm_stripes_buf[idx]/dm_stripes_total_buf[idx];
            }
            
            taskObj.dm_stripes.buf = dm_stripes_buf;
            
            /*
             taskObj.dm_stripes.update(idx, dm_stripes_buf);
            std::copy(std::begin(dm_stripes_buf),
                      std::end(dm_stripes_buf),
                      std::begin(taskObj.dm_stripes.buf) + );
             */
        }
    }
    
    Rcpp::Rcout << "unifracTT done\n";
}



std::vector<double> su::set_proportions_range(const su::BPTree & tree,
                                              uint32_t node,
                                              const su::Assay & table,
                                              unsigned int start,
                                              unsigned int end,
                                              PropMap & pm,
                                              bool normalize) {
    const unsigned int els = end-start;
    std::vector<double> props = std::vector(els, 0.0);
    if(tree.isleaf(node)) {
        props = table.get_obs_data_range(tree.names[node], start, end, normalize);
    } else {
        const unsigned int right = tree.rightchild(node);
        unsigned int current = tree.leftchild(node);
        
        while(current <= right && current != 0) {
            std::vector<double> vec = pm.get(current);  // pull from prop map
            pm.clear(current);  // remove from prop map, place back on stack
            
            for(unsigned int i = 0; i < els; i++)
                props[i] += vec[i];
            
            current = tree.rightsibling(current);
        }
    }
    pm.update(node, props);
    return props;
}





std::vector<double> su::stripes_to_condensed_form(su::StripeMap stripes,
                                   uint32_t n,
                                   unsigned int start,
                                   unsigned int stop) {
    // n must be >= 2, but that should be enforced upstream as that would imply
    // computing unifrac on a single sample.
    
    uint64_t comb_N = comb_2(n);
    std::vector<double> cf = std::vector<double>(comb_N, 0.0);
    
    for(unsigned int stripe = start; stripe < stop; stripe++) {
        //Does stripemap contain all the stripes or just one thread's stripes?
        std::vector dm_stripe = stripes.get(stripe);
        // compute the (i, j) position of each element in each stripe
        uint64_t i = 0;
        uint64_t j = stripe + 1;
        for(uint64_t k = 0; k < n; k++, i++, j++) {
            if(j == n) {
                i = 0;
                j = n - (stripe + 1);
            }
            // determine the position in the condensed form vector for a given (i, j)
            // based off of
            // https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
            uint64_t comb_N_minus_i = comb_2(n - i);
            cf[comb_N - comb_N_minus_i + (j - i - 1)] = dm_stripe[k];
        }
    }
    return cf;
}


