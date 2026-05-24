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

su::mat_t su::one_off(const su::Assay & table,
                      const su::BPTree & tree,
                      bool weighted,
                      bool bypass_tips) {
    
    //Number of stripes to be used, basically half of samples
    const unsigned int stripe_stop = (table.n_samples + 1) / 2;
    
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
    

    su::unifrac(std::ref(table),
                std::ref(tree),
                std::ref(dm_stripes),
                std::ref(dm_stripes_total),
                weighted,
                task);
        
    
    su::mat_t result;
    result.n_samples = table.n_samples;
    result.cf_size = su::comb_2(table.n_samples);
    result.sample_ids = table.sample_ids;
    result.is_upper_triangle = true;
    result.condensed_form =  su::stripes_to_condensed_form(dm_stripes,
                                                            table.n_samples,
                                                            task.start,
                                                            task.stop);
    
    return result;
}




void su::unifrac(const su::Assay &table,
                 const su::BPTree &tree,
                 su::StripeMap & dm_stripes,
                 su::StripeMap & dm_stripes_total,
                 bool weighted,
                 const su::task_parameters task_p)
{
    //unweighted
    if (weighted == false) {
        unifracTT<su::UnifracUnweightedTask>(
            table, tree, true, dm_stripes, dm_stripes_total,
            task_p );
    }
    //weighted unnormalized
    else {
        unifracTT<su::UnifracUnnormalizedWeightedTask>(
            table, tree, false, dm_stripes, dm_stripes_total,
            task_p );
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
    const unsigned int n_samples = task_p.n_samples;
    const uint64_t  n_samples_r = ((n_samples + UNIFRAC_BLOCK-1) /
                                   UNIFRAC_BLOCK)*UNIFRAC_BLOCK; // round up
    
    su::PropMapMulti propmap_multi(table.n_samples);
    
    const unsigned int max_emb =  TaskT::RECOMMENDED_MAX_EMBS;
    
    TaskT taskObj(dm_stripes, dm_stripes_total, max_emb, task_p);
    
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
    
    
    
    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens / 2) - 1;
    
    const unsigned int num_prop_chunks = propmap_multi.get_num_stacks();
    // num_prop_chunks = 1
    
    while (k<max_k)
    {
        const unsigned int k_start = k;
        unsigned int filled_emb = 0;
        
        // chunk the progress to maximize cache reuse
        for (unsigned int ck=0; ck<num_prop_chunks; ck++) {
            
            su::PropMap & propmap = propmap_multi.get_prop_map(ck);
            const unsigned int tstart = propmap_multi.get_start(ck);
            const unsigned int tend = propmap_multi.get_end(ck);

            unsigned int my_filled_emb = 0;
            unsigned int my_k=k_start;
            
            while ((my_filled_emb<max_emb) && (my_k<max_k)) {
                const uint32_t node = tree.postorderselect(my_k);
                my_k++;
                
                //calculate proportions range for given node
                std::vector<double> node_proportions = su::set_proportions_range(
                                                                        tree,
                                                                        node,
                                                                        table,
                                                                        tstart,
                                                                        tend,
                                                                        propmap);
                
                if(task_p.bypass_tips && tree.isleaf(node))
                    continue;
                
                
                if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                    lengths[filled_emb] = tree.lengths[node];
                    filled_emb++;
                }
                
                taskObj.embed_proportions_range(node_proportions,
                                                tstart,
                                                tend,
                                                my_filled_emb);
                my_filled_emb++;
            }
            
            if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                k=my_k;
            }
        }
        
        taskObj._run(filled_emb,lengths);
        filled_emb=0;
        
    }
    
    //want_total is used if you want the results as a percentage of the total?
    if(want_total) {
        const uint64_t start_idx = task_p.start;
        const uint64_t stop_idx = task_p.stop;
        
        for(uint64_t i = start_idx; i < stop_idx; i++){
            for(uint64_t j = 0; j < n_samples; j++) {
                uint64_t idx = ((i-start_idx)*n_samples_r)+j;
                taskObj.dm_stripes.buf[idx] = taskObj.dm_stripes.buf[idx] /
                                            taskObj.dm_stripes_total.buf[idx];
            }
        }
    }
}

std::vector<double> su::stripes_to_condensed_form(su::StripeMap & stripes,
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
            // determine the position in the condensed form vector for a given
            // (i, j)
            // based off of
            // https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
            uint64_t comb_N_minus_i = comb_2(n - i);
            cf[comb_N - comb_N_minus_i + (j - i - 1)] = dm_stripe[k];
        }
    }
    return cf;
}