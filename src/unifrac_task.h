/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC_TASK_H
#define __UNIFRAC_TASK_H 1

#include "stripemap.h"

// CPUs don't need such a big alignment
#define UNIFRAC_BLOCK 16

namespace su {

/* task specific compute parameters
*
* n_samples <int> the number of samples being processed
* start <uint> the first stripe to process
* stop <uint> the last stripe to process
* tid <uint> the thread identifier
* bypass_tips <bool> ignore tips on compute, reduces compute by ~50%
* g_unifrac_alpha <double> an alpha value for generalized unifrac
*/

struct task_parameters {
    uint32_t n_samples;          // number of samples
    unsigned int start;          // starting stripe
    unsigned int stop;           // stopping stripe
    unsigned int tid;            // thread ID
    bool bypass_tips;    // avoid compute at tips 
};

// Helper class that manages stripes
class UnifracTaskVector {
    public:
        su::StripeMap & dm_stripes;
        const unsigned int start_idx;
        const unsigned int n_samples;
        const uint64_t  n_samples_r;
        std::vector<double> buf;
        
        UnifracTaskVector(su::StripeMap & _dm_stripes,
                            const su::task_parameters _task_p);
        
        //Destructor copies the buffer values back into dm_stripes
        ~UnifracTaskVector();

    private:
        const su::task_parameters task_p;
};

// Base task class to be shared by all tasks
// Templated to allow proportions to be embedded as either doubles (weighted) or
// packed bools (unweighted)
template<class TEmb>
class UnifracTaskBase {
    public:
        UnifracTaskVector dm_stripes;
        UnifracTaskVector dm_stripes_total;
        
        su::task_parameters task_p;
        
        const unsigned int max_embs;
        //Continuous vector - each stripe has n_samples_r elements
        //Has at most max_embs stripes - when filled, results stored by the
        //task's _run() function and embeds are cleared for the next batch
        std::vector<TEmb> embedded_proportions; 
        
        UnifracTaskBase(su::StripeMap & _dm_stripes,
                        su::StripeMap & _dm_stripes_total,
                        unsigned int _max_embs,
                        su::task_parameters _task_p);
        
        virtual ~UnifracTaskBase();
        
        // Templated function used when initializing embeds
        static unsigned int get_emb_els(unsigned int max_embs);
        
        static std::vector<TEmb> initialize_embedded(
                                                const uint64_t  n_samples_r,
                                                unsigned int max_embs);
        
        // Store proportions from in into embedded_proportions
        void embed_proportions(const std::vector<double> & in,
                                unsigned int emb);
        
        void embed_proportions_range(const std::vector<double> & in,
                                        unsigned int start,
                                        unsigned int end,
                                        unsigned int emb);
        
    protected:
        void embed_proportions_range_straight(std::vector<double> & out,
                                                const std::vector<double> & in,
                                                unsigned int start,
                                                unsigned int end,
                                                unsigned int emb) const {
            const unsigned int n_samples  = dm_stripes.n_samples;
            const uint64_t n_samples_r  = dm_stripes.n_samples_r;
            const uint64_t offset = emb * n_samples_r;
            
            //Copy to stripe indicated by emb
            //Stripes are all contained in in/out in one mass
            //Start/end aren't necessarily the whole stripe?
            for(unsigned int i = start; i < end; i++){
                out[offset + i] = in[i-start];
            }
            
            if (end==n_samples){
                // avoid NaNs
                for(unsigned int i = n_samples; i < n_samples_r; i++){
                    out[offset + i] = 0.0;
                }
            }
        }
        
        // packed bool
        // Compute (in[:]>0) on each element, and store only the boolean bit.
        // The output values are stored in a multi-byte format, one bit per emb
        // index, so it will likely take multiple passes to store all the values
        // Note: assumes we are processing emb in increasing order, starting
        // from 0
        std::vector<uint64_t> embed_proportions_range_bool(
                                                std::vector<uint64_t>  out,
                                                const std::vector<double> & in,
                                                unsigned int start,
                                                unsigned int end,
                                                unsigned int emb) const {
            const unsigned int n_packed = sizeof(uint64_t)*8;
            const unsigned int n_samples  = dm_stripes.n_samples;
            const uint64_t n_samples_r  = dm_stripes.n_samples_r;
            
            // The output values are stored in a multi-byte format, one bit per
            // emb index
            // Compute the element to store the bit into, as well as which bit
            // in that element 
            unsigned int emb_block = emb/n_packed; // beginning of block
            unsigned int emb_bit = emb%n_packed;   // bit inside the elements
            const uint64_t offset = emb_block * n_samples_r;
            
            if (emb_bit == 0){
                // assign for emb_bit==0, so it clears the other bits
                // assumes we processing emb in increasing order starting from 0
                for(unsigned int i = start; i < end; i++){            
                    out[offset + i] = (in[i - start] > 0);
                }
                
                if (end == n_samples){
                    // avoid NaNs
                    for(unsigned int i = n_samples; i < n_samples_r; i++) {
                        out[offset + i] = 0;
                    }
                }
            }
            else {
                // just update my bit
                for(unsigned int i = start; i < end; i++){
                    out[offset + i] |= (uint64_t(in[i-start] > 0) << emb_bit);
                }
                
                // the rest of the els are already OK
            }
            return out;
        }
};


template<> inline void UnifracTaskBase<double>::embed_proportions_range(
        const std::vector<double> & in,
        unsigned int start,
        unsigned int end,
        unsigned int emb){
    embed_proportions_range_straight(embedded_proportions,in,start,end,emb);
}

template<> inline unsigned int UnifracTaskBase<double>::get_emb_els(
        unsigned int max_embs){
    return max_embs;
}

template<> inline void UnifracTaskBase<uint64_t>::embed_proportions_range(
        const std::vector<double> & in,
        unsigned int start,
        unsigned int end,
        unsigned int emb){
    embedded_proportions = embed_proportions_range_bool(embedded_proportions,
                                                        in,
                                                        start,
                                                        end,
                                                        emb);
}

template<> inline  unsigned int UnifracTaskBase<uint64_t>::get_emb_els(
        unsigned int max_embs){
    return (max_embs+63)/64;
}
    

/* void unifrac tasks
*
* all methods utilize the same function signature. that signature is as follows:
*
* dm_stripes vector<double> the stripes of the distance matrix being accumulated 
*      into for unique branch length
* dm_stripes vector<double> the stripes of the distance matrix being accumulated 
*      into for total branch length (e.g., to normalize unweighted unifrac)
* embedded_proportions <double*> the proportions vector for a sample, or rather
*      the counts vector normalized to 1. this vector is embedded as it is 
*      duplicated: if A, B and C are proportions for features A, B, and C, the
*      vector will look like [A B C A B C].
* length <double> the branch length of the current node to its parent.
* task_p <task_parameters*> task specific parameters.
*/

template<class TEmb>
class UnifracTask : public UnifracTaskBase<TEmb> {
    public:
        UnifracTask(su::StripeMap & _dm_stripes,
                    su::StripeMap & _dm_stripes_total,
                    unsigned int _max_embs,
                    su::task_parameters _task_p);
        
        virtual ~UnifracTask();
        
        virtual void run(unsigned int filled_embs,
                         const std::vector<double>  & lengths) = 0;
        
    protected:
        // Controls the size of inner loops in the calculation phase
        static const unsigned int step_size = 4;
        
        // Max embs are theoretically optimized for cache performance
        static const unsigned int RECOMMENDED_MAX_EMBS_STRAIGHT = 64-16;
        static const unsigned int RECOMMENDED_MAX_EMBS_BOOL = 64*32;
    
};
    
    

class UnifracUnweightedTask : public UnifracTask<uint64_t> {
    public:
        static const unsigned int RECOMMENDED_MAX_EMBS
            = UnifracTask<uint64_t>::RECOMMENDED_MAX_EMBS_BOOL;
        
        // Note: _max_emb MUST be multiple of 64
        UnifracUnweightedTask(su::StripeMap & _dm_stripes,
                                su::StripeMap & _dm_stripes_total,
                                unsigned int _max_embs,
                                su::task_parameters _task_p);
        
        virtual ~UnifracUnweightedTask();
        
        virtual void run(unsigned int filled_embs,
                            const std::vector<double> & lengths);
        
        void _run(unsigned int filled_embs,
                    const std::vector<double> & lengths);
    private:
        std::vector<double> sums; // temp buffer
};



class UnifracUnnormalizedWeightedTask : public UnifracTask<double> {
    public:
        static const unsigned int RECOMMENDED_MAX_EMBS
            = UnifracTask<double>::RECOMMENDED_MAX_EMBS_STRAIGHT;
        
        UnifracUnnormalizedWeightedTask(su::StripeMap & _dm_stripes,
                                        su::StripeMap & _dm_stripes_total,
                                        unsigned int _max_embs,
                                        su::task_parameters _task_p);
          
        virtual ~UnifracUnnormalizedWeightedTask();
      
        virtual void run(unsigned int filled_embs,
                            const std::vector<double> & lengths);
      
        void _run(unsigned int filled_embs,
                    const std::vector<double> & lengths);
          
    protected:
        // temp buffers
        std::vector<bool> zcheck;
        std::vector<double> sums;
};

}

#endif /* __UNIFRAC_TASK_H */
