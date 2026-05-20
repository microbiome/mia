/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC_TASKS
#define __UNIFRAC_TASKS 1

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
        
        // task specific arguments below
        double g_unifrac_alpha;      // generalized unifrac alpha
    };

    // Note: This adds a copy, which is suboptimal
    //       But was the easiest way to get a contiguous buffer
    //       And it does allow for fp32 compute, when desired
    
    class UnifracTaskVector {
    private:
      const su::task_parameters task_p;

    public:
      su::StripeMap & dm_stripes;
      const unsigned int start_idx;
      const unsigned int n_samples;
      const uint64_t  n_samples_r;
      std::vector<double> buf;

      UnifracTaskVector(su::StripeMap & _dm_stripes,
                        const su::task_parameters _task_p)
      : task_p(_task_p), dm_stripes(_dm_stripes)
      , start_idx(task_p.start), n_samples(task_p.n_samples)
      , n_samples_r(((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK) // round up
        //buf is just a new array with as many stripes as called for in task_p
        //n_samples_r tells us how many unifrac_blocks are required for n_samples.
        //Originally this was a null comparison, we might need to check what it does specifically 
      , buf((dm_stripes.is_empty(start_idx)) ?
                std::vector<double>() :
                std::vector<double>(n_samples_r*(task_p.stop-start_idx), 0.0)) // dm_stripes could be null, in which case keep it null
      {
        if (!buf.empty()) {
          //Initialize buffer to dm_stripe values
          for(unsigned int stripe=start_idx; stripe < task_p.stop; stripe++) {
             std::vector<double> dm_stripe = dm_stripes.get(stripe);
             //copy stripe to appropriate segment of buffer
             //The stripes themselves have n_samples elements,
             //but in the buffer each stripe gets n_samples_r elements?
             std::copy(std::begin(dm_stripe), std::end(dm_stripe),
                       std::begin(buf) + ((stripe-start_idx)*n_samples_r) );
           }
        }
      }
      
      //Destructor copies the buffer values back into dm_stripe
      ~UnifracTaskVector()
      {
        if (!buf.empty()) {
          for(unsigned int stripe=start_idx; stripe < task_p.stop; stripe++) {
             std::vector<double> vec = dm_stripes.get(stripe);
             std::copy( std::begin(buf) + ((stripe-start_idx)*n_samples_r),
                        std::begin(buf) + ((stripe-start_idx)*n_samples_r) + n_samples,
                        std::begin(vec) );
             dm_stripes.update(stripe, vec);
          }
        }
      }

    private:
      UnifracTaskVector() = delete;
      UnifracTaskVector operator=(const UnifracTaskVector&other) const = delete;
    };
    
    
    
    
    /***********************************************/   
    
    
    // Base task class to be shared by all tasks
    template<class TEmb>
    class UnifracTaskBase {
    public:
        //Two taskvectors for stripes and total
        UnifracTaskVector dm_stripes;
        UnifracTaskVector dm_stripes_total;
        
        su::task_parameters task_p;
        
        const unsigned int max_embs;
        std::vector<TEmb> embedded_proportions; //Continuous vector - each stripe has n_samples_r elements, for complex reasons?
        //Has at most max_embs stripes - when filled, results stored in task _run() and embeds cleared to continue
        
        UnifracTaskBase(su::StripeMap & _dm_stripes,
                        su::StripeMap & _dm_stripes_total,
                        unsigned int _max_embs,
                        su::task_parameters _task_p)
            : dm_stripes(_dm_stripes,_task_p),
              dm_stripes_total(_dm_stripes_total,_task_p),
              task_p(_task_p),
              max_embs(_max_embs),
              embedded_proportions(initialize_embedded(dm_stripes.n_samples_r,
                                                       _max_embs))
        {}
        
        virtual ~UnifracTaskBase() {}
        
        static unsigned int get_emb_els(unsigned int max_embs);
        
        static std::vector<TEmb> initialize_embedded(
                const uint64_t  n_samples_r,
                unsigned int max_embs )
        {
            uint64_t bsize = n_samples_r * get_emb_els(max_embs);
            return std::vector<TEmb>(bsize);
        }
        
        //Need to return a vector?
        void embed_proportions_range(
                const std::vector<double> & in,
                unsigned int start,
                unsigned int end,
                unsigned int emb);
        
        void embed_proportions(
                const std::vector<double> & in,
                unsigned int emb)
        {
            embed_proportions_range(in,0,dm_stripes.n_samples,emb);
        }
        
        
        
        //
        // ===== Internal, do not use directly =======
        //
        
        // Just copy from one buffer to another
        
        std::vector<double> embed_proportions_range_straight(
                                              std::vector<double> out,
                                              const std::vector<double> & in,
                                              unsigned int start,
                                              unsigned int end,
                                              unsigned int emb) const
        {
            const unsigned int n_samples  = dm_stripes.n_samples;
            const uint64_t n_samples_r  = dm_stripes.n_samples_r;
            const uint64_t offset = emb * n_samples_r;
            
            //Copy to stripe indicated by emb
            //Stripes are all contained in in/out in one mass
            //Start/end aren't necessarily the whole stripe?
            for(unsigned int i = start; i < end; i++) {
                out[offset + i] = in[i-start];
            }
            
            if (end==n_samples) {
                // avoid NaNs
                for(unsigned int i = n_samples; i < n_samples_r; i++) {
                    out[offset + i] = 0.0;
                }
            }
            return out;
        }
        
        
        
        // packed bool
        // Compute (in[:]>0) on each element, and store only the boolean bit.
        // The output values are stored in a multi-byte format, one bit per emb index,
        //    so it will likely take multiple passes to store all the values
        //
        // Note: assumes we are processing emb in increasing order, starting from 0
        
        //Only used with uint64_t
        std::vector<uint64_t> embed_proportions_range_bool(
                std::vector<uint64_t>  out,
                const std::vector<double> & in,
                unsigned int start,
                unsigned int end,
                unsigned int emb) const
        {
            
            const unsigned int n_packed = sizeof(uint64_t)*8;
            const unsigned int n_samples  = dm_stripes.n_samples;
            const uint64_t n_samples_r  = dm_stripes.n_samples_r;
            
            // The output values are stored in a multi-byte format, one bit per emb index
            // Compute the element to store the bit into, as well as whichbit in that element 
            unsigned int emb_block = emb/n_packed; // beginning of the element  block
            unsigned int emb_bit = emb%n_packed;   // bit inside the elements
            const uint64_t offset = emb_block * n_samples_r;
            
            if  (emb_bit == 0) {
                // assign for emb_bit==0, so it clears the other bits
                // assumes we processing emb in increasing order, starting from 0
                for(unsigned int i = start; i < end; i++) {            
                    out[offset + i] = (in[i - start] > 0);
                }
                
                if (end == n_samples) {
                    // avoid NaNs
                    for(unsigned int i = n_samples; i < n_samples_r; i++) {
                        out[offset + i] = 0;
                    }
                }
            } else {
                // just update my bit
                for(unsigned int i = start; i < end; i++) {
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
            unsigned int emb )
    {
        embedded_proportions = embed_proportions_range_straight(embedded_proportions,in,start,end,emb);
    }
    
    template<> inline unsigned int UnifracTaskBase<double>::get_emb_els(
            unsigned int max_embs )
    {
        return max_embs;
    }
    
    
    
    
    template<> inline void UnifracTaskBase<uint64_t>::embed_proportions_range(
            const std::vector<double> & in,
            unsigned int start,
            unsigned int end,
            unsigned int emb )
    {
        embedded_proportions = embed_proportions_range_bool(embedded_proportions,in,start,end,emb);
    }
    
    template<> inline  unsigned int UnifracTaskBase<uint64_t>::get_emb_els(
            unsigned int max_embs )
    {
        return (max_embs+63)/64;
    }
    
    
    
    
    
    
    /***********************************************/  
    
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
    protected:
        // Use one cache line on CPU
        // On GPU, sharing a cache line is actually a good thing
        static const unsigned int step_size = 16*4/sizeof(double);
        
    public:
        
        UnifracTask(su::StripeMap & _dm_stripes, su::StripeMap & _dm_stripes_total, unsigned int _max_embs, su::task_parameters _task_p)
            : UnifracTaskBase<TEmb>(_dm_stripes, _dm_stripes_total, _max_embs, _task_p) {}
        
        virtual ~UnifracTask() {}
        
        //Probably should return a vector?
        virtual void run(unsigned int filled_embs, std::vector<double> lengths) = 0;
        
    protected:
        static const unsigned int RECOMMENDED_MAX_EMBS_STRAIGHT = 128-16; // a little less to leave a bit of space of maxed-out L1
        // packed uses 32x less memory,so this should be 32x larger than straight... but there are additional structures, so use half of that
        static const unsigned int RECOMMENDED_MAX_EMBS_BOOL = 64*32;
        
    };
    
    /***********************************************/ 
    
    
    //Simplify all template stuff into doubles
    
    class UnifracUnweightedTask : public UnifracTask<uint64_t> {
    public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<uint64_t>::RECOMMENDED_MAX_EMBS_BOOL;
        
        // Note: _max_emb MUST be multiple of 64
        UnifracUnweightedTask(su::StripeMap & _dm_stripes, su::StripeMap & _dm_stripes_total, unsigned int _max_embs, su::task_parameters _task_p)
            : UnifracTask<uint64_t>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) 
            {
                const unsigned int bsize = _max_embs*32;
                sums = std::vector<double>(bsize, 0.0);
            }
        
        virtual ~UnifracUnweightedTask() {}
        
        virtual void run(unsigned int filled_embs, std::vector<double> lengths) {_run(filled_embs, lengths);}
        
        void _run(unsigned int filled_embs, std::vector<double> lengths);
    private:
        std::vector<double> sums; // temp buffer
    };


    /***********************************************/

    class UnifracNormalizedWeightedTask : public UnifracTask<double> {
      public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<double>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracNormalizedWeightedTask(su::StripeMap & _dm_stripes, su::StripeMap & _dm_stripes_total, unsigned int _max_embs, su::task_parameters _task_p)
        : UnifracTask<double>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p)
        {
          const unsigned int n_samples = this->task_p.n_samples;

          zcheck = std::vector<bool>(n_samples, 0);
          sums = std::vector<double>(n_samples, 0.0);
        }

        virtual ~UnifracNormalizedWeightedTask()
        {
        }

        virtual void run(unsigned int filled_embs, std::vector<double> lengths) {_run(filled_embs, lengths);}

        void _run(unsigned int filled_embs, std::vector<double> lengths);
      protected:
        // temp buffers
        std::vector<bool> zcheck;
        std::vector<double> sums;
    };

}

#endif
