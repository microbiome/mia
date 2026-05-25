/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <cstdlib>

#include "tree.h"
#include "unifrac_task.h"

using namespace su;

UnifracTaskVector::UnifracTaskVector(su::StripeMap & _dm_stripes,
                  const su::task_parameters _task_p)
    : dm_stripes(_dm_stripes)
    , start_idx(_task_p.start)
    , n_samples(_task_p.n_samples)
    , n_samples_r(((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK)
    , task_p(_task_p)
{
    
    // The buffer is only needed if the stripes are non-empty
    buf = (dm_stripes.is_empty(start_idx)) ? std::vector<double>() :
                                            std::vector<double>(n_samples_r*
                                                (task_p.stop-start_idx), 0.0);
    // Copy stripe values to buffer
    if (!buf.empty()){
        for(unsigned int stripe=start_idx; stripe < task_p.stop; stripe++) {
            std::vector<double> dm_stripe = dm_stripes.get(stripe);
            for(unsigned int k=0; k < dm_stripe.size(); k++) {
                buf[ (stripe-start_idx)*n_samples_r + k ] = dm_stripe[k];
            }
        }
    }
}

UnifracTaskVector::~UnifracTaskVector(){
    // If the buffer isn't empty, copy its values back to the stripes
    if (!buf.empty()){
        for(unsigned int stripe=start_idx; stripe < task_p.stop; stripe++) {
            std::vector<double> dm_stripe = dm_stripes.get(stripe);
            for(unsigned int k=0; k < dm_stripe.size(); k++) {
                dm_stripe[k] = buf[ (stripe-start_idx)*n_samples_r + k ];
            }
            dm_stripes.update(stripe, dm_stripe);
        }
    }
}



template<class TEmb>
UnifracTaskBase<TEmb>::UnifracTaskBase(su::StripeMap & _dm_stripes,
                                    su::StripeMap & _dm_stripes_total,
                                    unsigned int _max_embs,
                                    su::task_parameters _task_p)
    : dm_stripes(_dm_stripes,_task_p)
    , dm_stripes_total(_dm_stripes_total,_task_p)
    , task_p(_task_p)
    , max_embs(_max_embs)
{
    embedded_proportions = initialize_embedded(dm_stripes.n_samples_r,
                                                _max_embs);
}

template<class TEmb>
UnifracTaskBase<TEmb>::~UnifracTaskBase(){}

template<class TEmb>
std::vector<TEmb> UnifracTaskBase<TEmb>::initialize_embedded(
                                                const uint64_t  n_samples_r,
                                                unsigned int max_embs){
    uint64_t bsize = n_samples_r * get_emb_els(max_embs);
    return std::vector<TEmb>(bsize);
}

template<class TEmb>
void UnifracTaskBase<TEmb>::embed_proportions(const std::vector<double> & in,
                                                unsigned int emb){
    embed_proportions_range(in,0,dm_stripes.n_samples,emb);
}


template<class TEmb>
UnifracTask<TEmb>::UnifracTask(su::StripeMap & _dm_stripes,
                            su::StripeMap & _dm_stripes_total,
                            unsigned int _max_embs,
                            su::task_parameters _task_p)
    : UnifracTaskBase<TEmb>(_dm_stripes,
                            _dm_stripes_total,
                            _max_embs,
                            _task_p){}

template<class TEmb>
UnifracTask<TEmb>::~UnifracTask(){}



UnifracUnweightedTask::UnifracUnweightedTask(su::StripeMap & _dm_stripes,
                                            su::StripeMap & _dm_stripes_total,
                                            unsigned int _max_embs,
                                            su::task_parameters _task_p)
    : UnifracTask<uint64_t>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) 
{
    const unsigned int bsize = _max_embs*32;
    sums = std::vector<double>(bsize, 0.0);
}

UnifracUnweightedTask::~UnifracUnweightedTask(){}

void UnifracUnweightedTask::run(unsigned int filled_embs,
                                const std::vector<double> & lengths){
    _run(filled_embs, lengths);
}

void UnifracUnweightedTask::_run(unsigned int filled_embs,
                                const std::vector<double> & lengths){
    const uint64_t start_idx = this->task_p.start;
    const uint64_t stop_idx = this->task_p.stop;
    const uint64_t n_samples = this->task_p.n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
    
    const uint64_t step_size = UnifracUnweightedTask::step_size;
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size;
    
    const uint64_t filled_embs_els = filled_embs/64;
    const uint64_t filled_embs_rem = filled_embs%64; 
    
    const uint64_t filled_embs_els_round = (filled_embs+63)/64;
    
    // pre-compute sums of length elements, since they are likely to be accessed
    // many times
    // We will use a 8-bit map, to keep it small enough to keep in L1 cache
    for(uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
        for(uint64_t sub8=0; sub8<8; sub8++){
            const uint64_t emb8 = emb_el*8+sub8;
            
            std::vector<double> pl = std::vector<double>(8);
            
                
            uint64_t len_off = emb8*8;
            
            // compute all the combinations for this block (8-bits total)
            // psum[0] = 0.0   // +0*pl[0]+0*pl[1]+0*pl[2]+...
            // psum[1] = pl[0] // +0*pl[1]+0*pl[2]+...
            // psum[2] = pl[1] // +0*pl[0]+0*pl[2]+
            // psum[2] = pl[0] + pl[1]
            // ...
            // psum[255] = pl[1] +.. + pl[7] // + 0*pl[0]
            // psum[255] = pl[0] +pl[1] +.. + pl[7]
            for(uint64_t b8_i=0; b8_i<0x100; b8_i++){
                sums[(emb8<<8) + b8_i] =
                        (((b8_i >> 0) & 1) * lengths[len_off + 0]) +
                        (((b8_i >> 1) & 1) * lengths[len_off + 1]) +
                        (((b8_i >> 2) & 1) * lengths[len_off + 2]) +
                        (((b8_i >> 3) & 1) * lengths[len_off + 3]) +
                        (((b8_i >> 4) & 1) * lengths[len_off + 4]) +
                        (((b8_i >> 5) & 1) * lengths[len_off + 5]) +
                        (((b8_i >> 6) & 1) * lengths[len_off + 6]) +
                        (((b8_i >> 7) & 1) * lengths[len_off + 7]);
            }
        }
    }
    
    if (filled_embs_rem>0){ // add also the overflow elements
        const uint64_t emb_el=filled_embs_els;
        for(uint64_t sub8=0; sub8<8; sub8++){
            // we are summing we have enough buffer in sums
            const uint64_t emb8 = emb_el*8+sub8;
            
            // compute all the combinations for this block, set to 0 any past
            // the limit as above
            for(uint64_t b8_i=0; b8_i<0x100; b8_i++){
                double val= 0;
                for(uint64_t li=(emb8*8); li<filled_embs; li++){
                    val += ((b8_i >>  (li-(emb8*8))) & 1) * lengths[li];
                }
                sums[(emb8<<8) + b8_i] = val;
            }
            
        }
    }
    
    // point of thread
    for(uint64_t sk = 0; sk < sample_steps ; sk++){
        for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++){
            for(uint64_t ik = 0; ik < step_size ; ik++){
                // within-stripe index (0:n_samples-1)
                const uint64_t k = sk*step_size + ik; 
                //buffer index
                const uint64_t idx = (stripe-start_idx) * n_samples_r; 
                
                if (k>=n_samples) continue; // past the limit
                
                const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
                
                bool did_update = false;
                double my_stripe = 0.0;
                double my_stripe_total = 0.0;
                
                //Main calculation phase
                for(uint64_t emb_el=0; emb_el<filled_embs_els_round; emb_el++){
                    
                    const uint64_t offset = n_samples_r * emb_el;
                    
                    uint64_t sums_off = emb_el * 2048;
                    
                    uint64_t u1 = embedded_proportions[offset + k];
                    uint64_t v1 = embedded_proportions[offset + l1];
                    uint64_t o1 = u1 | v1;
                    
                    if (o1!=0) {  // zeros are prevalent
                        did_update=true;
                        uint64_t x1 = u1 ^ v1;
                        
                        // Use the pre-computed sums
                        // Each range of 8 lengths has already been pre-computed
                        // and stored in psum
                        // Since embedded_proportions packed format is in 64-bit
                        // format for performance reasons we need to add the 8
                        // sums using the four 8-bits for addressing inside psum
                        
                        my_stripe       += sums[sums_off + (x1 & 0xff)] + 
                            sums[sums_off + 0x100+((x1 >>  8) & 0xff)] +
                            sums[sums_off + 0x200+((x1 >> 16) & 0xff)] +
                            sums[sums_off + 0x300+((x1 >> 24) & 0xff)] +
                            sums[sums_off + 0x400+((x1 >> 32) & 0xff)] +
                            sums[sums_off + 0x500+((x1 >> 40) & 0xff)] +
                            sums[sums_off + 0x600+((x1 >> 48) & 0xff)] +
                            sums[sums_off + 0x700+((x1 >> 56)       )];
                        my_stripe_total += sums[sums_off + (o1 & 0xff)] +
                            sums[sums_off + 0x100+((o1 >>  8) & 0xff)] +
                            sums[sums_off + 0x200+((o1 >> 16) & 0xff)] +
                            sums[sums_off + 0x300+((o1 >> 24) & 0xff)] +
                            sums[sums_off + 0x400+((o1 >> 32) & 0xff)] +
                            sums[sums_off + 0x500+((o1 >> 40) & 0xff)] +
                            sums[sums_off + 0x600+((o1 >> 48) & 0xff)] +
                            sums[sums_off + 0x700+((o1 >> 56)       )];
                    }
                }
                
                if (did_update){
                    dm_stripes.buf[idx + k] += my_stripe;
                    dm_stripes_total.buf[idx + k] += my_stripe_total;
                }
            }
        }
    }
}




UnifracUnnormalizedWeightedTask::UnifracUnnormalizedWeightedTask(
                                            su::StripeMap & _dm_stripes,
                                            su::StripeMap & _dm_stripes_total,
                                            unsigned int _max_embs,
                                            su::task_parameters _task_p)
    : UnifracTask<double>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) 
{
    const unsigned int n_samples = this->task_p.n_samples;
    zcheck = std::vector<bool>(n_samples, 0);
    sums = std::vector<double>(n_samples, 0.0);
}

UnifracUnnormalizedWeightedTask::~UnifracUnnormalizedWeightedTask(){}

void UnifracUnnormalizedWeightedTask::run(unsigned int filled_embs,
                                        const std::vector<double> & lengths){
    _run(filled_embs, lengths);
}

void UnifracUnnormalizedWeightedTask::_run(unsigned int filled_embs,
                                        const std::vector<double> & lengths){
    const uint64_t start_idx = this->task_p.start;
    const uint64_t stop_idx = this->task_p.stop;
    const uint64_t n_samples = this->task_p.n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
    
    const uint64_t step_size = UnifracUnnormalizedWeightedTask::step_size;
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size;
    
    // check for zero values and pre-compute single column sums
    for(uint64_t k=0; k<n_samples; k++){
        bool all_zeros=true;
        double my_sum = 0.0;
        
        for(uint64_t emb=0; emb<filled_embs; emb++){
            const uint64_t offset = n_samples_r * emb;
            
            double u1 = embedded_proportions[offset + k];
            my_sum += u1*lengths[emb];
            all_zeros = all_zeros && (u1==0.0);
        }
        
        sums[k]     = my_sum;
        zcheck[k] = all_zeros;
    }
    
    // Main calculation phase
    for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++){
        for(uint64_t sk = 0; sk < sample_steps ; sk++){
            for(uint64_t ik = 0; ik < step_size ; ik++){
                
                // within-stripe index (0:n_samples-1)
                const uint64_t k = sk*step_size + ik;
                
                if (k>=n_samples) continue; // past the limit
                
                const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
                
                const bool allzero_k = zcheck[k];
                const bool allzero_l1 = zcheck[l1];
                
                if (allzero_k && allzero_l1) {
                    // nothing to do, would have to add 0
                }
                else {
                    double my_stripe;
                    
                    if (allzero_k || allzero_l1){
                        // one side has all zeros
                        // we can use the distributed property, and use the
                        // pre-computed values
                        
                        const uint64_t ridx = (allzero_k) ? l1 : k;
                        // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
                        // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
                        
                        my_stripe = sums[ridx];
                        
                    }
                    else {
                        // both sides non zero, use the explicit but slow
                        // approach
                        my_stripe = 0.0;
                        for(uint64_t emb=0; emb<filled_embs; emb++){
                            const uint64_t offset = n_samples_r * emb;
                            double u1 = embedded_proportions[offset + k];
                            double v1 = embedded_proportions[offset + l1];
                            double diff1 = u1 - v1;
                            double length = lengths[emb];
                            my_stripe     += fabs(diff1) * length;
                        } // for emb
                        
                    }
                    const uint64_t idx = (stripe-start_idx)*n_samples_r;
                    dm_stripes.buf[idx + k] += my_stripe;
                }
            } // for ik
        } // for stripe
    } // for sk
}


