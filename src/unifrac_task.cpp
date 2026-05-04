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

//for sleep
#include <windows.h>
#include <unistd.h>

#include "tree.h"
#include "unifrac_task.h"

void su::UnifracUnweightedTask::_run(unsigned int filled_embs, std::vector<double> lengths) {
    
    Rcpp::Rcout << "start UnifracUnweightedTask::_run\n";
    sleep(2);
    
    //Parameter finding
    
    //Task parameters determine stuff
    const uint64_t start_idx = this->task_p.start;
    const uint64_t stop_idx = this->task_p.stop;
    const uint64_t n_samples = this->task_p.n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
    
    /*
    // openacc only works well with local variables
    const uint64_t * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
    TFloat * const __restrict__ sums = this->sums;
     */
    
    std::vector<uint64_t> embedded_proportions = this->embedded_proportions;
    std::vector<double> dm_stripes_buf = this->dm_stripes.buf;
    std::vector<double> dm_stripes_total_buf = this->dm_stripes_total.buf;
    std::vector<double> sums = this->sums;
    
    Rcpp::Rcout << "buffers done\n";
    sleep(2);
    
    const uint64_t step_size = su::UnifracUnweightedTask::step_size;
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
    
    const uint64_t filled_embs_els = filled_embs/64;
    const uint64_t filled_embs_rem = filled_embs%64; 
    
    const uint64_t filled_embs_els_round = (filled_embs+63)/64;
    
    
    Rcpp::Rcout << "params done\n";
    sleep(2);
    
    Rcpp::Rcout << "start pre-compute\n";
    sleep(2);
    
    // pre-compute sums of length elements, since they are likely to be accessed many times
    // We will use a 8-bit map, to keep it small enough to keep in L1 cache
    for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
        for (uint64_t sub8=0; sub8<8; sub8++) {
            const uint64_t emb8 = emb_el*8+sub8;
            
            //TFloat * __restrict__ psum = &(sums[emb8<<8]);
            //const TFloat * __restrict__ pl   = &(lengths[emb8*8]);
            
            std::vector<double> psum   = std::vector<double>(256);
            std::vector<double> pl   = std::vector<double>(8);
            
            std::copy(  std::begin(sums) + (emb8<<8),
                        std::begin(sums) + (emb8<<8) + 256,
                        std::begin(psum) );
                
            std::copy(  std::begin(lengths) + (emb8*8),
                        std::begin(lengths) + (emb8*8) + 8,
                        std::begin(pl) );
            
            // compute all the combinations for this block (8-bits total)
            // psum[0] = 0.0   // +0*pl[0]+0*pl[1]+0*pl[2]+...
            // psum[1] = pl[0] // +0*pl[1]+0*pl[2]+...
            // psum[2] = pl[1] // +0*pl[0]+0*pl[2]+
            // psum[2] = pl[0] + pl[1]
            // ...
            // psum[255] = pl[1] +.. + pl[7] // + 0*pl[0]
            // psum[255] = pl[0] +pl[1] +.. + pl[7]
            for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
                psum[b8_i] = (((b8_i >> 0) & 1) * pl[0]) + (((b8_i >> 1) & 1) * pl[1]) + 
                    (((b8_i >> 2) & 1) * pl[2]) + (((b8_i >> 3) & 1) * pl[3]) +
                    (((b8_i >> 4) & 1) * pl[4]) + (((b8_i >> 5) & 1) * pl[5]) +
                    (((b8_i >> 6) & 1) * pl[6]) + (((b8_i >> 7) & 1) * pl[7]);
            }
            
            std::copy(std::begin(psum), std::end(psum),
                      std::begin(sums) + (emb8<<8));
        }
    }
    
    Rcpp::Rcout << "pre-compute done\n";
    sleep(2);
    Rcpp::Rcout << "start overflow elements\n";
    sleep(2);
    
    if (filled_embs_rem>0) { // add also the overflow elements
        const uint64_t emb_el=filled_embs_els;
        for (uint64_t sub8=0; sub8<8; sub8++) {
            // we are summing we have enough buffer in sums
            const uint64_t emb8 = emb_el*8+sub8;
            
            //TFloat * __restrict__ psum = &(sums[emb8<<8]);
            
            std::vector<double> psum   = std::vector<double>(256);
            std::copy(  std::begin(sums) + (emb8<<8),
                        std::begin(sums) + (emb8<<8) + 256,
                        std::begin(psum) );
                
            
            // compute all the combinations for this block, set to 0 any past the limit
            // as above
            for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
                double val= 0;
                for (uint64_t li=(emb8*8); li<filled_embs; li++) {
                    val += ((b8_i >>  (li-(emb8*8))) & 1) * lengths[li];
                }
                psum[b8_i] = val;
            }
            
            std::copy(std::begin(psum), std::end(psum),
                      std::begin(sums) + (emb8<<8));
        }
    }
    Rcpp::Rcout << "overflow elements done\n";
    sleep(2);
    
    //problem occurs here
    Rcpp::Rcout << "start stripes\n";
    sleep(2);
    // point of thread
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
        Rcpp::Rcout << "sk: " << sk << "\n";
        //sleep(2);
        
        for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
            
            Rcpp::Rcout << "    stripe: " << stripe << "\n";
            //sleep(2);
            
            const uint64_t idx = stripe-start_idx;
            
            std::vector<double> dm_stripe = this->dm_stripes.dm_stripes.get(idx);
            std::vector<double> dm_stripe_total = this->dm_stripes_total.dm_stripes.get(idx);
            
            for(uint64_t ik = 0; ik < step_size ; ik++) {
                
                Rcpp::Rcout << "        ik: " << ik << "\n";
                //sleep(2);
                
                const uint64_t k = sk*step_size + ik; // within-stripe index (0:n_samples-1)
                //const uint64_t idx = (stripe-start_idx) * n_samples_r; //n_samples_r seems to relate to continuous buffer shenanigans
                
                //TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
                //TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
                ////TFloat *dm_stripe = dm_stripes[stripe];
                ////TFloat *dm_stripe_total = dm_stripes_total[stripe];
                
                if (k>=n_samples) continue; // past the limit
                
                const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
                
                bool did_update = false;
                double my_stripe = 0.0;
                double my_stripe_total = 0.0;
                
                
                //This is the main calculation phase
                
                for (uint64_t emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                    Rcpp::Rcout << "            emb_el: " << emb_el << "\n";
                    //sleep(2);
                    const uint64_t offset = n_samples_r * emb_el;
                    //const TFloat * __restrict__ psum = &(sums[emb_el*0x800]);
                    
                    std::vector<double> psum   = std::vector<double>(2048);
                    std::copy(  std::begin(sums) + (emb_el * 2048),
                                std::begin(sums) + (emb_el * 2048) + 2048,
                                std::begin(psum) );
                    
                    uint64_t u1 = embedded_proportions[offset + k];
                    uint64_t v1 = embedded_proportions[offset + l1];
                    uint64_t o1 = u1 | v1;
                    
                    if (o1!=0) {  // zeros are prevalent
                        Rcpp::Rcout << "                update\n";
                        //sleep(2);
                        did_update=true;
                        uint64_t x1 = u1 ^ v1;
                        
                        // Use the pre-computed sums
                        // Each range of 8 lengths has already been pre-computed and stored in psum
                        // Since embedded_proportions packed format is in 64-bit format for performance reasons
                        //    we need to add the 8 sums using the four 8-bits for addressing inside psum
                        
                        my_stripe       += psum[              (x1 & 0xff)] + 
                            psum[0x100+((x1 >>  8) & 0xff)] +
                            psum[0x200+((x1 >> 16) & 0xff)] +
                            psum[0x300+((x1 >> 24) & 0xff)] +
                            psum[0x400+((x1 >> 32) & 0xff)] +
                            psum[0x500+((x1 >> 40) & 0xff)] +
                            psum[0x600+((x1 >> 48) & 0xff)] +
                            psum[0x700+((x1 >> 56)       )];
                        my_stripe_total += psum[              (o1 & 0xff)] +
                            psum[0x100+((o1 >>  8) & 0xff)] +
                            psum[0x200+((o1 >> 16) & 0xff)] +
                            psum[0x300+((o1 >> 24) & 0xff)] +
                            psum[0x400+((o1 >> 32) & 0xff)] +
                            psum[0x500+((o1 >> 40) & 0xff)] +
                            psum[0x600+((o1 >> 48) & 0xff)] +
                            psum[0x700+((o1 >> 56)       )];
                    }
                }
                
                if (did_update) {
                    dm_stripe[k]       += my_stripe;
                    dm_stripe_total[k] += my_stripe_total;
                }
                
            }
            
            this->dm_stripes.dm_stripes.update(idx, dm_stripe);
            this->dm_stripes_total.dm_stripes.update(idx, dm_stripe_total);
            
        }
    }
    Rcpp::Rcout << "stripes done\n";
    sleep(2);
    
    Rcpp::Rcout << "UnifracUnweightedTask::_run done\n";
    sleep(2);
}

// 
// 
// template<class TFloat>
// void SUCMP_NM::UnifracUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
// 
//     bool * const __restrict__ zcheck = this->zcheck;
//     TFloat * const __restrict__ sums = this->sums;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracUnnormalizedWeightedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     // check for zero values and pre-compute single column sums
// #ifdef _OPENACC
// #pragma acc parallel loop present(embedded_proportions,lengths,zcheck,sums)
// #else
// #pragma omp parallel for default(shared)
// #endif
//     for(uint64_t k=0; k<n_samples; k++) {
//             bool all_zeros=true;
//             TFloat my_sum = 0.0;
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat u1 = embedded_proportions[offset + k];
//                 my_sum += u1*lengths[emb];
//                 all_zeros = all_zeros && (u1==0.0);
//             }
// 
//             sums[k]     = my_sum;
//             zcheck[k] = all_zeros;
//     }
// 
// 
//     // now do the real compute
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracUnnormalizedWeightedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,lengths,zcheck,sums) async
// #else
//     // use dynamic scheduling due to non-homogeneity in the loop
// #pragma omp parallel for default(shared) schedule(dynamic,1)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//       for(uint64_t ik = 0; ik < step_size ; ik++) {
//        const uint64_t k = sk*step_size + ik;
// 
//        if (k>=n_samples) continue; // past the limit
// 
//        const bool zcheck_k = zcheck[sk]; // due to loop collapse in ACC, must load in here
// 
//        const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//        const bool allzero_k = zcheck[k];
//        const bool allzero_l1 = zcheck[l1];
// 
//        if (allzero_k && allzero_l1) {
//          // nothing to do, would have to add 0
//        } else {
//           TFloat my_stripe;
// 
//           if (allzero_k || allzero_l1) {
//             // one side has all zeros
//             // we can use the distributed property, and use the pre-computed values
// 
//             const uint64_t ridx = (allzero_k) ? l1 : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
//                                                 k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
// 
//             // keep reads in the same place to maximize GPU warp performance
//             my_stripe = sums[ridx];
// 
//           } else {
//             // both sides non zero, use the explicit but slow approach
//             my_stripe = 0.0;
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat u1 = embedded_proportions[offset + k];
//                 TFloat v1 = embedded_proportions[offset + l1];
//                 TFloat diff1 = u1 - v1;
//                 TFloat length = lengths[emb];
// 
//                 my_stripe     += fabs(diff1) * length;
//             } // for emb
// 
//           }
// 
//           const uint64_t idx = (stripe-start_idx)*n_samples_r;
//           TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//           //TFloat *dm_stripe = dm_stripes[stripe];
// 
//           // keep all writes in a single place, to maximize GPU warp performance
//           dm_stripe[k] += my_stripe;
//        } 
// 
//       } // for ik
//      } // for stripe
//     } // for sk
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracVawUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
//     const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracVawUnnormalizedWeightedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracVawUnnormalizedWeightedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
// #else
// #pragma omp parallel for default(shared) schedule(dynamic,1)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//       for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//         for(uint64_t ik = 0; ik < step_size ; ik++) {
//             const uint64_t k = sk*step_size + ik;
//             const uint64_t idx = (stripe-start_idx) * n_samples_r;
//             TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//             //TFloat *dm_stripe = dm_stripes[stripe];
// 
//             if (k>=n_samples) continue; // past the limit
// 
//             const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//             const TFloat m = sample_total_counts[k] + sample_total_counts[l1];
// 
//             TFloat my_stripe = dm_stripe[k];
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r*emb;
// 
//                 TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
//                 TFloat vaw = sqrt(mi * (m - mi));
// 
//                 if(vaw > 0) {
//                   TFloat u1 = embedded_proportions[offset + k];
//                   TFloat v1 = embedded_proportions[offset + l1];
//                   TFloat diff1 = fabs(u1 - v1);
//                   TFloat length = lengths[emb];
// 
//                    my_stripe += (diff1 * length) / vaw;
//                 }
//             }
// 
//             dm_stripe[k]     = my_stripe;
//         }
// 
//       }
//     }
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
//     TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
// 
//     bool * const __restrict__ zcheck = this->zcheck;
//     TFloat * const __restrict__ sums = this->sums;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracNormalizedWeightedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     // check for zero values and pre-compute single column sums
// #ifdef _OPENACC
// #pragma acc parallel loop present(embedded_proportions,lengths,zcheck,sums)
// #else
// #pragma omp parallel for default(shared)
// #endif
//     for(uint64_t k=0; k<n_samples; k++) {
//             bool all_zeros=true;
//             TFloat my_sum = 0.0;
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat u1 = embedded_proportions[offset + k];
//                 my_sum += u1*lengths[emb];
//                 all_zeros = all_zeros && (u1==0.0);
//             }
// 
//             sums[k]     = my_sum;
//             zcheck[k] = all_zeros;
//     }
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracNormalizedWeightedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths,zcheck,sums) async
// #else
//     // use dynamic scheduling due to non-homogeneity in the loop
// #pragma omp parallel for schedule(dynamic,1) default(shared)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//       for(uint64_t ik = 0; ik < step_size ; ik++) {
//        const uint64_t k = sk*step_size + ik;
// 
//        if (k>=n_samples) continue; // past the limit
// 
//        const bool zcheck_k = zcheck[sk]; // due to loop collapse in ACC, must load in here
// 
//        const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//        const bool allzero_k = zcheck[k];
//        const bool allzero_l1 = zcheck[l1];
// 
//        if (allzero_k && allzero_l1) {
//          // nothing to do, would have to add 0
//        } else {
//           const uint64_t idx = (stripe-start_idx) * n_samples_r;
//           TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//           TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
//           //TFloat *dm_stripe = dm_stripes[stripe];
//           //TFloat *dm_stripe_total = dm_stripes_total[stripe];
// 
//           // the totals can always use the distributed property
//           dm_stripe_total[k] += sums[k] + sums[l1];
//    
//           TFloat my_stripe;
// 
//           if (allzero_k || allzero_l1) {
//             // one side has all zeros
//             // we can use the distributed property, and use the pre-computed values
// 
//             const uint64_t ridx = (allzero_k) ? l1 : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
//                                                 k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
// 
//             // keep reads in the same place to maximize GPU warp performance
//             my_stripe = sums[ridx];
//           
//           } else {
//             // both sides non zero, use the explicit but slow approach
// 
//             my_stripe = 0.0;
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat u1 = embedded_proportions[offset + k];
//                 TFloat v1 = embedded_proportions[offset + l1];
//                 TFloat diff1 = u1 - v1;
//                 TFloat length = lengths[emb];
// 
//                 my_stripe     += fabs(diff1) * length;
//             }
// 
//           }
// 
//           // keep all writes in a single place, to maximize GPU warp performance
//           dm_stripe[k]       += my_stripe;
//        }
// 
//       } // for ik
//      } // for stripe
//     } // for sk
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracVawNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
//     const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
//     TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracVawNormalizedWeightedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracVawNormalizedWeightedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
// #else
// #pragma omp parallel for schedule(dynamic,1) default(shared)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//       for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//         for(uint64_t ik = 0; ik < step_size ; ik++) {
//             const uint64_t k = sk*step_size + ik;
//             const uint64_t idx = (stripe-start_idx) * n_samples_r;
//             TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
//             //TFloat *dm_stripe = dm_stripes[stripe];
//             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
// 
//             if (k>=n_samples) continue; // past the limit
// 
//             const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//             const TFloat m = sample_total_counts[k] + sample_total_counts[l1];
// 
//             TFloat my_stripe = dm_stripe[k];
//             TFloat my_stripe_total = dm_stripe_total[k];
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
//                 TFloat vaw = sqrt(mi * (m - mi));
// 
//                 if(vaw > 0) {
//                   TFloat u1 = embedded_proportions[offset + k];
//                   TFloat v1 = embedded_proportions[offset + l1];
//                   TFloat diff1 = fabs(u1 - v1);
//                   TFloat length = lengths[emb];
// 
//                   my_stripe += (diff1 * length) / vaw;
//                   my_stripe_total += ((u1 + v1) * length) / vaw;
//                 }
//             }
// 
//             dm_stripe[k]     = my_stripe;
//             dm_stripe_total[k]     = my_stripe_total;
// 
//         }
// 
//       }
//     }
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
//     TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
// 
//     const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracGeneralizedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracGeneralizedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
// #else
// #pragma omp parallel for schedule(dynamic,1) default(shared)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//       for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//         for(uint64_t ik = 0; ik < step_size ; ik++) {
//             const uint64_t k = sk*step_size + ik;
//             const uint64_t idx = (stripe-start_idx) * n_samples_r;
//             TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
//             //TFloat *dm_stripe = dm_stripes[stripe];
//             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
// 
//             if (k>=n_samples) continue; // past the limit
// 
//             const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//             TFloat my_stripe = dm_stripe[k];
//             TFloat my_stripe_total = dm_stripe_total[k];
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat u1 = embedded_proportions[offset + k];
//                 TFloat v1 = embedded_proportions[offset + l1];
//                 TFloat sum1 = u1 + v1;
// 
//                 if(sum1 != 0.0) { 
//                    TFloat length = lengths[emb];
//                    TFloat diff1 = fabs(u1 - v1);
//                    TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length; 
// 
//                    my_stripe += sum_pow1 * (diff1 / sum1); 
//                    my_stripe_total += sum_pow1; 
//                 }
//             }
// 
//             dm_stripe[k]     = my_stripe;
//             dm_stripe_total[k]     = my_stripe_total;
// 
//         }
// 
//       }
//     }
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracVawGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;
// 
//     // openacc only works well with local variables
//     const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
//     const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
//     const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
//     TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracVawGeneralizedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
//     // quick hack, to be finished
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracVawGeneralizedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
// #else
// #pragma omp parallel for schedule(dynamic,1) default(shared)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//       for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//         for(uint64_t ik = 0; ik < step_size ; ik++) {
//             const uint64_t k = sk*step_size + ik;
//             const uint64_t idx = (stripe-start_idx) * n_samples_r;
//             TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
//             //TFloat *dm_stripe = dm_stripes[stripe];
//             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
// 
//             if (k>=n_samples) continue; // past the limit
// 
//             const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//             const TFloat m = sample_total_counts[k] + sample_total_counts[l1];
// 
//             TFloat my_stripe = dm_stripe[k];
//             TFloat my_stripe_total = dm_stripe_total[k];
// 
// #pragma acc loop seq
//             for (uint64_t emb=0; emb<filled_embs; emb++) {
//                 const uint64_t offset = n_samples_r * emb;
// 
//                 TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
//                 TFloat vaw = sqrt(mi * (m - mi));
// 
//                 if(vaw > 0) {
//                   TFloat u1 = embedded_proportions[offset + k];
//                   TFloat v1 = embedded_proportions[offset + l1];
//                   TFloat length = lengths[emb];
// 
//                   TFloat sum1 = (u1 + v1) / vaw;
//                   TFloat sub1 = fabs(u1 - v1) / vaw;
//                   TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length;
// 
//                   my_stripe += sum_pow1 * (sub1 / sum1);
//                   my_stripe_total += sum_pow1;
//                 }
//             }
// 
//             dm_stripe[k]     = my_stripe;
//             dm_stripe_total[k]     = my_stripe_total;
// 
//         }
//       }
//     }
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
// template<class TFloat>
// void SUCMP_NM::UnifracVawUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
//     const uint64_t start_idx = this->task_p->start;
//     const uint64_t stop_idx = this->task_p->stop;
//     const uint64_t n_samples = this->task_p->n_samples;
//     const uint64_t n_samples_r = this->dm_stripes.n_samples_r;
// 
//     // openacc only works well with local variables
//     const uint32_t * const __restrict__ embedded_proportions = this->embedded_proportions;
//     const TFloat  * const __restrict__ embedded_counts = this->embedded_counts;
//     const TFloat  * const __restrict__ sample_total_counts = this->sample_total_counts;
//     TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
//     TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;
// 
//     const uint64_t step_size = SUCMP_NM::UnifracVawUnweightedTask<TFloat>::step_size;
//     const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
// 
//     const uint64_t filled_embs_els = (filled_embs+31)/32; // round up
// 
//     // point of thread
// #ifdef _OPENACC
//     const unsigned int acc_vector_size = SUCMP_NM::UnifracVawUnweightedTask<TFloat>::acc_vector_size;
// #pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
// #else
// #pragma omp parallel for schedule(dynamic,1) default(shared)
// #endif
//     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
//       for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
//         for(uint64_t ik = 0; ik < step_size ; ik++) {
//             const uint64_t k = sk*step_size + ik;
//             const uint64_t idx = (stripe-start_idx) * n_samples_r;
//             TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
//             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
//             //TFloat *dm_stripe = dm_stripes[stripe];
//             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
// 
//             if (k>=n_samples) continue; // past the limit
// 
//             const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
// 
//             TFloat my_stripe = dm_stripe[k];
//             TFloat my_stripe_total = dm_stripe_total[k];
// 
//             const TFloat m = sample_total_counts[k] + sample_total_counts[l1];
// 
// #pragma acc loop seq
//             for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
//                 const uint64_t offset_p = n_samples_r * emb_el;
//                 uint32_t u1 = embedded_proportions[offset_p + k];
//                 uint32_t v1 = embedded_proportions[offset_p + l1];
//                 uint32_t x1 = u1 ^ v1;
//                 uint32_t o1 = u1 | v1;
// 
//                 // embedded_proporions is packed
// 
// #pragma acc loop seq
//                 for (uint64_t ei=0; ei<32; ei++) {
//                    uint64_t emb=emb_el*32+ei;
//                    if (emb<filled_embs) {
//                      const uint64_t offset_c = n_samples_r * emb;
//                      TFloat mi = embedded_counts[offset_c + k] + embedded_counts[offset_c + l1];
//                      TFloat vaw = sqrt(mi * (m - mi));
// 
//                      // embedded_counts is not packed
// 
//                      if(vaw > 0) {
//                        TFloat length = lengths[emb];
//                        TFloat lv1 = length / vaw;
// 
//                        my_stripe +=       ((x1 >> ei) & 1)*lv1;
//                        my_stripe_total += ((o1 >> ei) & 1)*lv1;
//                      }
//                    }
//                 }
//             }
// 
//             dm_stripe[k]     = my_stripe;
//             dm_stripe_total[k]     = my_stripe_total;
// 
//         }
// 
//       }
//     }
// 
// #ifdef _OPENACC
//    // next iteration will use the alternative space
//    std::swap(this->embedded_proportions,this->embedded_proportions_alt);
// #endif
// }
// 
